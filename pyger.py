import sys
import json
import time
from itertools import count
import cPickle as pickle

from Chandra.Time import DateTime
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime
import Ska.Numpy

sys.path.append('/proj/sot/ska/analysis/thermal/nmass')
import nmass

def CtoF(cs):
    try:
        return [c * 1.8 + 32 for c in cs]
    except TypeError:
        return cs * 1.8 + 32

def FtoC(cs):
    try:
        return [(c-32) / 1.8 for c in cs]
    except TypeError:
        return (cs-32) / 1.8
        

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("--model",
                      default='minusz',
                      help="Model to predict")
    parser.add_option("--start",
                      default='2011:001',
                      help="Start time")
    parser.add_option("--n-sim",
                      default=2000,
                      type='int',
                      help="Number of simulation points")
    parser.add_option("--dt",
                      default=1000,
                      type='float',
                      help="Delta time for sims (sec)")
    parser.add_option("--tephin",
                      default=123.0,
                      type='float',
                      help="TEPHIN planning limit")
    parser.add_option("--tcylaft6",
                      default=93.0,
                      type='float',
                      help="TCYLAFT6 planning limit")
    parser.add_option("--max-dwell-ksec",
                      default=170.,
                      type='float',
                      help="Max allowed dwell time (ksec)")
    parser.add_option("--sim-inputs",
                      default="sim_inputs.pkl",
                      help="Simulation inputs pickle file")
    parser.add_option("--make-plots",
                      action="store_true",
                      help="Make plots")

    (opt, args) = parser.parse_args()
    return (opt, args)

def calc_sim_inputs(pars, sim_inputs_file, start):
    """Calculate the starting temperature vectors for the ensemble of pitch
    profiles at the given ``start`` time.  Creates sim_inputs[]['T0']
    values."""
    
    sim_inputs = pickle.load(open(sim_inputs_file))
    time0 = time.time()
    print 'Calculating sim inputs ...'
    for sim_input in sim_inputs:
        states = sim_input['states']
        time_adj = start.secs - states[-1]['tstop']
        for state in states:
            state['tstart'] += time_adj
            state['tstop'] += time_adj

        np_states = np.rec.fromrecords([(state['tstart'], state['tstop'], state['pitch']) for state in states],
                                    names=('tstart', 'tstop', 'pitch'))

        times = np.array([sim_input['tstop'] + time_adj - 1])
        T0s = np.array([sim_input['T0s'][x] for x in msids])
        Ts = nmass.calc_model(pars, np_states, times, T0s, msids, cache=True, state_only=True)
        sim_input['T0'] = Ts[:,-1]
    print 'Done in {0} secs'.format(time.time() - time0)
    return sim_inputs

def calc_dwells1(pars, sim_inputs, n_sim, start, stop, msids, times):
    dwells = []
    pitches1 = np.random.uniform(45., 169., size=n_sim)
    i_sims = np.random.randint(len(sim_inputs), size=n_sim)
    n_times = len(times)

    time0 = time.time()
    print 'Simulating {0} dwells ...'.format(n_sim)
    for i_sim, pitch1 in zip(i_sims, pitches1):
        sim_input = sim_inputs[i_sim]
        states1 = np.rec.fromrecords(((start.secs, stop.secs, pitch1),),
                                     names=('tstart', 'tstop', 'pitch'))
        Ts = nmass.calc_model(pars, states1, times, sim_input['T0'], msids, dt=opt.dt,
                              cache=True, max_dwell_ksec=opt.max_dwell_ksec)
        ok = np.ones(len(times), np.bool)
        for j, msid in enumerate(msids):
            if msid in limits:
                ok &= Ts[j, :] < limits[msid]
        try:
            i_ok = np.flatnonzero(ok)[-1]
            dwell_dur = times[i_ok] - times[0]
        except IndexError:
            i_ok = 0
            dwell_dur = 0.
        dwells.append((dwell_dur, pitch1, Ts[:, i_ok]))

    print 'Finished in {0:.1f} secs'.format(time.time() - time0)
    return np.rec.fromrecords(dwells, dtype=[('dur', np.float64),
                                             ('pitch', np.float64),
                                             ('T1', np.float64, (len(msids),))
                                             ])

def calc_med_dwells(dwells, pitch_bins):
    med_dwells = []
    for p0, p1 in zip(pitch_bins[:-1], pitch_bins[1:]):
        ok = (dwells['pitch'] >= p0) & (dwells['pitch'] < p1)
        dwells_ok = dwells[ok]
        dwells_ok.sort(order='dur')
        md = dwells_ok[len(dwells[ok]) // 2]
        med_dwells.append((md['dur'], (p0 + p1) / 2, md['T1']))
    return np.rec.fromrecords(med_dwells, dtype=[('dur', np.float64),
                                                 ('pitch', np.float64),
                                                 ('T1', np.float64, (len(msids),))
                                                 ])

def calc_dwell2_dwell3(times2, dwells1, pitches2, pitches3, med_dwells1):
    dwell_pitches3 = Ska.Numpy.interpolate(med_dwells1['dur'], med_dwells1['pitch'], pitches3)
    sim_outputs = []
    n_times2 = len(times2)

    for dwell1, pitch2, pitch3, dwell_pitch3, i in zip(
        dwells1, pitches2, pitches3, dwell_pitches3, count()):

        # Get the starting temp vector 
        # i = Ska.Numpy.interpolate(np.arange(len(med_dwells1)), med_dwells1['pitch'], [pitch1])[0]
        # T1 = med_dwells[i]

        pitch1 = dwell1['pitch']
        T1 = dwell1['T1']
        states2 = np.rec.fromrecords(((start.secs, stop.secs, pitch2),),
                                     names=('tstart', 'tstop', 'pitch'))
        T2s = nmass.calc_model(pars, states2, times2, T1, msids, dt=opt.dt,
                               cache=True, max_dwell_ksec=opt.max_dwell_ksec+5)

        tstop3 = start.secs + dwell_pitch3
        states3 = np.rec.fromrecords(((start.secs, tstop3, pitch3),),
                                     names=('tstart', 'tstop', 'pitch'))
        for i2 in reversed(range(n_times2)):
            T2 = T2s[:, i2]
            times3 = np.array([tstop3])
            T3s = nmass.calc_model(pars, states3, times3, T2, msids, state_only=True,
                               cache=True, max_dwell_ksec=opt.max_dwell_ksec+5)

            ok = True
            for j, msid in enumerate(msids):
                if msid in limits:
                    ok &= T3s[j, -1] < limits[msid]
            if not ok:
                i2 += 1
                break

        dwell2 = (times2[i2] - start.secs  if i2 < n_times2 else 500000.)
        print '{0:04d} {1:4.0f} {2:4.0f} {3:4.0f}'.format(i, pitch1, pitch2, dwell2/1000)
        sim_outputs.append((pitch1, pitch2, pitch3, dwell2/1000.))

    return np.rec.fromrecords(sim_outputs, names=('pitch1', 'pitch2', 'pitch3', 'dwell2'))

def calc_med_dwell2(sim_outputs, pitch_bins):
    p1s = sim_outputs['pitch1']
    p2s = sim_outputs['pitch2']
    dwells2 = sim_outputs['dwell2']
    n_bins = len(pitch_bins)
    med_dwell2 = np.ndarray((n_bins-1, n_bins-1))
    for i1, p1_0, p1_1 in zip(count(), pitch_bins[:-1], pitch_bins[1:]):
        ok1 = (p1s >= p1_0) & (p1s < p1_1)
        for i2, p2_0, p2_1 in zip(count(), pitch_bins[:-1], pitch_bins[1:]):
            ok2 = (p2s >= p2_0) & (p2s < p2_1)
            med_dwell2[i1, i2] = np.median(dwells2[ok1 & ok2])
    return med_dwell2

def plot_dwells1(dwells1, med_dwells1):
    plt.figure(1, figsize=(6,4))
    plt.clf()
    plt.plot(dwells1['pitch'], dwells1['dur'], '.', markersize=1.0)
    plt.plot(med_dwells1['pitch'], med_dwells1['dur'])
    plt.show()

def plot_dwells2(sim_outputs):
    plt.figure(2, figsize=(6,4))
    plt.clf()
    plt.scatter(sim_outputs['pitch1'], sim_outputs['pitch2'], c=sim_outputs['dwell2'],
                s=np.where(sim_outputs['dwell2'] > 400, 2, 10),
                edgecolors='none', vmin=0, vmax=170
                )
    plt.xlabel('Observation pitch')
    plt.ylabel('Cooling pitch')
    plt.title('Cooling times (blue=0, red=170 ksec)')
    plt.xlim(40,180)
    plt.ylim(40,180)

opt, args = get_options()

start = DateTime(opt.start)
stop = DateTime(start.secs + opt.max_dwell_ksec * 1000)
times = np.arange(start.secs, stop.secs, opt.dt)
limits = dict(tephin=FtoC(opt.tephin), tcylaft6=FtoC(opt.tcylaft6))
pitch_bins = np.arange(45, 170.1, 5)

pars = json.load(open('{0}/pars_{0}.json'.format(opt.model)))
msids = sorted(pars)

sim_inputs = calc_sim_inputs(pars, opt.sim_inputs, start)
dwells1 = calc_dwells1(pars, sim_inputs, opt.n_sim, start, stop, msids, times)
med_dwells1 = calc_med_dwells(dwells1, pitch_bins)

pitches2 = np.random.uniform(45., 169., size=opt.n_sim)
pitches3 = dwells1['pitch']
sim_outputs = calc_dwell2_dwell3(times, dwells1, pitches2, pitches3, med_dwells1)

med_dwell2 = calc_med_dwell2(sim_outputs, pitch_bins)
plot_dwells2(sim_outputs)
