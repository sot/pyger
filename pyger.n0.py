import sys
import json
import time
from itertools import cycle
import cPickle as pickle

from Chandra.Time import DateTime
import matplotlib.pyplot as plt
import numpy as np
from Ska.Matplotlib import plot_cxctime

sys.path.append('/proj/sot/ska/analysis/thermal/nmass')
import nmass

def CtoF(cs):
    return [c * 1.8 + 32 for c in cs]

def FtoC(cs):
    return [(c-32) / 1.8 for c in cs]

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("--model",
                      default='minusz',
                      help="Model to predict")
    parser.add_option("--start",
                      default='2009:100',
                      help="Start time")
    parser.add_option("--n-sim",
                      default=2000,
                      type='int',
                      help="Number of simulation points")
    parser.add_option("--dt",
                      default=328,
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
    parser.add_option("--sim-inputs",
                      default="sim_inputs.pkl",
                      help="Simulation inputs pickle file")
    parser.add_option("--make-plots",
                      action="store_true",
                      help="Make plots")

    (opt, args) = parser.parse_args()
    return (opt, args)

opt, args = get_options()
n_sim = opt.n_sim
limits = dict(tephin=opt.tephin, tcylaft6=opt.tcylaft6)

sim_inputs = pickle.load(open(opt.sim_inputs))

start = DateTime(opt.start)
stop = DateTime(start.secs + 170000.)

pars = json.load(open('{0}/pars_{0}.json'.format(opt.model)))

print 'Calculating sim inputs ...',
for sim_input in sim_inputs:
    states = sim_input['states']
    dt = start.secs - states[-1]['tstop']
    for state in states:
        state['tstart'] += dt
        state['tstop'] += dt

    np_states = np.rec.fromrecords([(state['tstart'], state['tstop'], state['pitch']) for state in states],
                                names=('tstart', 'tstop', 'pitch'))

    times = np.arange(sim_input['tstart']+dt, sim_input['tstop']+dt, 328)
    msids = sorted(sim_input['msids'])
    T0s = np.array([sim_input['T0s'][x] for x in msids])
    Ts = nmass.calc_model(pars, np_states, times, T0s, msids, dt=328.0)
    sim_input['Tsim0'] = Ts[:,-1]
print 'Done'

########################################################

dwells = []
pitches = np.random.uniform(45., 169., size=n_sim)
i_sims = np.random.randint(len(sim_inputs), size=n_sim)
times = np.arange(start.secs, stop.secs, opt.dt)

time0 = time.time()
print 'Simulating {0} dwells ...'.format(n_sim),
for i_sim, pitch in zip(i_sims, pitches):
    sim_input = sim_inputs[i_sim]
    states = np.rec.fromrecords(((start.secs, stop.secs, pitch),),
                                names=('tstart', 'tstop', 'pitch'))
    Ts = nmass.calc_model(pars, states, times, sim_input['Tsim0'], msids, dt=opt.dt)
    ok = np.ones(len(times), np.bool)
    for j, msid in enumerate(msids):
        if msid in limits:
            ok &= Ts[j, :] < (limits[msid] - 32) / 1.8
    try:
        dwell_ksec = (times[ok][-1] - times[0]) / 1000.
        dwells.append(dwell_ksec)
    except IndexError:
        dwells.append(0)
print 'Finished in {0:.1f} secs'.format(time.time() - time0)

if 1:
    plt.figure(1, figsize=(6,4))
    plt.clf()
    plt.plot(pitches, dwells, '.', markersize=0.7)
    plt.show()

