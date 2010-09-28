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

# sys.path.append('/proj/sot/ska/share/psmc')
import characteristics
import twodof

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
    
    parser.add_option("--start",
                      default='2011:001',
                      help="Start time")
    parser.add_option("--n-sim",
                      default=500,
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
    parser.add_option("--1pdeaat",
                      default=52.5,
                      type='float',
                      help="1PDEAAT planning limit")
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


class ConstraintModel(object):
    def calc_dwell1_T0s(self, start):
        """Calculate the starting temperature vectors for the ensemble of pitch
        profiles at the given ``start`` time.  Creates sim_inputs[]['dwell1_T0s']
        values."""

        print 'Calculating start temps for {0} dwells ...'.format(len(self.sim_inputs))
        for sim_input in self.sim_inputs:
            states = sim_input['states']
            state_cols = states[0].keys()
            time_adj = start.secs - states[-1]['tstop']
            for state in states:
                state['tstart'] += time_adj
                state['tstop'] += time_adj

            np_states = np.rec.fromrecords([[state[col] for col in state_cols] for state in states],
                                           names=state_cols)

            times = np.array([sim_input['tstop'] + time_adj - 1])
            T0s = np.array([sim_input['T0s'][x] for x in self.msids])
            Ts = self.calc_model(np_states, times, T0s, state_only=True)

            sim_input['dwell1_T0s'] = Ts[:,-1]

    def calc_dwells1(self, start, stop, times, pitches1, i_sims):
        sim_inputs = self.sim_inputs
        dwells1 = []
        n_times = len(times)

        print 'Simulating {0} {1} dwells ...'.format(len(i_sims), len(pitches1))
        for i_sim, pitch1 in zip(i_sims, pitches1):
            sim_input = sim_inputs[i_sim]
            states1 = np.rec.fromrecords(((start.secs, stop.secs, pitch1),),
                                         names=('tstart', 'tstop', 'pitch'))

            #Ts = nmass.calc_model(self.pars, states1, times, sim_input['dwell1_T0s'], self.msids, 
            #                      cache=True, max_dwell_ksec=self.max_dwell_ksec)

            Ts = self.calc_model(states1, times, sim_input['dwell1_T0s'])

            ok = np.ones(len(times), np.bool)
            for j, msid in enumerate(self.msids):
                if msid in self.limits:
                    ok &= Ts[j, :] < self.limits[msid]
            try:
                i_ok = np.flatnonzero(ok)[-1]
                dwell_dur = times[i_ok] - times[0]
            except IndexError:
                i_ok = 0
                dwell_dur = 0.
            dwells1.append((dwell_dur, pitch1, Ts[:, i_ok]))

        self.dwells1 = np.rec.fromrecords(dwells1,
                                          dtype=[('dur', np.float64),
                                                 ('pitch', np.float64),
                                                 ('T1', np.float64, (len(self.msids),))
                                                 ])

    def calc_dwell1_stats(self, pitch_bins):
        dwell1_stats = []
        n_msids = len(self.msids)
        for p0, p1 in zip(pitch_bins[:-1], pitch_bins[1:]):
            ok = (self.dwells1['pitch'] >= p0) & (self.dwells1['pitch'] < p1)
            dwells_ok = self.dwells1[ok]
            dwells_ok.sort(order='dur')
            n_dwells_ok = len(dwells_ok)
            dwell50 = dwells_ok[int(n_dwells_ok * 0.5)]
            dwell90 = dwells_ok[int(n_dwells_ok * 0.9)]
            dwell1_stats.append(((p0 + p1) / 2,
                                 dwell50['dur'], dwell90['dur'],
                                 dwell50['T1'], dwell90['T1'],
                                 ))
        self.dwell1_stats = np.rec.fromrecords(dwell1_stats,
                                               dtype=[('pitch', np.float64),
                                                      ('dur50', np.float64),
                                                      ('dur90', np.float64),
                                                      ('T1_50', np.float64, (n_msids,)),
                                                      ('T1_90', np.float64, (n_msids,)),
                                                      ])

    def plot_dwells1(self):
        plt.figure(1, figsize=(6,4))
        plt.clf()
        plt.plot(self.dwells1['pitch'], self.dwells1['dur'], '.', markersize=1.0)
        plt.plot(self.dwell1_stats['pitch'], self.dwell1_stats['dur50'])
        plt.plot(self.dwell1_stats['pitch'], self.dwell1_stats['dur90'])


class ConstraintMinusZ(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        self.name = 'minus_z'
        self.pars = json.load(open('minusz/pars_minusz.json'))
        self.sim_inputs = sim_inputs[self.name]
        self.limits = limits
        self.max_dwell_ksec = max_dwell_ksec
        self.msids = sorted(self.sim_inputs[0]['msids'])
        
    def calc_model(self, states, times, T0s, state_only=False, cache=True):
        Ts = nmass.calc_model(self.pars, states, times, T0s, self.msids, cache=cache,
                              state_only=state_only, max_dwell_ksec=self.max_dwell_ksec)
        return Ts


if __name__ == '__main__':
    opt, args = get_options()

    start = DateTime(opt.start)
    stop = DateTime(start.secs + opt.max_dwell_ksec * 1000)
    times = np.arange(start.secs, stop.secs, opt.dt)
    pitch_bins = np.arange(45, 170.1, 2)
    sim_inputs = pickle.load(open(opt.sim_inputs))
    i_sims = np.random.randint(len(sim_inputs['minus_z']), size=opt.n_sim)
    pitches1 = np.random.uniform(45., 169., size=opt.n_sim)

    minus_z = ConstraintMinusZ(sim_inputs,
                               limits=dict(tephin=FtoC(opt.tephin),
                                           tcylaft6=FtoC(opt.tcylaft6)),
                               max_dwell_ksec=opt.max_dwell_ksec)
    minus_z.calc_dwell1_T0s(start)
    minus_z.calc_dwells1(start, stop, times, pitches1, i_sims)
    minus_z.calc_dwell1_stats(pitch_bins)
    minus_z.plot_dwells1()
