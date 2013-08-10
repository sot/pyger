import sys
import os
import json

import numpy as np
import Ska.Numpy
from Chandra.Time import DateTime

import xija

from . import clogging


pkg_dir = os.path.dirname(os.path.abspath(__file__))
constraint_models = json.load(open(os.path.join(pkg_dir, 'constraint_models.json')))
logger = clogging.config_logger('pyger')

class ConstraintModel(object):
    def __init__(self, name, sim_inputs, limits, max_dwell_ksec):
        self.name = name
        self.sim_inputs = sim_inputs.get(name)
        self.limits = limits
        if name in constraint_models:
            self.msids = constraint_models[name]['msids']
            self.state_cols = constraint_models[name]['state_cols']
        self.max_dwell_ksec = max_dwell_ksec

    @property
    def n_sim(self):
        try:
            return len(self.dwells1)
        except AttributeError:
            return None

    def calc_dwell1_T0s(self, start):
        """Calculate the starting temperature vectors for the ensemble of pitch
        profiles at the given ``start`` time.  Creates sim_inputs[]['dwell1_T0s']
        values."""

        logger.info('{0}: calculating start temps for {1} dwells'.format(
            self.name.upper(), len(self.sim_inputs)))
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

            sim_input['dwell1_T0s'] = Ts[:, -1]

    def calc_dwells1(self, start, stop, times, pitches1, i_sims):
        self.start = start
        self.calc_dwell1_T0s(start)
        sim_inputs = self.sim_inputs
        dwells1 = []

        logger.info('{0}: simulating {1} dwells'.format(self.name.upper(), len(i_sims)))
        for i_sim, pitch1 in zip(i_sims, pitches1):
            sim_input = sim_inputs[i_sim]
            states1 = self.get_states1(start, stop, pitch1)
            Ts = self.calc_model(states1, times, sim_input['dwell1_T0s'])

            bad_idx = None
            constraint_name = 'none'
            for j, msid in enumerate(self.msids):
                if msid in self.limits:
                    bad_idxs = np.flatnonzero(Ts[j, :] >= self.limits[msid])
                    if len(bad_idxs) > 0 and (bad_idx is None or bad_idxs[0] < bad_idx):
                        bad_idx = bad_idxs[0]
                        constraint_name = msid

            ok_idx = -1 if (bad_idx is None) else max(bad_idx - 1, 0)
            dwell_dur = times[ok_idx] - times[0]
            dwells1.append((dwell_dur, pitch1, constraint_name, Ts[:, ok_idx]))

        self.dwells1 = np.rec.fromrecords(dwells1,
                                          dtype=[('dur', np.float64),
                                                 ('pitch', np.float64),
                                                 ('constraint_name', '|S10'),
                                                 ('T1', np.float64, (len(self.msids),))
                                                 ])

    def calc_dwell1_stats(self, pitch_bins):
        dwells1 = self.dwells1
        dwell1_stats = []
        for p0, p1 in zip(pitch_bins[:-1], pitch_bins[1:]):
            ok = (dwells1['pitch'] >= p0) & (dwells1['pitch'] < p1)
            dwells_ok = dwells1[ok]
            dwells_ok.sort(order='dur')
            n_dwells_ok = len(dwells_ok)
            dwell50 = dwells_ok[int(n_dwells_ok * 0.5)]
            dwell90 = dwells_ok[int(n_dwells_ok * 0.9)]
            dwell1_stats.append((p0, p1, (p0 + p1) / 2,
                                 dwell50['dur'], dwell90['dur'],
                                 ))
        self.dwell1_stats = np.rec.fromrecords(dwell1_stats,
                                               dtype=[('pitch_bin0', np.float64),
                                                      ('pitch_bin1', np.float64),
                                                      ('pitch', np.float64),
                                                      ('dur50', np.float64),
                                                      ('dur90', np.float64),
                                                      ])

