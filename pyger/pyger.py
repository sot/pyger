"""
Calculate Chandra dwell times given thermal constraints
"""
import sys
import os
import json
import time
from itertools import count, cycle
import cPickle as pickle
import re

import matplotlib
import numpy as np

from Chandra.Time import DateTime
import asciitable

from . import clogging
from . import nmass
from . import twodof
from . import characteristics

pkg_dir = os.path.dirname(os.path.abspath(__file__))
constraint_models = json.load(open(os.path.join(pkg_dir, 'constraint_models.json')))
__version__ = open(os.path.join(pkg_dir, 'VERSION')).read().strip()
logger = clogging.config_logger('pyger')

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

            sim_input['dwell1_T0s'] = Ts[:,-1]

    def calc_dwells1(self, start, stop, times, pitches1, i_sims):
        self.start = start
        self.calc_dwell1_T0s(start)
        sim_inputs = self.sim_inputs
        dwells1 = []
        n_times = len(times)

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


def hour_min_to_sec(hm):
    """Convert string in format hh:mm to seconds"""
    h, m = re.match(r'(\d{1,2}):(\d{2})', hm).groups()
    return 3600 * (int(h) + float(m) / 60.)


class ConstraintPline(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        ConstraintModel.__init__(self, 'pline', sim_inputs, limits, max_dwell_ksec)
        guidelines = asciitable.read(os.path.join(pkg_dir, 'pline_guidelines.dat'))
        # Convert guidelines to a more useful structure
        limits = []
        for guideline in guidelines:
            for colname in guidelines.dtype.names[2:]:
                cold_pitch_min, cold_pitch_max = re.search(r'dwell_(\d+)_(\d+)', colname).groups()
                limits.append(dict(warm_dwell=guideline['warm_dwell'] * 3600.,
                                   warm_pitch_max=guideline['warm_pitch_max'],
                                   cold_pitch_min=float(cold_pitch_min),
                                   cold_pitch_max=float(cold_pitch_max),
                                   cold_pitch_dur=hour_min_to_sec(guideline[colname])))
        self.limits = limits


    def calc_dwells1(self, start, stop, times, pitches1, i_sims):
        self.start = start
        self.calc_dwell1_T0s(start)

        dwells1 = []
        logger.info('{0}: simulating {1} dwells'.format(self.name.upper(), len(i_sims)))
        for i_sim, pitch1 in zip(i_sims, pitches1):
            sim_input = self.sim_inputs[i_sim]
            warm_dwell, warm_pitch_max = sim_input['dwell1_T0s']

            for limit in self.limits:
                if (warm_dwell >= limit['warm_dwell'] and
                    warm_pitch_max <= limit['warm_pitch_max'] and
                    pitch1 >= limit['cold_pitch_min'] and
                    pitch1 < limit['cold_pitch_max']):
                    dwell_dur = limit['cold_pitch_dur']
                    constraint_name = 'pline'
                    break
            else:
                dwell_dur = self.max_dwell_ksec * 1000
                constraint_name = 'none'

            dwells1.append((dwell_dur, pitch1, constraint_name, [warm_dwell, warm_pitch_max]))

        self.dwells1 = np.rec.fromrecords(dwells1,
                                          dtype=[('dur', np.float64),
                                                 ('pitch', np.float64),
                                                 ('constraint_name', '|S10'),
                                                 ('T1', np.float64, (2,))
                                                 ])

    def calc_dwell1_T0s(self, start):
        """Calculate the warm dwells for the ensemble of pitch
        profiles at the given ``start`` time.  Creates sim_inputs[]['dwell1_T0s']
        and sim_inputs[]['warm_pitch_max'] values.  Compare to this guideline table
        (values in hours):

        warm_dwell warm_pitch_max dwell_156_162  dwell_156_166  dwell_156_170 dwell_170_180
        30         80             9:20           7:50           4:20          0:00
        12         90             7:00           5:30           3:00          0:00
        10         90             6:15           5:10           2:40          0:00
        8          90             5:15           4:20           2:15          0:00
        8          110            5:00           4:15           2:00          0:00
        0          180            0:00           0:00           0:00          0:00
        """

        logger.info('{0}: calculating warm time for {1} dwells'.format(
            self.name.upper(), len(self.sim_inputs)))
        for sim_input in self.sim_inputs:
            warm_pitch_max = 0
            warm_dwell = 0
            for state in reversed(sim_input['states']):
                # States less than 10 minutes don't count toward warm/cold time
                dwell = state['tstop'] - state['tstart']
                if dwell < 600:
                    continue
                if state['pitch'] > 110:
                    break
                warm_pitch_max = max(state['pitch'], warm_pitch_max)
                warm_dwell += dwell
            sim_input['dwell1_T0s'] = [warm_dwell, warm_pitch_max]

class ConstraintMinusZ(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        self.pars = json.load(open(os.path.join(pkg_dir, 'pars_minusz.json')))
        ConstraintModel.__init__(self, 'minus_z', sim_inputs, limits, max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False, cache=True):
        states_dwell_ksec = (states[-1][1] - states[0][0]) / 1000.0
        max_dwell_ksec = max(self.max_dwell_ksec, states_dwell_ksec * 1.05)
        Ts = nmass.calc_model(self.pars, states, times, T0s, self.msids, cache=cache,
                              state_only=state_only, max_dwell_ksec=max_dwell_ksec)
        return Ts

    def get_states1(self, start, stop, pitch1):
        return np.rec.fromrecords(((start.secs, stop.secs, pitch1),),
                                  names=('tstart', 'tstop', 'pitch'))


class ConstraintPSMC(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6):
        self.pars = characteristics.model_par
        self.n_ccd = n_ccd
        self.powers = dict((x[0:3], x[3]) for x in characteristics.psmc_power)
        ConstraintModel.__init__(self, 'psmc', sim_inputs, limits, max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False, cache=False):
        pin0 = T0s[0]
        dea0 = T0s[1]
        T_pin, T_dea = twodof.calc_twodof_model(states, pin0, dea0, times,
                                                par=self.pars, dt=1000)
        return np.vstack([T_pin, T_dea])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, self.powers[self.n_ccd, 1, 1], pitch1, 75000)]
        return np.rec.fromrecords(states,
                                  names=('tstart', 'tstop', 'power', 'pitch', 'simpos'))

    def calc_dwell1_T0s(self, start):
        """Calculate the starting temperature vectors for the ensemble of pitch
        profiles at the given ``start`` time.  Creates sim_inputs[]['dwell1_T0s']
        values."""

        logger.info('{0}: calculating start temps for {1} dwells'.format(
            self.name.upper(), len(self.sim_inputs)))
        for sim_input in self.sim_inputs:
            sim_input['dwell1_T0s'] = [sim_input['T1s']['1pin1at'],
                                       sim_input['T1s']['1pdeaat']]

def plot_dwells1(constraint, plot_title=None, plot_file=None, figure=1):
    """Make a simple plot of the dwells and dwell statistics for the given ``constraint``.

    :param constraint: ConstraintModel object (e.g. constraints['all'])
    :param plot_title: plot title
    :param plot_file: output file for plot
    :param figure: matplotlib figure ID (default=1)
    """
    import matplotlib.pyplot as plt
    plt.rc("axes", labelsize=10, titlesize=12)
    plt.rc("xtick", labelsize=10)
    plt.rc("ytick", labelsize=10)
    plt.rc("font", size=10)
    plt.rc("legend", fontsize=10)

    dwells1 = constraint.dwells1
    dwell1_stats = constraint.dwell1_stats
    plt.figure(figure, figsize=(6,4))
    plt.clf()
    names = ('none', '1pdeaat', 'tcylaft6', 'tephin', 'pline')
    colors = ('r', 'g', 'k', 'c', 'b')
    for name, color in zip(names, colors):
        ok = dwells1['constraint_name'] == name
        plt.plot(dwells1['pitch'][ok], dwells1['dur'][ok] / 1000., '.' + color, markersize=3, label=name)
    plt.plot(dwell1_stats['pitch'], dwell1_stats['dur50'] / 1000., '-r')
    plt.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000., '-m')
    plt.grid()
    plt.title(plot_title or '')
    plt.xlabel('Pitch (deg)')
    plt.ylabel('Dwell (ksec)')
    plt.legend(loc='upper center')
    plt.ylim(0, constraint.max_dwell_ksec * 1.05)
    plt.subplots_adjust(bottom=0.12)
    if plot_file:
        logger.info('Writing constraint plot file {0}'.format(plot_file))
        plt.savefig(plot_file)


def merge_dwells1(constraints):
    """Merge the dwells in the ``constraints`` list, finding the shortest from among the
    constraints.

    :param constraints: list of ModelConstraint objects
    :returns: NumPy recarray of dwells with pitch, dur, constraint_name columns
    """
    dwells1 = []
    for i in range(constraints[0].n_sim):
        constraint = min(constraints, key=lambda x: x.dwells1['dur'][i])
        dwells1.append(tuple(constraint.dwells1[x][i] for x in ('pitch', 'dur', 'constraint_name')))
    dwells1 = np.rec.fromrecords(dwells1, names=('pitch', 'dur', 'constraint_name'))
    return dwells1


def calc_constraints(start='2011:001',
                     n_sim=500,
                     dt=1000.,
                     max_tephin=128.0,
                     max_tcylaft6=93.0,
                     max_1pdeaat=52.5,
                     n_ccd=6,
                     sim_file='sim_inputs.pkl',
                     max_dwell_ksec=200.,
                     min_pitch=45,
                     max_pitch=169,
                     bin_pitch=2,
                     constraint_models=('minus_z', 'psmc', 'pline'),
                     ):
    """
    Calculate allowed dwell times coming out of perigee given a set of
    constraint models.

    :param start: date at which to perform the constraint simulations
    :param n_sim: number of Monte-Carlo simulations of (pitch, starting condition) (default=500)
    :param dt: step size used in thermal model computations (default=1000 sec)
    :param max_tephin: TEPHIN planning limit (default=128 degF)
    :param max_tcylaft6: TCYLAFT6 planning limit (default=93 degF)
    :param max_1pdeaat: 1PDEAAT planning limit (default=52.5 degC)
    :param n_ccd: number of ACIS CCDs being used
    :param max_dwell_ksec: maximum allowed dwell time (default=200 ksec)
    :param sim_file: simulation inputs file from "pyger make" (default=sim_inputs.pkl)
    :param min_pitch: minimum pitch in simulations (default=45)
    :param max_pitch: maximum pitch in simulations (default=169)
    :param bin_pitch: pitch bin size for calculating stats (default=2)
    :param constraint_models: name of applicable constraint models (default=('minus_z', 'psmc', 'pline'))

    :returns: dict of computed constraint model objects
    """
    start = DateTime(start)
    stop = DateTime(start.secs + max_dwell_ksec * 1000)
    times = np.arange(start.secs, stop.secs, dt)
    try:
        sim_inputs = pickle.load(open(sim_file))
    except IOError:
        logger.error('ERROR: simulation inputs file "{0}" not found.'
                     '  Run "pyger make" or "pyger make --help".'.format(sim_file))
        sys.exit(1)
    i_sims = np.random.randint(len(sim_inputs['minus_z']), size=n_sim)
    pitches1 = np.random.uniform(min_pitch, max_pitch, size=n_sim)
    constraints = {}
    if 'minus_z' in constraint_models:
        constraints['minus_z'] = ConstraintMinusZ(sim_inputs,
                                                  limits=dict(tephin=FtoC(max_tephin),
                                                              tcylaft6=FtoC(max_tcylaft6)),
                                                  max_dwell_ksec=max_dwell_ksec)
    if 'psmc' in constraint_models:
        constraints['psmc'] = ConstraintPSMC(sim_inputs,
                                             limits={'1pdeaat': max_1pdeaat},
                                             max_dwell_ksec=max_dwell_ksec,
                                             n_ccd=n_ccd)
    if 'pline' in constraint_models:
        constraints['pline'] = ConstraintPline(sim_inputs,
                                               limits=None,
                                               max_dwell_ksec=max_dwell_ksec)

    constraints_list = [constraints[x] for x in constraint_models]
    pitch_bins = np.arange(min_pitch, max_pitch, bin_pitch)
    for constraint in constraints_list:
        constraint.calc_dwells1(start, stop, times, pitches1, i_sims)
        constraint.calc_dwell1_stats(pitch_bins)

    dwells1 = merge_dwells1(constraints_list)
    constraints['all'] = ConstraintModel('all', sim_inputs, limits=None,
                                         max_dwell_ksec=max_dwell_ksec)
    constraints['all'].dwells1 = dwells1
    constraints['all'].calc_dwell1_stats(pitch_bins)
    constraints['all'].start = start

    return constraints
