"""
Calculate Chandra dwell times given thermal constraints
"""
import sys
import os
import json
import cPickle as pickle
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import Ska.Numpy
from Chandra.Time import DateTime
import asciitable

import xija

from . import clogging
from .base import ConstraintModel, pkg_dir, constraint_models, logger

# pkg_dir = os.path.dirname(os.path.abspath(__file__))
__version__ = open(os.path.join(pkg_dir, 'VERSION')).read().strip()


def CtoF(cs):
    try:
        return [c * 1.8 + 32 for c in cs]
    except TypeError:
        return cs * 1.8 + 32


def FtoC(cs):
    try:
        return [(c - 32) / 1.8 for c in cs]
    except TypeError:
        return (cs - 32) / 1.8


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
        model_spec = os.path.join(pkg_dir, 'minusz_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'minus_z', sim_inputs, limits,
                                 max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False):
        model = xija.ThermalModel('minus_z', start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        state_times = np.array([states['tstart'], states['tstop']])
        model.comp['pitch'].set_data(states['pitch'], state_times)
        model.comp['eclipse'].set_data(False)
        model.comp['tcylaft6'].set_data(T0s[0])
        model.comp['tcylfmzm'].set_data(T0s[1])
        model.comp['tephin'].set_data(T0s[2])
        model.comp['tfssbkt1'].set_data(T0s[3])
        model.comp['tmzp_my'].set_data(T0s[4])

        model.make()
        model.calc()

        T_tcylaft6 = Ska.Numpy.interpolate(model.comp['tcylaft6'].mvals,
                                           xin=model.times, xout=times, sorted=True)
        T_tcylfmzm = Ska.Numpy.interpolate(model.comp['tcylfmzm'].mvals,
                                           xin=model.times, xout=times, sorted=True)
        T_tephin = Ska.Numpy.interpolate(model.comp['tephin'].mvals,
                                         xin=model.times, xout=times, sorted=True)
        T_tfssbkt1 = Ska.Numpy.interpolate(model.comp['tfssbkt1'].mvals,
                                           xin=model.times, xout=times, sorted=True)
        T_tmzp_my = Ska.Numpy.interpolate(model.comp['tmzp_my'].mvals,
                                          xin=model.times, xout=times, sorted=True)

        return np.vstack([T_tcylaft6, T_tcylfmzm, T_tephin, T_tfssbkt1, T_tmzp_my])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintMinusYZ(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, 'minusyz_model.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'minus_yz', sim_inputs, limits,
                                 max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False):
        model = xija.ThermalModel('minus_yz', start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        state_times = np.array([states['tstart'], states['tstop']])
        model.comp['pitch'].set_data(states['pitch'], state_times)
        model.comp['eclipse'].set_data(False)
        model.comp['pmtank3t'].set_data(T0s[0])
        model.comp['tmzp_my'].set_data(T0s[1])
        model.comp['tephin'].set_data(T0s[2])
        model.comp['tcylaft6'].set_data(T0s[3])
        model.comp['pseudo_0'].set_data(T0s[3] - 4)
        model.comp['pseudo_1'].set_data(T0s[0])

        model.make()
        model.calc()

        T_pmtank3t = Ska.Numpy.interpolate(model.comp['pmtank3t'].mvals,
                                           xin=model.times, xout=times, sorted=True)
        T_tmzp_my = Ska.Numpy.interpolate(model.comp['tmzp_my'].mvals,
                                          xin=model.times, xout=times, sorted=True)
        T_tephin = Ska.Numpy.interpolate(model.comp['tephin'].mvals,
                                         xin=model.times, xout=times, sorted=True)
        T_tcylaft6 = Ska.Numpy.interpolate(model.comp['tcylaft6'].mvals,
                                           xin=model.times, xout=times, sorted=True)

        return np.vstack([T_pmtank3t, T_tmzp_my, T_tephin, T_tcylaft6])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintDPA(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6):
        self.n_ccd = n_ccd
        model_spec = os.path.join(pkg_dir, 'dpa_model_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'dpa', sim_inputs, limits,
                                 max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False, cache=True):
        model = xija.ThermalModel('dpa', start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        state_times = np.array([states['tstart'], states['tstop']])
        model.comp['sim_z'].set_data(states['simpos'], state_times)
        model.comp['eclipse'].set_data(False)
        model.comp['1dpamzt'].set_data(T0s[0])
        model.comp['dpa_power'].set_data(0.0)

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            model.comp[name].set_data(states[name], state_times)

        model.make()
        model.calc()
        T_dpa = Ska.Numpy.interpolate(model.comp['1dpamzt'].mvals,
                                      xin=model.times, xout=times, sorted=True)

        return np.vstack([T_dpa])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, self.n_ccd, self.n_ccd, 1, 1,
                   pitch1, 75000)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos')
        return np.rec.fromrecords(states, names=names)


class ConstraintTank(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, 'pftank2t_model_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'tank', sim_inputs, limits,
                                 max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False):
        model = xija.ThermalModel('tank', start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        state_times = np.array([states['tstart'], states['tstop']])
        model.comp['pitch'].set_data(states['pitch'], state_times)
        model.comp['eclipse'].set_data(False)
        # Empirical formula from settling values for pf0tank2t and pftank2t.
        # The two values converge at 22 C (pitch > 140), while at pitch = 120
        # pftank2t = 39 and pf0tank2t = 36.
        model.comp['pf0tank2t'].set_data(22 + 14. / 17. * (T0s[0] - 22.0))
        model.comp['pftank2t'].set_data(T0s[0])

        model.make()
        model.calc()
        T_tank = Ska.Numpy.interpolate(model.comp['pftank2t'].mvals,
                                       xin=model.times, xout=times, sorted=True)

        return np.vstack([T_tank])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintPSMC(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6):
        self.n_ccd = n_ccd
        model_spec = os.path.join(pkg_dir, 'psmc_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'psmc', sim_inputs, limits,
                                 max_dwell_ksec)

    def calc_model(self, states, times, T0s, state_only=False, cache=True):
        model = xija.ThermalModel('psmc', start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        state_times = np.array([states['tstart'], states['tstop']])
        model.comp['sim_z'].set_data(states['simpos'], state_times)
        model.comp['pin1at'].set_data(T0s[0] - 10.5)  # pin1at ~ 1pdeaat - 10.5 C
        model.comp['1pdeaat'].set_data(T0s[0])
        model.comp['dpa_power'].set_data(0)

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            model.comp[name].set_data(states[name], state_times)

        model.make()
        model.calc()
        T_psmc = Ska.Numpy.interpolate(model.comp['1pdeaat'].mvals,
                                       xin=model.times, xout=times, sorted=True)

        return np.vstack([T_psmc])

    def get_states1(self, start, stop, pitch1):
        states = [(start.secs, stop.secs, self.n_ccd, self.n_ccd, 1, 1,
                   pitch1, 75000)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos')
        return np.rec.fromrecords(states, names=names)


def plot_dwells1(constraint, plot_title=None, plot_file=None, figure=1):
    """Make a simple plot of the dwells and dwell statistics for the given
    ``constraint``.

    :param constraint: ConstraintModel object (e.g. constraints['all'])
    :param plot_title: plot title
    :param plot_file: output file for plot
    :param figure: matplotlib figure ID (default=1)
    """
    plt.rc("axes", labelsize=10, titlesize=12)
    plt.rc("xtick", labelsize=10)
    plt.rc("ytick", labelsize=10)
    plt.rc("font", size=10)
    plt.rc("legend", fontsize=10)

    dwells1 = constraint.dwells1
    dwell1_stats = constraint.dwell1_stats
    plt.figure(figure, figsize=(6, 4))
    plt.clf()
    names = ('none', '1pdeaat', 'tcylaft6', 'tephin', 'pline',
             '1dpamzt', 'pftank2t')
    colors = ('r', 'g', 'k', 'c', 'b', 'm', 'y')
    for name, color in zip(names, colors):
        ok = dwells1['constraint_name'] == name
        plt.plot(dwells1['pitch'][ok], dwells1['dur'][ok] / 1000., '.' + color,
                 markersize=3, label=name, mec=color)
    plt.plot(dwell1_stats['pitch'], dwell1_stats['dur50'] / 1000., '-r')
    plt.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000., '-m')
    plt.grid()
    plt.title(plot_title or '')
    plt.xlabel('Pitch (deg)')
    plt.ylabel('Dwell (ksec)')
    plt.legend(loc='upper center')
    plt.ylim(constraint.max_dwell_ksec * -0.05, constraint.max_dwell_ksec * 1.05)
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


def calc_constraints(start='2013:001',
                     n_sim=500,
                     dt=1000.,
                     max_tephin=147.0,
                     max_tcylaft6=102.0,
                     max_1pdeaat=52.5,
                     max_1dpamzt=32.5,
                     max_pftank2t=93.0,
                     n_ccd=5,
                     sim_file='sim_inputs.pkl',
                     max_dwell_ksec=200.,
                     min_pitch=45,
                     max_pitch=169,
                     bin_pitch=2,
                     constraint_models=('minus_yz', 'psmc', 'pline', 'dpa', 'tank'),
                     ):
    """
    Calculate allowed dwell times coming out of perigee given a set of
    constraint models.

    :param start: date at which to perform the constraint simulations
    :param n_sim: number of Monte-Carlo simulations of (pitch, starting condition) (default=500)
    :param dt: step size used in thermal model computations (default=1000 sec)
    :param max_tephin: TEPHIN planning limit (default=138 degF)
    :param max_tcylaft6: TCYLAFT6 planning limit (default=99 degF)
    :param max_1pdeaat: 1PDEAAT planning limit (default=52.5 degC)
    :param max_1dpamzt: 1DPAMZT planning limit (default=32.5 degC)
    :param max_pftank2t: PFTANK2T planning limit (default=93.0 degF)
    :param n_ccd: number of ACIS CCDs being used
    :param max_dwell_ksec: maximum allowed dwell time (default=200 ksec)
    :param sim_file: simulation inputs file from "pyger make" (default=sim_inputs.pkl)
    :param min_pitch: minimum pitch in simulations (default=45)
    :param max_pitch: maximum pitch in simulations (default=169)
    :param bin_pitch: pitch bin size for calculating stats (default=2)
    :param constraint_models: constraint models (default=('minus_yz', 'psmc', 'pline', 'dpa', 'tank'))

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

    n_sim_inputs = len(sim_inputs[sim_inputs.keys()[0]])
    i_sims = np.random.randint(n_sim_inputs, size=n_sim)
    pitches1 = np.random.uniform(min_pitch, max_pitch, size=n_sim)
    constraints = {}
    if 'minus_z' in constraint_models:
        constraints['minus_z'] = ConstraintMinusZ(sim_inputs,
                                                  limits=dict(tephin=FtoC(max_tephin),
                                                              tcylaft6=FtoC(max_tcylaft6)),
                                                  max_dwell_ksec=max_dwell_ksec)
    if 'minus_yz' in constraint_models:
        constraints['minus_yz'] = ConstraintMinusYZ(sim_inputs,
                                                    limits=dict(tephin=FtoC(max_tephin),
                                                                tcylaft6=FtoC(max_tcylaft6)),
                                                    max_dwell_ksec=max_dwell_ksec)
    if 'psmc' in constraint_models:
        constraints['psmc'] = ConstraintPSMC(sim_inputs,
                                             limits={'1pdeaat': max_1pdeaat},
                                             max_dwell_ksec=max_dwell_ksec,
                                             n_ccd=n_ccd)
    if 'dpa' in constraint_models:
        constraints['dpa'] = ConstraintDPA(sim_inputs,
                                           limits={'1dpamzt': max_1dpamzt},
                                           max_dwell_ksec=max_dwell_ksec,
                                           n_ccd=n_ccd)
    if 'tank' in constraint_models:
        constraints['tank'] = ConstraintTank(sim_inputs,
                                             limits={'pftank2t': FtoC(max_pftank2t)},
                                             max_dwell_ksec=max_dwell_ksec)
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
