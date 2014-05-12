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
from .pyger_cases import read_cases, run_pyger_cases, PostPyger

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


def save_pyger_pickle(constraints, filename):
    """ Save pyger data to pickle file

    :param constraints: This can be either the output from calc_constraints() or a
                        numpy.core.records.recarray object, which is output by calc_constraints2()
    :param filename: Name of pickle file to write data to

    """
    pickleable = ['dwell1_stats', 'dwells1', 'times', 'limits', 'max_dwell_ksec', 'model_spec', 
                  'msids', 'name', 'n_ccd', 'n_sim', 'sim_inputs', 'start', 'state_col']

    all_pickle_data = {}
    for name in constraints:
        pickle_data = {}
        if isinstance(constraints[name], np.core.records.recarray):
            all_pickle_data.update({name:constraints[name]})
        else:
            for key in constraints[name].__dict__.keys():
                if key in pickleable:
                    if key == 'start':
                        pickle_data.update({key:constraints[name].__dict__[key].secs})
                    else:
                        pickle_data.update({key:constraints[name].__dict__[key]})
            all_pickle_data.update({name:pickle_data})
    pickle.dump(all_pickle_data, open(filename,'w'), protocol=2)


def load_pyger_pickle(filename):
    """ Load pyger data from pickle file back into object compatible with pyger plotting methods

    :param filename: File name of pickled output from calc_constraints()

    This is only meant to be used to read in the initial constraints object produced by
    calc_constraints(), not the cooldown data produced by calc_constraints2(). The data prduced
    by calc_constraints2() should be able to be read in with a simple pickle.load() function.
    """
    class saved_pyger_data(object):
        def __init__(self, pickled_constraint):
            for key in pickled_constraint:
                self.__dict__.update({key:pickled_constraint[key]})

    rawdata = pickle.load(open(filename,'r'))
    pyger_compatible_data = {}
    for name in rawdata.keys():
        constraint = saved_pyger_data(rawdata[name])
        pyger_compatible_data.update({name:constraint})

    return pyger_compatible_data


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
                                          dtype=[('duration', np.float64),
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

    def _get_init_comps(self, T0s, states):
        # Initialize all

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'eclipse': False,
                      'tcylaft6': T0s[0],
                      'tcylfmzm': T0s[1],
                      'tephin': T0s[2],
                      'tfssbkt1': T0s[3],
                      'tmzp_my': T0s[4]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, **stateskw):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintMinusYZ(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, 'minusyz_model.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'minus_yz', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'eclipse': False,
                      'pmtank3t': T0s[0],
                      'tmzp_my': T0s[1],
                      'tephin': T0s[2],
                      'tcylaft6': T0s[3],
                      'pseudo_0': T0s[3] - 4.0,
                      'pseudo_1': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, **stateskw):
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

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'sim_z': (states['simpos'], state_times),
                      'eclipse': False,
                      '1dpamzt': T0s[0],
                      'dpa_power': 0.0}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos')
        return np.rec.fromrecords(states, names=names)


class ConstraintTank(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, 'pftank2t_model_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'tank', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        # Empirical formula from settling values for pf0tank2t and pftank2t.
        # The two values converge at 22 C (pitch > 140), while at pitch = 120
        # pftank2t = 39 and pf0tank2t = 36.

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'eclipse': False,
                      'pf0tank2t': 22 + 14. / 17. * (T0s[0] - 22.0),
                      'pftank2t': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, **stateskw):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintFwdblkhd(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, '4rt700t_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'fwdblkhd', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'eclipse': False,
                      '4rt700t_0': T0s[0],
                      '4rt700t': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, **stateskw):
        states = [(start.secs, stop.secs, pitch1)]
        names = ('tstart', 'tstop', 'pitch')
        return np.rec.fromrecords(states, names=names)


class ConstraintAca(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec):
        model_spec = os.path.join(pkg_dir, 'aca_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'aca', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):
        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'eclipse': False,
                      'aca0': T0s[0],
                      'aacccdpt': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, **stateskw):
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

    def _get_init_comps(self, T0s, states):

        # pin1at ~ 1pdeaat - 10.5 C

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'sim_z': (states['simpos'], state_times),
                      'pin1at': T0s[0] - 10.5,
                      '1pdeaat': T0s[0],
                      'dpa_power': 0.0}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos)]
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
             '1dpamzt', 'pftank2t', 'aacccdpt')
    colors = ('r', 'g', 'k', 'c', 'b', 'm', 'y', 'g')
    for name, color in zip(names, colors):
        ok = dwells1['constraint_name'] == name
        plt.plot(dwells1['pitch'][ok], dwells1['duration'][ok] / 1000., '.' + color,
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


def plot_cooldown(constraints2, coolstats, hotstats, model, msid, limit, filename, 
                  save_to_file=True, additional_title_text=None):
    colorpalate = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    lightblue = "#81CCf4"

    if additional_title_text is None:
        additional_title_text = ''

    # Create plot framework
    fig = plt.figure(figsize=[14,8], facecolor='w')
    rect = [0.06, 0.1, 0.88, 0.8]
    ax = fig.add_axes(rect)


    if len(np.array(hotstats).flatten()) > 0:

        # fill in NaNs in cool stats for hot regions, sort by pitch
        nans = np.array([np.nan] * len(hotstats.pitch))
        coolpitch = np.concatenate((coolstats.pitch, hotstats.pitch), axis=0)
        coolperc10 = np.concatenate((coolstats.perc10, nans), axis=0)
        coolperc50 = np.concatenate((coolstats.perc50, nans), axis=0)
        coolperc90 = np.concatenate((coolstats.perc90, nans), axis=0)
        ind = coolpitch.argsort()
        coolpitch = coolpitch[ind]
        coolperc10 = coolperc10[ind]
        coolperc50 = coolperc50[ind]
        coolperc90 = coolperc90[ind]
        

        # Plot data
        ax.plot(constraints2.dwell2_pitch_set, constraints2.dwell2_times, '.', 
                color=lightblue, alpha=0.1)
        ax.fill_between(coolpitch, coolperc10, coolperc90, facecolor=colorpalate[1], alpha=0.5)
        ax.plot(coolpitch, coolperc50, label='50 Perc Cooldown Time', linewidth=3, color=colorpalate[1])
        ax.plot(coolpitch, coolperc10, label='10 Perc Cooldown Time', linewidth=2, color=colorpalate[1])
        ax.plot(coolpitch, coolperc90, label='90 Perc Cooldown Time', linewidth=2, color=colorpalate[1])

        ax.plot(hotstats.pitch, hotstats.dwell1_duration,
                linewidth=2, color=[0.4, 0.4, 0.4], label='Max Dwell Time')

        dwell1pitch = constraints2.dwell1_pitch
        #duration = constraints2.dwell1_duration
        duration_delta = constraints2.dwell1_duration_delta

        ax.plot(dwell1pitch, duration_delta, '.', color=colorpalate[0], alpha=0.4)
        ax.plot(hotstats.pitch, hotstats.dwell1_duration_delta, color=colorpalate[0], label='Heatup Time From Tcooled to Limit', linewidth=3)
        ax.fill_between(hotstats.pitch, 0, hotstats.dwell1_duration_delta, facecolor=colorpalate[0], alpha=0.2)

        # Annotate and format the plot
        ax.legend(loc='best')
        ylim = ax.get_ylim()
        yticks = np.arange(ylim[0], ylim[1] + 1, 25000)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks/1000.)
        ax.set_ylim(ylim)
        ax.set_xticks(range(45,175,5))
        ax.set_xlim(45, 170)
        ax.grid()
        ax.set_ylabel('Dwell Time in Kiloseconds')
        ax.set_xlabel('Pitch')
    else:
        ax.text(0.5, 0.5, 'Condition Not Limiting', ha='center', va='center', fontsize = 16)

    ax.set_title('{0}: Date:{1}, Limit={2}{3}'.format(msid.upper(),
                                                      DateTime(constraints2[0].dwell1_start).date[:8],
                                                      str(float(limit)),
                                                      additional_title_text))

    # Save plot to file
    if save_to_file:
        fig.savefig(filename)


def merge_dwells1(constraints):
    """Merge the dwells in the ``constraints`` list, finding the shortest from among the
    constraints.

    :param constraints: list of ModelConstraint objects
    :returns: NumPy recarray of dwells with pitch, duration, constraint_name columns
    """
    dwells1 = []
    for i in range(constraints[0].n_sim):
        constraint = min(constraints, key=lambda x: x.dwells1['duration'][i])
        dwells1.append(tuple(constraint.dwells1[x][
                       i] for x in ('pitch', 'duration', 'constraint_name')))
    dwells1 = np.rec.fromrecords(dwells1, names=('pitch', 'duration', 'constraint_name'))
    return dwells1


def calc_constraints(start='2013:001',
                     n_sim=500,
                     dt=1000.,
                     max_tephin=147.0,
                     max_tcylaft6=102.0,
                     max_1pdeaat=52.5,
                     max_1dpamzt=31.5,
                     max_pftank2t=93.0,
                     max_aacccdpt=-15.0,
                     max_4rt700t=79.0,
                     n_ccd=6,
                     sim_file='sim_inputs.pkl',
                     max_dwell_ksec=200.,
                     min_pitch=45,
                     max_pitch=169,
                     bin_pitch=2,
                     constraint_models=('minus_yz', 'psmc', 'pline', 'dpa', 'tank', 'aca',
                      'fwdblkhd')):
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
    :param max_aacccdpt: ACA CCD planning limit (default=-15.0 degC)
    :param max_4rt700t: OBA forward bulkhead planning limit (default=79.0 degF)
    :param n_ccd: number of ACIS CCDs being used
    :param max_dwell_ksec: maximum allowed dwell time (default=200 ksec)
    :param sim_file: simulation inputs file from "pyger make" (default=sim_inputs.pkl)
    :param min_pitch: minimum pitch in simulations (default=45)
    :param max_pitch: maximum pitch in simulations (default=169)
    :param bin_pitch: pitch bin size for calculating stats (default=2)
    :param constraint_models: constraint models, default=('minus_yz', 'psmc', 'pline', 'dpa', 'tank')

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
    if 'aca' in constraint_models:
        constraints['aca'] = ConstraintAca(sim_inputs,
                                           limits={'aacccdpt': max_aacccdpt},
                                           max_dwell_ksec=max_dwell_ksec)


    if '4rt700t' in constraint_models:
        constraints['4rt700t'] = ConstraintFwdblkhd(sim_inputs,
                                           limits={'4rt700t': max_4rt700t},
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




def calc_constraints2(constraints,
                      start='2013:001',
                      n_ccd=None,
                      max_dwell_ksec=400.,
                      pitch_num=50,
                      pitch_range=None,
                      hot_dwell_temp_ratio=0.9,
                      T_cool_ratio=0.9,
                      constraint_models=('minus_yz', 'psmc', 'dpa', 'tank', 'aca', 'fwdblkhd'),
                      msids=('tephin', 'tcylaft6', '1pdeaat', '1dpamzt', 'pftank2t', 'aacccdpt', '4rt700t')
                      ):
    """
    Calculate allowed dwell times coming out of perigee given a set of
    constraint models.

    :param start: date at which to perform the constraint simulations
    :param n_ccd: number of ACIS CCDs being used
    :param max_dwell_ksec: maximum allowed dwell time (default=200 ksec)
    :param pitch_num: Number of cool pitch values to simulate for each hot dwell
    :param pitch_range: Two element list or tuple defining target cooling pitch range. This
                        allows one to override the default cool pitch values.
    :param hot_dwell_temp_ratio: Time ratio during hot dwell for which to extract cooldown
                                 dwell starting conditions. This allows one to initiate a
                                 cooldown simulation at any point during a hot dwell, not
                                 just when the hot dwell reaches a thermal limit.
    :param T_cool_ratio: Temperature ratio with respect to the temperature increase during
                         each hot dwell, used to calculate the point at which a simulated
                         cooldown dwell has "cooled". Ultimately, the reported cooling time is
                         the time it takes for the MSID to reach the "cooled" temperature,
                         starting at the given hot conditions. A ratio of 0.9 means the MSID
                         has to have cooled back down 90% to the original hot dwell starting
                         temperature.
    :param constraint_models: constraint models included

    :returns: dict of computed constraint model objects
    """


    start = DateTime(start)
    stop = DateTime(start.secs + max_dwell_ksec * 1000)
    norm_profile = np.linspace(0, 1, 100)**10
    times = start.secs + norm_profile * (stop.secs - start.secs)


    constraints_list = [constraints[x] for x in constraint_models]
    cooldown = {}
    coolstats = {}
    hotstats = {}
    for constraint in constraints_list:
        for msid in msids:
            if msid in constraint.msids:

                if n_ccd is None:
                    if n_ccd in constraint.__dict__.keys():
                        n_ccd = constraint.n_ccd
                    else:
                        n_ccd = 6

                cooldown[msid] = constraint.calc_dwells2(msid,
                                                         start,
                                                         stop,
                                                         times,
                                                         pitch_num=pitch_num, 
                                                         pitch_range=pitch_range,
                                                         hot_dwell_temp_ratio=hot_dwell_temp_ratio, 
                                                         T_cool_ratio=T_cool_ratio,
                                                         ccd_count=n_ccd)
                coolstats[msid] = calc_dwell2_cool_stats(cooldown[msid])
                hotstats[msid] = calc_dwell2_hot_stats(cooldown[msid])

    return cooldown, coolstats, hotstats



def calc_dwell2_cool_stats(dwell2_case):
    """ Calculate relevant statistics for "cooldown" dwell simulations.
  
    :param dwell2_case: This is a Numpy recarray representing the output of calc_constraints2() for
                        a single MSID. If running calc_constraints2() for multiple MSIDs, run
                        calc_dwell2_stats() for each MSID individually.
  
    :returns: Numpy recarray of relevant statistical data
  
    """

    dwell2_pitches = dwell2_case['dwell2_pitch_set'][0] # Each set is identical
    dwell2_times = dwell2_case['dwell2_times'].swapaxes(0,1)
    coolstats = []
    for timeset, pitch in zip(dwell2_times, dwell2_pitches):
        t = np.sort(timeset)
        coolstats.append((pitch, 
                          t[int(len(t) * 0.1)],
                          t[int(len(t) * 0.5)],
                          t[int(len(t) * 0.9)]))

    dtype = np.dtype([('pitch', np.float64),
                      ('perc10', np.float64),
                      ('perc50', np.float64),
                      ('perc90', np.float64)])
    
    coolstats = np.rec.fromrecords(coolstats, dtype)
    coolstats.sort(order='pitch')

    return coolstats


def calc_dwell2_hot_stats(dwell2_case):
    """ Calculate relevant statistics for hot dwells used to seed cooldown simulations.
  
    :param dwell2_case: This is a Numpy recarray representing the output of calc_constraints2() for
                        a single MSID. If running calc_constraints2() for multiple MSIDs, run
                        calc_dwell2_stats() for each MSID individually.
  
    :returns: Numpy recarray of relevant statistical data
  
    Although the hot dwell data was already calculated and is present in the dwells1 datastructure,
    the full hot dwell durations are not 100% comparable to the cooldown dwells due to the
    cooldown temperature ratio (T_cool_ratio) used to determine when a cooldown dwell has reached
    "close enough" to the original starting temperature. The comparable hot dwell time has a
    portion of the initial dwell time "clipped" from the total duration. The amount of time clipped is equal to
    the amount of time it would take to reach the "close enough" cooldown temperature during the
    hot dwell starting at the dwell1 T0.

    """

    dwell2_case.sort(order='dwell1_pitch')
    
    pitch_set = np.arange(min(dwell2_case.dwell1_pitch), max(dwell2_case.dwell1_pitch) + 1, 2)
    pitch_pairs = [(pitch_set[n], pitch_set[n + 1]) for n in range(len(pitch_set) - 1)]
    
    hotstats = []
    for p in pitch_pairs:
        ind1 = dwell2_case.dwell1_pitch >= p[0]
        ind2 = dwell2_case.dwell1_pitch < p[1]
        ind = ind1 & ind2
        if any(ind):
            hotstats.append((np.mean(p), np.mean(dwell2_case.dwell1_duration[ind]), np.mean(dwell2_case.dwell1_duration_delta[ind])))

    dtype = np.dtype([('pitch', np.float64),
                      ('dwell1_duration', np.float64),
                      ('dwell1_duration_delta', np.float64)])

    
    # This is a bit of a hack, but prevents the recarray conversion from failing
    if len(hotstats) > 0:
        hotstats = np.rec.fromrecords(hotstats, dtype)
    else:
        hotstats = np.rec.fromrecords([[],[],[]], dtype)

    return hotstats



