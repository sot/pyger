# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Calculate Chandra dwell times given thermal constraints
"""
import sys
import os
import json
import pickle
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# import Ska.Numpy
from Chandra.Time import DateTime
# import asciitable
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
                  'msids', 'name', 'n_ccd', 'n_sim', 'sim_inputs', 'start', 'state_col', 'roll']

    all_pickle_data = {}
    for name in constraints:
        pickle_data = {}
        if isinstance(constraints[name], np.core.records.recarray):
            all_pickle_data.update({name:constraints[name]})
        else:
            for key in list(constraints[name].__dict__.keys()):
                if key in pickleable:
                    if key == 'start':
                        pickle_data.update({key:constraints[name].__dict__[key].secs})
                    else:
                        pickle_data.update({key:constraints[name].__dict__[key]})
            all_pickle_data.update({name:pickle_data})
    pickle.dump(all_pickle_data, open(filename,'wb'), protocol=2)


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

    rawdata = pickle.load(open(filename,'rb'))
    pyger_compatible_data = {}
    for name in list(rawdata.keys()):
        constraint = saved_pyger_data(rawdata[name])
        pyger_compatible_data.update({name:constraint})

    return pyger_compatible_data


class ConstraintDPA(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6, roll=0.0):
        self.n_ccd = n_ccd
        self.roll = roll
        model_spec = os.path.join(pkg_dir, 'dpa_model_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'dpa', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'sim_z': (states['simpos'], state_times),
                      'roll': (states['roll'], state_times),
                      'eclipse': False,
                      '1dpamzt': T0s[0],
                      'dpa0': T0s[0],
                      'dpa_power': 0.0}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch', 'roll'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, roll=None, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos, roll)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos', 'roll')
        return np.rec.fromrecords(states, names=names)


class ConstraintDEA(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6, roll=0.0):
        self.n_ccd = n_ccd
        self.roll = roll
        model_spec = os.path.join(pkg_dir, 'dea_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'dea', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'sim_z': (states['simpos'], state_times),
                      'roll': (states['roll'], state_times),
                      'eclipse': False,
                      '1deamzt': T0s[0],
                      'dpa_power': 0.0}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch', 'roll'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, roll=None, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos, roll)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos', 'roll')
        return np.rec.fromrecords(states, names=names)


class ConstraintTank(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, roll=0.0):
        model_spec = os.path.join(pkg_dir, 'pftank2t_model_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        self.roll = roll
        ConstraintModel.__init__(self, 'tank', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        # Empirical formula from settling values for pf0tank2t and pftank2t.
        # The two values converge at 22 C (pitch > 140), while at pitch = 120
        # pftank2t = 39 and pf0tank2t = 36.

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'roll': (states['roll'], state_times),
                      'eclipse': False,
                      'pf0tank2t': 22 + 14. / 17. * (T0s[0] - 22.0),
                      'pftank2t': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, roll=None, **stateskw):

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, pitch1, roll)]
        names = ('tstart', 'tstop', 'pitch', 'roll')
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


class ConstraintTcylaft6(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, roll=0.0):
        model_spec = os.path.join(pkg_dir, 'tcylaft6_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        self.roll = roll
        ConstraintModel.__init__(self, 'tcylaft6', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'pitch': (states['pitch'], state_times),
                      'roll': (states['roll'], state_times),
                      'eclipse': False,
                      'tcylaft6_0': T0s[0],
                      'tcylaft6': T0s[0]}

        return init_comps

    def _get_states1(self, start, stop, pitch1, roll=None, **stateskw):

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, pitch1, roll)]
        names = ('tstart', 'tstop', 'pitch', 'roll')
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
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6, dh_heater=True, roll=0.0):
        self.n_ccd = n_ccd
        self.roll = roll
        self.dh_heater = dh_heater
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
                      'dpa_power': 0.0,
                      'eclipse': False,
                      'dh_heater':self.dh_heater, 
                      'roll': (states['roll'], state_times)}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, roll=None, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos, roll)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos', 'roll')
        return np.rec.fromrecords(states, names=names)




class ConstraintACISFP(ConstraintModel):
    def __init__(self, sim_inputs, limits, max_dwell_ksec, n_ccd=6, dh_heater=True, roll=0.0):
        self.n_ccd = n_ccd
        self.roll = roll
        self.dh_heater = dh_heater
        model_spec = os.path.join(pkg_dir, 'acisfp_spec.json')
        self.model_spec = json.load(open(model_spec, 'r'))
        ConstraintModel.__init__(self, 'acisfp', sim_inputs, limits,
                                 max_dwell_ksec)

    def _get_init_comps(self, T0s, states):

        state_times = np.array([states['tstart'], states['tstop']])
        init_comps = {'sim_z': (states['simpos'], state_times),
                      'roll': (states['roll'], state_times),
                      'eclipse': False,
                      'fptemp_11': T0s[0],
                      '1cbat':-55.0,
                      'sim_px':-110.0,
                      'dpa_power': 0.0,
                      'orbitephem0_x': 25000e3,
                      'orbitephem0_y': 25000e3,
                      'orbitephem0_z': 25000e3,
                      'aoattqt1': 0.0,
                      'aoattqt2': 0.0,
                      'aoattqt3': 0.0,
                      'aoattqt4': 1.0,
                      'dh_heater':self.dh_heater}

        for name in ('ccd_count', 'fep_count', 'vid_board', 'clocking', 'pitch'):
            init_comps[name] = (states[name], state_times)

        return init_comps

    def _get_states1(self, start, stop, pitch1, ccd_count=None, fep_count=None, vid_board=1,
                     clocking=1, simpos=75000, roll=None, **stateskw):

        if ccd_count is None:
            ccd_count = self.n_ccd

        if fep_count is None:
            fep_count = ccd_count

        if roll is None:
            roll = self.roll

        states = [(start.secs, stop.secs, ccd_count, fep_count, vid_board, clocking,
                   pitch1, simpos, roll)]
        names = ('tstart', 'tstop', 'ccd_count', 'fep_count', 'vid_board',
                 'clocking', 'pitch', 'simpos', 'roll')
        return np.rec.fromrecords(states, names=names)


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
                     max_tcylaft6=999.0,
                     max_1pdeaat=999.0,
                     max_1dpamzt=999.0,
                     max_1deamzt=999.0,
                     max_pftank2t=999.0,
                     max_aacccdpt=999.0,
                     max_4rt700t=999.0,
                     max_fptemp_11=999.0,
                     n_ccd=6,
                     roll=0.0,
                     dh_heater=True,
                     sim_file='sim_inputs.pkl',
                     max_dwell_ksec=200.,
                     min_pitch=45,
                     max_pitch=169,
                     bin_pitch=2,
                     constraint_models=('psmc', 'dpa', 'dea', 'tank', 'aca', 'fwdblkhd',
                      'tcylaft6', 'acisfp')):
    """
    Calculate allowed dwell times coming out of perigee given a set of
    constraint models.

    :param start: date at which to perform the constraint simulations
    :param n_sim: number of Monte-Carlo simulations of (pitch, starting condition) (default=500)
    :param dt: step size used in thermal model computations (default=1000 sec)
    :param max_tcylaft6: TCYLAFT6 planning
    :param max_1pdeaat: 1PDEAAT planning limit
    :param max_1dpamzt: 1DPAMZT planning limit
    :param max_1deamzt: 1DEAMZT planning limit
    :param max_pftank2t: PFTANK2T planning limit
    :param max_aacccdpt: ACA CCD planning limit
    :param max_4rt700t: OBA forward bulkhead planning limit
    :param max_acisfp: ACIS Focal Plane limit
    :param n_ccd: number of ACIS CCDs being used
    :param roll: roll used for simulations
    :param max_dwell_ksec: maximum allowed dwell time (default=200 ksec)
    :param sim_file: simulation inputs file from "pyger make" (default=sim_inputs.pkl)
    :param min_pitch: minimum pitch in simulations (default=45)
    :param max_pitch: maximum pitch in simulations (default=169)
    :param bin_pitch: pitch bin size for calculating stats (default=2)
    :param constraint_models: constraint models, default=('psmc', 'dpa', 'dea', 'tank', 'aca',
        'fwdblkhd', 'tcylaft6', 'acisfp')

    :returns: dict of computed constraint model objects
    """
    start = DateTime(start)
    stop = DateTime(start.secs + max_dwell_ksec * 1000)
    times = np.arange(start.secs, stop.secs, dt)
    try:
        sim_inputs = pickle.load(open(sim_file, 'rb'))
    except IOError:
        logger.error('ERROR: simulation inputs file "{0}" not found.'
                     '  Run "pyger make" or "pyger make --help".'.format(sim_file))
        sys.exit(1)

    n_sim_inputs = len(sim_inputs[list(sim_inputs.keys())[0]])
    i_sims = np.random.randint(n_sim_inputs, size=n_sim)
    pitches1 = np.random.uniform(min_pitch, max_pitch, size=n_sim)
    constraints = {}
    if 'psmc' in constraint_models:
        constraints['psmc'] = ConstraintPSMC(sim_inputs,
                                             limits={'1pdeaat': max_1pdeaat},
                                             max_dwell_ksec=max_dwell_ksec,
                                             n_ccd=n_ccd,
                                             dh_heater=dh_heater,
                                             roll=roll)
    if 'dpa' in constraint_models:
        constraints['dpa'] = ConstraintDPA(sim_inputs,
                                           limits={'1dpamzt': max_1dpamzt},
                                           max_dwell_ksec=max_dwell_ksec,
                                           n_ccd=n_ccd,
                                           roll=roll)
    if 'dea' in constraint_models:
        constraints['dea'] = ConstraintDEA(sim_inputs,
                                           limits={'1deamzt': max_1deamzt},
                                           max_dwell_ksec=max_dwell_ksec,
                                           n_ccd=n_ccd,
                                           roll=roll)
    if 'tank' in constraint_models:
        constraints['tank'] = ConstraintTank(sim_inputs,
                                             limits={'pftank2t': FtoC(max_pftank2t)},
                                             max_dwell_ksec=max_dwell_ksec,
                                             roll=roll)
    if 'aca' in constraint_models:
        constraints['aca'] = ConstraintAca(sim_inputs,
                                           limits={'aacccdpt': max_aacccdpt},
                                           max_dwell_ksec=max_dwell_ksec)

    if 'fwdblkhd' in constraint_models:
        constraints['fwdblkhd'] = ConstraintFwdblkhd(sim_inputs,
                                                     limits={'4rt700t': FtoC(max_4rt700t)},
                                                     max_dwell_ksec=max_dwell_ksec)

    if 'tcylaft6' in constraint_models:
        constraints['tcylaft6'] = ConstraintTcylaft6(sim_inputs, 
                                                     limits={'tcylaft6': FtoC(max_tcylaft6)},
                                                     max_dwell_ksec=max_dwell_ksec,
                                                     roll=roll)              

    if 'acisfp' in constraint_models:
        constraints['acisfp'] = ConstraintACISFP(sim_inputs,
                                           limits={'fptemp_11': max_fptemp_11},
                                           max_dwell_ksec=max_dwell_ksec,
                                           n_ccd=n_ccd,
                                           dh_heater=dh_heater,
                                           roll=roll)


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
                      dh_heater=True,
                      max_dwell_ksec=400.,
                      pitch_num=50,
                      pitch_range=None,
                      hot_dwell_temp_ratio=0.9,
                      T_cool_ratio=0.9,
                      constraint_models=('psmc', 'dpa', 'dea', 'tank', 'aca', 'fwdblkhd',
                        'tcylaft6', 'acisfp'),
                      msids=('tcylaft6', '1pdeaat', '1dpamzt', '1deamzt', 'pftank2t',
                        'aacccdpt', '4rt700t', 'acisfp_11')
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
                    if n_ccd in list(constraint.__dict__.keys()):
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

    if (len(dwell2_times) > 0) and ('none' not in dwell2_case['msid']):
      coolstats = []
      for timeset, pitch in zip(dwell2_times, dwell2_pitches):
          t = np.sort(timeset)
          coolstats.append((pitch, 
                            t[int(len(t) * 0.1)],
                            t[int(len(t) * 0.5)],
                            t[int(len(t) * 0.9)]))
    else:
      logger.info(('No cooldown calculations performed, tstart={}').format(dwell2_case['dwell1_start']))
      # print('No cooldown calculations performed, tstart={}'.format(dwell2_case['dwell1_start']))

      coolstats = [[None, None, None, None], [None, None, None, None]]


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



