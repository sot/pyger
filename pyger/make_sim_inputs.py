''' Generate a set of simulation profiles from which starting conditions can be calculated.
'''
import os

from itertools import count
import pickle as pickle
import json
import numpy as np

import Ska.engarchive.fetch as fetch
from Chandra.Time import DateTime
import Chandra.cmd_states as cmd_states

from . import clogging

logger = clogging.config_logger('make_sim_inputs')


def make_sim_inputs(start=None, stop=None, sim_file='sim_inputs.pkl', n_days=3, min_dwell_sec=1000,
                    max_dwell_num=300):
    """
    Generate a set of simulation profiles from which starting conditions can be calculated.

    :param start: Start date for time period to find simulation dwells
    :param stop: Stop date for time period used to find simulation dwells
    :param sim_file: Filename to save simulation profiles (pickle file)
    :param n_days: Number of days of dwell history to save for propagation
    :param min_dwell_sec: Minimum number of seconds to use as a propagation ending dwell
    :param max_dwell_num: Maximum number of sample dwells to save

    Generates a dictionary, "sim_inputs" of simulation profiles for each model in the following
    form:

    >>> sim_inputs
    Returns a dictionary of simulation profiles for each model (too long to show).

    >>> sim_inputs[modelname]
    Returns a list of simulation profiles (too long to show). "modelname" is the name of the
    model queried.

    >>> sim_inputs[modelname][sim_number]
    {'msids', 'states', 'T0s', 'tstart', 'T1s' 'tstop'}
    """

    stop = DateTime(stop or DateTime().secs - 86400. * 10)
    start = DateTime(start or DateTime().secs - 86400. * 375)

    (sim_start_times, sim_stop_times) = get_sim_times(start, stop, n_days, min_dwell_sec,
                                                      max_dwell_num)

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    constraint_models = json.load(open(os.path.join(pkg_dir, 'constraint_models.json')))

    sim_inputs = {}
    for name, model in list(constraint_models.items()):
        logger.info('  Assembling simulation propagation data for model: {0}'.format(name))

        sim_data = get_sim_data(start, stop, n_days, sim_start_times, sim_stop_times, model)
        sim_inputs[name] = sim_data

    if sim_file:
        write_sim_inputs(sim_file, sim_inputs)
    else:
        return sim_inputs


def K2F(k):
    return (k - 273.15) * 1.8 + 32.


def get_states(start, stop, n_days, state_vals=['pitch']):
    """
    Fetch states for the specified values.

    :param state_vals: List of state types to query

    Pitch is used as a default in the case where this function is called without a state
    value.

    An extra day before the start time and an extra hour after the stop time are fetched to
    prevent missing data problems.
    """

    state_vals = [x.encode('ascii') for x in state_vals]
    start_time = start.secs - 86400 * (n_days + 1)
    stop_time = stop.secs + 3600
    states = cmd_states.fetch_states(start_time, stop_time, state_vals)

    return states


def get_telemetry(start, stop, n_days, sim_stop_times, msids):
    """
    Fetch telemetry for the specified msids.

    :param msids: List of msids to fetch

    An extra day before the start time and an extra hour after the stop time are fetched to
    prevent missing data problems.
    """

    start_time = start.secs - 86400 * (n_days + 1)
    stop_time = np.max((stop.secs, np.max(sim_stop_times))) + 3600
    dats = fetch.MSIDset(msids, start_time, stop_time, stat='5min')
    dats.interpolate()

    return dats


def get_sim_times(start, stop, n_days, min_dwell_sec, max_dwell_num):
    """
    Generate a list of simulation start and stop times

    Simulation stop times are generated by finding all dwells with minimum durations
    specified by "min_dwell_sec". A minimum dwell time is specified to filter out
    intermediate states during maneuvers. A default value of 1Ks is used since this is the
    smallest dwell mission planning will generally schedule.

    Simulation start times are calculated by subtracting "n_days" in seconds from the
    simulation stop times.

    These start and stop times are used to propagate starting temperatures for the desired
    date in the pyger simulation (not calculated in sim_inputs() ).

    sim_start_times and sim_stop_times are lists of time values in units of seconds
    """

    logger.info('  Finding suitable propagation ending dwells between: \n' +
                '    Start:{0}\n'.format(start.date) +
                '    Stop:{0}'.format(stop.date))
    logger.info('  Limiting number of sample dwells to {0}'.format(max_dwell_num))

    states = get_states(start, stop, n_days)
    duration = states['tstop'] - states['tstart']
    sim_ind_all = np.argwhere(duration > min_dwell_sec).flatten()
    rand_ind = np.random.randint(len(sim_ind_all), size=max_dwell_num)
    sim_ind = sim_ind_all[rand_ind]

    sim_stop_times = states[sim_ind]['tstop']
    sim_start_times = sim_stop_times - 86400 * n_days

    logger.info('  Found {0} suitable propagation ending dwells at least {1} seconds long.'
                .format(len(sim_stop_times), int(min_dwell_sec)))

    return (sim_start_times, sim_stop_times)


def get_sim_data(start, stop, n_days, sim_start_times, sim_stop_times, model):
    """
    Gather the state and telemetry data required to propagate simulation starting values

    :param model: Dictionary of model information from "constraint_models"

    This gathers the state and telemetry data for one model for N dwells, where N is the
    number of suitable dwells found in telemetry in the specified time span (usually most
    recent year of data).
    """

    msids = model['msids']
    state_cols = model['state_cols']

    logger.info('    Fetching state values: {0}'.format(', '.join(state_cols)))
    states = get_states(start, stop, n_days, state_vals=state_cols)

    logger.info('    Fetching telemetry for: {0}'.format(', '.join(msids)))
    dats = get_telemetry(start, stop, n_days, sim_stop_times, msids)

    start_ind_telem = np.searchsorted(dats.times, sim_start_times)
    stop_ind_telem = np.searchsorted(dats.times, sim_stop_times)

    sim_model_inputs = []

    for i, sim_start, sim_stop in zip(count(), sim_start_times, sim_stop_times):
        sim_states_mask = (states['tstop'] > sim_start) & (states['tstart'] < sim_stop)

        T0s = dict((x, dats[x].vals[start_ind_telem[i]] - 273.15) for x in msids)
        T1s = dict((x, dats[x].vals[stop_ind_telem[i]] - 273.15) for x in msids)

        out_states = []

        for state in states[sim_states_mask]:

            out_state = dict(tstart=np.max((state['tstart'], sim_start)),
                             tstop=np.min((state['tstop'], sim_stop)))

            for state_col in state_cols:
                out_state[state_col] = state[state_col]

            out_states.append(out_state)

        sim_model_inputs.append(dict(msids=msids, tstart=sim_start, tstop=sim_stop, T0s=T0s,
                                     T1s=T1s, states=out_states))

    return sim_model_inputs


def write_sim_inputs(sim_file, sim_inputs):
    """
    Write sim input data to a pickle file.
    """

    logger.info('  Writing simulation inputs to {0}'.format(sim_file))
    with open(sim_file, 'w') as f:
        pickle.dump(sim_inputs, f, protocol=-1)
