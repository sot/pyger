''' Generate a set of simulation profiles from which starting conditions can be calculated.
'''
import os

from itertools import count
import cPickle as pickle
import json
import numpy as np

import Ska.engarchive.fetch as fetch
from Chandra.Time import DateTime
import Chandra.cmd_states as cmd_states

from . import clogging

logger = clogging.config_logger('make_sim_inputs')


class sim_inputs(object):
    def __init__(self, start=None, stop=None, sim_file='sim_inputs.pkl', n_days=3, min_dwell_sec=1000):
        ''' Generate a set of simulation profiles from which starting conditions can be calculated.

        :param start: Start date for time period to find simulation dwells
        :param stop: Stop date for time period used to find simulation dwells
        :param sim_file: Filename to save simulation profiles (pickle file)
        :param n_days: Number of days of dwell history to save for propagation

        Generates a dictionary, "sim_inputs" of simulation profiles for each model in the following
        form:

        >>> sim_inputs returns:
        Returns a dictionary of simulation profiles for each model (too long to show).

        >>> sim_inputs[modelname]
        Returns a list of simulation profiles (too long to show). "modelname" is the name of the
        model queried.

        >>> sim_inputs[modelname][sim_number]
        {'msids', 'states', 'T0s', 'tstart', 'T1s' 'tstop'}

        '''

        self.stop = DateTime(stop or DateTime().secs - 86400. * 10)
        self.start = DateTime(start or DateTime().secs - 86400. * 375)
        self.sim_file = sim_file
        self.n_days = n_days
        self.min_dwell_sec = min_dwell_sec

        (sim_start_times, sim_stop_times) = self._get_sim_times()
        self.sim_start_times = sim_start_times
        self.sim_stop_times = sim_stop_times

        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        self.constraint_models = json.load(open(os.path.join(pkg_dir, 'constraint_models.json')))

        self.sim_inputs = {}
        for name, model in self.constraint_models.items():
            logger.info('  Assembling simulation propagation data for model: {0}'.format(name))

            sim_data = self._get_sim_data(model)
            self.sim_inputs[name] = sim_data

        if self.sim_file:
            self._write_sim_inputs()

    def K2F(k):
        return (k - 273.15) * 1.8 + 32.

    def _get_states(self, state_vals=['pitch']):
        ''' Fetch states for the specified values.

        :param state_vals: List of state types to query

        Pitch is used as a default in the case where this function is called without a state
        value.

        An extra day before the start time and an extra hour after the stop time are fetched to
        prevent missing data problems.

        '''

        state_vals = [x.encode('ascii') for x in state_vals]
        start_time = self.start.secs - 86400 * (self.n_days + 1)
        stop_time = self.stop.secs + 3600
        states = cmd_states.fetch_states(start_time, stop_time, state_vals)

        return states

    def _get_telemetry(self, msids):
        ''' Fetch telemetry for the specified msids.

        :param msids: List of msids to fetch

        An extra day before the start time and an extra hour after the stop time are fetched to
        prevent missing data problems.

        '''

        start_time = self.start.secs - 86400 * (self.n_days + 1)
        stop_time = np.max((self.stop.secs, self.sim_stop_times[-1])) + 3600
        dats = fetch.MSIDset(msids, start_time, stop_time, stat='5min')
        dats.interpolate()

        return dats

    def _get_sim_times(self):
        ''' Generate a list of simulation start and stop times

        Simulation stop times are generated by finding all dwells with minimum durations
        specified by "self.min_dwell_sec". A minimum dwell time is specified to filter out
        intermediate states during maneuvers. A default value of 1Ks is used since this is the
        smallest dwell mission planning will generally schedule.

        Simulation start times are calculated by subtracting "n_days" in seconds from the
        simulation stop times.

        These start and stop times are used to propagate starting temperatures for the desired
        date in the pyger simulation (not calculated in sim_inputs() ).

        sim_start_times and sim_stop_times are lists of time values in units of seconds

        '''

        logger.info('  Finding suitable simulation dwells between: \n' +
                    '    Start:{0}\n'.format(self.start.date) +
                    '    Stop:{0}'.format(self.stop.date))

        states = self._get_states()
        duration = states['tstop'] - states['tstart']
        sim_ind = np.argwhere(duration > self.min_dwell_sec).flatten()

        sim_stop_times = states[sim_ind]['tstop']
        sim_start_times = sim_stop_times - 86400 * self.n_days

        logger.info('  Found {0} suitable dwells at least {1} seconds long.'
                    .format(len(sim_stop_times), int(self.min_dwell_sec)))

        return (sim_start_times, sim_stop_times)

    def _get_sim_data(self, model):
        ''' Gather the state and telemetry data required to propagate simulation starting values

        :param model: Dictionary of model information from "constraint_models"

        This gathers the state and telemetry data for one model for N dwells, where N is the
        number of suitable dwells found in telemetry in the specified time span (usually most
        recent year of data).

        '''

        msids = model['msids']
        state_cols = model['state_cols']

        logger.info('    Fetching state values: {0}'.format(', '.join(state_cols)))
        states = self._get_states(state_vals=state_cols)

        logger.info('    Fetching telemetry for: {0}'.format(', '.join(msids)))
        dats = self._get_telemetry(msids)

        start_ind_telem = np.searchsorted(dats.times, self.sim_start_times)
        stop_ind_telem = np.searchsorted(dats.times, self.sim_stop_times)

        sim_model_inputs = []

        for i, sim_start, sim_stop in zip(count(), self.sim_start_times, self.sim_stop_times):
            sim_states_mask = (states['tstop'] > sim_start) & (states['tstart'] < sim_stop)

            # This currently assumes all fetched values are temperatures.
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

    def _write_sim_inputs(self):
        ''' Write sim input data to a pickle file.

        '''
        logger.info('  Writing simulation inputs to {0}'.format(self.sim_file))
        with open(self.sim_file, 'w') as f:
            pickle.dump(self.sim_inputs, f, protocol=-1)
