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

    def _set_init_comps(self, init_comps, model):
        """ Set initial values based on states defined in specialized model classes

        :param init_comps: dict of node names and initial values
        :param model: xija.model instance

        This function is called by calc_model()

        """

        for key, val in init_comps.items():
            if key in self.state_cols or key == 'sim_z':
                model.comp[key].set_data(val[0], val[1])
            else:
                model.comp[key].set_data(val)

    def _get_model_temps(self, model, times):
        """ Retrieve calculated temperatures after xija model run

        :param model: xija.model instance after calculating temperatures
        :param times: list or array of times in seconds onto which model temperatures are
                      interpolated
        :returns: list of temperatures, corresponding to each msid (by order), interpolated to the
                  "times" list/array passed to this function
        """

        Ts = []
        for msid in self.msids:
            Ts.append(Ska.Numpy.interpolate(model.comp[msid].mvals, xin=model.times, xout=times,
                                            sorted=True))
        return Ts

    def _calc_model(self, states, times, T0s):
        """ Model calculation routine

        :param states: recarray of initial states
        :param times: list or array of times in seconds to which temperatures will be interpolated
        :param T0s: list of initial temperatures, in order of msid as defined in self.msids

        :returns: list of predicted temperatures for each msid, in order of msid as defined in
                  self.msids
        """

        model = xija.XijaModel(self.name, start=states['tstart'][0],
                                  stop=states['tstop'][-1],
                                  model_spec=self.model_spec)

        init_comps = self._get_init_comps(T0s, states)
        self._set_init_comps(init_comps, model)
        model.make()
        model.calc()

        return np.vstack(self._get_model_temps(model, times))

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
            Ts = self._calc_model(np_states, times, T0s)

            sim_input['dwell1_T0s'] = Ts[:, -1]

    def calc_dwells1(self, start, stop, times, pitches1, i_sims):
        """ Calculate initial dwell temperatures.

        :param start: Chandra.DateTime object for simulation start time
        :param stop: Chandra.DateTime object for simulation stop time
        :param times: Numpy array of times for which to calculate temperatures
        :param pitches1: Numpy array of pitch values to be simulated
        :param i_sims: List of indices into sim_inputs to use for each simulation. Each i_sim
                       index is paired with a pitch value so len(pitches1) == len(i_sims).

        These are the initial dwell model runs that are used to generate the traditional pyger
        plots (in conjunction with find_limit_crossings function).

        Allowed duration is calculated by a separate function "find_limit_crossings". This function
        is run at the end of calc_dwells1, however since full dwell data is also saved, one can
        recalculate allowed dwell time (initial dwell) for different limit sets.

        """

        self.start = start
        self.times = times
        self.calc_dwell1_T0s(start)

        logger.info('{0}: simulating {1} dwells'.format(self.name.upper(), len(i_sims)))

        # Calculate temperatures for each simulation.
        dwells = []
        for i_sim, pitch1 in zip(i_sims, pitches1):
            sim_input = self.sim_inputs[i_sim]
            states1 = self._get_states1(start, stop, pitch1)
            Ts = self._calc_model(states1, times, sim_input['dwell1_T0s'])

            dwells.append((pitch1, sim_input['dwell1_T0s'], times, Ts, 'none',
                           self.max_dwell_ksec * 1000))

        # Copy simulation info to dwells1 datastructure
        self.dwells1 = np.rec.fromrecords(dwells,
                                          dtype=[('pitch', np.float64),
                                                 ('T0s', np.float64, (len(self.msids),)),
                                                 ('times', np.float64, (len(times),)),
                                                 ('Ts', np.float64, (len(self.msids), len(times))),
                                                 ('constraint_name', '|S10'),
                                                 ('duration', np.float64)
                                                 ])

        self.find_limit_crossings()


    def find_limit_crossings(self, limits=None):
        """ Determine which dwells reach specified MSID limits.

        :param limits: dict of limits, keys are msids

        This function modifies the dwells1 datastructure by filling in constraint information
        """

        def interp_dwell_time(msid, msid_num, dwell, bad_idxs):
            ind1 = bad_idxs[0]
            ind0 = ind1 - 1
            dT = limits[msid] - dwell.Ts[msid_num, ind0]
            DT = dwell.Ts[msid_num, ind1] - dwell.Ts[msid_num, ind0]
            Dt = dwell.times[ind1] - dwell.times[ind0]
            return dwell.times[ind0] + Dt * dT / DT

        if not limits:
            limits = self.limits
        else:
            self.limits = limits

        for dwell1_num, dwell in enumerate(self.dwells1):

            bad_idx = None
            constraint_name = 'none'
            duration = self.max_dwell_ksec * 1000

            for n, msid in enumerate(self.msids):
                if msid in limits:

                    # Find the indexes where temperatures are higher than the limit
                    bad_idxs = np.flatnonzero(dwell.Ts[n, :] >= limits[msid])

                    # Explanation of subsequent "if" statement:
                    # If there are more than zero instances of a temperature higher than the limit
                    # and if either a constraint hasn't already been reached, or if the point in
                    # time when the prior limit was reached happened later than when the current
                    # constraint was reached, then a more constraining condition been found.
                    if len(bad_idxs) > 0 and (bad_idx is None or bad_idxs[0] < bad_idx):
                        bad_idx = bad_idxs[0]
                        constraint_name = msid

                        if bad_idxs[0] > 0:
                            # If the first temperature is below the limit (normal dwell)
                            t_interp = interp_dwell_time(msid, n, dwell, bad_idxs)
                        else:
                            # If the first temperature is above the limit (e.g. safe mode)
                            t_interp = dwell.times[0]

                        duration = t_interp - dwell.times[0]

                ok_idx = -1 if (bad_idx is None) else max(bad_idx - 1, 0)

            self.dwells1[dwell1_num]['constraint_name'] = constraint_name
            self.dwells1[dwell1_num]['duration'] = duration


    def calc_dwell1_stats(self, pitch_bins):
        dwells1 = self.dwells1
        dwell1_stats = []
        for p0, p1 in zip(pitch_bins[:-1], pitch_bins[1:]):
            ok = (dwells1['pitch'] >= p0) & (dwells1['pitch'] < p1)
            dwells_ok = dwells1[ok]
            dwells_ok.sort(order='duration')
            n_dwells_ok = len(dwells_ok)
            dwell50 = dwells_ok[int(n_dwells_ok * 0.5)]
            dwell90 = dwells_ok[int(n_dwells_ok * 0.9)]
            dwell1_stats.append((p0, p1, (p0 + p1) / 2,
                                 dwell50['duration'], dwell90['duration'],
                                 ))
        self.dwell1_stats = np.rec.fromrecords(dwell1_stats,
                                               dtype=[('pitch_bin0', np.float64),
                                                      ('pitch_bin1', np.float64),
                                                      ('pitch', np.float64),
                                                      ('dur50', np.float64),
                                                      ('dur90', np.float64),
                                                      ])


    def find_long_hot_pitch_dwells(self, msid, perc=0.2):
        """ Select the initial hot pitch sims that start with a low temperature (bottom 20%) 

        :param msid: MSID used to filter hot pitch dwells
        :param perc: fraction of starting temperatures (0.2 = 20% lowest temperatures)

        :returns: index into self.dwells1 datastructure
        
        This returns the index into the self.dwells1 datastructure, indicating which of those sims
        that reached a limit, started at the specified percentile/fraction of coldest starting
        temperatures for the MSID indicated. 

        Note about the data returned:
        Invariably, most of these will sims will have startedwithin a narrow pitch range before
        "maneuvering" to the simulated "hot" pitch. This isimportant to keep in mind since these
        lowest temperatures will often not be able to be reached during the second dwell set of
        simulations (think of these as the "cooldown" simulations). This is the purpose of adding a
        pad to the "cooled" value when calculating cooldown times.
 
        """

        msid_ind = self.msids.index(msid)

        # Eliminate sims that have zero duration 
        good_sim = self.dwells1['duration'] > 0

        # Generate a list of hot pitch values from which to choose the best dwells 
        hot_pitch_ind_all = np.flatnonzero((self.dwells1['constraint_name'] != 'none') & good_sim)

        # You want to look at the starting temperatures by msid, not by dwell 
        start_temps = self.dwells1['T0s'][hot_pitch_ind_all]
        start_temps_transpose = start_temps.transpose()

        # Assemble a list of indices for the sorted array for each msid 
        sort_ind = start_temps_transpose[msid_ind].argsort()

        # Build a list of indices from 0 to the index representing "perc" of the total 
        # length of the array, i.e. 0..19 for an array 100 elements long 
        ind = range(int(len(start_temps_transpose[0]) * perc))

        # Get the first N indices in the sorted arrays of each msid. This represents the
        # lowest N starting values for each msid. N is the length of 'ind' 
        lowest_temps_ind = sort_ind[ind]

        # Check to see if enough hot pitch values were found. If too few are found then 
        # either this constraint is not limited in this configuration, or there is not enough
        # data to work with. If enough hot pitch values were found, then return the simulations
        # that started with the coldest temperatures.
        hot_pitch_len = len(hot_pitch_ind_all)

        if hot_pitch_len > 20:
            # Return the union of each of these index arrays so the lowest N values for each
            # msid are included, return the indices in the context of the original dwells 
            # array
            hot_pitch_ind = hot_pitch_ind_all[lowest_temps_ind]

        else:
            sentence_fragment = ('hot pitch dwells were found to seed the cooldown ' +
                                 'simulations. This constraint either is not sufficiently ' +
                                 'limited in this configuration or there are too few ' +
                                 'simulations to pick from. Try increasing the number of ' +
                                 'simulations used (n-sim keyword for calc_constraints())')

            logger.info('{0}: {1} {2}'.format(self.name.upper(), hot_pitch_len, 
                                                   sentence_fragment))
            hot_pitch_ind = []

        return hot_pitch_ind


    def calc_dwells2(self, msid, start, stop, times, pitch_num, hot_dwell_temp_ratio, 
                     T_cool_ratio, pitch_range=None, match_cool_time=True, **statekw):
        """ Calculate model temperature at "non-hot" pitch ranges, starting at hot conditions

        :param msid: MSID for which to calculate cooling times
        :param start: Chandra.DateTime object for cooldown simulation start time
        :param stop: Chandra.DateTime object for cooldown simulation stop time
        :param times: Numpy array of times for which to calculate temperatures, not necessarily
                      the same times as used in the initial dwell
        :param pitch_num: Number of cool pitch values to simulate for each hot dwell
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
        :param pitch_range: Two element list or tuple defining target cooling pitch range. This
                            allows one to override the default cool pitch values.
        :param match_cool_time: Flag used to indicate whether or not the user wants to clip the
                                time from t=0 to the time where the temperature reaches the
                                T_cool_ratio from the heatup dwell time. This facilitates an equal
                                comparison of dwell times (i.e. delta temp for heatup equals delta
                                temp for cooldown so time comparisons are relevant).
        :param **statekw: Model state information; allows one to specify different states for the
                          cooldown dwell. For example one could simulate the initial dwell with 6
                          ccds and simulate the cooldown dwell with 4 ccds.

        :returns: Numpy recarray containing cooldown simulation information


        Regarding pitch values used:
        1) By default pitch values previously used are re-used as the sample set for the cooldown 
           dwells. This makes dividing up the hot vs. cold pitch values much easier. Dividing up 
           the pitch values between hot and "not hot" (i.e. cold) pitch values reduces the number
           of model runs since there is no need to run a "cooldown" dwell for at hot pitch value.
           If new pitch values were used, the current constraint profile would have to be mapped
           in some fashion, complicating the process of finding the required balance of hot and
           cold time. Also, it doesn't make sense to run a cooldown dwell for an already cold 
           dwell.

        2) Since some constraints are strongly dependant on other factors besides pitch (e.g. dpa),
           the ability to override the default behavior described above is built in. If this 
           override is triggered, the "pitch_range" keyword should return a two element list or 
           tuple including the minimum and maximum pitch to use for the cooldown simulations.
           Pitch values are chosen from those used in the original simulations.

        This second set of dwells is sometimes referred to "cooldown" dwells or sims since the 
        original intent was to determine the time to cool from hot conditions.

        Regarding time span of second set of dwells ("cooldown dwells"):
        Cooldown dwells are run from "start" to "stop"; both are inputs to this function. The 
        total cooling time (defined in start, stop, times) should be equal to or longer than the 
        maximum dwell time used in the initial dwell so cooling times are adequately mapped.

        Regarding filtering:
        Dwells that have zero duration (as specifed in the self.dwells1 datastructure) are 
        filtered out. These will likely have final temperatures that are higher than the specified
        limit.
        
        """


        def find_dwell2_pitch_ind(pitch_range, pitch_num):
            """ Select a specified number of cool pitch values randomly from the set of sims in 
            self.dwells1

            """

            # Eliminate sims that have zero duration 
            good_sim = self.dwells1['duration'] > 0

            if not pitch_range:
                # Find the indices to all sims that do not reach the limit for any msid, and which
                # do not have a zero duration (may be redundant) 
                dwell2_pitch_ind_all = np.flatnonzero((self.dwells1['constraint_name'] != 
                                                      msid.lower()) & good_sim)

            else:
                # If a pitch range is specified, then select only pitch values within this range
                ind1 = np.flatnonzero(self.dwells1['pitch'] >= pitch_range[0])
                ind2 = np.flatnonzero(self.dwells1['pitch'] <= pitch_range[1])
                dwell2_pitch_ind_all = ind1 & ind2 & good_sim

            # Select N of these indices randomly, where N is specified by 'pitch_num' 
            ind = np.random.randint(len(dwell2_pitch_ind_all), size=pitch_num)

            return dwell2_pitch_ind_all[ind]


        def calculate_init_data(self, msid_ind, i_hot, ratio):
            """ Determine starting temperatures for the second dwell set

            Starting temperatures are determined by using a ratio of delta temperature between
            the hot dwell starting temperature and the hot dwell final temperature at or just
            before reaching the limit for the constraining msid. The time at this interpolated
            temperature is then used to interpolate the starting temperatures for the rest of the
            nodes in the model.

            """
            dwells1 = self.dwells1
            T = dwells1['Ts'][i_hot]
            times = dwells1['times'][i_hot]

            Tlim = self.limits[self.msids[msid_ind]]

            # Note that if the ratio is small enough, it is possible the resulting temp appears
            # twice in the dwell since temperatures will dip sometimes before rising
            T_ratio = T[msid_ind][0] + ratio * (Tlim - T[msid_ind][0])
            t_ratio = np.interp(T_ratio, T[msid_ind], times)

            T_dwell2_0 = np.array([np.interp(t_ratio, times, Ts_hot) for Ts_hot in T])
            t_dwell2 = t_ratio - times[0]

            return T_dwell2_0, t_dwell2


        def calc_cooldown_times(self, i_hot, times, Ts_dwell2, msid_ind, T_cool_ratio):
            times = times - times[0]
            Ts_dwell1 = self.dwells1['Ts'][i_hot][msid_ind]
            Ts_dwell2 = Ts_dwell2[msid_ind]
            T_threshold = Ts_dwell1[0] + (1 - T_cool_ratio) * (Ts_dwell1[-1] - Ts_dwell1[0])
            t_interp = np.interp(T_threshold, Ts_dwell2[::-1], times[::-1])
            
            return t_interp, T_threshold

        def clip_dwell1_time(self, msid_ind, i_hot, T_cool_ratio, dur_delta):
            T = self.dwells1['Ts'][i_hot][msid_ind]
            times = self.dwells1['times'][i_hot]
            times = times - times[0]

            Tlim = self.limits[self.msids[msid_ind]]

            # If the user does run with more than one limit active (e.g. tephin and tcylaft6),
            # then the temperature may not actually reach the limit for the msid specified in the
            # calc_dwells2() function call. More specfically, the dwell will be registered in the
            # self.dwells1 datastructure due to reaching a different limit, however this particular
            # dwell may be unconstraining for the msid specified. In this case, just calculate the
            # T_ratio using the first and last temperatures (specified in the code under "else")
            if np.any(T >= Tlim):
                T_ratio = T[0] + (1 - T_cool_ratio) * (Tlim - T[0])
            else:
                T_ratio = T[0] + (1 - T_cool_ratio) * (T[-1] - T[0])

            # It is possible that the "T_ratio" temperature was reached twice in the heatup dwell
            # due to the varying time constants of the nodes in the model. Make sure the first
            # instance is used.
            ind = np.arange(len(T))[T > T_ratio]

            t_clip =  np.interp(T_ratio, T[:(ind[0] + 1)], times[:(ind[0] + 1)])

            dur_delta = dur_delta - t_clip

            return dur_delta



        # Find the indices to the appropriate inital hot pitch values, and random cool pitch 
        # values 
        hot_pitch_ind = self.find_long_hot_pitch_dwells(msid)
        dwell2_pitch_ind = find_dwell2_pitch_ind(pitch_range, pitch_num=pitch_num)

        msid_ind = self.msids.index(msid)

        logger.info('{0} - {1}: simulating {3} cooldown dwells for each of the {2} hot dwells'\
                     .format(msid, self.name.upper(), len(hot_pitch_ind), len(dwell2_pitch_ind)))
        logger.info(('Using a hot dwell temperature ratio of {0} to determine cooldown dwell' + 
                       'starting conditions').format(hot_dwell_temp_ratio))

        dwells2 = []
        for j, i_hot in enumerate(hot_pitch_ind):
            logger.info('{0}: simulating cooldown dwells for hot dwell {1} of {2}'\
                        .format(self.name.upper(), j + 1, len(hot_pitch_ind)))

            T_dwell2_0, dur_delta = calculate_init_data(self, msid_ind, i_hot, hot_dwell_temp_ratio)

            if match_cool_time is True:
                dur_delta = clip_dwell1_time(self, msid_ind, i_hot, T_cool_ratio, dur_delta)

            dwell2_pitch_set = self.dwells1['pitch'][dwell2_pitch_ind]

            t_cool_set = []
            for cool_pitch in dwell2_pitch_set:
                states2 = self._get_states1(start, stop, cool_pitch, **statekw)
                Ts = self._calc_model(states2, times, T_dwell2_0)
                t_cool, T_cool = calc_cooldown_times(self, i_hot, times, Ts, msid_ind, T_cool_ratio)
                t_cool_set.append(t_cool)

            dwell1duration = self.dwells1['duration'][i_hot]
            dwell2_case = (i_hot, self.dwells1['pitch'][i_hot], msid, dur_delta, dwell1duration,
                           self.dwells1['T0s'][i_hot], T_dwell2_0, dwell2_pitch_set, 
                           np.array(t_cool_set), T_cool, self.start.date, DateTime(start).date)

            dwells2.append(dwell2_case)
        
        # Bogus data intended to show zero required cooling time (because no limiting condition
        # was found). This is a bit of a hack, and better solutions will be considered in the
        # future. This helps to avoid errors when automating pyger data generation.
        limits = np.array([1000,] * len(self.msids))
        dwell2_pitch_set = np.array(self.dwells1['pitch'][dwell2_pitch_ind])
        zerostart = np.array([np.nan,] * len(self.msids))
        allnantime = np.array((np.nan,) * len(dwell2_pitch_ind))
        if len(dwells2) == 0:
            dwells2 = [(-1, 45, 'none', self.max_dwell_ksec * 1000, self.max_dwell_ksec * 1000, 
                         zerostart, limits, dwell2_pitch_set, allnantime, 0, self.start.date, DateTime(start).date),
                        (-2, 169, 'none', self.max_dwell_ksec * 1000, self.max_dwell_ksec * 1000, 
                         zerostart, limits, dwell2_pitch_set, allnantime, 0, self.start.date, DateTime(start).date)]

        dtype = np.dtype([('dwell1_ind', np.int32), 
                          ('dwell1_pitch', np.float64),
                          ('msid', '|S10'),
                          ('dwell1_duration_delta', np.float64),
                          ('dwell1_duration', np.float64),
                          ('T_dwell1_0', np.float64, len(self.msids) ),
                          ('T_dwell2_0', np.float64, len(self.msids) ),
                          ('dwell2_pitch_set', np.float64, len(dwell2_pitch_ind) ),
                          ('dwell2_times', np.float64, len(dwell2_pitch_ind) ),
                          ('dwell2_cool_temp', np.float64),
                          ('dwell1_start', '|S21'),
                          ('dwell2_start', '|S21')
                          ])

        return  np.rec.fromrecords(dwells2, dtype=dtype)

