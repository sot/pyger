"""
N-DOF thermal model of connected nodes
"""
import re
import json
import sys

import numpy as np
import scipy.interpolate

import Ska.Numpy
from Chandra.Time import DateTime
import scipy.weave

EPOCH_START = DateTime('2005:001:00:00:00')
EPOCH_END = DateTime('2010:001:00:00:00')
EPOCH_DELTA = EPOCH_END.mjd - EPOCH_START.mjd
TWOPI365 = 2 * np.pi / 365.25
cache_eig = {}
cache_exp = {}
cache_pwrs = {}

class Power(object):
    def __init__(self, pars):
        # Get the power for the state
        Pi_pars = sorted(x for x in pars if x.startswith('pi_'))
        Pf_pars = sorted(x for x in pars if x.startswith('pf_'))
        P_pitches = np.array([float(re.search(r'_(\d+)', x).group(1)) for x in Pi_pars])
        Pi_vals = np.array([pars[x] for x in Pi_pars])
        Pf_vals = np.array([pars[x] for x in Pf_pars])
        self.Pi_interp = scipy.interpolate.interp1d(P_pitches, Pi_vals, kind='cubic')
        self.Pf_interp = scipy.interpolate.interp1d(P_pitches, Pf_vals, kind='cubic')
        self.tau_sc = pars['tau_sc']
        self.p_ampl = pars['p_ampl']
        self.exp_factor = (1 - np.exp(-EPOCH_DELTA / self.tau_sc))

    def calc(self, t_days, pitchs):
        Pi = self.Pi_interp([pitchs])
        Pf = self.Pf_interp([pitchs])
        M = (Pi - Pf) / self.exp_factor
        B = Pi - M
        P = M * np.exp(-t_days / self.tau_sc) + B + self.p_ampl * np.cos(t_days * TWOPI365)
        return P.reshape(-1)

def calc_model(pars, states, times, T0s, msids, dt=328.,
               max_dwell_ksec=200, cache=False, state_only=False):
    """Calculate temperatures for given ``states`` and initial temperatures
    ``T0s``, and model parameters ``U``, ``Ue``, ``P``.

    The ``states`` input must be an iterable list with the following keys::
      tstart  tstop  pitch

    Solve for dT/dt = A T + B where T is a vector, A is a tri-diagonal matrix
    and B is vector.

    For reference the 2-dof model calculates:

      dT0_dt = U01 * T1(t) + P0(t) + U0e * Te - (U01 + U0e) * T0(t)
             = U01 * (T1 - T0) + P0 + U0e * (Te - T0)
             = U01 * (T1 - T0) + (P0 + U0e * Te) - U0e * T0
      where T1 is the coupled node temperature and U0e is the external coupling.

    Looking at the last equation there is a degeneracy with P0(pitch) and
    U0e * Te.  So T_e is actually always frozen at 0.0 and U0e * Te gets slurped
    into P0.

      dT0_dt = U01 * (T1 - T0) + P0 - U0e * T0

    To add an internal power term (e.g. for ACIS modeling):

      dT0_dt = U01 * (T1 - T0) + P0(pitch) + C0 * P_int(state)  - U0e * T0
    """

    # Populate the A matrix.  Depends on tau_*  (pars)
    N = len(msids)
    id_pars = tuple(val for msid_pars in pars.values() for val in msid_pars.values())
    if cache and id_pars in cache_eig:
        eigvals, eigvecs, eigvecinvs, Ue, Te = cache_eig[id_pars]
    else:
        A = np.zeros((N, N))
        idx = dict((x, i) for (i, x) in enumerate(msids))
        for msid0, pars0 in pars.items():
            i0 = idx[msid0]
            A[i0, i0] = -1.0 / pars0['tau_ext']
            for msid1 in msids:
                tau_msid1 = 'tau_' + msid1
                if tau_msid1 in pars0:
                    U1 = 1.0 / pars0[tau_msid1]
                    i1 = idx[msid1]
                    A[i0, i0] += -U1
                    A[i0, i1] += U1

        eigvals, eigvecs = np.linalg.eig(A)
        eigvecinvs = np.linalg.inv(eigvecs)

        Ue = np.array([1.0/pars[msid]['tau_ext'] for msid in msids])
        Te = np.array([pars[msid]['T_e'] for msid in msids])

        cache_eig[id_pars] = eigvals, eigvecs, eigvecinvs, Ue, Te

    # Depends on pars, dt, N
    key = (id_pars, dt, N, max_dwell_ksec)
    if cache and key in cache_exp:
        all_state_times, T1, T2 = cache_exp[key]
    else:
        all_state_times = np.arange(0., max_dwell_ksec, dt/1000.)
        n_all_state_times = len(all_state_times)
        exp_l_t = []
        for i in range(N):
            exp_l_t.append(np.exp(eigvals[i] * all_state_times))

        M1 = np.zeros((n_all_state_times, N, N))
        for i in range(N):
            M1[:, i, i] = (exp_l_t[i] - 1) / eigvals[i]
        T1 = np.dot(eigvecs, np.dot(M1, eigvecinvs))

        for i in range(N):
            M1[:, i, i] = exp_l_t[i]
        T2 = np.dot(eigvecs, np.dot(M1, eigvecinvs))

        cache_exp[key] = all_state_times, T1, T2

    # Depends on states and pitch values (to nearest int) and pars
    if cache:
        key = (tuple(np.int32(states['tstart'])) +
               tuple(np.int32(states['tstart'])) +
               tuple(np.int32(states['pitch']+0.5)) +
               (id_pars,))
    if cache and key in cache_pwrs:
        pwrs = cache_pwrs[key]
    else:
        t_days = ((states['tstart'] + states['tstop']) / 2. - EPOCH_START.secs) / 86400.
        power = dict((msid, Power(pars[msid])) for msid in msids)
        pwrs = np.array([power[msid].calc(t_days, states['pitch']) for msid in msids])
        cache_pwrs[key] = pwrs

    pred_Ts = []
    pred_times = []

    # Depends on:
    # pwrs, Ue, Te, states, all_state_times, T1, B, T2, T0s
    # 
    for i, state in enumerate(states):
        t0 = state['tstart']
        t1 = state['tstop']

        B = pwrs[:,i] + Ue * Te 

        n_t = int((t1 - t0) / dt) + 2
        if state_only:
            i_times = np.array([0, n_t-2, n_t-1] if n_t > 2 else [n_t-2, n_t-1])
        else:
            i_times = slice(0, n_t)

        # Calculate predicted temperatures for this state
        T = np.dot(T1[:, i_times, :], B.reshape(N,1)) + np.dot(T2[:, i_times, :], T0s.reshape(N,1))
        T = T.reshape(N, -1)
        
        state_times = all_state_times[i_times] * 1000
        dt1 = t1 - t0            # duration of state in sec
        tm0 = state_times[-2]    # second to last time (sec) in calculated temps
        tm1 = state_times[-1]    # last time in calculated temps (extends beyond state)

        # If dt1==tm0 (state ends at second-to-last point) then T_state_end = T[:,-2].
        # If dt1==tm1 (state ends at last point) then T_state_end = T[:,-1].
        # Linearly interpolate temps between these two extremes.
        frac0 = (tm1 - dt1) / (tm1 - tm0)
        frac1 = (dt1 - tm0) / (tm1 - tm0)
        T[:,-1] = T[:,-2] * frac0 + T[:,-1] * frac1
        state_times[-1] = dt1    # Set last "calculated" time to exact end of state

        pred_Ts.append(T)
        pred_times.append(state_times + t0)
        T0s = T[:,-1]

    # Stack the list of np.arrays into a single np.array
    pred_T = np.hstack(pred_Ts)
    pred_time = np.hstack(pred_times)

    # Interpolate predicted temperatures in degC at desired output times
    Ts = [Ska.Numpy.interpolate(pred_T[i,:], pred_time, times) for i in range(N)]

    return np.vstack(Ts)
