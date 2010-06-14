"""
N-mass thermal model of connected nodes.

These routines make use of a thermal model that is specified by the following
set of parameters.  Each node (or "mass") in the model has a corresponding
parameter set in a file "pars_<MSID>.json"::
   
  Input power to the node for a particular sun pitch at 2004:001

    pi_045: 0.487916 
    pi_060: 0.936244 
    pi_090: 1.17665 
    pi_120: 1.16191 
    pi_145: 0.800581 
    pi_170: 0.271429 

  Input power to the node for a particular sun pitch at 2010:001

    pf_045: 1.033705 
    pf_060: 1.496659 
    pf_090: 1.955456 
    pf_120: 1.958443 
    pf_145: 1.493900 
    pf_170: 0.663974 

  Long-term and annual variation is specified as follows:

    p_ampl: 0.06623    # Sinusoidal amplitude added to input power
    tau_sc: 1732.25294 # Expoenential decay timescale

  The external heat sink has temperature and coupling described by:
  
    T_e: -0.615463
    tau_ext: 37.737    # ksec

  Couplings to other nodes are given by tau_<MSID>:

    tau_tcylaft6: 118.118   # ksec
    tau_tmzp_my: 90.465

"""
import re
import json

import numpy as np
import scipy.interpolate

import Ska.Numpy
from Chandra.Time import DateTime

DEBUG = False
EPOCH_START = DateTime('2004:001:00:00:00')
EPOCH_END = DateTime('2010:001:00:00:00')
EPOCH_DELTA = EPOCH_END.mjd - EPOCH_START.mjd
TWOPI365 = 2 * np.pi / 365.25

def power(pars, t_day, pitch):
    """Calculate the effective input power term for the given model parameters
    ``pars`` (corresponding to a particular MSID) at the specified ``t_day``
    and ``pitch``.

    This routine first does a cubic spline interpolation in pitch to find the
    power at the given ``pitch`` at the EPOCH_START time using the
    corresponding "pi_<pitch>" parameters. Then the same interpolation is done
    using the EPOCH_END parameters "pf_<pitch>".

    Next the long-term exponential decay connecting the EPOCH_START and
    EPOCH_END power values is computed using the ``tau_sc`` parameter.  Finally
    the yearly sinusoidal variation is included with the ``p_ampl`` parameter::

    :param pars: dict of model parameters
    :param t_day: days since EPOCH_START
    :param pitch: sun pitch angle

    :returns: effective power term
    """

    Pi_pars = sorted(x for x in pars if x.startswith('pi_'))
    Pf_pars = sorted(x for x in pars if x.startswith('pf_'))
    P_pitches = np.array([float(re.search(r'_(\d+)', x).group(1)) for x in Pi_pars])
    Pi_vals = np.array([pars[x] for x in Pi_pars])
    Pf_vals = np.array([pars[x] for x in Pf_pars])
    Pi = scipy.interpolate.interp1d(P_pitches, Pi_vals, kind='cubic')([pitch])[0]
    Pf = scipy.interpolate.interp1d(P_pitches, Pf_vals, kind='cubic')([pitch])[0]

    tau_sc = pars['tau_sc']
    p_ampl = pars['p_ampl']

    M = (Pi - Pf) / (1 - np.exp(-EPOCH_DELTA / tau_sc))
    B = Pi - M
    P = M * np.exp(-t_day / tau_sc) + B + p_ampl * np.cos(t_day * TWOPI365)

    return P

def calc_model(pars, states, times, T0s, dt=328.):
    """Calculate temperatures for given ``states`` and initial temperatures
    ``T0s``.  Within each state calculate the model at ``dt`` second intervals,
    then interpolate back onto the requested ``times``.  The ``states`` input
    must be an iterable list of dicts with the following keys: tstart tstop
    pitch.

    Solve for dT/dt = A T + B where T is a vector of temperatures, A is a
    fixed matrix defining the node couplings, and B is a vector defining the
    pitch-dependent power input to each node.

    :param pars: dict of parameters
    :param states: iterable list of states (must be contiguous)
    :param times: array of times at which to return the model temperatures
    :param T0s: N-vector of initial temperatures at each node
    :param dt: time step used for model propagation within each state

    :returns: predicted temperature arrays with shape [len(msids), len(times)]
    """

    # Populate the A matrix to define node couplings
    msids = sorted(pars.keys())
    N = len(msids)
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

    # Define couplings to the external constant-temperature sinks
    Ue = np.array([1.0 / pars[msid]['tau_ext'] for msid in msids])
    Te = np.array([pars[msid]['T_e'] for msid in msids])

    pred_Ts = []
    pred_times = []
    for state in states:
        t0 = state['tstart']

        t_day = ((state['tstart'] + state['tstop']) / 2. - EPOCH_START.secs) / 86400.
        pwrs = np.array([power(pars[msid], t_day, state['pitch']) for msid in msids])
        B = pwrs + Ue * Te 

        if DEBUG and state == states[0]:
            print 'T0s =', T0s
            print 'A = \n', A
            print 'B =', B
            print 'dT/dt =', np.dot(A, T0s) + B
            print 'Ue =', Ue
            print 'Te =', Te
            print 'pwrs =', pwrs
        
        # Make array of times within state 
        n_t = int((state['tstop'] - state['tstart']) / dt)
        state_times = np.linspace(state['tstart'], state['tstop'], n_t+2)

        # Calculated predicted temperatures for this state For now use an
        # explicit loop.  Vectorizing would improve perfomance significantly.
        T = np.ones((N, len(state_times)))
        for i, t in enumerate((state_times - t0) / 1000):
            exp_eigvals = np.exp(eigvals * t)
            L1 = np.dot(np.dot(np.diag((exp_eigvals - 1) / eigvals), eigvecinvs), B)
            L2 = np.dot(np.dot(np.diag(exp_eigvals), eigvecinvs), T0s)
            T[:,i] = np.dot(eigvecs, L1 + L2)

        pred_Ts.append(T)
        pred_times.append(state_times)
        T0s = T[:,-1]

    # Stack the list of np.arrays into a single np.array
    pred_T = np.hstack(pred_Ts)
    pred_time = np.hstack(pred_times)

    # Interpolate predicted temperatures in degC at desired output times
    Ts = []
    for i in range(N):
        Ts.append(Ska.Numpy.interpolate(pred_T[i,:], pred_time, times))

    return np.vstack(Ts)
