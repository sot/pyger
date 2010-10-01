import os
import sys
from itertools import count
import cPickle as pickle
import json

import numpy

import Ska.engarchive.fetch as fetch
from Chandra.Time import DateTime
import Chandra.cmd_states as cmd_states
import Ska.DBI

from . import characteristics
psmc_powers = dict((x[0:3], x[3]) for x in characteristics.psmc_power)

pkg_dir = os.path.dirname(os.path.abspath(__file__))
constraint_models = json.load(open(os.path.join(pkg_dir, 'constraint_models.json')))

def K2F(k):
    return (k-273.15) * 1.8 + 32.

def get_states(datestart, datestop, state_vals):
    db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read', database='aca')
    states = db.fetchall("""select * from cmd_states
                            where datestop > '%s' and
                            datestart < '%s'
                            order by datestart""" %
                         (datestart, datestop))
    db.conn.close()
    states = cmd_states.reduce_states(states, cols=state_vals, allow_identical=False)

    bad = states['pitch'] > 169.9
    states['pitch'][bad] = 169.9

    return states

def make_sim_inputs(start=None, stop=None, outfile='sim_inputs.pkl', n_days=3):
    stop = DateTime(stop or DateTime().secs - 86400. * 10)
    start = DateTime(start or DateTime().secs - 86400. * 365)
    cobsrqid = fetch.MSID('cobsrqid', start.secs, stop.secs)

    # Find transitions between engineering obsids and science obsids
    obsids = cobsrqid.vals
    start_eng = (obsids[:-1] < 50000) & (obsids[1:] > 50000)
    stop_eng = (obsids[:-1] > 50000) & (obsids[1:] < 50000)

    # Find the approximate end of perigee times as the end of a
    # continuous block of eng obsids at least 6 hours long
    start_times = cobsrqid.times[start_eng]
    stop_times = cobsrqid.times[stop_eng]
    if stop_times[0] < start_times[0]:
        stop_times = stop_times[1:]
    if len(start_times) > len(stop_times):
        start_times = start_times[:len(stop_times)]
    if numpy.any(stop_times - start_times <= 0.0):
        raise ValueError('Problem finding perigee times')

    perigee_dwells = stop_times - start_times > 6 * 3600
    eop_times = stop_times[perigee_dwells]
    dts = eop_times[1:] - eop_times[:-1]
    ok = (dts > 225000) & (dts < 235000)
    eop_times_ok = eop_times[:-1][ok]

    

    sim_inputs = {}

    for name, model in constraint_models.items():
        msids = model['msids']
        state_cols = model['state_cols']
        dats = {}
        idx_starts = {}
        idx_stops = {}
        sim_stop_times = eop_times_ok
        sim_start_times = sim_stop_times - 86400 * n_days
        for msid in msids:
            dats[msid] = fetch.MSID(msid, start.secs - 86400*(n_days+1), stop.secs, stat='5min')
            idx_starts[msid] = numpy.searchsorted(dats[msid].times, sim_start_times)
            idx_stops[msid] = numpy.searchsorted(dats[msid].times, sim_stop_times)

        states = get_states(start.date, stop.date, state_cols)

        sim_inputs[name] = []
        for i, tstart, tstop in zip(count(), sim_start_times, sim_stop_times):
            ok = (states['tstop'] > tstart) & (states['tstart'] < tstop)
            T0s = dict((x, dats[x].vals[idx_starts[x][i]] - 273.15) for x in msids)
            T1s = dict((x, dats[x].vals[idx_stops[x][i]] - 273.15) for x in msids)
            out_states = []
            for state in states[ok]:
                out_state = dict(tstart=max(state['tstart'], tstart), tstop=min(state['tstop'], tstop))
                for state_col in state_cols:
                    out_state[state_col] = state[state_col]
                if name == 'psmc':   # special case pseudo- state col
                    out_state['power'] = psmc_powers[state['ccd_count'], 1, 1]
                out_states.append(out_state)

            sim_inputs[name].append(dict(msids=msids, tstart=tstart, tstop=tstop, T0s=T0s, T1s=T1s, states=out_states))

    with open(outfile, 'w') as f:
        pickle.dump(sim_inputs, f, protocol=-1)

