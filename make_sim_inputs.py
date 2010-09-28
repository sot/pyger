import sys
from itertools import count
import cPickle as pickle
import json

import numpy

import Ska.Numpy
import Ska.engarchive.fetch as fetch
from Chandra.Time import DateTime
import Chandra.cmd_states as cmd_states
import Ska.DBI

import characteristics
psmc_powers = dict((x[0:3], x[3]) for x in characteristics.psmc_power)

n_days_prop = 3

def K2F(k):
    return (k-273.15) * 1.8 + 32.

def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("--out",
                      help="Output file (default = sim_inputs_<start>_<stop>.pkl)")
    parser.add_option("--start",
                      help="Sim input start time (default = Now - 1 year)")
    parser.add_option("--stop",
                      help="Sim input stop time (default = Now - 10 days)")
    (opt, args) = parser.parse_args()
    return (opt, args)


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


opt, args = get_options()

# Get obsids over the duration
stop = DateTime(opt.stop or DateTime().secs - 86400. * 10)
start = DateTime(opt.start or DateTime().secs - 86400. * 365)
print start.date, stop.date
cobsrqid = fetch.MSID('cobsrqid', start.secs, stop.secs)
print 'got cobsrqid'

# Find transitions between engineering obsids and science obsids
obsids = cobsrqid.vals
start_eng = (obsids[:-1] < 50000) & (obsids[1:] > 50000)
stop_eng = (obsids[:-1] > 50000) & (obsids[1:] < 50000)

# Find the approximate end of perigee times as the end of a
# continuous block of eng obsids at least 6 hours long
start_times = cobsrqid.times[start_eng]
stop_times = cobsrqid.times[stop_eng]
perigee_dwells = stop_times - start_times > 6 * 3600
eop_times = stop_times[perigee_dwells]
dts = eop_times[1:] - eop_times[:-1]
ok = (dts > 225000) & (dts < 235000)
eop_times_ok = eop_times[:-1][ok]

models = json.load(open('constraint_models.json'))

sim_inputs = {}

for name, model in models.items():
    msids = model['msids']
    state_cols = model['state_cols']
    dats = {}
    idx_starts = {}
    idx_stops = {}
    sim_stop_times = eop_times_ok
    sim_start_times = sim_stop_times - 86400 * n_days_prop
    for msid in msids:
        dats[msid] = fetch.MSID(msid, start.secs - 86400*(n_days_prop+1), stop.secs, stat='5min')
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

    outfile = opt.out or 'sim_inputs_{0}_{1}.pkl'.format(start.greta[:7], stop.greta[:7])

with open(outfile, 'w') as f:
    pickle.dump(sim_inputs, f, protocol=-1)
    
