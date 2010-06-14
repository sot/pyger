import sys
from itertools import count
import cPickle as pickle

import numpy

import Ska.Numpy
import Ska.engarchive.fetch as fetch
from Ska.Matplotlib import plot_cxctime
from Chandra.Time import DateTime
import Chandra.cmd_states as cmd_states
import Ska.DBI

n_days_prop = 3

def K2F(k):
    return (k-273.15) * 1.8 + 32.


def get_states(datestart, datestop):
    db = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read', database='aca')
    states = db.fetchall("""select * from cmd_states
                            where datestop > '%s' and
                            datestart < '%s'
                            order by datestart""" %
                         (datestart, datestop))
    db.conn.close()
    states = Ska.Numpy.compress(
        states,
        delta=dict(pitch=2.0),
        colnames=('pitch', 'tstart', 'tstop','datestart','datestop'),
        avg=dict(tstart=numpy.min, tstop=numpy.max,
                 datestart=min, datestop=max),
        diff=dict((key, lambda x,y: False) for key
                  in ('tstart', 'tstop','datestart','datestop')))

    bad = states['pitch'] > 169.9
    states['pitch'][bad] = 169.9

    return states

# Get the last 4 years of data
datestart = '2009:146'   # 2006:002 works also
# datestart = '2010:136'   # 2006:002 works also
datestop = '2010:146'
cobsrqid = fetch.MSID('cobsrqid', datestart, datestop)

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
# eop_times_ok = eop_times_ok[2:]  # chop the first couple

msids = sorted(['tcylaft6', 'tmzp_my', 'tephin', 'tfssbkt1', 'tcylfmzm'])
dats = {}
idx_starts = {}
idx_stops = {}
sim_stop_times = eop_times_ok
sim_start_times = sim_stop_times - 86400 * n_days_prop
for msid in msids:
    dats[msid] = fetch.MSID(msid, DateTime(datestart).secs - 86400*(n_days_prop+1), datestop, stat='5min')
    idx_starts[msid] = numpy.searchsorted(dats[msid].times, sim_start_times)
    idx_stops[msid] = numpy.searchsorted(dats[msid].times, sim_stop_times)

states = get_states(datestart, datestop)

sim_inputs = []
for i, tstart, tstop in zip(count(), sim_start_times, sim_stop_times):
    ok = (states['tstop'] > tstart) & (states['tstart'] < tstop)
    T0s = dict((x, dats[x].vals[idx_starts[x][i]]-273.15) for x in msids)
    T1s = dict((x, dats[x].vals[idx_stops[x][i]]-273.15) for x in msids)
    out_states = []
    for state in states[ok]:
        state_tstart = max(state['tstart'], tstart)
        state_tstop = min(state['tstop'], tstop)
        out_states.append(dict(tstart=state_tstart,
                               tstop=state_tstop,
                               # datestart=DateTime(state_tstart).date,
                               # datestop=DateTime(state_tstop).date,
                               pitch=state['pitch']))
    sim_inputs.append(dict(msids=msids, tstart=tstart, tstop=tstop, T0s=T0s, T1s=T1s, states=out_states))

with open('sim_inputs.pkl', 'w') as f:
    pickle.dump(sim_inputs, f, protocol=-1)
    
