import asciitable
from Chandra.Time import DateTime

dat = asciitable.read('start_science_times.dat')

times = dat['ss_time']
dt = times[1:] - times[:-1]
bad = (dt < 215000) | (dt > 240000)
ibads = np.flatnonzero(bad)
dates = DateTime(times).date

for ibad in ibads:
    try:
        print dates[ibad-1], '{0:6.1f}'.format(dt[ibad-1]/1000), \
              dates[ibad], '{0:6.1f}'.format(dt[ibad]/1000), \
              dates[ibad+1], '{0:6.1f}'.format(dt[ibad+1]/1000), \
              dates[ibad+2]
    except:
        pass
