# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Plot best-case constraints for 2012 in 2-month intervals"""
import matplotlib.pyplot as plt

import pyger

# Set up a few values that might be varied for different cases
max_tephin = 147.0
max_tcylaft6 = 102.0
max_1pdeaat = 52.5
max_1dpamzt = 32.5
max_pftank2t = 93.0
n_ccd = 5
max_dwell_ksec = 200
n_sim = 500   # Use 5000 for final plots but 500 for testing

# Set up the font sizes and initialize figure
plt.rc("axes", labelsize=10, titlesize=12)
plt.rc("legend", fontsize=10)
plt.figure(1, figsize=(6,4.5))
plt.clf()

# Step through 2012 in 2-month intervals
for month in range(0, 12, 2):
    start = '2012-%02d-01T00:00:00' % (month+1)
    print '*** Calculating constraints for %s ***' % start
    constraints = pyger.calc_constraints(start=start,
                                         max_tephin=max_tephin,
                                         max_tcylaft6=max_tcylaft6,
                                         max_1pdeaat=max_1pdeaat,
                                         max_1dpamzt=max_1dpamzt,
                                         max_pftank2t=max_pftank2t,
                                         n_ccd=n_ccd,
                                         n_sim=n_sim,
                                         max_dwell_ksec=max_dwell_ksec)
    dwell1_stats = constraints['all'].dwell1_stats
    plt.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000, label=start[:7])

# Finish making the plot with labels and then save to a file
plt.title('Best-case dwells for 2012 by month')
plt.legend(loc='upper center')
plt.xlabel('Pitch (deg)')
plt.ylabel('Dwell (ksec)')
plt.grid()
plt.ylim(0, max_dwell_ksec * 1.05)
plt.savefig('dwells_2012_month.png')


