"""Plot best-case constraints for 2011 in 2-month intervals"""
import pyger

# Set up a few values that might be varied for different cases
max_tephin = 128.0
max_tcylaft6 = 93.0
n_ccd = 6
max_dwell_ksec = 200
n_sim = 500   # Use 5000 for final plots but 500 for testing

# Set up the font sizes and initialize figure
rc("axes", labelsize=10, titlesize=12)
rc("legend", fontsize=10)
figure(1, figsize=(6,4.5))
clf()

# Step through 2011 in 2-month intervals
for month in range(0, 12, 2):
    start = '2011-%02d-01T00:00:00' % (month+1)
    print '*** Calculating constraints for %s ***' % start
    constraints = pyger.calc_constraints(start=start,
                                         max_tephin=max_tephin,
                                         max_tcylaft6=max_tcylaft6,
                                         n_ccd=n_ccd,
                                         n_sim=n_sim,
                                         max_dwell_ksec=max_dwell_ksec)
    dwell1_stats = constraints['all'].dwell1_stats
    plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000, label=start[:7])

# Finish making the plot with labels and then save to a file
title('Best-case dwells for 2011 by month')
legend(loc='upper center')
xlabel('Pitch (deg)')
ylabel('Dwell (ksec)')
grid()
ylim(0, max_dwell_ksec * 1.05)
savefig('dwells_2011_month.png')


