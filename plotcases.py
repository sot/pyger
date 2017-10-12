# Licensed under a 3-clause BSD style license - see LICENSE.rst
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle as pickle
from matplotlib.ticker import AutoMinorLocator
import numpy as np

from Chandra.Time import DateTime

import pyger


mpl.rcParams['xtick.major.pad'] = 10
mpl.rcParams['ytick.major.pad'] = 10
# mpl.rcParams['savefig.edgecolor'] = [1,0.8,.8]
mpl.rc('font', family='sans-serif') 
mpl.rc('font', serif='Helvetica Neue') 
mpl.rc('font', weight='light')
helveticaneue = mpl.font_manager.FontProperties(family='Helvetica Neue', style='normal', variant='normal', weight='light')
lightfont = mpl.font_manager.FontProperties(weight='light')
minorLocator   = AutoMinorLocator(2)



def plot_cooldown(constraints2, coolstats, hotstats, model, msid, limit, filename, 
                  save_to_file=True, additional_title_text=None):
    colorpalate = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    lightblue = "#81CCf4"

    if additional_title_text is None:
        additional_title_text = ''

    # Create plot framework
    fig = plt.figure(figsize=[14,8], facecolor='w')
    rect = [0.1, 0.1, 0.88, 0.8]
    ax = fig.add_axes(rect)


    if len(np.array(hotstats).flatten()) > 0:

        # fill in NaNs in cool stats for hot regions, sort by pitch
        nans = np.array([np.nan] * len(hotstats.pitch))
        coolpitch = np.concatenate((coolstats.pitch, hotstats.pitch), axis=0)
        coolperc10 = np.concatenate((coolstats.perc10, nans), axis=0)
        coolperc50 = np.concatenate((coolstats.perc50, nans), axis=0)
        coolperc90 = np.concatenate((coolstats.perc90, nans), axis=0)
        ind = coolpitch.argsort()
        coolpitch = coolpitch[ind]
        coolperc10 = coolperc10[ind]
        coolperc50 = coolperc50[ind]
        coolperc90 = coolperc90[ind]
        

        # Plot data
        ax.plot(constraints2.dwell2_pitch_set, constraints2.dwell2_times, '.', 
                color=lightblue, alpha=0.1)
        ax.fill_between(coolpitch, coolperc10, coolperc90, facecolor=colorpalate[1], alpha=0.5)
        ax.plot(coolpitch, coolperc50, label='50 Perc Cooldown Time', linewidth=3, color=colorpalate[1])
        ax.plot(coolpitch, coolperc10, label='10 Perc Cooldown Time', linewidth=2, color=colorpalate[1])
        ax.plot(coolpitch, coolperc90, label='90 Perc Cooldown Time', linewidth=2, color=colorpalate[1])

        ax.plot(hotstats.pitch, hotstats.dwell1_duration,
                linewidth=4, color=[0.4, 0.4, 0.4], label='Max Dwell Time')

        dwell1pitch = constraints2.dwell1_pitch
        #duration = constraints2.dwell1_duration
        duration_delta = constraints2.dwell1_duration_delta

        ax.plot(dwell1pitch, duration_delta, '.', color=colorpalate[0], alpha=0.4)
        ax.plot(hotstats.pitch, hotstats.dwell1_duration_delta, color=colorpalate[0], label='Heatup Time From Tcooled to Limit', linewidth=3)
        ax.fill_between(hotstats.pitch, 0, hotstats.dwell1_duration_delta, facecolor=colorpalate[0], alpha=0.2)

        # Annotate and format the plot
        ax.legend(loc='best', fontsize=20, framealpha=0.5)
        # ylim = ax.get_ylim()
        ylim = ax.set_ylim(0, 400000)
        yticks = np.arange(ylim[0], ylim[1] + 1, 50000)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks/1000., fontsize=20)
        ax.set_ylim(ylim)
        ax.yaxis.set_minor_locator(minorLocator)

        ax.set_xticks(list(range(45,175,10)))
        ax.set_xticklabels(list(range(45,175,10)), fontsize=20)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.set_xlim(45, 170)
        ax.grid(which='both')
        ax.set_ylabel('Dwell Time in Kiloseconds', fontsize=20, fontproperties=lightfont)
        ax.set_xlabel('Pitch', fontsize=20, fontproperties=lightfont)
    else:
        ax.text(0.5, 0.5, 'Condition Not Limiting', ha='center', va='center', fontsize = 30)


    if len(additional_title_text) > 1:
        titlefontsize = 22
    else:
        titlefontsize = 26
    ax.set_title('{0}: Date:{1}, Limit={2}{3}'.format(msid.upper(),
                                                      DateTime(constraints2[0].dwell1_start).date[:8],
                                                      str(float(limit)),
                                                      additional_title_text), fontsize=titlefontsize)

    # Save plot to file
    if save_to_file:
        fig.savefig(filename)


def plotset(cases):
    for n in range(len(cases)):

        filename = './pygerdata/' + cases[n]['filename'] + '_dwell2.pkl'
        constraints2, coolstats, hotstats = pickle.load(open(filename, 'rb'))

        model = cases[n]['constraint_model']
        msid = cases[n]['msid'].lower()
        limit = cases[n]['max_' + msid]
        if 'true' in cases[n]['dh_heater'].lower():
            dh = 'ON'
        else:
            dh = 'OFF'
        n_ccd_dwell1 = cases[n]['n_ccd_dwell1']
        n_ccd_dwell2 = cases[n]['n_ccd_dwell2']
        constraints2 = constraints2[msid]
        picfile = cases[n]['filename'] + '.png'

        coolstats = coolstats[msid]
        hotstats = hotstats[msid]

        if ('psmc' in model.lower()) or ('dpa' in model.lower()) or ('dea' in model.lower()):
            additionaltitletext = ', {}CCDs Hot, {}CCDs Cool, DH={}'.format(n_ccd_dwell1, n_ccd_dwell2, dh)
        else:
            additionaltitletext = ''

        plot_cooldown(constraints2, coolstats, hotstats, model, msid, limit, picfile, save_to_file=True, additional_title_text=additionaltitletext)

        plt.close(plt.gcf())



if __name__ == "__main__":

    cases = pyger.read_cases('cases_2016_001_090.csv')
    plotset(cases)

    cases = pyger.read_cases('cases_2016_182_274.csv')
    plotset(cases)

    cases = pyger.read_cases('cases_2017_001_090.csv')
    plotset(cases)

    cases = pyger.read_cases('cases_2017_182_274.csv')
    plotset(cases)

    cases = pyger.read_cases('cases_2018_001_090.csv')
    plotset(cases)

    cases = pyger.read_cases('cases_2018_182_274.csv')
    plotset(cases)

