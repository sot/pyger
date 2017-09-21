#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import matplotlib.pyplot as plt
import cPickle as pickle

import pyger


def plot_constraints(constraints):
    """Make a plot of the dwells and dwell statistics for each limited msid

    :param constraint: Main obect containing data for all constraint models

    """

    plotinfo = {'bgcolor': [1, 1, 1],
                'fgcolor': [0.0, 0.0, 0.0],
                'top': 0.92,
                'hspace': 0.06,
                'width': 8.5,
                'height': 11,
                'left': 0.1,
                'bottom': 0.05,
                'right': 0.94}

    plt.rc('axes', edgecolor=plotinfo['fgcolor'])

    fig = plt.figure(figsize=(plotinfo['width'], plotinfo['height']),
                     facecolor=plotinfo['bgcolor'])

    # Add title and subtitles for test plot page
    fig.text(0.5, 0.98, 'Pyger Test Plots', ha="center",
             va="center", size=14, color=plotinfo['fgcolor'])

    fig.text(0.5, 0.96, 'Pyger Version %s' % pyger.__version__,
             ha="center", va="center", size=12,
             color=plotinfo['fgcolor'])

    # Determine individual plot height
    numplots = len(constraints.keys())
    plotheight = ((1 - plotinfo['bottom'] -
                  (1 - plotinfo['top']) - plotinfo['hspace'] * (numplots - 1))
                  / float(numplots))

    # Plot the dwell limitations for each model individually
    for n, key in enumerate(constraints.keys()):

        # Determine plot location and add axis to figure
        bottom = plotinfo['bottom'] + n * (plotheight + plotinfo['hspace'])
        loc = [plotinfo['left'],
               bottom,
               plotinfo['right'] - plotinfo['left'],
               plotheight]  # left, bottom, width, height
        ax = fig.add_axes(loc, axisbg=plotinfo['bgcolor'])

        # Get individual limit items and dwell data
        # limitedlist = list_individual_constraints(constraints[key])
        dwells1 = constraints[key].dwells1
        dwell1_stats = constraints[key].dwell1_stats
        limitedlist = list(set(dwells1['constraint_name']))

        # Plot dwell times for each individual simulation, for each limited
        # item separately
        for name in limitedlist:

            # If name is not None, then plot each limited item separately. This
            # is only an issue because None is
            if name:
                ok = dwells1['constraint_name'] == name
                ax.plot(dwells1['pitch'][ok], dwells1['duration'][ok] / 1000., '.',
                        markersize=3, label=name)

        # If name is None, then it is the 'pline' model, plot all points.
        if not name:
            ax.plot(dwells1['pitch'], dwells1['duration'] / 1000., '.',
                    markersize=3, label='all')

        ax.plot(dwell1_stats['pitch'], dwell1_stats['dur50'] / 1000., '-r')
        ax.plot(dwell1_stats['pitch'], dwell1_stats['dur90'] / 1000., '-m')
        ax.grid()
        ax.legend(loc='best', fontsize=6, labelspacing=0.05)
        ax.set_title('%s' % key)
        ax.set_xlabel('Pitch (deg)')
        ax.set_ylabel('Dwell (ksec)')
        ax.set_ylim(constraints[key].max_dwell_ksec * -0.05,
                    constraints[key].max_dwell_ksec * 1.05)

    outfile = 'test_plot_constraints_{}.png'.format(pyger.__version__)
    print('\nSaving plots to {}'.format(outfile))
    fig.savefig(outfile)


def runtest():
    sim_inputs = pickle.load(open('sim_inputs.pkl'))
    allmodels = sim_inputs.keys()
    constraints = pyger.calc_constraints(start='2013:050', max_dwell_ksec=400,
                                         constraint_models=tuple(allmodels))
    plot_constraints(constraints)


if __name__ == '__main__':
    print('******* Running test_plot_constraints *********')
    print('Pyger module file is {}\n'.format(pyger.__file__))

    runtest()
