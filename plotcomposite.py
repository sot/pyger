import numpy as np
import matplotlib.pyplot as plt

from Chandra.Time import DateTime

import pyger
from functools import reduce


def assemble_cases(options, n_sim, max_dwell_ksec, dwell_2_pitch_num, dates):

    template = {'max_1dpamzt': '999',
                'max_1deamzt': '999',
                'max_1pdeaat': '999',
                'max_aacccdpt': '999',
                'max_pftank2t': '999',
                'max_tcylaft6': '999',
                'max_tephin': '999',
                'max_4rt700t': '999',
                'max_fptemp_11': '999'}

    newcases = []
    for date in dates:
        for msid in list(options.keys()):
            lims = np.arange(
                options[msid]['startlim'],
                options[msid]['stoplim'] +
                options[msid]['limstep'],
                options[msid]['limstep'])
            for lim in lims:
                if options[msid]['sameccd']:
                    nccds = list(zip(options[msid]['nccd1'], options[msid]['nccd1']))
                else:
                    nccds = list(
                        itertools.product(
                            options[msid]['nccd1'],
                            options[msid]['nccd2']))

                if 'cool_pitch_min' in list(options[msid].keys()):
                    coolpitchmin = options[msid]['cool_pitch_min']
                    coolpitchmax = options[msid]['cool_pitch_max']
                else:
                    coolpitchmin = None
                    coolpitchmax = None

                for nccd1, nccd2 in nccds:
                    for dh in options[msid]['dh']:
                        s = template.copy()
                        s['msid'] = msid
                        s['constraint_model'] = options[msid.upper()]['model']
                        s['filename'] = 'pyger_single_msid_{}{}_{}_{}_{}ccd-{}ccd_DH-{}'.format(
                            date[:4], date[-3:], msid.lower(), lim, nccd1, nccd2, dh)
                        s['max_' + msid.lower()] = str(lim)
                        s['title'] = msid.upper() + ': ' + str(lim) + \
                            options[msid]['units'] + ' ' + date
                        s['n_ccd_dwell1'] = str(nccd1)
                        s['n_ccd_dwell2'] = str(nccd2)
                        s['dh_heater'] = str(dh == 'ON')
                        s['start'] = date
                        s['n_sim'] = str(n_sim)
                        s['max_dwell_ksec'] = max_dwell_ksec
                        s['dwell_2_pitch_num'] = str(dwell_2_pitch_num)
                        s['cool_pitch_min'] = coolpitchmin
                        s['cool_pitch_max'] = coolpitchmax
                        newcases.append(s)

    return newcases


def write_cases(newcases, filename):
    fid = open(filename, 'w')

    header = list(newcases[0].keys())

    # Write Header
    for name in header[:-1]:
        fid.write('{},'.format(name))
    fid.write('{}\n'.format(header[-1]))

    # Write Cases
    for case in newcases:
        for name in header[:-1]:
            fid.write('{},'.format(case[name]))
        fid.write('{}\n'.format(case[header[-1]]))

    fid.close()


def chunks(listofstuff, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(listofstuff), n):
        yield listofstuff[i:i + n]


def factors(n):
    return set(reduce(list.__add__, ([i, n // i]
                                     for i in range(1, int(n**0.5) + 1) if n % i == 0)))


def add_pline(postdata):
    dwell_156_162 = 60 * (9 * 60 + 20)
    dwell_156_166 = 60 * (7 * 60 + 50)
    dwell_156_170 = 60 * (4 * 60 + 20)

    postdata.max_hot_time.max_time[156] = dwell_156_162
    postdata.max_hot_time.max_time[157] = dwell_156_162
    postdata.max_hot_time.max_time[158] = dwell_156_162
    postdata.max_hot_time.max_time[159] = dwell_156_162
    postdata.max_hot_time.max_time[160] = dwell_156_162
    postdata.max_hot_time.max_time[161] = dwell_156_162
    postdata.max_hot_time.max_time[162] = dwell_156_162
    postdata.max_hot_time.max_time[163] = dwell_156_166
    postdata.max_hot_time.max_time[164] = dwell_156_166
    postdata.max_hot_time.max_time[165] = dwell_156_166
    postdata.max_hot_time.max_time[166] = dwell_156_166
    postdata.max_hot_time.max_time[167] = dwell_156_170
    postdata.max_hot_time.max_time[168] = dwell_156_170
    postdata.max_hot_time.max_time[169] = dwell_156_170
    postdata.max_hot_time.max_time[170] = dwell_156_170

    for n in range(156, 171, 1):
        postdata.max_hot_time.limiting_msid[n] = 'PLINE'

    return postdata


def get_max_times(postdata, maxsec=600000):
    # Use a maximum dwell time in place of Nans.
    indnan = np.isnan(postdata.max_hot_time.max_time)
    maxtimes = postdata.max_hot_time.max_time
    maxtimes[indnan] = maxsec
    return maxtimes


def get_msid_plot_info(cases):
    msids = [
        '1PDEAAT',
        'TCYLAFT6',
        'PFTANK2T',
        'AACCCDPT',
        '4RT700T',
        '1DPAMZT',
        '1DEAMZT',
        'FPTEMP_11',
        'PLINE']
    units = {
        '1PDEAAT': 'C',
        'TCYLAFT6': 'F',
        'PFTANK2T': 'F',
        'AACCCDPT': 'C',
        '4RT700T': 'F',
        '1DPAMZT': 'C',
        '1DEAMZT': 'C',
        'FPTEMP_11': 'C',
        'PLINE': ''}
    colors = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
              '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3', '#444444']
    limits = dict(list(zip(msids, [999] * len(msids))))
    for msid in msids:
        if msid == 'PLINE':
            limits[msid] = 'See Table'
        else:
            for case in cases:
                limits[msid] = np.min(
                    (np.float(case['max_' + msid.lower()]), limits[msid]))
    return msids, units, colors, limits


def gen_composite_plot(
        post_aphelion,
        post_perihelion,
        extra_text='',
        maxy=500000,
        ap_annot_pitch=70,
        peri_annot_pitch=60):

    post_perihelion = add_pline(post_perihelion)
    post_aphelion = add_pline(post_aphelion)

    maxtimes_perihelion = get_max_times(post_perihelion)
    maxtimes_aphelion = get_max_times(post_aphelion)

    msids, units, colors, limits = get_msid_plot_info(cases_aphelion)

    annot = '\n'.join(['{} ({}{})'.format(
        msid, limits[msid], units[msid]) for msid in msids])

    fig = plt.figure(figsize=[14, 7])
    ax = fig.add_axes([0.1, 0.15, 0.8, 0.75])
    ax.fill_between(post_aphelion.pitch_set, maxtimes_aphelion, maxtimes_perihelion,
                    facecolor=[.7, .7, .7], hatch='/', interpolate=True)
    ax.plot(post_aphelion.pitch_set, maxtimes_aphelion,
            color=[.2, .2, .2], linewidth=3)
    ax.plot(post_aphelion.pitch_set, maxtimes_perihelion,
            color=[.2, .2, .2], linewidth=3)

    line_handles = []
    line_labels = []
    for msid, color in zip(msids, colors):
        ind = np.array(post_aphelion.max_hot_time.limiting_msid == msid)
        ax.plot(
            np.array(post_aphelion.pitch_set)[ind],
            maxtimes_aphelion[ind],
            'o',
            color=color,
            markersize=8)
        ind = np.array(post_perihelion.max_hot_time.limiting_msid == msid)
        ax.plot(
            np.array(post_perihelion.pitch_set)[ind],
            maxtimes_perihelion[ind],
            'o',
            color=color,
            markersize=8)
        leg = plt.Line2D((0, 0), (1, 1),
                         color=color,
                         marker='o',
                         linestyle=None)
        line_handles.append(leg)
        line_labels.append('{} ({}{})'.format(msid, limits[msid], units[msid]))

    ax.set_yticks(list(range(0, maxy + 50000, 50000)))
    ax.set_yticklabels(list(range(0, maxy / 1000 + 50, 50)), fontsize=20)
    ax.set_ylabel('Dwell Time (Kiloseconds)', fontsize=24)
    ax.set_ylim(0, maxy)

    ax.set_xticks(list(range(40, 180, 10)))
    ax.set_xticklabels(list(range(40, 180, 10)), fontsize=20)
    ax.set_xlabel('Dwell Pitch', fontsize=24)
    ax.set_xlim(47, 170)

    ap_year = post_aphelion.dates[0][:4]
    ap_arrow_tip = np.interp(
        ap_annot_pitch,
        post_aphelion.pitch_set,
        maxtimes_aphelion) + 3000
    peri_year = post_perihelion.dates[0][:4]
    peri_arrow_tip = np.interp(peri_annot_pitch, post_perihelion.pitch_set,
                               maxtimes_perihelion) - 3000

    ax.annotate('Aphelion {}'.format(ap_year),
                xy=(ap_annot_pitch,
                    ap_arrow_tip),
                xycoords='data',
                xytext=(ap_annot_pitch,
                        ap_arrow_tip + 50000),
                fontsize=16,
                textcoords='data',
                ha='center',
                va='center',
                arrowprops=dict(arrowstyle="->",
                                linewidth=2))
    ax.annotate('Perihelion {}'.format(peri_year),
                xy=(peri_annot_pitch,
                    peri_arrow_tip),
                xycoords='data',
                xytext=(peri_annot_pitch,
                        30000),
                fontsize=16,
                textcoords='data',
                ha='center',
                va='center',
                arrowprops=dict(arrowstyle="->",
                                linewidth=2))

    ax.legend(
        line_handles,
        line_labels,
        loc='upper center',
        fontsize=16,
        numpoints=1)
    ax.grid(True)

    for case in post_perihelion.cases:
        if '1pdeaat' in case['msid'].lower():
            nccd1 = case['n_ccd_dwell1']
            nccd2 = case['n_ccd_dwell2']
            if 'true' in str(case['dh_heater']).lower():
                dh = 'ON'
            else:
                dh = 'OFF'
    if len(extra_text) > 0:
        title_extra_text = ' ({})'.format(extra_text)
        filename_extra_text = '_{}'.format(
            '_'.join(extra_text.split(' ')).lower())
    else:
        title_extra_text = ''
        filename_extra_text = ''

    title = '{} {}CCD DH:{} Composite Dwell Capability{}'.format(
        ap_year, nccd1, dh, title_extra_text)
    ax.set_title(title, fontsize=26, y=1.02)

    filename = 'Composite_{}_{}ccd_dh_{}_max_dwell_range{}.png'.format(
        ap_year, nccd1, dh, filename_extra_text)
    fig.savefig(filename)


if __name__ == '__main__':

    n_sim = 2000
    max_dwell_ksec = 700
    dwell_2_pitch_num = 200

    options = {'TCYLAFT6': {'units': 'F', 'startlim': 108, 'stoplim': 108, 'limstep': 1,
                            'nccd1': [4, ], 'dh': ['ON', ], 'model': 'tcylaft6', 'sameccd': True},
               'PFTANK2T': {'units': 'F', 'startlim': 93, 'stoplim': 93, 'limstep': 1,
                            'nccd1': [4, ], 'dh': ['ON', ], 'model': 'tank', 'sameccd': True},
               'AACCCDPT': {'units': 'C', 'startlim': -12.5, 'stoplim': -12.5, 'limstep': 1,
                            'nccd1': [4, ], 'dh': ['ON', ], 'model': 'aca', 'sameccd': True},
               '4RT700T': {'units': 'F', 'startlim': 84, 'stoplim': 84, 'limstep': 1,
                           'nccd1': [4, ], 'dh': ['ON', ], 'model': 'fwdblkhd', 'sameccd': True},
               '1PDEAAT': {'units': 'C', 'startlim': 52.5, 'stoplim': 52.5, 'limstep': 0.5,
                           'nccd1': [4, ], 'dh': ['ON', ], 'model': 'psmc', 'sameccd': True},
               '1DPAMZT': {'units': 'C', 'startlim': 35.5, 'stoplim': 35.5, 'limstep': 0.5,
                           'nccd1': [4, ], 'dh': ['ON', ], 'model': 'dpa', 'sameccd': True},
               '1DEAMZT': {'units': 'C', 'startlim': 35.5, 'stoplim': 35.5, 'limstep': 0.5,
                           'nccd1': [4, ], 'dh': ['ON', ], 'model': 'dea', 'sameccd': True},
               'FPTEMP_11': {'units': 'C', 'startlim': -114, 'stoplim': -114, 'limstep': 1,
                             'nccd1': [4, ], 'dh': ['ON', ], 'model': 'acisfp', 'sameccd': True}}

    dates = [year + ':' + day for day in ['001', ] for year in ['2016', ]]
    cases_perihelion = assemble_cases(
        options, n_sim, max_dwell_ksec, dwell_2_pitch_num, dates)

    dates = [year + ':' + day for day in ['182', ] for year in ['2016', ]]
    cases_aphelion = assemble_cases(
        options,
        n_sim,
        max_dwell_ksec,
        dwell_2_pitch_num,
        dates)

    post_perihelion = pyger.PostPyger(
        cases_perihelion, '/Users/wolf809/AXAFAUTO/Pyger_Prod/pygerdata/')
    post_aphelion = pyger.PostPyger(
        cases_aphelion,
        '/Users/wolf809/AXAFAUTO/Pyger_Prod/pygerdata/')

    gen_composite_plot(
        post_aphelion,
        post_perihelion,
        extra_text='July 2016 Limits',
        maxy=500000)
