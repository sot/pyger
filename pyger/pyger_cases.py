
import numpy as np
import pickle
import csv
import pandas as pd
import os

from Chandra.Time import DateTime
import pyger


def read_cases(csv_casefile):
    cases = []
    with open(csv_casefile,'r') as f:
        reader = csv.reader(f)
        headers = reader.next()
        for row in reader:
            cases.append(dict(zip(headers,row)))
    return cases


def run_pyger_cases(cases, savedwell1=False):
    
    for case in cases:

        models = (case['constraint_model'],)
        msids = (case['msid'].lower(),)
        constraints1 = pyger.calc_constraints(start=case['start'], 
                                              max_tephin=float(case['max_tephin']),
                                              max_tcylaft6=float(case['max_tcylaft6']),
                                              max_1pdeaat=float(case['max_1pdeaat']),
                                              max_1dpamzt=float(case['max_1dpamzt']),
                                              max_pftank2t=float(case['max_pftank2t']),
                                              max_aacccdpt=float(case['max_aacccdpt']),
                                              n_ccd=int(case['n_ccd_dwell1']),
                                              n_sim=int(case['n_sim']),
                                              max_dwell_ksec=float(case['max_dwell_ksec']),
                                              constraint_models=models)
        if savedwell1:
            pyger.save_pyger_pickle(constraints1, case['filename'] + '_dwell1.pkl')
            print('Saving to {0}'.format(case['filename'] + '_dwell1.pkl'))

        constraints2, coolstats, hotstats = pyger.calc_constraints2(
                                                    constraints1,
                                                    start=case['start'],
                                                    max_dwell_ksec=float(case['max_dwell_ksec']),
                                                    pitch_num=int(case['dwell_2_pitch_num']),
                                                    hot_dwell_temp_ratio=1.0,
                                                    T_cool_ratio=0.9,
                                                    constraint_models=models,
                                                    msids=msids,
                                                    n_ccd=int(case['n_ccd_dwell2']))
        filename = case['filename'] + '_dwell2.pkl'
        pickle.dump((constraints2, coolstats, hotstats), open(filename,'w'), protocol=2)
        print('Saving to {0}'.format(filename))


class PostPyger(object):
    """ Consolidate results from a set of models

    """

    def __init__(self, cases, pickledir):
        """
        :param cases: list of cases to include in post processing run.

        There should be one case for each MSID, with no conflicting states (e.g. 5 ACIS ccds for
        one msid, and 4 for another)
        """
        self.cases = cases
        self.pickledir = pickledir
        self.dates = np.array([case['start'] for case in cases])
        self.date_set = set(self.dates)
        self.pitch_set = range(46,171,1)
        
        panel_dict= {}
        for case in self.cases:
            frame = self.get_case_frame(case)
            panel_dict[case['msid']] = frame
            print('Imported info for {0}'.format(case['msid']))
                  
        self.constraint_panel = pd.Panel(panel_dict)
        self.calc_stats()
        self.get_constraint_info()
        

    def get_constraint_info(self):
        self.info = {'date':self.cases[0]['start'],
                     'models':[]}
        for case in self.cases:
            msid = case['msid'].lower()
            self.info[msid] = case['max_' + msid]
            if msid == '1dpamzt' or msid == '1pdeaat':
                self.info['dwell1_ccds'] = case['n_ccd_dwell1']
                self.info['dwell2_ccds'] = case['n_ccd_dwell2']                          
            self.info['models'].append(case['constraint_model'])
     

    def get_case_frame(self, case):
        filename = os.path.join(self.pickledir, (case['filename'] + '_dwell1.pkl'))
        dwell1data = pickle.load(open(filename,'r'))
        filename = os.path.join(self.pickledir, (case['filename'] + '_dwell2.pkl'))
        data, coolstats, hotstats = pickle.load(open(filename,'r'))

        msid = case['msid'].lower()
        model = case['constraint_model']
        max_dwell_sec = dwell1data[model]['max_dwell_ksec'] * 1000

        if len(data) > 0:
            
            dwell1_stats = dwell1data[model]['dwell1_stats']
            dur90 = np.interp(self.pitch_set, dwell1_stats.pitch, dwell1_stats.dur90)
            dur50 = np.interp(self.pitch_set, dwell1_stats.pitch, dwell1_stats.dur50)
            hot_delta = np.interp(self.pitch_set, hotstats[msid].pitch, 
                                  hotstats[msid].dwell1_duration_delta, left=np.NaN, right=np.NaN)
            dur90[dur90 > max_dwell_sec*0.95] = np.NaN
            dur50[dur90 > max_dwell_sec*0.95] = np.NaN

            # For anyone reading this code, dur90 and cool90 are not apples to apples comparable.
            # dur90 represents the "best case" capability generated by the original Pyger methodology.
            # cool90 represents the time it takes to cool down 90% from the max time calculated
            # separately in the cooldown portion of the pyger code.
            cool90 = np.interp(self.pitch_set, coolstats[msid].pitch, coolstats[msid].perc90, 
                               left=np.NaN, right=np.NaN)
            cool90[np.isnan(cool90)] = max_dwell_sec
        else:
            # these conventions may change, keeping some nans and other max dwell sec doesn't make sense
            dur50 = np.array([np.nan] * len(self.pitch_set))
            dur90 = np.array([np.nan] * len(self.pitch_set))
            hot_delta = np.array([np.nan] * len(self.pitch_set))
            cool90 = np.array([max_dwell_sec] * len(self.pitch_set))

        return pd.DataFrame({'avg_max_hot_time':dur50,
                             'max_hot_time':dur90,
                             'hot_dwell':hot_delta,
                             'cool_dwell':cool90},
                             index=self.pitch_set)
    

    def calc_stats(self):
        avg_max_hot_time = self.constraint_panel.minor_xs('avg_max_hot_time').min(axis=1)
        min_max_hot_time = self.constraint_panel.minor_xs('max_hot_time').min(axis=1)
        min_max_hot_time_msid = self.constraint_panel.minor_xs('max_hot_time').idxmin(axis=1)
        self.max_hot_time = pd.DataFrame({'max_time':min_max_hot_time, 
                                          'avg_max_time':avg_max_hot_time,
                                          'limiting_msid':min_max_hot_time_msid})
        
        min_hot_dwell = self.constraint_panel.minor_xs('hot_dwell').min(axis=1)
        min_hot_dwell_msid = self.constraint_panel.minor_xs('hot_dwell').idxmin(axis=1)
        self.hot_dwell = pd.DataFrame({'heating_time':min_hot_dwell,
                                       'limiting_msid':min_hot_dwell_msid})

        self.cool_dwells = self.constraint_panel.minor_xs('cool_dwell').max(axis=1)


