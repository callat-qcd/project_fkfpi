import numpy as np
import gvar as gv
import lsqfit
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import time

sys.path.append("../")
from fitter import data_loader as dl
from fitter import bootstrapper as bs

p_dict = {
    'order' : {
        'fit' : 'nnlo', # 'nlo', 'nnlo', or 'nnnlo'
        'vol' : 10, # max 10
        'include_log' : True, # alphaS
        'include_log2' : True
    },
    'bias_correct' : True,

    'abbrs' : [u'a09m220', u'a09m310', u'a09m350', u'a09m400', u'a12m130',
           u'a12m220',  u'a12m220S', u'a12m220L', u'a12m310', u'a12m350',
           u'a12m400',  u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'
}

for fit_type in ['ma-ratio']:
    for F2 in ['FKFK', 'FKFpi', 'FpiFpi']:
        for vol in [0, 10]:
            for include_alphaS in [False, True]:
                for include_logSq in [False, True]:

                    t0 = time.time()

                    # Loop through combinations
                    p_dict['order']['vol'] = vol
                    p_dict['order']['include_log'] = include_alphaS
                    p_dict['order']['include_log2'] = include_logSq
                    p_dict['F2'] = F2
                    p_dict['fit_type'] = fit_type

                    if vol > 0:
                        include_FV = True
                    else:
                        include_FV = False

                    name = fit_type +"_"+F2
                    if include_FV:
                        name = name + '_FV'
                    if include_alphaS:
                        name = name + '_alphaS'
                    if include_logSq:
                        name = name + '_logSq'

                    print name
                    print "(model, F2, Vol, alphaS, logSq):", (p_dict['fit_type'], p_dict['F2'], include_FV, include_alphaS, include_logSq)

                    # Load data
                    data_loader = dl.data_loader()
                    fit_data = data_loader.get_fit_data()

                    temp_prior = data_loader.get_prior(fit_type=p_dict['fit_type'], F2=p_dict['F2'],
                                                       include_FV=include_FV, include_alphaS=include_alphaS, include_logSq=include_logSq)

                    if temp_prior is None:
                        bootstrapper = bs.bootstrapper(fit_data, prior=None, order=p_dict['order'], F2=p_dict['F2'],
                                                       fit_type=p_dict['fit_type'], abbrs=p_dict['abbrs'], bias_correct=p_dict['bias_correct'])

                        new_prior = bootstrapper.create_prior_from_fit()
                        data_loader.save_prior(new_prior, fit_type=p_dict['fit_type'], F2=p_dict['F2'],
                                               include_FV=include_FV, include_alphaS=include_alphaS, include_logSq=include_logSq)

                    t1 = time.time()


                    print "\nTotal time (s): ", t1 - t0, "\n"
