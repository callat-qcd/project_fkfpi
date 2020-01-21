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
from fitter import fitter as fit
from fitter import special_functions as sf

for j in range(10): # Sometimes this needs to be loaded twice...
    matplotlib.rcParams['figure.figsize'] = [10, 10]

p_dict = {
    'order' : {
        'fit' : 'nnlo', # 'nlo', 'nnlo', or 'nnnlo'
        'exclude' : [], # put LECs here
    },
    'bs_N' : 1,  # if 0, use full list
    'chained' : False,
    'bias_correct' : True,

    'use_prior' : True,
    'abbrs' : [u'a06m310L', u'a09m220', u'a09m310', u'a09m350', u'a09m400', u'a12m130',
           u'a12m220',  u'a12m220S', u'a12m220L', u'a12m310', u'a12m350',
           u'a12m400',  u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'

    'make_plots' : True,
    'show_error_ellipses' : False,
    'show_bs_histograms' : False,

    'save_results' : False,
}

for fit_type in ['ma', 'ma-ratio', 'xpt', 'xpt-ratio']:
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

                    print(name)
                    print("(model, F2, Vol, alphaS, logSq):", (p_dict['fit_type'], p_dict['F2'], include_FV, include_alphaS, include_logSq))

                    # Load data
                    data_loader = dl.data_loader()
                    fit_data = data_loader.get_fit_data()

                    # Get prior
                    prior = None
                    if p_dict['use_prior']:
                        prior = data_loader.get_prior(fit_type=p_dict['fit_type'], F2=p_dict['F2'], include_FV=(p_dict['order']['vol']>0),
                                                     include_alphaS=p_dict['order']['include_log'], include_logSq=p_dict['order']['include_log2'])

                    # Make bootstrapper
                    bootstrapper = bs.bootstrapper(fit_data, prior=prior, order=p_dict['order'], F2=p_dict['F2'], chain_fits=p_dict['chained'],
                                                   fit_type=p_dict['fit_type'], bs_N=p_dict['bs_N'], abbrs=p_dict['abbrs'], bias_correct=p_dict['bias_correct'],)

                    if p_dict['make_plots']:
                        data_loader.save_plots(
                            bootstrapper.make_plots(
                                show_error_ellipses=p_dict['show_error_ellipses'],
                                show_bootstrap_histograms=p_dict['show_bs_histograms']),
                            output_filename='fits/'+name
                        )
                    else:
                        print(bootstrapper)

                    if p_dict['save_results']:
                        data_loader.save_fit_info(bootstrapper.get_fit_info())


                    t1 = time.time()

                    print("\nTotal time (s): ", t1 - t0, "\n")
