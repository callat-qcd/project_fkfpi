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
import itertools

sys.path.append("../")
import fitter.data_loader as dl
import fitter.fit_manager as fm
import fitter.fitter as fit
import fitter.special_functions as sf

for j in range(10): # Sometimes this needs to be loaded twice...
    matplotlib.rcParams['figure.figsize'] = [10, 10]

params = {
    'use_prior' : True,
    'bias_correct' : True,
    'fast_sunset' : False,

    'abbrs' : [u'a06m310L',
               u'a09m135', u'a09m220', u'a09m310', u'a09m350', u'a09m400',
               u'a12m130', u'a12m220', u'a12m220S', u'a12m220L', u'a12m310', u'a12m350', u'a12m400',
               u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'

    'save_results' : True, # If fast_sunset = True, this should be set to False
    'save_pickles' : True,
    'replace_fits' : False,
    'make_plots' : False,
}

if params['save_results']:
    if len(sys.argv) > 1:
        params['collection_name'] = sys.argv[1]
    else:
        params['collection_name'] = input('Name for fit collection: ')

t0_all = time.time()

data_loader = dl.data_loader()
model_list = data_loader.get_model_names(params['collection_name'])
finished_models = []
if data_loader.get_fit_results(params['collection_name']) is not None:
    finished_models = list(data_loader.get_fit_results(params['collection_name']))

# Get all enumerations of these choices
for j, model in enumerate(model_list):
    print('\n---\nMaking model', j+1, 'of', len(model_list))

    t0 = time.time()

    # Load data
    fit_data = data_loader.get_fit_data()
    phys_point_data = data_loader.get_phys_point_data()
    params['model_info'] = data_loader.get_model_info_from_name(model)

    # Get prior
    prior = None
    if params['use_prior']:
        prior = data_loader.get_prior(
            collection_name=params['collection_name'],
            model_info=params['model_info']
        )

    # Make fit_manager
    fit_manager = fm.fit_manager(fit_data, phys_point_data, prior=prior, **params)

    if (params['replace_fits']
            or finished_models is None
            or not (fit_manager.model in finished_models)
        ):

        print(fit_manager)

        # Save results
        data_loader.save_fit_info(fit_manager.fit_info,
                                  collection_name=params['collection_name'], save_pickles=params['save_pickles'])

    if params['make_plots']:
        data_loader.save_plots(
            fit_manager.make_plots(
                show_error_ellipses=params['show_error_ellipses'],
                show_bootstrap_histograms=params['show_bs_histograms']),
            output_filename='fits/'+name
        )


    t1 = time.time()

    print("\nTotal time (s): ", t1 - t0, "\n")


t1_all = time.time()
print("\nTotal time for all fits (s): ", t1_all - t0_all, "\n")
