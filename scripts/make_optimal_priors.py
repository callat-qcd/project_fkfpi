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



params = {
    'bias_correct' : True,
    'fast_sunset' : True,

    'abbrs' : [u'a06m310L', u'a09m135', u'a09m220', u'a09m310', u'a09m350', u'a09m400', u'a12m130',
           u'a12m220',  u'a12m220S', u'a12m220L', u'a12m310', u'a12m350',
           u'a12m400',  u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'

    'save_results' : True, # If fast_sunset = True, this should be set to False
    'save_pickles' : False,
    'replace_entries' : False,
}

if len(sys.argv) > 1:
    params['collection_name'] = sys.argv[1]
else:
    params['collection_name'] = input('Name for fit collection: ')


t0_all = time.time()

data_loader = dl.data_loader()
model_list = data_loader.get_model_names(params['collection_name'])

for model in model_list:
    t0 = time.time()

    # Load data
    fit_data = data_loader.get_fit_data()
    phys_point_data = data_loader.get_phys_point_data()
    params['model_info'] = data_loader.get_model_info_from_name(model)

    # See if a prior has been generated; if not, create one
    temp_prior = data_loader.get_prior(collection_name=params['collection_name'], model_info=params['model_info'])

    if temp_prior is None or params['replace_entries']:

        # Make fit_manager
        fit_manager = fm.fit_manager(fit_data, phys_point_data, prior=temp_prior, **params)

        new_prior = fit_manager.create_prior_from_fit()
        data_loader.save_prior(
            collection_name=params['collection_name'],
            prior=new_prior,
            name=fit_manager.model
        )

        if params['save_results']:
            data_loader.save_fit_info(
                fit_manager.fit_info,
                collection_name=params['collection_name'],
                save_pickles=params['save_pickles']
            )

    t1 = time.time()

    print("\nTotal time (s): ", t1 - t0, "\n")

t1_all = time.time()
print("\nTotal time for all fits (s): ", t1_all - t0_all, "\n")
