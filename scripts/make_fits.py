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
    'bias_correct' : True,
    'include_su2_isospin_corrrection' : False,

    'use_prior' : False,
    'abbrs' : [u'a06m310L', u'a09m220', u'a09m310', u'a09m350', u'a09m400', u'a12m130',
               u'a12m220',  u'a12m220S', u'a12m220L', u'a12m310', u'a12m350',
               u'a12m400',  u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'

    'save_results' : True,
    'save_pickles' : False,
    'replace_fits' : True,
}

if p_dict['save_results']:
    if len(sys.argv) > 1:
        p_dict['output_name'] = sys.argv[1]
    else:
        p_dict['output_name'] = input('Name for fit collection: ')


choices = {
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],
    'vol' : [0, 10],

    # semi-nnlo corrections
    'include_alpha_s' : [False, True],
    'semi-nnlo_corrections' : ['nnlo-full_bijnens', 'nnlo-full', 'nnlo-logSq_ct', 'nnlo-ct'],

    # nnnlo corrections
    'include_latt_n3lo' : [False, True],
}

t0_all = time.time()

# Get all enumerations of these choices
for j, choice in enumerate(dict(zip(choices, x)) for x in itertools.product(*choices.values())):
    print(choice)
    print(j+1, 'of', len(list(itertools.product(*[choices[key] for key in list(choices)]))))

    p_dict['fit_type'] = choice['fit_type']
    p_dict['F2'] = choice['F2']

    p_dict['order']['vol'] = choice['vol']
    p_dict['order']['include_alpha_s'] = choice['include_alpha_s']
    p_dict['order']['include_latt_n3lo'] = choice['include_latt_n3lo']

    if choice['semi-nnlo_corrections'] == 'nnlo-full_bijnens':
        p_dict['order']['include_log'] = True
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = True

    elif choice['semi-nnlo_corrections'] == 'nnlo-full':
        p_dict['order']['include_log'] = True
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-logSq_ct':
        p_dict['order']['include_log'] = False
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-ct':
        p_dict['order']['include_log'] = False
        p_dict['order']['include_log2'] = False
        p_dict['order']['include_sunset'] = False
        p_dict['use_bijnens_central_value'] = False


    t0 = time.time()

    # Load data
    data_loader = dl.data_loader()
    fit_data = data_loader.get_fit_data()
    phys_point_data = data_loader.get_phys_point_data()

    # Get prior
    prior = None
    if p_dict['use_prior']:
        prior = data_loader.get_prior(fit_type=p_dict['fit_type'], order=p_dict['order']['fit'], F2=p_dict['F2'],
                      include_log=p_dict['order']['include_log'], include_log2=p_dict['order']['include_log2'],
                      include_sunset=p_dict['order']['include_sunset'], include_alpha_s=p_dict['order']['include_alpha_s'],
                      include_latt_n3lo=p_dict['order']['include_latt_n3lo'], include_FV=(p_dict['order']['vol'] > 6),
                      use_bijnens_central_value=p_dict['use_bijnens_central_value']
                 )

    # Make bootstrapper
    bootstrapper = bs.bootstrapper(
        fit_data, phys_point_data, prior=prior, order=p_dict['order'], F2=p_dict['F2'],
        include_su2_isospin_corrrection=p_dict['include_su2_isospin_corrrection'], use_bijnens_central_value=p_dict['use_bijnens_central_value'],
        fit_type=p_dict['fit_type'], abbrs=p_dict['abbrs'], bias_correct=p_dict['bias_correct']
    )

    if p_dict['replace_fits'] or not (bootstrapper.get_name() in data_loader.get_models()):

        print(bootstrapper)

        # Save results
        data_loader.save_fit_info(bootstrapper.get_fit_info(),
                                  output_name=p_dict['output_name'], save_pickles=p_dict['save_pickles'])


    t1 = time.time()

    print("\nTotal time (s): ", t1 - t0, "\n")


t1_all = time.time()
print("\nTotal time for all fits (s): ", t1_all - t0_all, "\n")
