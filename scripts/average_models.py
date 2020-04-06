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
import fitter.model_average as md
import fitter.misc as ms
import fitter.data_loader as dl
import fitter.fit_manager as fm

for j in range(10): # Sometimes this needs to be loaded twice...
    matplotlib.rcParams['figure.figsize'] = [10, 10]

# Promp for output
if len(sys.argv) > 1:
    collection_name = sys.argv[1]
else:
    collection_name = input('Specify collection to average: ')

data_loader = dl.data_loader()

fit_results = data_loader.get_fit_results(collection_name)
other_results = data_loader.get_fit_results('other_collabs')
model_average = md.model_average(fit_results)

# Save fit info to /results/{collection_name}/README.md
data_loader.save_model_average(collection_name, str(model_average))

# Make histograms
for vary_choice in ['fit_type', 'F2', 'include_alpha_s', 'latt_ct', 'include_FV', 'semi-n2lo_corrections']:
    fig = model_average.plot_histogram('FK/Fpi', vary_choice=vary_choice)
    data_loader.save_plots(fig, output_filename=collection_name+'/figs/histogram_fit_'+vary_choice)

fig = model_average.plot_histogram('delta_su2')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/histogram_delta_su2')


fig = model_average.plot_comparison(other_results=other_results, param='FK/Fpi_pm')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/comparison_fits')

fig = model_average.plot_comparison(other_results=other_results, param='delta_su2')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/comparison_delta_su2')


fig = model_average.plot_fits('a')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/all_fits_vs_latt_spacing')

fig = model_average.plot_fits('mpi')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/all_fits_vs_mpi')

fig = model_average.plot_fits('volume')
data_loader.save_plots(fig, output_filename=collection_name+'/figs/all_fits_vs_volume')

for j, model in enumerate(model_average.get_model_names(by_weight=True)[:5]):
    print('Making fits for model:', model)
    params = {
        'use_prior' : True,
        'bias_correct' : True,
        'fast_sunset' : False,

        'abbrs' : [u'a06m310L',
                   u'a09m220', u'a09m310', u'a09m350', u'a09m400', #u'a09m135',
                   u'a12m130', u'a12m220', u'a12m220S', u'a12m220L', u'a12m310', u'a12m350', u'a12m400',
                   u'a15m135XL', u'a15m220', u'a15m310', u'a15m350', 'a15m400'], # u'a15m130'
    }

    data_loader = dl.data_loader()
    fit_data = data_loader.get_fit_data()
    phys_point_data = data_loader.get_phys_point_data()
    params['model_info'] = data_loader.get_model_info_from_name(model)
    prior = data_loader.get_prior(
        collection_name=collection_name,
        model_info=params['model_info']
    )
    fit_manager = fm.fit_manager(fit_data, phys_point_data, prior=prior, **params)

    fig = fit_manager.plot_fit('a')
    data_loader.save_plots(fig, output_filename=collection_name+'/figs/fits/'+str(j+1)+'-'+fit_manager.model+'-vs_a')

    fig = fit_manager.plot_fit('mpi')
    data_loader.save_plots(fig, output_filename=collection_name+'/figs/fits/'+str(j+1)+'-'+fit_manager.model+'-vs_mpi')
