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

# Make histograms
for vary_choice in ['fit_type', 'F2', 'include_alpha_s', 'include_latt_n3lo', 'include_FV', 'semi-nnlo_corrections']:
    fig = model_average.plot_histogram('FK/Fpi', vary_choice=vary_choice)
    data_loader.save_plots(fig, output_filename=output_name+'/histogram_fit_'+vary_choice)

fig = model_average.plot_histogram('delta_su2')
data_loader.save_plots(fig, output_filename=output_name+'/histogram_delta_su2')


fig = model_average.plot_comparison(other_results=other_results, param='FK/Fpi_pm')
data_loader.save_plots(fig, output_filename=output_name+'/comparison_fits')

fig = model_average.plot_comparison(other_results=other_results, param='delta_su2')
data_loader.save_plots(fig, output_filename=output_name+'/comparison_delta_su2')


fig = model_average.plot_fits('a')
data_loader.save_plots(fig, output_filename=output_name+'/all_fits_vs_latt_spacing')

fig = model_average.plot_fits('mpi')
data_loader.save_plots(fig, output_filename=output_name+'/all_fits_vs_mpi')

fig = model_average.plot_fits('volume')
data_loader.save_plots(fig, output_filename=output_name+'/all_fits_vs_volume')
