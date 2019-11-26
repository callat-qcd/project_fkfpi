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
from fitter import misc as ms
from fitter import model_average as md

for j in range(10): # Sometimes this needs to be loaded twice...
    matplotlib.rcParams['figure.figsize'] = [10, 10]

data_loader = dl.data_loader()
fit_results = data_loader.get_fit_info('fit_results.csv')
other_results = data_loader.get_fit_info('other_results.csv')

model_average = md.model_average(fit_results)


fig = model_average.plot_histogram('FK/Fpi')
data_loader.save_plots(fig, output_filename='histogram_fit')

fig = model_average.plot_histogram('delta_su2')
data_loader.save_plots(fig, output_filename='histogram_delta_su2')


fig = model_average.plot_comparison(other_results=other_results, param='FK/Fpi_pm')
data_loader.save_plots(fig, output_filename='comparison_fits')

fig = model_average.plot_comparison(other_results=other_results, param='delta_su2')
data_loader.save_plots(fig, output_filename='comparison_delta_su2')


fig = model_average.plot_fits('a')
data_loader.save_plots(fig, output_filename='all_fits_vs_latt_spacing')

fig = model_average.plot_fits('mpi')
data_loader.save_plots(fig, output_filename='all_fits_vs_mpi')

fig = model_average.plot_fits('volume')
data_loader.save_plots(fig, output_filename='all_fits_vs_volume')
