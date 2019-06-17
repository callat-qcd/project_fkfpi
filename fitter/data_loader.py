import pandas as pd
import numpy as np
import gvar as gv
import sys
import time
import os
import h5py
from matplotlib.backends.backend_pdf import PdfPages

class data_loader(object):

    def __init__(self):
        self.project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
        self.file_h5 = os.path.normpath(self.project_path+'/FK_Fpi_data.h5')

    def get_fit_data(self):
        data = {}
        with h5py.File(self.file_h5, "r") as f:
            for ensemble in f.keys():
                data[ensemble] = {}
                for key in f[ensemble].keys():
                    dset = f[ensemble][key]
                    shape = dset.shape
                    if shape != ():
                        data[ensemble][key] = dset[:]
                    else:
                        data[ensemble][key] = dset[()]
        return data

    def get_ensembles(self):
        with h5py.File(self.file_h5, "r") as f:
            ensembles = f.keys()
        return sorted(ensembles)

    def get_prior(self, fit_type):
        filepath = self.project_path + '/priors/'+fit_type+'.csv'
        df_read = pd.read_csv(filepath, index_col=0)
        return gv.gvar({key : df_read.to_dict("index")[key]["0"]
                for key in df_read.to_dict("index").keys()})

    def get_variable_names(self):
        names = []
        with h5py.File(self.file_h5, "r") as f:
            for ensemble in f.keys():
                for key in f[ensemble].keys():
                    names.append(key)
        return sorted(np.unique([names]))

    def save_plots(self, figs=None, bootstrapper=None, output_file=None):
        if figs is None:
            figs = bootstrapper.make_plots()

        if output_file is None:
            if not os.path.exists(os.path.normpath(self.project_path+'/tmp/')):
                os.makedirs(os.path.normpath(self.project_path+'/tmp/'))
            output_file = os.path.normpath(self.project_path+'/tmp/temp.pdf')

        output_pdf = PdfPages(output_file)
        try:
            for fig in figs:
                output_pdf.savefig(fig)
        except TypeError: # save figs directly if figs is just a figure, not an array
            output_pdf.savefig(figs)

        output_pdf.close()
        print "Done."

    def save_prior(self, prior, fit_type):
        filepath = self.project_path + '/priors/'+fit_type+'.csv'
        out_prior = {}
        for key in sorted(prior.keys()):
            out_prior[key] = [str(prior[key])]

        df = pd.DataFrame.from_dict(out_prior).T
        df.to_csv(filepath)
