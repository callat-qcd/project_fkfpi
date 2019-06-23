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

    def get_ensembles(self):
        with h5py.File(self.file_h5, "r") as f:
            ensembles = f.keys()
        return sorted(ensembles)

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

    # whose: 'mine', 'others'
    def get_fit_info(self, whose=None):
        filepath = None
        if whose=='mine' or whose is None:
            filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')
        elif whose=='others':
            filepath = os.path.normpath(self.project_path + '/results/other_results.csv')

        if not os.path.isfile(filepath):
            return None


        df_fit = pd.read_csv(filepath, header=0)
        cols = ['fit', 'logGBF', 'chi2/df', 'Q', 'vol corr', 'latt corr']
        cols = np.intersect1d(cols, df_fit.columns.values)

        fit_types = df_fit['name'].values
        output_dict = {}
        for name in fit_types:
            index = np.argwhere(df_fit['name'].values == name)

            output_dict[name]= {key: np.asscalar(df_fit[key].values[index]) for key in cols}

        return output_dict

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

    def save_fit_info(self, bootstrapper):
        print "Saving..."

        if not os.path.exists(self.project_path + '/results/'):
            os.makedirs(self.project_path + '/results/')
        filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')

        fit_type = bootstrapper.fit_type+'_'+bootstrapper.order['fit']
        fit_info = {
            'name' : bootstrapper.fit_type+'_'+bootstrapper.order['fit'],
            'fit' : bootstrapper.extrapolate_to_phys_point(),
            'logGBF' : bootstrapper.fits[0].logGBF,
            'chi2/df' : bootstrapper.fits[0].chi2 / bootstrapper.fits[0].dof,
            'Q' : bootstrapper.fits[0].Q,
            'vol corr' : bootstrapper.order['vol'],
            'latt corr' : bootstrapper.order['latt_spacing']
        }
        cols = ['name', 'fit', 'logGBF', 'chi2/df', 'Q', 'vol corr', 'latt corr']

        #return fit_info

        # Add result to file if file exists
        if os.path.isfile(filepath):
            df_best_fits = pd.read_csv(filepath, index_col=0).to_dict()

            output_dict = {}
            for key in df_best_fits.keys():
                output_dict[key] = np.array(df_best_fits[key].values())

            if fit_type in output_dict['name']:
                index = np.asscalar(np.argwhere(output_dict['name'] == fit_type))
                for key in fit_info.keys():
                    output_dict[key][index] = fit_info[key]
            else:
                # We want the abbr to be listed in alphabetical order
                index = np.asscalar(np.argwhere(np.sort(np.append(output_dict['name'], fit_type), axis=None)
                                  ==  fit_type))
                for key in output_dict.keys():
                    if index > len(output_dict[key]) - 1:
                        output_dict[key]= np.append(output_dict[key], fit_info[key])
                    else:
                        output_dict[key] = np.insert(output_dict[key], index, fit_info[key])

            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols]
            df.to_csv(filepath)

        # Create new file if it doesn't exist
        else:
            output_dict = {key : [fit_info[key]] for key in fit_info.keys()}
            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols] # rearrange in logical order
            df.to_csv(filepath)

        return None

    def save_plots(self, figs=None, bootstrapper=None, output_filename=None):
        if figs is None:
            figs = bootstrapper.make_plots()

        if not os.path.exists(os.path.normpath(self.project_path+'/tmp/')):
            os.makedirs(os.path.normpath(self.project_path+'/tmp/'))

        if output_filename is None:
            output_file = os.path.normpath(self.project_path+'/tmp/temp.pdf')
        else:
            output_file = os.path.normpath(self.project_path+'/tmp/'+output_filename+'.pdf')

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
