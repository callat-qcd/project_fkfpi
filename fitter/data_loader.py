import pandas as pd
import numpy as np
import gvar as gv
import sys
import time
import os
import h5py
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict

class data_loader(object):

    def __init__(self):
        self.project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
        self.file_h5 = os.path.normpath(self.project_path+'/FK_Fpi_data.h5')

    def _pickle_fit_parameters(self, fit_info):
        name = fit_info['name']
        filename = self.project_path+'/pickles/'+name+'.p'

        output = {}
        for key in fit_info['prior'].keys():
            output[key] = [fit_info['prior'][key], fit_info['posterior'][key]]

        gv.dump(output, filename)
        return None

    def _unpickle_fit_parameters(self, name):
        filepath = self.project_path+'/pickles/'+name+'.p'
        if os.path.isfile(filepath):
            return gv.load(filepath)
        else:
            return None

    def get_ensembles(self):
        with h5py.File(self.file_h5, "r") as f:
            ensembles = list(f.keys())

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
    def get_fit_info(self, filename=None):
        if filename is None:
            filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')
        else:
            filepath = os.path.normpath(self.project_path + '/results/'+filename)

        if not os.path.isfile(filepath):
            return None


        df_fit = pd.read_csv(filepath, header=0)
        cols = df_fit.columns.values
        #cols = ['fit', 'logGBF', 'chi2/df', 'Q', 'vol corr', 'latt corr']
        #cols = np.intersect1d(cols, df_fit.columns.values)

        fit_types = df_fit['name'].values
        output_dict = OrderedDict()
        for name in sorted(fit_types):
            index = np.argwhere(df_fit['name'].values == name)

            output_dict[name]= {key: np.asscalar(df_fit[key].values[index]) for key in cols}

            fit_params = self._unpickle_fit_parameters(name)
            output_dict[name]['prior'] = {key : fit_params[key][0] for key in fit_params.keys()}
            output_dict[name]['posterior'] = {key : fit_params[key][1] for key in fit_params.keys()}

        return output_dict

    def get_prior(self, fit_type, F2, include_FV, include_alphaS, include_logSq):
        filepath = os.path.normpath(self.project_path + '/priors/'+fit_type+'.csv')

        print(filepath)

        if not os.path.isfile(filepath):
            return None

        name = F2
        if include_FV:
            name = name + '_FV'
        if include_alphaS:
            name = name + '_alphaS'
        if include_logSq:
            name = name + '_logSq'

        df_prior = pd.read_csv(filepath, header=0)
        if name not in df_prior['name'].values:
            return None

        index = np.argwhere(df_prior['name'].values == name)
        prior = gv.BufferDict()
        for key in ['L_4', 'L_5', 'A_k', 'A_p', 'A_a', 'A_loga']:
            value = np.asscalar(df_prior[key].values[index])
            if value is not np.nan:
                prior[key] = gv.gvar(value)

        return prior

    def get_variable_names(self):
        names = []
        with h5py.File(self.file_h5, "r") as f:
            for ensemble in f.keys():
                for key in f[ensemble].keys():
                    names.append(key)
        return sorted(np.unique([names]))

    # pickle correlated prior/posterior,
    # save rest to csv file
    def save_fit_info(self, fit_info):
        print("Saving...")

        self._pickle_fit_parameters(fit_info)

        if not os.path.exists(self.project_path + '/results/'):
            os.makedirs(self.project_path + '/results/')
        filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')

        # get fit info
        cols = np.array(['name', 'FK/Fpi', 'delta_su2', 'logGBF', 'chi2/df', 'Q', 'vol'])
        lecs_cols = ['L_4', 'L_5', # nlo terms
                     'L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9', #semi-nnlo terms
                     'A_a', 'A_k', 'A_p', 'A_loga'] # nnlo terms
                     #'A_aa', 'A_ak', 'A_ap', # nnnlo terms
                     #'A_kk', 'A_kp', 'A_pp'] # more nnnlo terms

        posterior_params = {}
        for key in lecs_cols:
            if key in fit_info['posterior'].keys():
                posterior_params[key] = fit_info['posterior'][key]
            else:
                posterior_params[key] = np.nan


        # Add result to file if file exists
        if os.path.isfile(filepath):
            df_best_fits = pd.read_csv(filepath, index_col=0).to_dict()

            output_dict = {}
            for key in df_best_fits.keys():
                output_dict[key] = np.array(list(df_best_fits[key].values()), dtype="object")

            if fit_info['name'] in output_dict['name']:
                index = np.argwhere(output_dict['name'] == fit_info['name']).item()

                for key in cols:
                    output_dict[key][index] = fit_info[key]
                for key in lecs_cols:
                    output_dict[key][index] = posterior_params[key]

            else:
                for key in cols:
                    output_dict[key] = np.append(output_dict[key], fit_info[key])
                for key in lecs_cols:
                    output_dict[key] = np.append(output_dict[key], posterior_params[key])

            output_cols = np.concatenate((cols,lecs_cols), axis=0)
            df = pd.DataFrame.from_dict(output_dict)
            df = df[output_cols]
            df.sort_values('name')
            df.to_csv(filepath)

        # Create new file if it doesn't exist
        else:
            output_dict = {}
            for key in cols:
                output_dict[key] = [fit_info[key]]
            for key in lecs_cols:
                output_dict[key] = [posterior_params[key]]

            output_cols = np.concatenate((cols,lecs_cols), axis=0)
            df = pd.DataFrame.from_dict(output_dict)
            df = df[output_cols] # rearrange in logical order
            df.to_csv(filepath)

        print("Done.")
        return None

    def save_plots(self, figs=None, output_filename=None):

        if not os.path.exists(os.path.normpath(self.project_path+'/tmp/')):
            os.makedirs(os.path.normpath(self.project_path+'/tmp/'))

        if output_filename is None:
            output_file = os.path.normpath(self.project_path+'/tmp/temp.pdf')
        else:
            output_file = os.path.normpath(self.project_path+'/figs/generated/'+output_filename+'.pdf')

        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))

        output_pdf = PdfPages(output_file)
        try:
            for fig in figs:
                output_pdf.savefig(fig)
        except TypeError: # save figs directly if figs is just a figure, not an array
            output_pdf.savefig(figs)

        output_pdf.close()
        print("Done.")

        return None

    def save_prior(self, prior, fit_type, F2, include_FV, include_alphaS, include_logSq):
        print("Saving...")

        if not os.path.exists(self.project_path + '/priors/'):
            os.makedirs(self.project_path + '/priors/')
        filepath = os.path.normpath(self.project_path + '/priors/'+fit_type+'.csv')

        # get fit info
        cols = np.array(['name', 'L_4', 'L_5', 'A_k', 'A_p', 'A_a', 'A_loga'])

        name = F2
        if include_FV:
            name = name + '_FV'
        if include_alphaS:
            name = name + '_alphaS'
        if include_logSq:
            name = name + '_logSq'


        # fit_info keys not in cols -> create key in fit_info
        for key in cols:
            if key not in prior.keys():
                prior[key] = np.nan

        # Add result to file if file exists
        if os.path.isfile(filepath):
            df_prior = pd.read_csv(filepath, index_col=0).to_dict()

            output_dict = {}
            for key in df_prior.keys():
                output_dict[key] = np.array(list(df_prior[key].values()), dtype="object")

            if name in output_dict['name']:
                index = (np.argwhere(output_dict['name'] == name)).item()
                for key in prior.keys():
                    output_dict[key][index] = prior[key]
                output_dict['name'][index] = name
            else:
                for key in output_dict.keys():
                    if key == 'name':
                        output_dict['name'] = np.append(output_dict['name'], name)
                    else:
                        output_dict[key] = np.append(output_dict[key], prior[key])

            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols]
            df.sort_values('name')
            df.to_csv(filepath)

        # Create new file if it doesn't exist
        else:
            output_dict = {key : [prior[key]] for key in prior.keys()}
            output_dict['name'] = name
            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols] # rearrange in logical order
            df.to_csv(filepath)

        print("Done.")
        return None
