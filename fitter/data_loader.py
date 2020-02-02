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
        self.models = None

    def _pickle_fit_parameters(self, fit_info, output_name):
        model = fit_info['name']
        filename = self.project_path+'/pickles/'+ output_name +'/'+ model +'.p'
        if not os.path.exists(self.project_path +'/pickles/'+ output_name):
            os.makedirs(self.project_path +'/pickles/'+ output_name)

        output = {
            'FK/Fpi' : fit_info['FK/Fpi'],
            'delta_su2' : fit_info['delta_su2'],
            'logGBF' : gv.gvar(fit_info['logGBF']),
            'chi2/df' : gv.gvar(fit_info['chi2/df']),
            'Q' : gv.gvar(fit_info['Q']),
        }

        for key in fit_info['prior'].keys():
            output['prior:'+key] = fit_info['prior'][key]

            output[key] = [fit_info['prior'][key], fit_info['posterior'][key]]

        for key in fit_info['posterior'].keys():
            output['posterior:'+key] = fit_info['posterior'][key]

        for key in fit_info['phys_point'].keys():
            print(key)
            # gvar can't handle integers -- entries not in correlation matrix
            if key not in ['a/w0', 'a', 'L', 'alpha_s', 'a2DI']:
                output['phys_point:'+key] = fit_info['phys_point'][key]

        gv.dump(output, filename)
        return None

    def _unpickle_fit_parameters(self, model, output_name):
        filepath = self.project_path+'/pickles/'+ output_name +'/'+ model +'.p'
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

    def get_fit_info(self, output_name):
        if os.path.exists(self.project_path +'/pickles/'+ output_name):
            fit_info = {}
            for model in self.get_models(output_name):
                fit_info_model = self._unpickle_fit_parameters(model=model, output_name=output_name)
                fit_info[model] = {}
                fit_info[model]['name'] = model
                fit_info[model]['FK/Fpi'] = fit_info_model['FK/Fpi']
                fit_info[model]['delta_su2'] = fit_info_model['delta_su2']
                fit_info[model]['logGBF'] = fit_info_model['logGBF'].mean
                fit_info[model]['chi2/df'] = fit_info_model['chi2/df'].mean
                fit_info[model]['Q'] = fit_info_model['Q'].mean
                fit_info[model]['prior'] = {}
                fit_info[model]['posterior'] = {}
                fit_info[model]['phys_point'] = {}

                for key in fit_info_model.keys():
                    if key.startswith('prior'):
                        fit_info[model]['prior'][key.split(':')[-1]] = fit_info_model[key]
                    elif key.startswith('posterior'):
                        fit_info[model]['posterior'][key.split(':')[-1]] = fit_info_model[key]
                    elif key.startswith('phys_point'):
                        fit_info[model]['phys_point'][key.split(':')[-1]] = fit_info_model[key]

                for key in ['a/w0', 'a', 'L', 'alpha_s', 'a2DI']:
                    fit_info[model]['phys_point'][key] = self.get_phys_point_data(key)


            return fit_info

        # Alternatively, read info from csv file.
        # Used for getting fit results from other collabs
        else:
            filepath = os.path.normpath(self.project_path + '/results/'+output_name+'.csv')

            df_fit = pd.read_csv(filepath, header=0)
            cols = df_fit.columns.values

            models = df_fit['name'].values
            output_dict = OrderedDict()
            for model in sorted(models):
                index = np.argwhere(df_fit['name'].values == model)

                output_dict[model] = {}
                for key in cols:
                    if key in df_fit.keys():
                        output_dict[model][key] = (df_fit[key].values[index]).item()

            return output_dict

        if not os.path.isfile(filepath):
            return None


    def get_phys_point_data(self, parameter=None):
        phys_point_data = {
            'a/w0' : 0,
            'a' : 0,
            'L' : np.infty,
            'alpha_s' : 0, # Not sure, but not important since it comes with a^2

            'mpi' : gv.gvar('134.8(3)'), # '138.05638(37)'
            'mk' : gv.gvar('494.2(3)'), # '495.6479(92)'
            'mss' : gv.gvar('688.5(2.2)'), # Taken from arxiv/1303.1670

            'a2DI' : 0,
            'Fpi' : gv.gvar(130.2/np.sqrt(2), 0.8/np.sqrt(2)), # PDG
            'FK' : gv.gvar(155.5/np.sqrt(2), 0.7/np.sqrt(2)), # PDG
            'w0' : gv.gvar('0.175(10)'),

            'FK/Fpi_pm' : gv.gvar('1.1932(19)'), # FLAG, SU(2) isospin corrected value (arxiv/1902.08191, eqn 80)
        }


        # Or get mss, mrs with Gell-Mann-Oakes-Renner relations: arxiv/0505265 (3.45)
        mpi = phys_point_data['mpi']
        mk = phys_point_data['mk']
        phys_point_data['mss'] = np.sqrt(2 *(mk)**2 - (mpi)**2) *1.000000001 # prevents division by 0

        # ma pion
        phys_point_data['mju'] = phys_point_data['mpi']

        # ma kaon
        phys_point_data['mru'] = phys_point_data['mk']
        phys_point_data['mjs'] = phys_point_data['mk']

        # ma eta_s
        phys_point_data['mrs'] = phys_point_data['mss']

        FK = phys_point_data['FK']
        Fpi = phys_point_data['Fpi']
        phys_point_data['lam2_chi_kpi'] = (4*np.pi)**2 *FK *Fpi
        phys_point_data['lam2_chi_pipi'] = (4*np.pi)**2 *Fpi *Fpi
        phys_point_data['lam2_chi_kk'] = (4*np.pi)**2 *FK *FK
        phys_point_data['lam2_chi_00'] = (4*np.pi)**2 *(gv.gvar('131.5(0.1)') / np.sqrt(2))**2
        
        if parameter is None:
            return phys_point_data
        else:
            return phys_point_data[parameter]

    def get_prior(self, fit_type, order, F2,
                  include_log, include_log2, include_sunset,
                  include_alpha_s, include_latt_n3lo, include_FV, use_bijnens_central_value):
        filepath = os.path.normpath(self.project_path + '/priors/'+fit_type+'.csv')

        print(filepath)

        if not os.path.isfile(filepath):
            return None

        name = fit_type +'_'+ F2+'_'+order
        if include_log:
            name = name + '_log'
        if include_log2:
            name = name + '_logSq'
        if include_sunset:
            name = name + '_sunset'
        if include_alpha_s:
            name = name + '_alphaS'
        if include_latt_n3lo:
            name = name + '_a4'
        if include_FV:
            name = name + '_FV'
        if use_bijnens_central_value:
            name = name + '_bijnens'

        df_prior = pd.read_csv(filepath, header=0)
        if name not in df_prior['name'].values:
            return None

        index = np.argwhere(df_prior['name'].values == name)
        prior = gv.BufferDict()
        for key in ['L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8',
                    'L_4', 'L_5', 'A_k', 'A_p', 'A_a', 'A_loga', 'A_aa']:
            value = np.asscalar(df_prior[key].values[index])
            if value is not np.nan:
                prior[key] = gv.gvar(value)

        return prior

    def get_models(self, output_name):
        if self.models is None:
            models = []
            for file in os.listdir(self.project_path +'/pickles/'+ output_name):
                if(file.endswith('.p')):
                    models.append(file.split('.')[0])
            self.models = models
        return self.models


    # pickle correlated prior/posterior,
    # save rest to csv file
    def save_fit_info(self, fit_info, output_name=None, save_pickles=True):
        print("Saving...")

        if output_name is None:
            output_name = 'fit_results'
        if save_pickles:
            self._pickle_fit_parameters(fit_info, output_name)

        if not os.path.exists(self.project_path + '/results/'):
            os.makedirs(self.project_path + '/results/')
        filepath = os.path.normpath(self.project_path + '/results/'+output_name+'.csv')

        # get fit info
        cols = np.array(['name', 'FK/Fpi', 'delta_su2', 'logGBF', 'chi2/df', 'Q', 'vol'])
        lecs_cols = ['L_4', 'L_5', # nlo terms
                     'L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9', #semi-nnlo terms
                     'A_a', 'A_k', 'A_p', 'A_loga', 'A_aa'] # nnlo terms
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

    def save_prior(self, prior, name):
        print("Saving...")

        fit_type = name.split('_')[0]

        if not os.path.exists(self.project_path + '/priors/'):
            os.makedirs(self.project_path + '/priors/')
        filepath = os.path.normpath(self.project_path + '/priors/'+fit_type+'.csv')

        # get fit info
        cols = np.array(['name', 'L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8',
                         'L_4', 'L_5', 'A_k', 'A_p', 'A_a', 'A_loga', 'A_aa'])

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
