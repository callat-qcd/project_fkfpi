import pandas as pd
import numpy as np
import gvar as gv
import sys
import re
import time
import os
import h5py
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict

class data_loader(object):

    def __init__(self):
        self.project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
        self.file_h5 = os.path.normpath(self.project_path+'/FK_Fpi_data.h5')

    def _pickle_fit_parameters(self, fit_info, collection_name):
        model = fit_info['name']
        filename = self.project_path +'/results/'+ collection_name +'/pickles/'+ model +'.p'
        if not os.path.exists(self.project_path +'/results/'+ collection_name +'/pickles/'):
            os.makedirs(self.project_path +'/results/'+ collection_name +'/pickles/')

        output = {
            'FK/Fpi' : fit_info['FK/Fpi'],
            'delta_su2' : fit_info['delta_su2'],
            'logGBF' : gv.gvar(fit_info['logGBF']),
            'chi2/df' : gv.gvar(fit_info['chi2/df']),
            'Q' : gv.gvar(fit_info['Q']),
        }

        for key in fit_info['prior'].keys():
            output['prior:'+key] = fit_info['prior'][key]
            #output[key] = [fit_info['prior'][key], fit_info['posterior'][key]]

        for key in fit_info['posterior'].keys():
            output['posterior:'+key] = fit_info['posterior'][key]

        for key in fit_info['phys_point'].keys():
            # gvar can't handle integers -- entries not in correlation matrix
            output['phys_point:'+key] = fit_info['phys_point'][key]

        for key in fit_info['error_budget']:
            output['error_budget:'+key] = gv.gvar(fit_info['error_budget'][key])

        gv.dump(output, filename)
        return None

    def _unpickle_fit_parameters(self, model, collection_name):
        filepath = self.project_path +'/results/'+ collection_name +'/pickles/'+ model +'.p'
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

    def get_fit_results(self, collection_name):
        if os.path.exists(self.project_path +'/results/'+ collection_name +'/pickles/'):
            fit_info = {}

            models = []
            for file in os.listdir(self.project_path +'/results/'+ collection_name +'/pickles/'):
                if(file.endswith('.p')):
                    models.append(file.split('.')[0])

            for model in models:
                fit_info_model = self._unpickle_fit_parameters(model=model, collection_name=collection_name)
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
                fit_info[model]['error_budget'] = {}

                for key in fit_info_model.keys():
                    if key.startswith('prior'):
                        fit_info[model]['prior'][key.split(':')[-1]] = fit_info_model[key]
                    elif key.startswith('posterior'):
                        fit_info[model]['posterior'][key.split(':')[-1]] = fit_info_model[key]
                    elif key.startswith('phys_point'):
                        fit_info[model]['phys_point'][key.split(':')[-1]] = fit_info_model[key]
                    elif key.startswith('error_budget'):
                        fit_info[model]['error_budget'][key.split(':')[-1]] = fit_info_model[key].mean

            return fit_info

        # Alternatively, read info from csv file.
        # Used for getting fit results from other collabs
        elif os.path.exists(self.project_path + '/results/' +collection_name+ '/results.csv'):
            filepath = os.path.normpath(self.project_path + '/results/' +collection_name+ '/results.csv')
            if not os.path.exists(filepath):
                return None

            df_fit = pd.read_csv(filepath, header=0)
            cols = df_fit.columns.values

            models = df_fit['name'].values
            output_dict = OrderedDict()
            for model in sorted(models):
                index = np.argwhere(df_fit['name'].values == model)

                output_dict[model] = {}
                output_dict[model]['posterior'] = {}
                output_dict[model]['error_budget'] = {}
                for key in cols:
                    if key.startswith('po:'):
                        output_dict[model]['posterior'][key.split(':')[-1]] = (df_fit[key].values[index]).item()
                    elif key.startswith('eb:'):
                        output_dict[model]['error_budget'][key.split(':')[-1]] = (df_fit[key].values[index]).item()
                    elif key.startswith('Unnamed:'):
                        pass
                    elif key in df_fit.keys():
                        output_dict[model][key] = (df_fit[key].values[index]).item()

            return output_dict

        else:
            return None

    def get_model_info_from_name(self, name):
        model_info = {}
        model_info['name'] = name
        model_info['fit_type'] = name.split('_')[0] # eg, 'ma-ratio'
        model_info['F2'] = name.split('_')[1] # eg, 'FKFPi'
        model_info['order'] = name.split('_')[2] # eg, 'n2lo'

        model_info['include_FV'] = False
        model_info['exclude'] = []
        model_info['use_bijnens_central_value'] = False

        model_info['include_alpha_s'] = False
        model_info['include_log'] = False
        model_info['include_log2'] = False
        model_info['include_sunset'] = False

        model_info['latt_ct'] = 'n2lo'


        if '_FV' in name:
            model_info['include_FV'] = True
        if '_alphaS' in name:
            model_info['include_alpha_s'] = True
        if '_log_' in name or name.endswith('_log'): # (consider '_logSq')
            model_info['include_log'] = True
        if '_logSq' in name:
            model_info['include_log2'] = True
        if '_sunset' in name:
            model_info['include_sunset'] = True
        if '_bijnens' in name:
            model_info['use_bijnens_central_value'] = True

        if '_a6' in name:
            model_info['latt_ct'] = 'n4lo'
        elif '_a4' in name:
            model_info['latt_ct'] = 'n3lo'

        if (model_info['include_log'] == True
                and model_info['include_log2'] == True
                and model_info['include_sunset'] == True
                and model_info['use_bijnens_central_value'] == True):
            model_info['semi-n2lo_corrections'] = 'n2lo-full_bijnens'

        elif (model_info['include_log'] == True
                and model_info['include_log2'] == True
                and model_info['include_sunset'] == True
                and model_info['use_bijnens_central_value'] == False):
            model_info['semi-n2lo_corrections'] = 'n2lo-full'

        elif (model_info['include_log'] == False
                and model_info['include_log2'] == True
                and model_info['include_sunset'] == True
                and model_info['use_bijnens_central_value'] == False):
            model_info['semi-n2lo_corrections'] = 'n2lo-logSq_ct'

        elif (model_info['include_log'] == False
                and model_info['include_log2'] == False
                and model_info['include_sunset'] == False
                and model_info['use_bijnens_central_value'] == False):
            model_info['semi-n2lo_corrections'] = 'n2lo-ct'

        return model_info

    def get_model_name_from_model_info(self, model_info):
        name = model_info['fit_type'] +'_'+ model_info['F2'] +'_'+ model_info['order']
        if model_info['include_log']:
            name = name + '_log'
        if model_info['include_log2']:
            name = name + '_logSq'
        if model_info['include_sunset']:
            name = name + '_sunset'
        if model_info['include_alpha_s']:
            name = name + '_alphaS'
        if (model_info['latt_ct'] in ['n3lo', 'n4lo']
                or model_info['order'] in ['n3lo']):
            name = name + '_a4'
        if model_info['latt_ct'] in ['n4lo']:
            name = name + '_a6'
        if model_info['include_FV']:
            name = name + '_FV'
        if model_info['use_bijnens_central_value']:
            name = name + '_bijnens'
        return name


    def get_model_names(self, collection_name):
        if not os.path.exists(os.path.normpath(self.project_path+'/results/'+collection_name+'/model_list.txt')):
            print('No such collection!')
            return None

        output = []
        with open(self.project_path+'/results/'+collection_name+'/model_list.txt', 'r') as file:
            output = [current_line.rstrip() for current_line in file.readlines()]

        return output



    def get_phys_point_data(self, parameter=None):
        phys_point_data = {
            'a/w0' : gv.gvar(0),
            'a' : gv.gvar(0),
            'L' : gv.gvar(np.infty),
            'alpha_s' : gv.gvar(0), # Not sure, but not important since it comes with a^2

            'mpi' : gv.gvar('134.8(3)'), # '138.05638(37)'
            'mk' : gv.gvar('494.2(3)'), # '495.6479(92)'
            'mss' : gv.gvar('688.5(2.2)'), # Taken from arxiv/1303.1670

            'a2DI' : gv.gvar(0),
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

    def get_prior(self, collection_name, model_info, **kwargs):
        filepath = os.path.normpath(self.project_path + '/results/'+ collection_name +'/prior.csv')

        if not os.path.isfile(filepath):
            return None

        name = self.get_model_name_from_model_info(self, model_info)

        df_prior = pd.read_csv(filepath, header=0)
        if name not in df_prior['name'].values:
            return None

        index = np.argwhere(df_prior['name'].values == name)
        prior = gv.BufferDict()
        lecs_keys =  ['A_x', 'L_4', 'L_5', # nlo terms
                      'L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9', #semi-n2lo terms
                      'A_a', 'A_k', 'A_p', 'A_loga', 'A_aa', # n2lo terms
                      'A_aa', 'A_ak', 'A_ap', 'A_kk', 'A_kp', 'A_pp'] # n3lo terms
        for key in lecs_keys:
            if key in df_prior.keys():
                value = np.asscalar(df_prior[key].values[index])
                if value is not np.nan:
                    prior[key] = gv.gvar(value)

        return prior

    # pickle correlated prior/posterior,
    # save rest to csv file
    def save_fit_info(self, fit_info, collection_name=None, save_pickles=True):
        print("Saving...")

        if collection_name is None:
            collection_name = 'fit_results'
        if save_pickles:
            self._pickle_fit_parameters(fit_info, collection_name)

        if not os.path.exists(self.project_path +'/results/'+ collection_name):
            os.makedirs(self.project_path +'/results/'+ collection_name)
        filepath = os.path.normpath(self.project_path +'/results/'+ collection_name +'/results.csv')

        # get fit info
        cols = np.array(['name', 'FK/Fpi', 'delta_su2', 'logGBF', 'chi2/df', 'Q', 'vol'])
        lecs_cols =  ['A_x', 'L_4', 'L_5', # nlo terms
                      'L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9', #semi-n2lo terms
                      'A_a', 'A_k', 'A_p', 'A_loga', # n2lo terms
                      'A_aa', 'A_ak', 'A_ap', 'A_kk', 'A_kp', 'A_pp'] # n3lo terms
        eb_cols = ['disc', 'chiral', 'pp_input', 'stat']
        output_cols = []
        for key in cols:
            output_cols.append(key)
        for key in lecs_cols:
            output_cols.append('po:'+key)
        for key in eb_cols:
            output_cols.append('eb:'+key)

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
                    output_dict['po:'+key][index] = posterior_params[key]
                for key in eb_cols:
                    output_dict['eb:'+key][index] = fit_info['error_budget'][key]

            else:
                for key in cols:
                    output_dict[key] = np.append(output_dict[key], fit_info[key])
                for key in lecs_cols:
                    output_dict['po:'+key] = np.append(output_dict['po:'+key], posterior_params[key])
                for key in eb_cols:
                    output_dict['eb:'+key] = np.append(output_dict['eb:'+key], fit_info['error_budget'][key])

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
                output_dict['po:'+key] = [posterior_params[key]]
            for key in eb_cols:
                output_dict['eb:'+key] = fit_info['error_budget'][key]

            df = pd.DataFrame.from_dict(output_dict)
            df = df[output_cols] # rearrange in logical order
            df.to_csv(filepath)

        print("Done.")
        return None

    def save_model_average(self, collection_name, model_avg_str):
        filename = self.project_path +'/results/'+ collection_name + '/README.md'

        if os.path.exists(filename):
            with open(filename, 'r') as file:
                file_content = file.read()
            if '## Model Average' in file_content:
                file_content = re.sub(r'## Model Average\s```(.*\s*)*```',
                                      '## Model Average\n```yaml\n'+model_avg_str+'```',
                                      file_content)
                with open(filename, 'w') as file:
                    file.write(file_content)

        else:
            with open(filename, 'a+') as file:
                file_content = '\n## Model Average'
                file_content += '\n```yaml\n'
                file_content += model_avg_str
                file_content += '```'
                file.write(file_content)

        return None

    def save_model_names(self, collection_name, model_list):
        if not os.path.exists(os.path.normpath(self.project_path+'/results/'+collection_name)):
            os.makedirs(os.path.normpath(self.project_path+'/results/'+collection_name))

        with open(self.project_path+'/results/'+collection_name+'/model_list.txt', 'w') as file:
            file.writelines('%s\n'%model for model in model_list)

        return None

    def save_plots(self, figs=None, output_filename=None):

        if figs is None:
            print('Nothing here!')
            return None

        if not os.path.exists(os.path.normpath(self.project_path+'/tmp/')):
            os.makedirs(os.path.normpath(self.project_path+'/tmp/'))

        if output_filename is None:
            output_file = os.path.normpath(self.project_path+'/tmp/temp.pdf')
        else:
            output_file = os.path.normpath(self.project_path+'/results/'+output_filename+'.pdf')

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

    def save_prior(self, collection_name, prior, name):

        fit_type = name.split('_')[0]

        if not os.path.exists(self.project_path +'/results/'+ collection_name):
            os.makedirs(self.project_path +'/results/'+ collection_name)
        filepath = os.path.normpath(self.project_path +'/results/'+ collection_name +'/prior.csv')

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

        print("Saving prior.")
        return None
