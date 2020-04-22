#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import os, sys, shutil
import yaml

import gvar as gv


class BayesModelAvg:
    ''' Under Bayes Theorem, models which are used to fit the same dataset
        can easily be averaged together with the exp(logGBF) factor providing
        a weight, such that the relative weigth of any given model with respect
        to the most likely model (of those selected) is given by

        w_i = exp( logGBF_i - logGBF_max )

        The final uncertainty under such model averaging is given simply by a
        combination of the weighted variance from each fit, plus the variance
        over the models:

        mu  = sum_i w_i mu_i
        Var = sum_i w_i * var_i   + [ sum_j mu_j**2 * w_j - mu**2 ]

        This class performs this weighted average, and creates histograms of
        the final answer including splitting the historgram in various ways to
        understand the dominant contributions to the model average uncertainty
    '''
    def __init__(self, results):
        ''' results is a dictionary where the key is the model name and the
            value is an lsqfit object
        '''
        self.results = results
        # create a list of the keys so we gaurantee the same order
        self.r_list  = [a_res for a_res in self.results]
        self.weights = self.get_weights()

    def get_weights(self):
        weights = []
        for a_res in self.r_list:
            weights.append(np.exp(self.results[a_res].logGBF))
        weights = np.array(weights)
        weights = weights / np.max(weights)
        weights = weights / weights.sum()
        return weights

    def print_weighted_models(self):
        i_weights = np.argsort(self.weights)[::-1]
        print(r"model & $\chi^2/{\rm dof}$ & $Q$& logGBF& weight& $F_K/F_\pi$\\")
        print(r'\hline')
        for a_i, a_model in enumerate(np.array(self.r_list)[i_weights]):
            chi2   = self.results[a_model].chi2
            dof    = self.results[a_model].dof
            Q      = self.results[a_model].Q
            logGBF = self.results[a_model].logGBF
            w_i    = self.weights[i_weights][a_i]
            phys   = self.results[a_model].phys
            print(r'%33s &  %.3f&  %.3f&  %.3f&  %.3f&  %s\\' %(a_model, chi2/dof, Q, logGBF, w_i, phys))



def check_for_duplicates(list_of_elems):
    ''' Check if given list contains any duplicates '''
    if len(list_of_elems) == len(set(list_of_elems)):
        return False
    else:
        return True


def sys_models(switches):
    def check_model(sys_val,models,nnlo=False,nnnlo=False):
        for model in models:
            if 'taylor' not in model:
                new_model = model+sys_val
                if nnlo:
                    if (sys_val not in model) and (new_model not in models):
                        if nnnlo:
                            if 'nnlo' in model and 'nnnlo' not in model:
                                models.append(new_model)
                        else:
                            if 'nnlo' in model:
                                models.append(new_model)
                else:
                    if (sys_val not in model) and (new_model not in models):
                        models.append(new_model)

    models = switches['ansatz']['models'].copy()
    if switches['sys']['FV']:
        check_model('_FV',models)
    if switches['sys']['alphaS']:
        check_model('_alphaS',models,nnlo=True)
    if switches['sys']['nnlo_ct']:
        check_model('_ct',models,nnlo=True)
    if switches['sys']['logSq']:
        check_model('_logSq',models,nnlo=True)
    if switches['sys']['a4']:
        check_model('_a4',models,nnlo=True,nnnlo=True)
    if switches['sys']['ratio']:
        for model in models:
            model_ratio = model.replace('xpt','xpt-ratio').replace('ma','ma-ratio')
            if '-ratio' not in model and model_ratio not in models:
                models.append(model_ratio)
    models_FPK = []
    for model in models:
        if switches['sys']['Lam_chi']:
            for FPK in switches['scales']:
                models_FPK.append(model+'_'+FPK)
        else:
            models_FPK.append(model+'_'+switches['scale'])
    if switches['debug_models']:
        for model in models_FPK:
            print(model)
    print(len(models_FPK),'models')
    print('Duplicate models?',check_for_duplicates(models_FPK))
    return models_FPK

def gather_model_elements(model):
    eft    = model.split('_')[0]
    order  = model.split('_')[1]
    FF     = model.split('_')[-1]
    fv     = 'FV'     in model
    alphaS = 'alphaS' in model
    ct     = 'ct'     in model
    a4     = 'a4'     in model

    if FF not in ['PP','PK','KK']:
        sys.exit('unrecognized FF choice [PP, PK, KK]: '+FF)

    if 'ratio' in eft:
        model_elements = [eft.replace('-','_')+'_nlo']
    else:
        model_elements = [eft+'_nlo']
    if eft == 'taylor':
        if order in ['nnlo', 'nnnlo']:
            model_elements += ['nnlo_ct']
            if alphaS:
                model_elements += ['nnlo_alphaS']
        if order in ['nnnlo']:
            model_elements += ['nnnlo_ct']
    else:
        if order in ['nnlo','nnnlo']:
            if alphaS:
                model_elements += ['nnlo_alphaS']
            if ct:
                model_elements += ['nnlo_ct']
            else:
                model_elements += ['nnlo_ct','xpt_nnlo_logSq','xpt_nnlo_log','xpt_nnlo_FF_'+FF]
                if 'ratio' in eft:
                    model_elements += ['xpt_nnlo_ratio']


            if a4 and order == 'nnlo':
                model_elements += ['nnnlo_a4']
            if order == 'nnnlo':
                model_elements += ['nnnlo_ct']

    return model_elements, FF, fv

def prior_width_scan(model, fitEnv, fit_model, priors, switches):
    new_priors = dict(priors)
    p2_vals = [2, 5, 10]
    p3_vals = [2, 5, 10]

    logGBF_array = []
    if os.path.exists('data/saved_prior_search.yaml'):
        with open('data/saved_prior_search.yaml','r') as fin:
            prior_grid = yaml.safe_load(fin.read())
            if model in prior_grid:
                for k1 in prior_grid[model]:
                    for k2 in prior_grid[model][k1]:
                        for k3 in prior_grid[model][k1][k2]:
                            for k4 in prior_grid[model][k1][k2][k3]:
                                logGBF_array.append([k1, k2, k3, k4, float(prior_grid[model][k1][k2][k3][k4])])
        shutil.copyfile('data/saved_prior_search.yaml','data/saved_prior_search.yaml.bak')
    else:
        prior_grid = dict()

    if model not in prior_grid:
        prior_grid[model] = dict()
    nnlo_x_lst = ['p_4','k_4']
    nnlo_a_lst = ['s_4','saS_4']
    n3lo_x_lst = ['kp_6','k_6','p_6']
    n3lo_a_lst = ['s_6','sk_6','sp_6']

    n_p2 = len(p2_vals)
    n_p3 = len(p3_vals)
    i_search = 1
    i_tot = n_p2**2 * n_p3**2
    for i2x,p2x in enumerate(p2_vals):
        if p2x not in prior_grid[model]:
            prior_grid[model][p2x] = dict()
        for k2x in nnlo_x_lst:
            new_priors[k2x] = gv.gvar(0,p2x)
        for i2a,p2a in enumerate(p2_vals):
            if p2a not in prior_grid[model][p2x]:
                prior_grid[model][p2x][p2a] = dict()
            for k2a in nnlo_a_lst:
                new_priors[k2a] = gv.gvar(0,p2a)
            for i3x,p3x in enumerate(p3_vals):
                if p3x not in prior_grid[model][p2x][p2a]:
                    prior_grid[model][p2x][p2a][p3x] = dict()
                for k3x in n3lo_x_lst:
                    new_priors[k3x] = gv.gvar(0,p3x)
                for i3a,p3a in enumerate(p3_vals):
                    for k3a in n3lo_a_lst:
                        new_priors[k3a] = gv.gvar(0,p3a)
                    if p3a not in prior_grid[model][p2x][p2a][p3x]:
                        tmp_result = fitEnv.fit_data(new_priors)
                        lgbf = float(tmp_result.logGBF)
                        prior_grid[model][p2x][p2a][p3x][p3a] = lgbf
                        logGBF_array.append([ p2x, p2a, p3x, p3a, lgbf] )
                    else:
                        lgbf = float(prior_grid[model][p2x][p2a][p3x][p3a])
                        logGBF_array.append([ p2x, p2a, p3x, p3a, lgbf] )
                    if switches['prior_verbose']:
                        print('%4d / %4d: nnlo_x = %.2f,   nnlo_a = %.2f,   n3lo_x = %.2f,   n3lo_a = %.2f:   logGBF=%f' \
                            %(i_search, i_tot, p2x, p2a, p3x, p3a, lgbf))
                    else:
                        sys.stdout.write("%4d / %4d\r" %(i_search, i_tot))
                        sys.stdout.flush()
                    i_search += 1
    logGBF_array = np.array(logGBF_array)
    i_l_max = np.argmax(logGBF_array[:,4])
    vals = ''
    vals = vals.join("& %5s   " %str(k) for k in logGBF_array[i_l_max][:-1])
    vals = vals + "& %f" %logGBF_array[i_l_max][-1]
    print('%33s' %model, vals)
    prior_file = open('data/saved_prior_search.yaml', 'w')
    yaml.dump(prior_grid, prior_file)
    prior_file.close()
