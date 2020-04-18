#!/usr/bin/env python3
# python libraries
import os, sys
# Peter Lepage's libraries
import gvar as gv
import lsqfit

# FK / Fpi libraries
sys.path.append('utils')
import io_utils
import chipt
import analysis
import plotting

''' Write description
'''

def main():
    print("python     version:", sys.version)
    import input_params as ip
    # Load input params
    switches   = ip.switches
    priors     = ip.priors
    phys_point = ip.phys_point
    check_fit  = ip.check_fit

    # if check_fit: - add support

    # load data
    gv_data = io_utils.format_h5_data('data/FK_Fpi_data.h5',switches)
    #print('y',gv_data['y'])

    models = sys_models(switches)
    for model in models:
        print('===============================================================')
        print(model)
        print('---------------------------------------------------------------')
        model_list, FF, fv = gather_model_elements(model)
        fit_model  = chipt.FitModel(model_list, _fv=fv, _FF=FF)
        fitEnv     = FitEnv(gv_data, fit_model, switches)
        fit_result = fitEnv.fit_data(priors)
        #print(fit_result.format(maxline=True))
        report_phys_point(fit_result, phys_point, model_list, FF)
        #gv.dump(fit_result, 'pickled_fits/'+model+'.p', add_dependencies=True)

'''
    This is the main class that runs a given fit specified by a model
'''
class FitEnv:
    # xyp_dict is a dictionary with keys 'x', 'y', 'p'
    # 'y' value is a dict{ensemble : yval},
    # 'x' value is a dict{ensemble : {<all the x variables which are not priors>}
    # 'p' value is a dict{(ensemble, pKey) : aPal} for all the
    #       ensemble-specific priors like meson masses as well as teh LECs
    def __init__(self, xyp_dict, model, switches):
        self.switches   = switches
        self.ensembles  = switches['ensembles_fit']
        self.pruned_x   = {ens : { k : v for k, v in xyp_dict['x'][ens].items()}
                                for ens in self.ensembles}
        self.pruned_y   = {ens : xyp_dict['y'][ens] for ens in self.ensembles}
        required_params = model.get_required_parameters()
        self.pruned_p   = {(ens, k) : v for (ens, k), v in xyp_dict['p'].items()
                                if k in required_params and ens in self.ensembles}
        self.model      = model

    def fit_function(self, x, p):
        r = dict()
        for ens in x.keys():
            p_ens = dict()
            for k, v in p.items():
                if type(k) == tuple and k[0] == ens:
                    p_ens[k[1]] = v # the x-params which are priors
                else:
                    p_ens[k] = v    # the LECs of the fit
            r[ens] = self.model(x[ens], p_ens)
        return r

    def fit_data(self, lec_priors):
        required_params = self.model.get_required_parameters()
        # add the LEC priors to our list of priors for the fit
        self.pruned_p.update({ k:v for k,v in lec_priors.items() if k in required_params})
        x = self.pruned_x
        y = self.pruned_y
        p = self.pruned_p
        if self.switches['scipy']:
            fitter='scipy_least_squares'
        else:
            fitter='gsl_multifit'
        #print('p',p)
        return lsqfit.nonlinear_fit(data=(x,y), prior=p, fcn=self.fit_function, fitter=fitter, debug=True)


def check_for_duplicates(list_of_elems):
    ''' Check if given list contains any duplicates '''
    if len(list_of_elems) == len(set(list_of_elems)):
        return False
    else:
        return True

def sys_models(switches):
    def check_model(sys_val,models,nnlo=False,nnnlo=False):
        for model in models:
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
        #return models
    models = switches['ansatz']['models'].copy()
    #print(len(models),'models')
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

    model_elements = [eft+'_nlo']
    if eft == 'taylor':
        if order in ['nnlo', 'nnnlo']:
            model_elements += ['nnlo_ct']
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

def report_phys_point(fit_result, phys_point, model_list, FF):
    phys_data = {'x':dict(phys_point['x']), 'y':dict(phys_point['y']), 'p':dict(phys_point['p'])}
    fit_model = chipt.FitModel(model_list, _fv=False, _FF=FF)
    fitEnv    = FitEnv(phys_data, fit_model, {'ensembles_fit':['phys']})
    for k in fit_result.p:
        if isinstance(k,str):
            phys_data['p'][k] = fit_result.p[k]
    result    = fitEnv.fit_function(phys_data['x'], phys_data['p'])
    print('FK/Fpi = %s' %result['phys'])


if __name__ == "__main__":
    main()
