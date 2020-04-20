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
    if switches['check_fit']:
        models = sys_models(switches)
        print('DEBUGGING FIT FUNCTION')
        print('p')
        for k in check_fit['p']:
            print("%7s" %k,check_fit['p'][k])
        print('x')
        for k in check_fit['x']:
            print("%7s" %k,check_fit['x'][k])
        for model in models:
            print('===============================================================')
            print('DEBUGGING Terms in ',model)
            print('---------------------------------------------------------------')
            model_list, FF, fv = gather_model_elements(model)
            debug_fit_function(check_fit, model_list, FF, fv, switches)
        sys.exit()

    # load data
    gv_data = io_utils.format_h5_data('data/FK_Fpi_data.h5',switches)

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

    # create a callable function that acts on a single x and p (not all ensembles)
    def _fit_function(a_model, x, p):
        return a_model(x,p)

    def fit_function(self, x, p):
        a_result = dict()
        for ens in x.keys():
            p_ens = dict()
            for k, v in p.items():
                if type(k) == tuple and k[0] == ens:
                    p_ens[k[1]] = v # the x-params which are priors
                else:
                    p_ens[k] = v    # the LECs of the fit
            model = self.model
            a_result[ens] = FitEnv._fit_function(model, x[ens], p_ens)
        return a_result

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
    phys_data = dict(phys_point)
    fit_model = chipt.FitModel(model_list, _fv=False, _FF=FF)
    for k in fit_result.p:
        if isinstance(k,str):
            phys_data['p'][k] = fit_result.p[k]
    result    = FitEnv._fit_function(fit_model, phys_data['x'], phys_data['p'])
    print('FK/Fpi = %s' %result)

def debug_fit_function(check_fit, model_list, FF, fv, switches):
    switches['ensembles_fit'] = ['']
    x = check_fit['x']
    p = check_fit['p']
    fit_model = chipt.FitModel(model_list, _fv=False, _FF=FF)
    cP        = chipt.ConvenienceDict(fit_model, x, p)
    if fv:
        fit_model_fv = chipt.FitModel(model_list, _fv=True, _FF=FF)
        cP_FV        = chipt.ConvenienceDict(fit_model_fv, x, p)
        for term in model_list:
            if '_nlo' in term:
                print('%16s ' %(term+'_FV'), getattr(chipt.FitModel, term)(fit_model_fv, x, p, cP_FV))
                print('%16s ' %(term), getattr(chipt.FitModel, term)(fit_model, x, p, cP))
            else:
                print('%16s ' %(term), getattr(chipt.FitModel, term)(fit_model, x, p, cP))
    else:
        for term in model_list:
            print('%16s ' %(term), getattr(chipt.FitModel, term)(fit_model, x, p, cP))


if __name__ == "__main__":
    main()
