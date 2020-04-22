#!/usr/bin/env python3
# python libraries
import os, sys, shutil, copy
import matplotlib.pyplot as plt
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

''' useful for changing matplotlib.plt.show behavior and accessing results
    from main() if run interactively in iPython
'''
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

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
        models = analysis.sys_models(switches)
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
            model_list, FF, fv = analysis.gather_model_elements(model)
            debug_fit_function(check_fit, model_list, FF, fv)
        sys.exit()

    # load data
    gv_data = io_utils.format_h5_data('data/FK_Fpi_data.h5',switches)

    models = analysis.sys_models(switches)
    if switches['prior_search']:
        print('Performing Prior Width Scan')
        print(r'%33s &  nnlo_x &  nnlo_a &  n3lo_x &  n3lo_a & logGBF_max' %'model')
        print('======================================================================================')
        for model in models:
            model_list, FF, fv = analysis.gather_model_elements(model)
            fit_model  = chipt.FitModel(model_list, _fv=fv, _FF=FF)
            fitEnv     = FitEnv(gv_data, fit_model, switches)
            analysis.prior_width_scan(model, fitEnv, fit_model, priors, switches)
    else:
        fit_results = dict()
        plt.ion()
        for model in models:
            print('===============================================================')
            print(model)

            do_fit = False
            if switches['save_fits'] or switches['debug_save_fit']:
                if os.path.exists('pickled_fits/'+model+'.p'):
                    print('reading pickled_fits/'+model+'.p')
                    fit_result = gv.load('pickled_fits/'+model+'.p')
                else:
                    do_fit = True
            else:
                do_fit = True
            if do_fit or switches['debug_save_fit']:
                model_list, FF, fv = analysis.gather_model_elements(model)
                fit_model  = chipt.FitModel(model_list, _fv=fv, _FF=FF)
                fitEnv     = FitEnv(gv_data, fit_model, switches)
                tmp_result = fitEnv.fit_data(priors)
                tmp_result.phys = report_phys_point(tmp_result, phys_point, model_list, FF)
                if switches['print_fit']:
                    print(tmp_result.format(maxline=True))
                if not os.path.exists('pickled_fits/'+model+'.p') and switches['save_fits']:
                    gv.dump(tmp_result, 'pickled_fits/'+model+'.p', add_dependencies=True)
                    fit_result = gv.load('pickled_fits/'+model+'.p')
                if switches['debug_save_fit']:
                    print('live fit')
                    print('dF/dy =', tmp_result.phys.partialsdev(tmp_result.y))
                    print('dF/dprior =', tmp_result.phys.partialsdev(tmp_result.prior))
                    print('pickled fit')
                    print('dF/dy =', fit_result.phys.partialsdev(fit_result.y))
                    print('dF/dprior =', fit_result.phys.partialsdev(fit_result.prior))
                if do_fit:
                    fit_result = tmp_result

            fit_results[model] = fit_result
            if switches['make_plots']:
                plots = plotting.ExtrapolationPlots(model, model_list, fitEnv, fit_result, switches)
                if 'alphaS' not in model and 'ma' not in model:
                    plots.plot_vs_eps_asq(phys_point)
                if 'ma' not in model:
                    plots.plot_vs_eps_pi(phys_point)

        model_avg = analysis.BayesModelAvg(fit_results)
        model_avg.print_weighted_models()

        plt.ioff()
        if run_from_ipython():
            plt.show(block=False)
            return fit_results
        else:
            plt.show()


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
        self.x          = xyp_dict['x']
        self.pruned_x   = {ens : { k : v for k, v in xyp_dict['x'][ens].items()}
                                for ens in self.ensembles}
        self.y          = xyp_dict['y']
        self.pruned_y   = {ens : xyp_dict['y'][ens] for ens in self.ensembles}
        required_params = model.get_required_parameters()
        self.p          = xyp_dict['p']
        self.pruned_p   = {(ens, k) : v for (ens, k), v in xyp_dict['p'].items()
                                if k in required_params and ens in self.ensembles}
        self.model      = model

    # create a callable function that acts on a single x and p (not all ensembles)
    @classmethod
    def _fit_function(cls, a_model, x, p):
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

def report_phys_point(fit_result, phys_point_params, model_list, FF):
    phys_data = copy.deepcopy(phys_point_params)
    if 'ma' in model_list[0]:
        model_list[0] = model_list[0].replace('ma','xpt')
    fit_model = chipt.FitModel(model_list, _fv=False, _FF=FF)
    for k in fit_result.p:
        if isinstance(k,str):
            phys_data['p'][k] = fit_result.p[k]
    result    = FitEnv._fit_function(fit_model, phys_data['x'], phys_data['p'])
    print('  chi2/dof [dof] = %.2f [%d]   Q=%.3f   logGBF = %.3f' \
        %(fit_result.chi2/fit_result.dof, fit_result.dof, fit_result.Q, fit_result.logGBF))
    print('  FK/Fpi = %s' %result)
    return result

def debug_fit_function(check_fit, model_list, FF, fv):
    x = check_fit['x']
    p = check_fit['p']
    fit_model = chipt.FitModel(model_list, _fv=False, _FF=FF)
    cP        = chipt.ConvenienceDict(fit_model, x, p)
    result    = 0.
    result_FV = 0.
    if fv:
        fit_model_fv = chipt.FitModel(model_list, _fv=True, _FF=FF)
        cP_FV        = chipt.ConvenienceDict(fit_model_fv, x, p)
        for term in model_list:
            if '_nlo' in term:
                t_FV = getattr(chipt.FitModel, term)(fit_model_fv, x, p, cP_FV)
                t    = getattr(chipt.FitModel, term)(fit_model, x, p, cP)
                result    += t
                result_FV += t_FV
                print('%16s   ' %(term+'_FV'), t_FV)
                print('%16s   ' %(term), t)
            else:
                t = getattr(chipt.FitModel, term)(fit_model, x, p, cP)
                result_FV += t
                result    += t
                print('%16s   ' %(term), t)
    else:
        for term in model_list:
            t = getattr(chipt.FitModel, term)(fit_model, x, p, cP)
            result += t
            print('%16s   ' %(term), t)
    print('---------------------------------------')
    if fv:
        print('%16s   %f' %('total_FV', result_FV))
    print('%16s   %f' %('total', result))

if __name__ == "__main__":
    main()
