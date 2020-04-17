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

    for model in ['xpt_nlo', 'ma_nlo']:
        print(model)
        fit_model  = chipt.FitModel([model], _fv=False, _FF='PP')
        fitEnv = FitEnv(gv_data, fit_model, switches)
        fit_result = fitEnv.fit_data(priors)
        print(fit_result.format(maxline=True))

        fit_model  = chipt.FitModel([model], _fv=True, _FF='PP')
        fitEnv = FitEnv(gv_data, fit_model, switches)
        fit_result = fitEnv.fit_data(priors)
        print(fit_result.format(maxline=True))

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



if __name__ == "__main__":
    main()
