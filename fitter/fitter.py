import lsqfit
import numpy as np
import gvar as gv
import sys
import os

#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fitter.special_functions as sf

class fitter(object):

    def __init__(self,
            fit_type, order, latt_ct, F2,
            include_FV, exclude,
            include_alpha_s, include_log, include_log2, include_sunset,
            fit_data=None, prior=None, fast_sunset=False, **kwargs
        ):
        self.prior = prior
        self.fit_data = fit_data

        # Model info
        self.model_info = {
            'fit_type' : fit_type,
            'order' : order,
            'latt_ct' : latt_ct,
            'F2' : F2,

            'include_FV' : include_FV,
            'exclude' : exclude,

            'include_alpha_s' : include_alpha_s,
            'include_log' : include_log,
            'include_log2' : include_log2,
            'include_sunset' : include_sunset,

        }

        # attributes of fitter object
        self.counter = {'iters' : 0, 'evals' : 0} # To force empbayes_fit to converge?
        self.fit = None
        self.empbayes_fit = None

    def _make_fitargs(self, z):
        y_data = self._make_y_data()
        prior = self._make_prior()

        # Ideally:
            # Don't bother with more than the hundredth place
            # Don't let z=0 (=> null GBF)
            # Don't bother with negative values (meaningless)
        # But for some reason, these restrictions (other than the last) cause empbayes_fit not to converge
        multiplicity = {}
        for key in z:
            multiplicity[key] = 0
            z[key] = np.abs(z[key])



        # Helps with convergence (minimizer doesn't use extra digits -- bug in lsqfit?)
        sig_fig = lambda x : np.around(x, int(np.floor(-np.log10(x))+3)) # Round to 3 sig figs
        capped = lambda x, x_min, x_max : np.max([np.min([x, x_max]), x_min])

        # Min/max values for unc
        z_min = {}
        z_max = {}

        z_min['chiral_nlo'] = 1e-2
        z_max['chiral_nlo'] = 1e2

        z_min['chiral_n2lo'] = 1e-2
        z_max['chiral_n2lo'] = 1e2

        z_min['chiral_n3lo'] = 1e-2
        z_max['chiral_n3lo'] = 1e3

        z_min['latt_n2lo'] = 1e0
        z_max['latt_n2lo'] = 1e3

        z_min['latt_n3lo'] = 1e0
        z_max['latt_n3lo'] = 1e3

        z_min['latt_n4lo'] = 1e0
        z_max['latt_n4lo'] = 1e4


        for key in prior.keys():
            # polynomial fit
            if key in ['A_x']:
                multiplicity['chiral_nlo'] += 1
                z['chiral_nlo'] = sig_fig(capped(z['chiral_nlo'], z_min['chiral_nlo'], z_max['chiral_nlo']))
                prior[key] = gv.gvar(0, 1) *z['chiral_nlo']

            # chiral_n2lo
            if key in ['A_p', 'A_k']:
                multiplicity['chiral_n2lo'] += 1
                z['chiral_n2lo'] = sig_fig(capped(z['chiral_n2lo'], z_min['chiral_n2lo'], z_max['chiral_n2lo']))
                prior[key] = gv.gvar(0, 1) *z['chiral_n2lo']

            # chiral_n3lo
            if key in ['A_ak', 'A_ap', 'A_kk', 'A_kp', 'A_pp']:
                multiplicity['chiral_n3lo'] += 1
                z['chiral_n3lo'] = sig_fig(capped(z['chiral_n3lo'], z_min['chiral_n3lo'], z_max['chiral_n3lo']))
                prior[key] = gv.gvar(0, 1) *z['chiral_n3lo']

            # latt_n2lo
            elif key in ['A_loga', 'A_a']:
                multiplicity['latt_n2lo'] += 1
                z['latt_n2lo'] = sig_fig(capped(z['latt_n2lo'], z_min['latt_n2lo'], z_max['latt_n2lo']))
                prior[key] = gv.gvar(0, 1) *z['latt_n2lo']

            # latt_n3lo
            elif key in ['A_aa']:
                multiplicity['latt_n3lo'] += 1
                z['latt_n3lo'] = sig_fig(capped(z['latt_n3lo'], z_min['latt_n3lo'], z_max['latt_n3lo']))
                prior[key] = gv.gvar(0, 1) *z['latt_n3lo']

            # latt_n4lo
            elif key in ['A_aaa']:
                multiplicity['latt_n4lo'] += 1
                z['latt_n4lo'] = sig_fig(capped(z['latt_n4lo'], z_min['latt_n4lo'], z_max['latt_n4lo']))
                prior[key] = gv.gvar(0, 1) *z['latt_n4lo']


        self.counter['iters'] += 1
        fitfcn = self._make_models()[-1].fitfcn
        print(self.counter['iters'], ' ', z)#{key : np.round(1. / z[key], 8) for key in z.keys()})

        # Jeffrey's prior
        def plausibility(s):
            plaus = 0
            for key in s:
                k = 1 / np.log(z_max[key]/z_min[key])
                plaus -= np.log(k/s[key]) *multiplicity[key]

            return plaus

        plaus = 0# plausibility(z)
        #print(plaus)

        return (dict(data=y_data, fcn=fitfcn, prior=prior), plaus)

    def _make_empbayes_fit(self):

        z0 = gv.BufferDict()

        # chiral terms
        if self.model_info['fit_type'] == 'poly':
            z0['chiral_nlo'] = 10.0
        if self.model_info['order'] in ['n2lo', 'n3lo']:
            z0['chiral_n2lo'] = 10.0
            z0['latt_n2lo'] = 10.0
        if self.model_info['order'] in ['n3lo']:
            z0['chiral_n3lo'] = 10.0
            z0['latt_n3lo'] = 10.0

        # latt terms
        if self.model_info['latt_ct'] in ['n2lo', 'n3lo', 'n4lo']:
            z0['latt_n2lo'] = 10.0
        if self.model_info['latt_ct'] in ['n3lo', 'n4lo']:
            z0['latt_n3lo'] = 10.0
        if self.model_info['latt_ct'] in ['n4lo']:
            z0['latt_n4lo'] = 10.0


        # Might need to change minargs default values for empbayes_fit to converge:
        # tol=1e-8, svdcut=1e-12, debug=False, maxit=1000, add_svdnoise=False, add_priornoise=False
        # Note: maxit != maxfev. See https://github.com/scipy/scipy/issues/3334
        # For Nelder-Mead algorithm, maxfev < maxit < 3 maxfev?

        # For debugging. Same as 'callback':
        # https://github.com/scipy/scipy/blob/c0dc7fccc53d8a8569cde5d55673fca284bca191/scipy/optimize/optimize.py#L651
        def analyzer(arg):
            self.counter['evals'] += 1
            print("\nEvals: ", self.counter['evals'], arg,"\n")
            print(type(arg[0]))
            return None




        fit, z = lsqfit.empbayes_fit(z0, fitargs=self._make_fitargs,
                                     tol=0.01, maxit=100, analyzer=None)

        print(z)
        self.empbayes_fit = fit
        return fit

    def _make_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()

        if self.model_info['fit_type'] == 'simultaneous':
            fitter = lsqfit.MultiFitter(models=models)
            fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
        else:
            fitter = lsqfit.MultiFitter(models=models)
            fit = fitter.lsqfit(data=y_data, prior=prior, fast=False, mopt=True)

        self.fit = fit
        print("Done!")
        return fit

    def _make_models(self):
        models = np.array([])

        # simultaneous fits -- don't chain nlo -> n2lo
        if self.model_info['fit_type'] == 'simultaneous':
            model_info = self.model_info.copy()
            for fit_type in ['ma', 'xpt']:
                model_info = self.model_info.copy()
                model_info['fit_type'] = fit_type
                models = np.append(models, fk_fpi_model(datatag=fit_type+'_'+model_info['order'],
                            **model_info))
            return models


        model_info = self.model_info.copy()
        models = np.append(models, fk_fpi_model(datatag=self.model_info['fit_type']+'_'+model_info['order'],
                    **self.model_info))

        return models

    def _make_prior(self, fit_data=None):
        if fit_data is None:
            fit_data = self.fit_data

        prior = self.prior

        newprior = gv.BufferDict()
        #model_info = self.model_info

        # Move fit_data into prior
        for key in fit_data:
            if key != 'y':
                newprior[key] = fit_data[key]

        if self.model_info['include_alpha_s'] :
            newprior['alpha_s'] = fit_data['alpha_s']
            newprior['A_loga'] = prior['A_loga']

        # Fit parameters, depending on fit type
        if self.model_info['fit_type'] in ['ma', 'ma-ratio', 'xpt', 'xpt-ratio']:
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']
        elif self.model_info['fit_type'] in ['poly']:
            newprior['A_x'] = prior['A_x']


        if self.model_info['include_log']:
            newprior['L_1'] = prior['L_1']
            newprior['L_2'] = prior['L_2']
            newprior['L_3'] = prior['L_3']
            newprior['L_4'] = prior['L_4']
            newprior['L_6'] = prior['L_6']
            newprior['L_7'] = prior['L_7']
            newprior['L_8'] = prior['L_8']

        if self.model_info['order'] in ['n2lo', 'n3lo']:
            newprior['A_k'] = prior['A_k']
            newprior['A_p'] = prior['A_p']
            newprior['A_a'] = prior['A_a']

        if self.model_info['order'] in ['n3lo']:
            newprior['A_aa'] = prior['A_aa']
            newprior['A_ak'] = prior['A_ak']
            newprior['A_ap'] = prior['A_ap']
            newprior['A_kk'] = prior['A_kk']
            newprior['A_kp'] = prior['A_kp']
            newprior['A_pp'] = prior['A_pp']


        # Lattice spacing terms
        if self.model_info['latt_ct'] in ['n2lo', 'n3lo', 'n4lo']:
            newprior['A_a'] = prior['A_a']
        if self.model_info['latt_ct'] in ['n3lo', 'n4lo']:
            newprior['A_aa'] = prior['A_aa']
        if self.model_info['latt_ct'] in ['n4lo']:
            newprior['A_aaa'] = prior['A_aaa']

        for key in self.model_info['exclude']:
            if key in newprior.keys():
                del(newprior[key])

        '''
        output = gv.BufferDict()
        for key in newprior:
            if key not in fit_data:
                sdev = np.log(gv.sdev(newprior[key]))
                output['log('+key+')'] = gv.gvar(sdev, sdev)
            else:
                output[key] = newprior[key]

        return output
        '''

        return newprior

    def _make_y_data(self, fit_data=None):
        if fit_data is None:
            fit_data = self.fit_data

        if 'y' not in fit_data.keys():
            return None
        else:
            return fit_data['y']

    def get_fit(self):
        if self.fit is None:
            return self._make_fit()
        else:
            return self.fit

    def get_empbayes_fit(self):
        if self.empbayes_fit is None:
            return self._make_empbayes_fit()
        else:
            return self.empbayes_fit


class fk_fpi_model(lsqfit.MultiFitterModel):

    def __init__(self, datatag,
            fit_type, order, latt_ct, F2,
            include_FV, exclude,
            include_alpha_s, include_log, include_log2, include_sunset,
            fast_sunset=False, **kwargs
        ):
        super(fk_fpi_model, self).__init__(datatag)

        # Model info
        self.model_info = {
            'fit_type' : fit_type,
            'order' : order,
            'latt_ct' : latt_ct,
            'F2' : F2,

            'include_FV' : include_FV,
            'exclude' : exclude,

            'include_alpha_s' : include_alpha_s,
            'include_log' : include_log,
            'include_log2' : include_log2,
            'include_sunset' : include_sunset,
        }

        self.fast_sunset = fast_sunset # Use correlated gvar instead
        self.debug = False


    def fitfcn(self, p, fit_data=None, debug=None):

        if debug:
            self.debug = debug
            self.debug_table = {}

        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        for key in self.model_info['exclude']:
            p[key] = 0

        # nlo fits
        if self.model_info['order'] in ['nlo', 'n2lo', 'n3lo']:
            # mixed-action/xpt fits
            if self.model_info['fit_type'] == 'ma-ratio':
                output = self.fitfcn_nlo_ma_ratio(p)
            elif self.model_info['fit_type'] == 'ma':
                output = self.fitfcn_nlo_ma(p)
            elif self.model_info['fit_type'] == 'xpt-ratio':
                output = self.fitfcn_nlo_xpt_ratio(p)
            elif self.model_info['fit_type'] == 'xpt':
                output = self.fitfcn_nlo_xpt(p)
            elif self.model_info['fit_type'] == 'poly':
                output = self.fitfcn_nlo_polynomial(p)

        if debug:
            temp = output

        # n2lo corrections
        if self.model_info['order'] in ['n2lo', 'n3lo']:
            output = output + self.fitfcn_n2lo_pure_ct(p)
            if self.model_info['fit_type'] != 'poly':
                output = output + self.fitfcn_n2lo_ratio(p)
                output = output + self.fitfcn_n2lo_renormalization_ct(p)
        elif self.model_info['latt_ct'] in ['n2lo', 'n3lo', 'n4lo']:
            output = output + self.fitfcn_n2lo_latt_spacing_ct(p)

        # semi-n2lo corrections
        if self.model_info['include_log']:
            output = output + self.fitfcn_semi_n2lo_log_ct(p)
        if self.model_info['include_log2']:
            output = output + self.fitfcn_semi_n2lo_log_squared_ct(p)
        if self.model_info['include_sunset']:
            output = output + self.fitfcn_semi_n2lo_sunset_ct(p)
        if self.model_info['include_alpha_s']:
            output = output + self.fitfcn_semi_n2lo_alpha_s_ct(p)


        # n3lo corrections
        if self.model_info['order'] in ['n3lo']:
            output = output + self.fitfcn_n3lo_pure_ct(p)
        elif self.model_info['latt_ct'] in ['n3lo', 'n4lo']:
            output = output + self.fitfcn_n3lo_latt_spacing_ct(p)

        # n4lo corrections
        if self.model_info['latt_ct'] in ['n4lo']:
            output = output + self.fitfcn_n4lo_latt_spacing_ct(p)

        for key in self.model_info['exclude']:
            del(p[key])

        if debug:
            #print(gv.tabulate(self.debug_table))
            temp_string = ''
            for key in self.debug_table:
                temp_string +='  % .6f:  %s\n' %(self.debug_table[key], key)
            temp_string +='   -----\n'
            temp_string +='  % .6f:  %s\n' %(output, 'total')
            print(temp_string)

        return output

    def fitfcn_nlo_ma(self, p):
        # Constants
        if self.model_info['include_FV']:
            order_vol = 10
        else:
            order_vol = 0

        # Independent variables
        mju = p['mju']
        mpi = p['mpi']
        mk  = p['mk']
        mru = p['mru']
        msj = p['mjs'] #msj = mjs
        mss = p['mss']
        mrs = p['mrs']
        mx = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2) + p['a2DI'])

        del2_ju = p['a2DI']
        del2_rs = del2_ju
        L = p['L']

        lam2_chi = p['lam2_chi']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        output = (
            1

            + (1.0/2.0) *sf.fcn_I_m(mju, L, mu, order_vol) / F2
            + (1.0/8.0) *sf.fcn_I_m(mpi, L, mu, order_vol) / F2
            + (1.0/4.0) *sf.fcn_I_m(mru, L, mu, order_vol) / F2
            - (1.0/2.0) *sf.fcn_I_m(msj, L, mu, order_vol) / F2
            + (1.0/4.0) *sf.fcn_I_m(mss, L, mu, order_vol) / F2
            - (1.0/4.0) *sf.fcn_I_m(mrs, L, mu, order_vol) / F2
            - (3.0/8.0) *sf.fcn_I_m(mx, L, mu, order_vol) / F2

            + del2_ju *(
                - sf.fcn_dI_m(mpi, L, mu, order_vol) / (8 *F2)
                + sf.fcn_K_mM((mpi, mx), L, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *(
                + sf.fcn_K21_mM((mpi, mx), L, mu, order_vol) / (24 *F2)
            )
            + del2_ju *del2_rs *(
                - sf.fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2)
                + sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + sf.fcn_K_mM((mss, mx), L, mu, order_vol) / (4 *F2)
                - (mk)**2 *sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
                + (mpi)**2 *sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
            )

            + 4 *(4 *np.pi)**2 *(mk**2 - mpi**2) / lam2_chi *p['L_5']
        )

        if self.debug:
            self.debug_table['nlo_ma'] = output


        return output


    def fitfcn_nlo_ma_ratio(self, p):
        # Constants
        if self.model_info['include_FV']:
            order_vol = 10
        else:
            order_vol = 0

        # Independent variables
        mju = p['mju']
        mpi = p['mpi']
        mk = p['mk']
        mru = p['mru']
        msj = p['mjs'] #msj = mjs
        mss = p['mss']
        mrs = p['mrs']
        mx = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2) + p['a2DI'])

        del2_ju = p['a2DI']
        del2_rs = del2_ju
        L = p['L']

        lam2_chi = p['lam2_chi']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        eps2_pi = mpi**2 / lam2_chi
        eps2_k = mk**2 / lam2_chi

        Fpi_nlo_per_F0 = (
            + 1

            - sf.fcn_I_m(mju, L, mu, order_vol) / F2
            - (1.0/2.0) *sf.fcn_I_m(mru, L, mu, order_vol) / F2

            + 4 *eps2_pi *(4 *np.pi)**2 *(p['L_4'] + p['L_5'])
            + 8 *eps2_k *(4 *np.pi)**2 *p['L_4']
        )

        FK_nlo_per_F0 = (
            + 1

            - (1.0/2.0) *sf.fcn_I_m(mju, L, mu, order_vol) / F2
            + (1.0/8.0) *sf.fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/4.0) *sf.fcn_I_m(mru, L, mu, order_vol) / F2
            - (1.0/2.0) *sf.fcn_I_m(msj, L, mu, order_vol) / F2
            - (1.0/4.0) *sf.fcn_I_m(mrs, L, mu, order_vol) / F2
            + (1.0/4.0) *sf.fcn_I_m(mss, L, mu, order_vol) / F2
            - (3.0/8.0) *sf.fcn_I_m(mx, L, mu, order_vol) / F2

            + 4 *eps2_pi *(4 *np.pi)**2 *p['L_4']
            + 4 *eps2_k *(4 *np.pi)**2 *(2 *p['L_4'] + p['L_5'])

            + del2_ju *(
                - sf.fcn_dI_m(mpi, L, mu, order_vol) / (8 *F2)
                + sf.fcn_K_mM((mpi, mx), L, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *(
                + sf.fcn_K21_mM((mpi, mx), L, mu, order_vol) / (24 *F2)
            )
            + del2_ju *del2_rs *(
                - sf.fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2) # May need to change sign
                + sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + sf.fcn_K_mM((mss, mx), L, mu, order_vol) / (4 *F2)
                - (mk)**2 *sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
                + (mpi)**2 *sf.fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
            )
        )

        output = FK_nlo_per_F0 / Fpi_nlo_per_F0

        if self.debug:
            self.debug_table['nlo_ma-ratio'] = output

        return output


    def fitfcn_nlo_polynomial(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = 1 + (eps2_k - eps2_pi) *p['A_x']

        if self.debug:
            self.debug_table['nlo_poly'] = output

        return output


    def fitfcn_nlo_xpt(self, p):

        # Constants
        if self.model_info['include_FV']:
            order_vol = 10
        else:
            order_vol = 0

        # Independent variables
        mpi = p['mpi']
        mk = p['mk']
        meta = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2))

        lam2_chi = p['lam2_chi']
        eps2_pi = mpi**2 / lam2_chi
        eps2_k = mk**2 / lam2_chi

        L = p['L']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        output = (
            1
            + 1*(
            + (5.0/8.0) *sf.fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/4.0) *sf.fcn_I_m(mk, L, mu, order_vol) / F2
            - (3.0/8.0) *sf.fcn_I_m(meta, L, mu, order_vol) / F2)
            + 4 *(eps2_k - eps2_pi) *(4 *np.pi)**2 *p['L_5']
        )

        if self.debug:
            self.debug_table['nlo_xpt'] = output

        return output


    def fitfcn_nlo_xpt_ratio(self, p):
       # Constants
        if self.model_info['include_FV']:
            order_vol = 10
        else:
            order_vol = 0

        # Independent variables
        mpi = p['mpi']
        mk = p['mk']
        meta = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2))

        lam2_chi = p['lam2_chi']
        eps2_pi = mpi**2 / lam2_chi
        eps2_k = mk**2 / lam2_chi
        #eps2_eta = meta**2 / lam2_chi

        L = p['L']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        Fpi_nlo_per_F0 = (1
            - sf.fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/2.0) *sf.fcn_I_m(mk, L, mu, order_vol) / F2
            + 4 *eps2_pi *(4 *np.pi)**2 *(p['L_4'] + p['L_5'])
            + 8 *eps2_k *(4 *np.pi)**2 *p['L_4']
        )

        FK_nlo_per_F0 = (1
            - (3.0/8.0) *sf.fcn_I_m(mpi, L, mu, order_vol) / F2
            - (3.0/4.0) *sf.fcn_I_m(mk, L, mu, order_vol) / F2
            - (3.0/8.0) *sf.fcn_I_m(meta, L, mu, order_vol) / F2
            + 4 *eps2_pi *(4 *np.pi)**2 *p['L_4']
            + 4 *eps2_k *(4 *np.pi)**2 *(2 *p['L_4'] + p['L_5'])
        )

        output = FK_nlo_per_F0 / Fpi_nlo_per_F0

        if self.debug:
            self.debug_table['nlo_xpt'] = output

        return output


    def fitfcn_n2lo_pure_ct(self, p, include_log=None):
        if include_log is None:
            include_log = self.model_info['include_log']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + (eps2_a) *p['A_a']
            + (eps2_k) *p['A_k']
            + (eps2_pi) *p['A_p']
        ) *(eps2_k - eps2_pi)

        if include_log:
            output += (4 *np.pi)**2 *(eps2_k - eps2_pi) *(
                + eps2_k *(
                    + 8 *(4 *np.pi)**2 *p['L_5']*(
                        + 8 *p['L_4'] + 3 *p['L_5'] - 16 *p['L_6'] - 8 *p['L_8']
                    )
                    - 2 *p['L_1'] - p['L_2'] - 1./18 *p['L_3'] +4./3 *p['L_5'] - 16 *p['L_7'] - 8 *p['L_8']
                )
                + eps2_pi *(
                    + 8 *(4 *np.pi)**2 *p['L_5']*(
                        + 4 *p['L_4'] + 5 *p['L_5'] - 8 *p['L_6'] - 8 *p['L_8']
                    )
                    - 2 *p['L_1'] - p['L_2'] - 5./18 *p['L_3'] - 4./3 *p['L_5'] + 16 *p['L_7'] + 8 *p['L_8']
                )
            )

        if self.debug:
            self.debug_table['n2lo_ct'] = output

        return output

    def fitfcn_n2lo_ratio(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        fcn_l = lambda x : x *np.log(x)

        if self.model_info['fit_type'] in ['ma-ratio', 'xpt-ratio']:
            dFK = -(
                    + 3./8 *fcn_l(eps2_pi) + 3./4 *fcn_l(eps2_k) + 3./8 *fcn_l(eps2_eta)
                    - 4 *(4*np.pi)**2 *(eps2_pi *p['L_4'] + eps2_k *(2 *p['L_4'] + p['L_5']))
                )

            dFpi = -(
                    + fcn_l(eps2_pi) + 1./2 *fcn_l(eps2_k)
                    - 4 *(4*np.pi)**2 *(eps2_pi *(p['L_4'] + p['L_5']) + 2 *eps2_k *p['L_4'])
                )

            output = dFK *dFpi - dFpi**2
        else:
            output = 0

        if self.debug:
            self.debug_table['n2lo_ratio'] = output

        return output

    def fitfcn_n2lo_renormalization_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        fcn_l = lambda x : x *np.log(x)

        dFK = -(
                + 3./8 *fcn_l(eps2_pi) + 3./4 *fcn_l(eps2_k) + 3./8 *fcn_l(eps2_eta)
                - 4 *(4*np.pi)**2 *(eps2_pi *p['L_4'] + eps2_k *(2 *p['L_4'] + p['L_5']))
            )

        dFpi = -(
                + fcn_l(eps2_pi) + 1./2 *fcn_l(eps2_k)
                - 4 *(4*np.pi)**2 *(eps2_pi *(p['L_4'] + p['L_5']) + 2 *eps2_k *p['L_4'])
            )

        if self.model_info['F2'] == 'FpiFpi':
            output = -3 / 2. *(eps2_k - eps2_pi) *dFpi

        elif self.model_info['F2'] == 'FKFpi':
            output = -3 / 4. *(eps2_k - eps2_pi) *(dFK + dFpi)
            output += (dFK - dFpi)**2

        elif self.model_info['F2'] == 'FKFK':
            output = -3 / 2. *(eps2_k - eps2_pi) *dFK
            output += 2 *(dFK - dFpi)**2

        if self.debug:
            self.debug_table['n2lo_scale'] = output

        return output

    def fitfcn_semi_n2lo_alpha_s_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        alpha_s = p['alpha_s']

        output = (
            alpha_s *eps2_a *(eps2_k - eps2_pi) *p['A_loga']
        )

        if self.debug:
            self.debug_table['n2lo_alphaS'] = output

        return output


    def fitfcn_semi_n2lo_log_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        output = (
            + sf.fcn_Cr_j(1, eps2_pi, eps2_k, p) *np.log(eps2_pi)
            + sf.fcn_Cr_j(2, eps2_pi, eps2_k, p) *np.log(eps2_k)
            + sf.fcn_Cr_j(3, eps2_pi, eps2_k, p) *np.log(eps2_eta)
        )

        if self.debug:
            self.debug_table['n2lo_log'] = output

        return output

    def fitfcn_semi_n2lo_log_squared_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        output = (
            + sf.fcn_Kr_j(1, eps2_pi, eps2_k) *np.log(eps2_pi) *np.log(eps2_pi)
            + sf.fcn_Kr_j(2, eps2_pi, eps2_k) *np.log(eps2_pi) *np.log(eps2_k)
            + sf.fcn_Kr_j(3, eps2_pi, eps2_k) *np.log(eps2_pi) *np.log(eps2_eta)
            + sf.fcn_Kr_j(4, eps2_pi, eps2_k) *np.log(eps2_k) *np.log(eps2_k)
            + sf.fcn_Kr_j(5, eps2_pi, eps2_k) *np.log(eps2_k) *np.log(eps2_eta)
            + sf.fcn_Kr_j(6, eps2_pi, eps2_k) *np.log(eps2_eta) *np.log(eps2_eta)
        )

        if self.debug:
            self.debug_table['n2lo_logSq'] = output

        return output

    def fitfcn_semi_n2lo_sunset_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        if self.fast_sunset and 'sunset' in p:
            output = eps2_k**2 *p['sunset']
        else:
            output = eps2_k**2 *sf.fcn_FF(eps2_pi/eps2_k)

        if self.debug:
            self.debug_table['n2lo_sunset'] = output

        return output

    def fitfcn_n2lo_latt_spacing_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + (eps2_a) *p['A_a']
        ) *(eps2_k - eps2_pi)

        if self.debug:
            self.debug_table['n2lo_a2'] = output

        return output

    def fitfcn_n3lo_latt_spacing_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + (eps2_a)**2 *p['A_aa']
        ) *(eps2_k - eps2_pi)

        if self.debug:
            self.debug_table['n3lo_a4'] = output

        return output

    def fitfcn_n3lo_pure_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + eps2_a *(
                + eps2_a  *p['A_aa']
                + eps2_k  *p['A_ak']
                + eps2_pi *p['A_ap']
            )
            + eps2_k *(
                + eps2_k  *p['A_kk']
                + eps2_pi *p['A_kp']
            )
            + eps2_pi *(
                + eps2_pi *p['A_pp']
            )
        ) *(eps2_k - eps2_pi)

        if self.debug:
            self.debug_table['n3lo_ct'] = output

        return output

    def fitfcn_n4lo_latt_spacing_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / 4
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + (eps2_a)**4 *p['A_aaa']
        ) *(eps2_k - eps2_pi)

        if self.debug:
            self.debug_table['n4lo_a6'] = output

        return output

    def buildprior(self, prior, mopt=None, extend=False):
        return prior

    def builddata(self, data):
        return data
