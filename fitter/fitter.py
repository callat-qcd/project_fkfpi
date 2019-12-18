import lsqfit
import numpy as np
import gvar as gv
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import special_functions as sf

class fitter(object):

    def __init__(self, order, fit_type, F2, fit_data=None, prior=None, fast_sunset=False):
        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.empbayes_fit = None
        self.order = order
        self.fit_type = fit_type
        self.F2 = F2
        self.fast_sunset=fast_sunset

        # To force empbayes_fit to converge?
        self.counter = 0
        #self.z = {}

    def _make_fitargs(self, z):
        y_data = self._make_y_data()
        prior = self._make_prior()

        #Force convergence
        #plausibility = 0
        #if self.counter > 100:
        #    plausibility = -10 *np.random.random()
        #    z['chiral'] = self.z['chiral'] #np.abs(z['chiral'])
        #    z['spacing_n2lo'] = self.z['spacing_n2lo']
        #    if self.order['include_latt_n3lo']:
        #        z['spacing_n3lo'] = self.z['spacing_n3lo']
        #else:


        # Ideally:
            # Don't bother with more than the hundredth place
            # Don't let z=0 (=> null GBF)
            # Don't bother with negative values (meaningless)
        # But for some reason, these restrictions (other than the last) cause empbayes_fit not to converge
        z['chiral'] = np.abs(z['chiral']) #np.max([np.abs(np.around(z['chiral'], 2)), 0.01])
        z['spacing_n2lo'] = np.abs(z['spacing_n2lo']) #np.max([np.abs(np.around(z['spacing_n2lo'], 2)), 0.01]) #
        if self.order['include_latt_n3lo']:
            z['spacing_n3lo'] = np.abs(z['spacing_n3lo']) #np.max([np.abs(np.around(z['spacing_n3lo'], 2)), 0.01]) #

            # Force convergence
            #self.z ={}
            #for key in z:
            #    self.z[key] = z[key]

        for key in prior.keys():
            if key in ['A_p', 'A_k']:
                prior[key] = gv.gvar(0, 1) *z['chiral']
            if key in ['A_loga', 'A_a']:
                prior[key] = gv.gvar(0, 1) *z['spacing_n2lo']
            if key in ['A_aa']:
                prior[key] = gv.gvar(0, 1) *z['spacing_n3lo']

        self.counter += 1
        fitfcn = self._make_models()[-1].fitfcn
        print(self.counter, ' ', z)

        return dict(data=y_data, fcn=fitfcn, prior=prior)#, plausibility

    def _make_empbayes_fit(self):
        #models = self._make_models(fast_sunset=True)
        #y_data = self._make_y_data()
        #prior = self._make_prior()

        z0 = gv.BufferDict()
        z0['chiral'] = 1.0
        z0['spacing_n2lo'] = 1.0
        if self.order['include_latt_n3lo']:
            z0['spacing_n3lo'] = 1.0

        fit, z = lsqfit.empbayes_fit(z0, fitargs = self._make_fitargs, tol=0.01)
        self.empbayes_fit = fit
        return fit

    def _make_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()

        if self.fit_type == 'simultaneous':
            fitter = lsqfit.MultiFitter(models=models)
            fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
        else:
            fitter = lsqfit.MultiFitter(models=models)
            fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)

        self.fit = fit
        print("Done!")
        return fit

    def _make_models(self):
        models = np.array([])

        # simultaneous fits -- don't chain nlo -> nnlo
        if self.fit_type == 'simultaneous':
            order = self.order.copy()
            for fit_type in ['ma', 'xpt']:
                models = np.append(models, fk_fpi_model(datatag=fit_type+'_'+order['fit'],
                            F2=self.F2, order=self.order, fit_type=fit_type, fast_sunset=self.fast_sunset))
            return models


        order = self.order.copy()
        models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                    F2=self.F2, order=order, fit_type=self.fit_type, fast_sunset=self.fast_sunset))

        return models

    def _make_prior(self, fit_data=None):
        if fit_data is None:
            fit_data = self.fit_data

        prior = self.prior

        newprior = gv.BufferDict()
        order = self.order

        # Move fit_data into prior
        for key in fit_data:
            if key != 'y':
                newprior[key] = fit_data[key]

        if order['include_alpha_s'] :
            newprior['alpha_s'] = fit_data['alpha_s']
            newprior['A_loga'] = prior['A_loga']

        # Fit parameters, depending on fit type
        newprior['L_4'] = prior['L_4']
        newprior['L_5'] = prior['L_5']

        if order['include_log']:
            newprior['L_1'] = prior['L_1']
            newprior['L_2'] = prior['L_2']
            newprior['L_3'] = prior['L_3']
            newprior['L_4'] = prior['L_4']
            newprior['L_6'] = prior['L_6']
            newprior['L_7'] = prior['L_7']
            newprior['L_8'] = prior['L_8']

        if order['fit'] in ['nnlo', 'nnnlo']:
            newprior['A_k'] = prior['A_k']
            newprior['A_p'] = prior['A_p']
            newprior['A_a'] = prior['A_a']

        if order['include_latt_n3lo']:
            newprior['A_aa'] = prior['A_aa']

        for key in self.order['exclude']:
            if key in newprior.keys():
                del(newprior[key])

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

    def __init__(self, datatag, fit_type, order, F2, fast_sunset=False):
        super(fk_fpi_model, self).__init__(datatag)

        self.fast_sunset = fast_sunset # Use correlated gvar instead
        self.fit_type = fit_type
        self.order = order
        self.F2 = F2
        self.debug = False


    def fitfcn(self, p, fit_data=None, debug=None):
        if debug is not None:
            self.debug = debug

        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        for key in self.order['exclude']:
            p[key] = 0

        # NLO fits
        if self.order['fit'] in ['nlo', 'nnlo', 'nnnlo']:
            # mixed-action/xpt fits
            if self.fit_type == 'ma-ratio':
                output = self.fitfcn_nlo_ma_ratio(p)
            elif self.fit_type == 'ma':
                output = self.fitfcn_nlo_ma(p)
            elif self.fit_type == 'xpt-ratio':
                output = self.fitfcn_nlo_xpt_ratio(p)
            elif self.fit_type == 'xpt':
                output = self.fitfcn_nlo_xpt(p)

        if debug:
            temp = output

        if self.order['fit'] in ['nnlo', 'nnnlo']:
            output = output + self.fitfcn_nnlo_pure_ct(p)
            output = output + self.fitfcn_nnlo_renormalization_ct(p)

        # semi-nnlo corrections
        if self.order['include_log']:
            output = output + self.fitfcn_seminnlo_log_ct(p)
        if self.order['include_log2']:
            output = output + self.fitfcn_seminnlo_log_squared_ct(p)
        if self.order['include_sunset']:
            output = output + self.fitfcn_seminnlo_sunset_ct(p)
        if self.order['include_alpha_s']:
            output = output + self.fitfcn_seminnlo_alpha_s_ct(p)

        if self.order['include_latt_n3lo']:
            output = output + self.fitfcn_nnnlo_latt_spacing_ct(p)

        for key in self.order['exclude']:
            del(p[key])

        if self.debug:
            print('nnlo', output - temp)
            print('p', p)

        return output

    def fitfcn_nlo_ma(self, p):
        # Constants
        order_vol = self.order['vol']

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

        return output


    def fitfcn_nlo_ma_ratio(self, p):
        # Constants
        order_vol = self.order['vol']

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

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_nlo_xpt(self, p):

        # Constants
        order_vol = self.order['vol']

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

        return output


    def fitfcn_nlo_xpt_ratio(self, p):
       # Constants
        order_vol = self.order['vol']

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

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_nnlo_pure_ct(self, p, include_log=None):
        if include_log is None:
            include_log = self.order['include_log']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
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

        return output

    def fitfcn_nnlo_renormalization_ct(self, p):
        lam2_chi = p['lam2_chi']
        #eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        fcn_l = lambda x : x *np.log(x)

        xi = (
            + (5./8) *fcn_l(eps2_pi)
            - (1./4) *fcn_l(eps2_k)
            - (3./8) *fcn_l(eps2_eta)

            + 4 *(4 *np.pi)**2 *(
                + eps2_pi *p['L_5']
            )
        )

        if self.F2 == 'FKFK':
            output = (3./2) *(eps2_k - eps2_pi) *(
                + (3./8) *fcn_l(eps2_pi)
                + (3./4) *fcn_l(eps2_k)
                + (3./8) *fcn_l(eps2_eta)

                - 4 *(4 *np.pi)**2 *(
                    + eps2_pi *p['L_4']
                    + eps2_k *(2 *p['L_4'] + p['L_5'])
                )

                + 2 *xi**2
            )

        elif self.F2 == 'FKFpi':
            output = (3./2) *(eps2_k - eps2_pi) *(
                + (11./8) *fcn_l(eps2_pi)
                + (5./4) *fcn_l(eps2_k)
                + (3./8) *fcn_l(eps2_eta)

                - 4 *(4 *np.pi)**2 *(
                    + eps2_pi *(2 *p['L_4'] + p['L_5'])
                    + eps2_k *(4 *p['L_4'] + p['L_5'])
                )

                + xi**2
            )

        elif self.F2 == 'FpiFpi':
            output = (3./2) *(eps2_k - eps2_pi) *(
                + fcn_l(eps2_pi)
                + (1./2) * fcn_l(eps2_k)

                - 4 *(4 *np.pi)**2 *(
                    + eps2_pi *(p['L_4'] + p['L_5'])
                    + eps2_k *(2 *p['L_4'])
                )
            )

        if self.fit_type in ['ma-ratio', 'xpt-ratio']:
            output = output + (
                (
                    + fcn_l(eps2_pi) + 1./2 *fcn_l(eps2_k)
                    - 4 *(4*np.pi)**2 *(eps2_pi *(p['L_4'] + p['L_5']) + 2 *eps2_k *p['L_4'])
                ) *(
                    + 3./8 *fcn_l(eps2_pi) + 3./4 *fcn_l(eps2_k) + 3./8 *fcn_l(eps2_eta)
                    + 4 *(4*np.pi)**2 *(eps2_pi *p['L_4'] + eps2_k *(2 *p['L_4'] + p['L_5']))
                )

                + 1./2 *(
                    + fcn_l(eps2_pi) + 1./2 *fcn_l(eps2_k)
                    - 4 *(4*np.pi)**2 *(eps2_pi *(p['L_4'] + p['L_5']) + 2 *eps2_k *p['L_4'])
                )**2
            )

        if self.debug:
            print('mu fix', output)
        return output

    def fitfcn_seminnlo_alpha_s_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        alpha_s = p['alpha_s']

        output = (
            alpha_s *eps2_a *(eps2_k - eps2_pi) *p['A_loga']
        )
        return output


    def fitfcn_seminnlo_log_ct(self, p):
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
            print('log', output)
        return output

    def fitfcn_seminnlo_log_squared_ct(self, p):
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
            print('logSq', output)
        return output

    def fitfcn_seminnlo_sunset_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        if self.fast_sunset:
            print('yes')
            output = eps2_k**2 *p['sunset']
        else:
            output = eps2_k**2 *sf.fcn_FF(eps2_pi/eps2_k)

        if self.debug:
            print('FF', output)
        return output


    def fitfcn_nnnlo_latt_spacing_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            + (eps2_a)**2 *p['A_aa']
        ) *(eps2_k - eps2_pi)

        return output

    def buildprior(self, prior, mopt=None, extend=False):
        return prior

    def builddata(self, data):
        return data
