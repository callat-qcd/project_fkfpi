import lsqfit
import numpy as np
import gvar as gv

# import fcn_T_m, fcn_dI_m, etc
from special_functions import *

class fitter(object):

    def __init__(self, order, fit_type, fit_data=None, prior=None, chain_fits=True):
        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.empbayes_fit = None
        self.order = order
        self.fit_type = fit_type
        self.chain_fits = chain_fits

    def _make_fitargs(self, z):
        y_data = self._make_y_data()
        prior = self._make_prior()

        if self.order['fit'] == 'nlo':
            for key in prior.keys():
                if key in ['L_5', 'L_4']:
                    prior[key] = prior[key] *z

        if self.order['fit'] == 'nnlo' and not self.order['include_log']:
            for key in prior.keys():
                if key in ['A_p', 'A_k']:
                    prior[key] = prior[key] *z

        if self.order['fit'] == 'nnlo' and self.order['include_log']:
            for key in prior.keys():
                if key in ['A_loga', 'A_a']:
                    prior[key] = prior[key] *z

        #for key in prior.keys():
        #    if key in ['A_a', 'A_p', 'A_k']:
        #        prior[key] = prior[key] *z

        #if self.order['include_log']:
        #    prior['A_loga'] = prior['A_loga'] *z

        fitfcn = self._make_models()[-1].fitfcn

        return dict(data=y_data, fcn=fitfcn, prior=prior)

    def _make_empbayes_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()
        for j, model in enumerate(models):
            # Make model up till nlo/nnlo
            if j < len(models) - 1:
                fitter = lsqfit.MultiFitter(models=model)
                fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
                for key in fit.p:
                    if key in ['L_4', 'L_5']:
                        prior[key] = gv.gvar(fit.pmean[key], 3*fit.psdev[key])
                    elif key in ['A_a', 'A_p', 'A_k']:
                        prior[key] = fit.p[key]
                    elif (key in ['A_loga']) and (self.order['include_log']):
                        prior[key] = fit.p[key]

            # For nnlo/nnnlo, determine best parameter using empircal Bayes criterion
            else:
                fit, z = lsqfit.empbayes_fit(1.0, self._make_fitargs)
                #print z

        return fit

    def _make_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()

        if self.fit_type == 'simultaneous':
            fitter = lsqfit.MultiFitter(models=models)
            fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
        else:
            for model in models:
                fitter = lsqfit.MultiFitter(models=model)
                fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
                for key in fit.p:
                    if key in ['L_4', 'L_5']:
                        #pass
                        prior[key] = gv.gvar(fit.pmean[key], 3*fit.psdev[key])
                    elif key in ['A_a', 'A_p', 'A_k']:
                        prior[key] = fit.p[key]
                        #pass

        self.fit = fit
        print "Done!"
        return fit

    def _make_models(self):
        models = np.array([])

        # simultaneous fits -- don't chain nlo -> nnlo
        if self.fit_type == 'simultaneous':
            order = self.order.copy()
            for fit_type in ['ma', 'xpt']:
                models = np.append(models, fk_fpi_model(datatag=fit_type+'_'+order['fit'],
                            order=self.order, fit_type=fit_type))
            return models

        if not self.chain_fits:
            order = self.order.copy()
            models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                        order=order, fit_type=self.fit_type))
            return models

        if self.order['fit'] in ['nlo', 'nnlo', 'nnnlo']:
            order = self.order.copy()
            order['fit'] = 'nlo'
            models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                        order=order, fit_type=self.fit_type))
        if self.order['fit'] in ['nnlo', 'nnnlo']:
            order = self.order.copy()
            order['fit'] = 'nnlo'
            models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                        order=order, fit_type=self.fit_type))
        if self.order['fit'] in ['nnnlo']:
            order = self.order.copy()
            order['fit'] = 'nnnlo'
            models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                        order=order, fit_type=self.fit_type))
        return models

    def _make_prior(self, fit_data=None):
        if fit_data is None:
            fit_data = self.fit_data

        prior = self.prior

        newprior = gv.BufferDict()
        order = self.order

        # Move fit_data into prior
        #newprior['Fpi'] = fit_data['Fpi']
        #newprior['FK'] = fit_data['FK']
        newprior['mpi'] = fit_data['mpi']
        newprior['mk'] = fit_data['mk']
        newprior['mss'] = fit_data['mss']
        newprior['mju'] = fit_data['mju']
        newprior['mjs'] = fit_data['mjs']
        newprior['mru'] = fit_data['mru']
        newprior['mrs'] = fit_data['mrs']
        newprior['a2DI'] = fit_data['a2DI']
        newprior['lam2_chi'] = fit_data['lam2_chi']
        #newprior['a'] = fit_data['a']
        newprior['L'] = fit_data['L']
        #newprior['w0'] = fit_data['w0']
        newprior['a/w0'] = fit_data['a/w0']

        #for key in newprior.keys():
        #    newprior[key] = gv.gvar(gv.mean(newprior[key]), gv.sdev(newprior[key])/100)

        # Fit parameters, depending on fit type
        if self.fit_type == 'xpt-ratio':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'xpt':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma-ratio':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma-old':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'simultaneous':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']


        if order['fit'] in ['nnlo', 'nnnlo']:
            newprior['A_a'] = prior['A_a']
            newprior['A_k'] = prior['A_k']
            newprior['A_p'] = prior['A_p']

        if order['fit'] in ['nnnlo']:
            newprior['A_aa'] = prior['A_aa']
            newprior['A_ak'] = prior['A_ak']
            newprior['A_ap'] = prior['A_ap']
            newprior['A_kk'] = prior['A_kk']
            newprior['A_kp'] = prior['A_kp']
            newprior['A_pp'] = prior['A_pp']

        if order['include_log'] == True:
            newprior['A_loga'] = prior['A_loga']
            newprior['alpha_s'] = fit_data['alpha_s']

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

    def __init__(self, datatag, fit_type, order):
        super(fk_fpi_model, self).__init__(datatag)

        # Fit data
        self.order = order
        self.fit_type = fit_type

    def fitfcn(self, p, fit_data=None):
        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        if self.order['fit'] in ['nlo', 'nnlo', 'nnnlo']:
            # mixed-action/xpt fits
            if self.fit_type == 'ma-ratio':
                output = self.fitfcn_ma_ratio(p)
            elif self.fit_type == 'ma':
                output = self.fitfcn_ma(p)
            elif self.fit_type == 'xpt-ratio':
                output = self.fitfcn_xpt_ratio(p)
            elif self.fit_type == 'xpt':
                output = self.fitfcn_xpt(p)
            elif self.fit_type == 'ma-old':
                output = self.fitfcn_ma_old(p)

        #print "model: ", output

        if self.order['fit'] in ['nnlo', 'nnnlo']:
            output = output + self.fitfcn_nnlo_cts(p)

        if self.order['fit'] in ['nnnlo']:
            output = output + self.fitfcn_nnnlo_cts(p)

        if self.order['include_log']:
            output = output + self.fitfcn_nnlo_log_ct(p)

        if self.order['include_log2']:
            output = output + self.fitfcn_nnlo_log2_ct(p)

        return output

    def fitfcn_nnlo_cts(self, p):
        #w0 = p['w0']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        for key in self.order['exclude']:
            p[key] = 0

        output = (
            + (eps2_a) *p['A_a']
            + (eps2_k) *p['A_k']
            + (eps2_pi) *p['A_p']
        )

        for key in self.order['exclude']:
            del(p[key])
        #print "+ nnlo", output *(eps2_k - eps2_pi)

        #print "A_a: ", np.median(eps2_a)
        #print "A_k: ", np.median(eps2_k)
        #print "A_pi: ", np.median(eps2_pi)

        #print "\n\n"

        #try:
        #    print gv.corr(p['mpi'][0], p['mk'][0])
        #except ValueError:
        #    pass

        return output *(eps2_k - eps2_pi)

    def fitfcn_nnnlo_cts(self, p):

        #w0 = p['w0']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        for key in self.order['exclude']:
            p[key] = 0

        output = (
              + (eps2_a) *(
                + eps2_a *p['A_aa']
                #+ (eps2_k - eps2_pi) *p['A_ax']
                + eps2_k *p['A_ak']
                + eps2_pi *p['A_ap']
              )
              #+ (eps2_k - eps2_pi) *(
                #+ (eps2_k - eps2_pi) *p['A_xx']
                #+ eps2_k *p['A_xk']
                #+ eps2_pi *p['A_xp']
              #)
              + (eps2_k) *(
                + eps2_k *p['A_kk']
                + eps2_pi *p['A_kp']
              )
              + (eps2_pi) *(
                + eps2_pi *p['A_pp']
              )
        )

        for key in self.order['exclude']:
            del(p[key])

        return output *(eps2_k - eps2_pi)

    def fitfcn_nnlo_log_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        alpha_s = p['alpha_s']

        output = (
            alpha_s *eps2_a *(eps2_k - eps2_pi) *p['A_loga']
        )
        return output

    def fitfcn_nnlo_log2_ct(self, p):
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
            (1.0 / 96.0) *(eps2_k - eps2_pi) *(17 *eps2_k + 37 *eps2_pi) *(np.log(np.sqrt(eps2_pi *eps2_k)))**2
        )
        return output

    def fitfcn_ma_ratio(self, p):
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

            - fcn_I_m(mju, L, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(mru, L, mu, order_vol) / F2

            + 4 *eps2_pi *(4 *np.pi)**2 *(p['L_4'] + p['L_5'])
            + 8 *eps2_k *(4 *np.pi)**2 *p['L_4']
        )

        FK_nlo_per_F0 = (
            + 1

            - (1.0/2.0) *fcn_I_m(mju, L, mu, order_vol) / F2
            + (1.0/8.0) *fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mru, L, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(msj, L, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mrs, L, mu, order_vol) / F2
            + (1.0/4.0) *fcn_I_m(mss, L, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(mx, L, mu, order_vol) / F2

            + 4 *eps2_pi *(4 *np.pi)**2 *p['L_4']
            + 4 *eps2_k *(4 *np.pi)**2 *(2 *p['L_4'] + p['L_5'])

            + del2_ju *(
                - fcn_dI_m(mpi, L, mu, order_vol) / (8 *F2)
                + fcn_K_mM((mpi, mx), L, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *(
                + fcn_K21_mM((mpi, mx), L, mu, order_vol) / (24 *F2)
            )
            + del2_ju *del2_rs *(
                - fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2) # May need to change sign
                + fcn_K21_mM((mss, mx), L, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + fcn_K_mM((mss, mx), L, mu, order_vol) / (4 *F2)
                - (mk)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
                + (mpi)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
            )
        )

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_ma(self, p):
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

            + (1.0/2.0) *fcn_I_m(mju, L, mu, order_vol) / F2
            + (1.0/8.0) *fcn_I_m(mpi, L, mu, order_vol) / F2
            + (1.0/4.0) *fcn_I_m(mru, L, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(msj, L, mu, order_vol) / F2
            + (1.0/4.0) *fcn_I_m(mss, L, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mrs, L, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(mx, L, mu, order_vol) / F2

            + del2_ju *(
                - fcn_dI_m(mpi, L, mu, order_vol) / (8 *F2)
                + fcn_K_mM((mpi, mx), L, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *(
                + fcn_K21_mM((mpi, mx), L, mu, order_vol) / (24 *F2)
            )
            + del2_ju *del2_rs *(
                - fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2)
                + fcn_K21_mM((mss, mx), L, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + fcn_K_mM((mss, mx), L, mu, order_vol) / (4 *F2)
                - (mk)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
                + (mpi)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
            )

            + 4 *(4 *np.pi)**2 *(mk**2 - mpi**2) / lam2_chi *p['L_5']
        )

        return output


    def fitfcn_xpt_ratio(self, p):
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
            - fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(mk, L, mu, order_vol) / F2
            + 4 *eps2_pi *(4 *np.pi)**2 *(p['L_4'] + p['L_5'])
            + 8 *eps2_k *(4 *np.pi)**2 *p['L_4']
        )

        FK_nlo_per_F0 = (1
            - (3.0/8.0) *fcn_I_m(mpi, L, mu, order_vol) / F2
            - (3.0/4.0) *fcn_I_m(mk, L, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(meta, L, mu, order_vol) / F2
            + 4 *eps2_pi *(4 *np.pi)**2 *p['L_4']
            + 4 *eps2_k *(4 *np.pi)**2 *(2 *p['L_4'] + p['L_5'])
        )

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_xpt(self, p):

        # Constants
        order_vol = self.order['vol']

        # Independent variables
        mpi = p['mpi']
        mk = p['mk']
        meta = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2))

        lam2_chi = p['lam2_chi']
        eps2_pi = mpi**2 / lam2_chi
        eps2_k = mk**2 / lam2_chi
        #eps2_eta = mx**2 / lam2_chi

        L = p['L']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        output = (
            1
            + 1*(
            + (5.0/8.0) *fcn_I_m(mpi, L, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mk, L, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(meta, L, mu, order_vol) / F2)
            + 4 *(eps2_k - eps2_pi) *(4 *np.pi)**2 *p['L_5']
        )

        return output


    # Taken from arxiv/1701.07559, eqn (17);
    # finite volume correction from notes
    def fitfcn_ma_old(self, p):

        # Constants
        pi = np.pi

        # Independent variables
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_ss = p['mss']**2 / lam2_chi
        eps2_ju = p['mju']**2 / lam2_chi
        eps2_ru = p['mru']**2 / lam2_chi
        eps2_sj = p['mjs']**2 / lam2_chi
        eps2_rs = p['mrs']**2 / lam2_chi
        del2_pq = p['a2DI'] / lam2_chi
        #del2_pq = p['a2DI'] / lam2_chi
        eps2_x = (4.0/3.0) *eps2_k - (1.0/3.0) *eps2_pi + del2_pq
        L = p['L']

        # Log terms
        l_ju, l_pi, l_sj, l_ru, l_rs, l_ss, l_x = [np.log(eps2) for eps2
                                                     in [eps2_ju, eps2_pi, eps2_sj,
                                                         eps2_ru, eps2_rs, eps2_ss, eps2_x]]

        output = (1 +
                 + 1/2.0 *eps2_ju *l_ju
                 + 1/8.0 *l_pi *(
                    + eps2_pi
                    - del2_pq *(eps2_x + eps2_pi) / (eps2_x - eps2_pi)
                    + (del2_pq)**2 *eps2_x / (3.0 *(eps2_x - eps2_pi)**2)
                    - 4 *(del2_pq)**2 *eps2_pi / (3.0 *(eps2_x - eps2_pi) *(eps2_ss - eps2_pi))
                 )
                 - 1/2.0 *eps2_sj *l_sj
                 + 1/4.0 *eps2_ru *l_ru
                 - 1/4.0 *eps2_rs *l_rs
                 + 1/4.0 *l_ss *(
                    + eps2_ss
                    + del2_pq *(3 *(eps2_ss)**2 + 2 *(eps2_k - eps2_pi) *eps2_x - 3 *eps2_ss *eps2_x) / (3.0 *(eps2_x - eps2_ss)**2)
                    - (del2_pq)**2 *(2 *(eps2_ss)**2 - eps2_x *(eps2_ss + eps2_pi)) / (3.0 *(eps2_x - eps2_ss)**2 *(eps2_ss - eps2_pi))
                 )
                 - 3/8.0 *eps2_x *l_x *(
                    + 1
                    - 2 *del2_pq / (3.0 *(eps2_x -eps2_pi))
                    + del2_pq *(4 *(eps2_k - eps2_pi) + 6 *(eps2_ss - eps2_x)) / (9.0 *(eps2_x - eps2_ss)**2)
                    + (del2_pq)**2 / (9.0 *(eps2_x - eps2_pi)**2)
                    - 2 *(del2_pq)**2 *(2 *eps2_ss - eps2_pi - eps2_x) / (9.0 *(eps2_x - eps2_ss)**2 *(eps2_x - eps2_pi))
                 )
                  + (del2_pq *(eps2_k - eps2_pi)) / (6.0 *(eps2_x - eps2_ss))
                  + (del2_pq)**2 / (24.0 *(eps2_x - eps2_pi))
                  - (del2_pq)**2 / (12.0 *(eps2_x - eps2_ss))
                  - del2_pq / 8.0
                  + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['L_5']
        )

        # Next add finite volume corrections
        # Constants
        c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]
        order_vol = self.order['vol']

        # Variables
        L = L = p['L']
        lam2_chi = p['lam2_chi']
        to_eps2 = lambda x : x**2/lam2_chi

        mju = p['mju']
        mpi = p['mpi']
        mss = p['mss']
        mk = p['mk']
        mx = np.sqrt((4.0/3.0) *(mk**2) - (1.0/3.0) *(mpi**2) + p['a2DI'])

        del2_ju = p['a2DI']
        del2_rs = p['a2DI']

        eps2_ju = to_eps2(mju)
        eps2_pi = to_eps2(mpi)
        eps2_ss = to_eps2(mss)
        eps2_x = to_eps2(mx)

        eps2_D_ju = (del2_ju/lam2_chi)
        eps2_D_rs = (del2_rs/lam2_chi)

        for n in range(1, np.min((order_vol, 11))):
            xju = mju *L *n
            xpi = mpi *L *n

            output = output + c[n]*(
                + 2 *eps2_ju *fcn_Kn(1, xju) / xju
                + (1.0/2.0) *eps2_pi *fcn_Kn(1, xpi) / xpi
                - (1.0/8.0) *eps2_D_ju *(
                    + 2 *fcn_Kn(1, xpi) / xpi
                    - fcn_Kn(0, xpi)
                    - fcn_Kn(2, xpi)
                )
                - eps2_D_ju *eps2_pi / (eps2_x - eps2_pi) *fcn_Kn(1, xpi) / xpi
                + (1/24) *(eps2_D_ju)**2 *(
                    + 4 *eps2_pi / (eps2_x - eps2_pi)**2 *fcn_Kn(1, xpi) / xpi
                    + 1 / (eps2_x - eps2_pi) *(
                        + 2 *fcn_Kn(1, xpi) / xpi
                        - fcn_Kn(0, xpi)
                        - fcn_Kn(2, xpi)
                    )
                )
                + (2.0/3.0) *eps2_D_ju *eps2_D_rs *eps2_pi / ((eps2_pi - eps2_x) *(eps2_pi - eps2_ss)) *fcn_Kn(1, xpi) / xpi
            )

        return output

    def buildprior(self, prior, mopt=None, extend=False):
        return prior

    def builddata(self, data):
        return data
