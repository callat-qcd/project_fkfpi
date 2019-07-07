import lsqfit
import numpy as np
import gvar as gv

# import fcn_T_m, fcn_dI_m, etc
from special_functions import *

class fitter(object):

    def __init__(self, fit_data=None, prior=None, fit_type=None, order=None):

        # Renormalization momentum
        if order is None:
            order = {
                'fit' : 'nlo',
                'latt_spacing' : 2, # no order 1 term -- starts at 2
                'vol' : 1
            }

        if fit_type is None:
            fit_type = 'ma-taylor'

        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.order = order
        self.fit_type = fit_type

    def _make_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()

        fitter = lsqfit.MultiFitter(models=models)
        fit = fitter.lsqfit(data=y_data, prior=prior, fast=False)
        #fit = fitter.chained_lsqfit(data=y_data, prior=prior)


        self.fit = fit
        return fit

    def _make_models(self):
        models = np.array([])
        if self.fit_type == 'simultaneous':
            for fit_type in ['xpt', 'xpt-taylor']:
                models = np.append(models, fk_fpi_model(datatag=fit_type,
                            order=self.order, fit_type=fit_type))
        else:
            models = np.append(models, fk_fpi_model(datatag=self.fit_type,
                        order=self.order, fit_type=self.fit_type))
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
        newprior['a'] = fit_data['a']
        newprior['L'] = fit_data['L']
        newprior['w0'] = fit_data['w0']

        #for key in newprior.keys():
        #    newprior[key] = gv.gvar(gv.mean(newprior[key]), gv.sdev(newprior[key])/100)

        # Fit parameters, depending on fit type
        if self.fit_type == 'xpt':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'xpt-taylor':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma-taylor':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma-old':
            newprior['L_5'] = prior['L_5']
        elif self.fit_type == 'simultaneous':
            newprior['L_4'] = prior['L_4']
            newprior['L_5'] = prior['L_5']


        # Fit parameters, depending on order
        if order['fit'] in ['nnlo', 'nnnlo']:
            newprior['A_22'] = prior['A_22']
            newprior['M_400'] = prior['M_400']
            newprior['M_220'] = prior['M_220']
            newprior['M_202'] = prior['M_202']

        if order['fit'] in ['nnnlo']:
            newprior['A_42'] = prior['A_42']
            newprior['A_24'] = prior['A_24']
            newprior['A_2220'] = prior['A_2220']
            newprior['A_2202'] = prior['A_2202']


        # Lattice artifacts
        if order['latt_spacing'] >= 2:
            newprior['c_a2'] = prior['c_a2']

        if order['latt_spacing'] >= 3:
            newprior['c_a3'] = prior['c_a3']

        if order['latt_spacing'] >= 4:
            newprior['c_a4'] = prior['c_a4']

        newprior['c_mpia2'] = prior['c_mpia2']

        # Fudge factor
        #newprior['c_fudge'] = prior['c_fudge']

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


class fk_fpi_model(lsqfit.MultiFitterModel):

    def __init__(self, datatag, fit_type, order):
        super(fk_fpi_model, self).__init__(datatag)

        #if order is None:
        #    order = {
        #        'fit' : 'nlo',
        #        'latt_spacing' : 2, # no order 1 term -- starts at 2
        #        'vol' : 1
        #    }

        #if fit_type is None:
        #    fit_type = 'ma-taylor'

        # Fit data
        self.order = order
        self.fit_type = fit_type

    def fitfcn(self, p, fit_data=None):
        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        #print self.fit_type
        # Lattice artifact terms
        output = (self.fitfcn_latt_spacing_corrections(p)
                  #+ self.fitfcn_finite_vol_corrections(p) # Don't need this -- already in fcn_I_m definitions
                 + self.fitfcn_mpia_corrections(p)) # Doesn't seem to be doing anything

        if self.order['fit'] in ['nlo', 'nnlo', 'nnnlo']:
            # mixed-action/xpt fits
            if self.fit_type == 'ma':
                output = output + self.fitfcn_ma(p)
            elif self.fit_type == 'ma-taylor':
                output = output + self.fitfcn_ma_taylor(p)
            elif self.fit_type == 'xpt':
                output = output + self.fitfcn_xpt(p)
            elif self.fit_type == 'xpt-taylor':
                output = output + self.fitfcn_xpt_taylor(p)
            elif self.fit_type == 'ma-old':
                output = output + self.fitfcn_ma_old(p)

        if self.order['fit'] in ['nnlo', 'nnnlo']:
            output = output + self.fitfcn_nnlo_cts(p)

        if self.order['fit'] in ['nnnlo']:
            output = output + self.fitfcn_nnnlo_cts(p)

        return output

    def fitfcn_mpia_corrections(self, p):
        a = p['a']
        mpi = p['mpi']
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi


        output = (eps2_k - eps2_pi) *(a *mpi)**2 /lam2_chi *p['c_mpia2']
        return output

    def fitfcn_latt_spacing_corrections(self, p):
        order = self.order['latt_spacing']
        lam2_chi = p['lam2_chi']
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        a = p['a']
        output = 0

        if order < 2:
            return output

        if order >= 2:
            output = output + a**2 *p['c_a2']
        if order >= 3:
            output = output + a**3 *p['c_a3']
        if order >= 4:
            output = output + a**4 *p['c_a4']

        return output *(eps2_k - eps2_pi)


    # This function shouldn't be used since the finite volume corrections
    # are already embedded in the models
    def fitfcn_finite_vol_corrections(self, p):

        # Constants
        c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]
        order_vol = self.order['vol']
        order_fit = self.order['fit']

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

        output = 0
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

    def fitfcn_nnlo_cts(self, p):
        w0 = p['w0']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a'] / (w0 *np.sqrt(4 *np.pi)))**2
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
              eps2_a *(eps2_k - eps2_pi) *p['A_22']
            + (eps2_k - eps2_pi)**2 *p['M_400']
            + eps2_k *(eps2_k - eps2_pi) *p['M_220']
            + eps2_pi *(eps2_k - eps2_pi) *p['M_202']
        )
        return output

    def fitfcn_nnnlo_cts(self, p):

        w0 = p['w0']

        lam2_chi = p['lam2_chi']
        eps2_a = (p['a'] / (w0 *np.sqrt(4 *np.pi)))**2
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi

        output = (
              (eps2_a)**2 *(eps2_k - eps2_pi) *p['A_42']
            + eps2_a *(eps2_k - eps2_pi)**2 *p['A_24']
            + eps2_a *eps2_k *(eps2_k - eps2_pi) *p['A_2220']
            + eps2_a *eps2_pi *(eps2_k - eps2_pi) *p['A_2202']
        )
        return output




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
            + (1.0/2.0) *fcn_I_m(mss, L, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(mx, L, mu, order_vol) / F2

            + 4 *eps2_pi *p['L_4']
            + 4 *eps2_k *(2 *p['L_4'] + p['L_5'])

            + del2_ju *(
                - fcn_dI_m(mpi, L, mu, order_vol) / (8 *F2)
                + fcn_K_mM((mpi, mx), L, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *(
                + fcn_K21_mM((mpi, mx), L, mu, order_vol) / (24 *F2)
            )
            + del2_ju *del2_rs *(
                + fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2)
                + fcn_K21_mM((mss, mx), L, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + fcn_K_mM((mss, mx), L, mu, order_vol) / (4 *F2)
                - (mk)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
                + (mpi)**2 *fcn_K21_mM((mss, mx), L, mu, order_vol) / (6 *F2)
            )
        )

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_ma_taylor(self, p):
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
                + fcn_K_m1m2m3((mpi, mss, mx), L, mu, order_vol) / (6 *F2)
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
            + 4 *eps2_pi *p['L_4']
            + 4 *eps2_k *(2 *p['L_4'] + p['L_5'])
        )

        return FK_nlo_per_F0 / Fpi_nlo_per_F0


    def fitfcn_xpt_taylor(self, p):
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
        eps2_sj = eps2_ru
        eps2_rs = p['mrs']**2 / lam2_chi
        del2_pq = p['a2DI'] / lam2_chi
        #del2_pq = p['a2DI'] / lam2_chi
        eps2_x = (4.0/3.0) *eps2_k - (1.0/3.0) *eps2_pi + del2_pq
        L = p['L']

        # Log terms
        l_ju, l_pi, l_sj, l_ru, l_rs, l_ss, l_x = [np.log(eps2) for eps2
                                                     in [eps2_ju, eps2_pi, eps2_sj,
                                                         eps2_ru, eps2_rs, eps2_ss, eps2_x]]

        # Order = 0
        output = (1
                 + (del2_pq *(eps2_k - eps2_pi)) / (6.0 *(eps2_x - eps2_ss))
                 + (del2_pq)**2 / (24.0 *(eps2_x - eps2_pi))
                 - (del2_pq)**2 / (12.0 *(eps2_x - eps2_ss))
                 - del2_pq / 8.0
                 + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['L_5']
                 )

        # Order = 1
        output = (output +
                 + 1/2.0 *eps2_ju *l_ju
                 + 1/8.0 *l_pi *(
                      eps2_pi
                    - del2_pq *(eps2_x + eps2_pi) / (eps2_x - eps2_pi)
                    + (del2_pq)**2 *eps2_x / (3.0 *(eps2_x - eps2_pi)**2)
                    - 4 *(del2_pq)**2 *eps2_pi / (3.0 *(eps2_x - eps2_pi) *(eps2_ss - eps2_pi))
                 )
                 - 1/2.0 *eps2_sj *l_sj
                 + 1/4.0 *eps2_ru *l_ru
                 - 1/4.0 *eps2_rs *l_rs
                 + 1/4.0 *l_ss *(
                      eps2_ss
                    + del2_pq *(3 *(eps2_ss)**2 +2 *(eps2_k -eps2_pi) *eps2_x - 3 *eps2_ss *eps2_x) / (3.0 *(eps2_x - eps2_ss)**2)
                    - (del2_pq)**2 *(2 *(eps2_ss)**2 - eps2_x *(eps2_ss + eps2_pi)) / (3.0 *(eps2_x - eps2_ss)**2 *(eps2_ss - eps2_pi))
                 )
                 - 3/8.0 *eps2_x *l_x *(
                      1
                    - 2 *del2_pq / (3.0 *(eps2_x -eps2_pi))
                    + del2_pq *(4 *(eps2_k - eps2_pi) + 6 *(eps2_ss - eps2_x)) / (9.0 *(eps2_x - eps2_ss)**2)
                    + (del2_pq)**2 / (9.0 *(eps2_x - eps2_pi)**2)
                    - 2 *(del2_pq)**2 *(2 *eps2_ss - eps2_pi - eps2_x) / (9.0 *(eps2_x - eps2_ss)**2 *(eps2_x - eps2_pi))
                 )
        )

        return output

    def buildprior(self, prior, mopt=None, extend=False):
        order = self.order

        mprior = gv.BufferDict()
        mprior['mpi'] = prior['mpi']
        mprior['mk'] = prior['mk']
        mprior['lam2_chi'] = prior['lam2_chi']
        mprior['a'] = prior['a']
        mprior['L'] = prior['L']
        mprior['w0'] = prior['w0']

        # Fit parameters, depending on fit type
        if self.fit_type == 'xpt':
            mprior['L_4'] = prior['L_4']
            mprior['L_5'] = prior['L_5']
        elif self.fit_type == 'xpt-taylor':
            mprior['L_5'] = prior['L_5']
        elif self.fit_type == 'ma':
            mprior['L_4'] = prior['L_4']
            mprior['L_5'] = prior['L_5']
            mprior['a2DI'] = prior['a2DI']
            mprior['mss'] = prior['mss']
            mprior['mju'] = prior['mju']
            mprior['mjs'] = prior['mjs']
            mprior['mru'] = prior['mru']
            mprior['mrs'] = prior['mrs']
        elif self.fit_type == 'ma-taylor':
            mprior['L_5'] = prior['L_5']
            mprior['a2DI'] = prior['a2DI']
            mprior['mss'] = prior['mss']
            mprior['mju'] = prior['mju']
            mprior['mjs'] = prior['mjs']
            mprior['mru'] = prior['mru']
            mprior['mrs'] = prior['mrs']
        elif self.fit_type == 'ma-old':
            mprior['L_5'] = prior['L_5']


        # Fit parameters, depending on order
        if order['fit'] in ['nnlo', 'nnnlo']:
            mprior['A_22'] = prior['A_22']
            mprior['M_400'] = prior['M_400']
            mprior['M_220'] = prior['M_220']
            mprior['M_202'] = prior['M_202']

        if order['fit'] in ['nnnlo']:
            mprior['A_42'] = prior['A_42']
            mprior['A_24'] = prior['A_24']
            mprior['A_2220'] = prior['A_2220']
            mprior['A_2202'] = prior['A_2202']


        # Lattice artifacts
        if order['latt_spacing'] >= 2:
            mprior['c_a2'] = prior['c_a2']

        if order['latt_spacing'] >= 3:
            mprior['c_a3'] = prior['c_a3']

        if order['latt_spacing'] >= 4:
            mprior['c_a4'] = prior['c_a4']

        mprior['c_mpia2'] = prior['c_mpia2']

        return mprior

    def builddata(self, data):
        return data
