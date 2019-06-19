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
            fit_type = 'mix'

        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.order = order
        self.fit_type = fit_type

    def get_fit(self):
        if self.fit is None:
            return self._make_fit()
        else:
            return self.fit

    def _make_fit(self):
        models = self._make_models()
        fitter = lsqfit.MultiFitter(models=models)
        fit = fitter.lsqfit(data=self.fit_data, prior=self.prior)
        self.fit = fit
        return fit

    def _make_models(self, fit_data=None):
        if fit_data is None:
            fit_data=self.fit_data

        models = np.array([])
        models = np.append(models, fk_fpi_model(datatag='fk_fpi', fit_data=fit_data,
                    order=self.order, fit_type=self.fit_type))
        return models



class fk_fpi_model(lsqfit.MultiFitterModel):

    def __init__(self, datatag, fit_data, fit_type=None, order=None):
        super(fk_fpi_model, self).__init__(datatag)

        if order is None:
            order = {
                'fit' : 'nlo',
                'latt_spacing' : 2, # no order 1 term -- starts at 2
                'vol' : 1
            }

        if fit_type is None:
            fit_type = 'mix'

        # Fit data
        self.fit_data = fit_data
        self.order = order
        self.fit_type = fit_type

    def fitfcn(self, p, fit_data=None):
        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        #print self.fit_type
        # Lattice artifact terms
        output = (self.fitfcn_latt_spacing_corections(p)
                 + self.fitfcn_mpia_corrections(p))

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

        if self.order['fit'] in ['nnlo', 'nnnlo']:
            output = output + self.fitfcn_nnlo_cts(p)

        if self.order['fit'] in ['nnnlo']:
            output = output + self.fitfcn_nnnlo_cts(p)

        return output

    def fitfcn_mpia_corrections(self, p):
        a = p['a']
        mpi = p['mpi']
        lam2_chi = p['lam2_chi']
        output = (a *mpi)**2 /lam2_chi *p['c_mpia2']
        return output

    def fitfcn_latt_spacing_corections(self, p):
        order = self.order['latt_spacing']
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

        return output

    def fitfcn_nnlo_cts(self, p):
        # Will need to adjust this later!
        w0 = 5.81743

        lam2_chi = p['lam2_chi']
        eps2_a = p['a'] / (w0 *np.sqrt(4 *np.pi))
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


        # Will need to adjust this later!
        w0 = 5.81743

        lam2_chi = p['lam2_chi']
        eps2_a = p['a'] / (w0 *np.sqrt(4 *np.pi))
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
        pi = np.pi
        order = self.order['fit']

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
        mpil = p['MpiL']

        # Log terms
        l_ju, l_pi, l_sj, l_ru, l_rs, l_ss, l_x = [np.log(eps2) for eps2
                                                     in [eps2_ju, eps2_pi, eps2_sj,
                                                         eps2_ru, eps2_rs, eps2_ss, eps2_x]]

        # Force output array to have the correct shape
        output = 0 *eps2_pi



        # Order = 0
        output = (output + 1
                 + (del2_pq *(eps2_k - eps2_pi)) / (6.0 *(eps2_x - eps2_ss))
                 + (del2_pq)**2 / (24.0 *(eps2_x - eps2_pi))
                 - (del2_pq)**2 / (12.0 *(eps2_x - eps2_ss))
                 - del2_pq / 8.0
                 + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['L_5']
                 )
        if order == 0:
            return output

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

        if order == 1:
            return output

        return output

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
        MpiL = p['MpiL']

        lam2_chi = p['lam2_chi']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        output = (
            1

            + (1.0/2.0) *fcn_I_m(mju, MpiL, mu, order_vol) / F2
            + (1.0/8.0) *fcn_I_m(mpi, MpiL, mu, order_vol) / F2
            + (1.0/4.0) *fcn_I_m(mru, MpiL, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(msj, MpiL, mu, order_vol) / F2
            + (1.0/4.0) *fcn_I_m(mss, MpiL, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mrs, MpiL, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(mx, MpiL, mu, order_vol) / F2

            + del2_ju *(
                - fcn_dI_m(mpi, MpiL, mu, order_vol) / (8 *F2)
                + fcn_K_mM((mpi, mx), MpiL, mu, order_vol) / (4 *F2)
            )
            - (del2_ju)**2 *fcn_K21_mM((mpi, mx), MpiL, mu, order_vol) / (24 *F2)
            + del2_ju *del2_rs *(
                + fcn_K_m1m2m3((mpi, mss, mx), MpiL, mu, order_vol) / (6 *F2)
                + fcn_K21_mM((mss, mx), MpiL, mu, order_vol) / (12 *F2)
            )
            + del2_rs *(
                + fcn_K_mM((mss, mx), MpiL, mu, order_vol) / (4 *F2)
                - (mk)**2 *fcn_K21_mM((mss, mx), MpiL, mu, order_vol) / (6 *F2)
                + (mpi)**2 *fcn_K21_mM((mss, mx), MpiL, mu, order_vol) / (6 *F2)
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

        MpiL = p['MpiL']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        Fpi_nlo_per_F0 = (1
            - fcn_I_m(mpi, MpiL, mu, order_vol) / F2
            - (1.0/2.0) *fcn_I_m(mk, MpiL, mu, order_vol) / F2
            + 4 *eps2_pi *(4 *np.pi)**2 *(p['L_4'] + p['L_5'])
            + 8 *eps2_k *(4 *np.pi)**2 *p['L_4']
        )

        FK_nlo_per_F0 = (1
            - (3.0/8.0) *fcn_I_m(mpi, MpiL, mu, order_vol) / F2
            - (3.0/4.0) *fcn_I_m(mk, MpiL, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(meta, MpiL, mu, order_vol) / F2
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

        MpiL = p['MpiL']
        mu = np.sqrt(lam2_chi)
        F2 = lam2_chi /(4*np.pi)**2

        output = (
            1
            + 1*(
            + (5.0/8.0) *fcn_I_m(mpi, MpiL, mu, order_vol) / F2
            - (1.0/4.0) *fcn_I_m(mk, MpiL, mu, order_vol) / F2
            - (3.0/8.0) *fcn_I_m(meta, MpiL, mu, order_vol) / F2)
            + 4 *(eps2_k - eps2_pi) *(4 *np.pi)**2 *p['L_5']
        )
        return output



    def buildprior(self, prior, mopt=None, extend=False):
        newprior = gv.BufferDict()
        order = self.order

        # Move fit_data into prior
        fit_data = self.fit_data
        newprior['Fpi'] = fit_data['Fpi']
        newprior['FK'] = fit_data['FK']
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
        newprior['MpiL'] = fit_data['MpiL']

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

    def builddata(self, data):
        return data['y']
