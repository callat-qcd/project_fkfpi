import lsqfit
import numpy as np
import gvar as gv


class fitter(object):

    def __init__(self, fit_data=None, prior=None, order=None):

        # Renormalization momentum
        if order is None:
            order = 0
        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.order = order

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
        models = np.append(models, fk_fpi_model(datatag='fk_fpi', fit_data=fit_data, order=self.order))
        return models



class fk_fpi_model(lsqfit.MultiFitterModel):

    def __init__(self, datatag, fit_data, order=None):
        super(fk_fpi_model, self).__init__(datatag)

        if order is None:
            order = 0

        # Fit data
        self.fit_data = fit_data
        self.order = order

    def fitfcn(self, p, fit_data=None):
        if fit_data is None:
            fit_data = self.fit_data

        # Constants
        pi = np.pi
        order = self.order

        # Independent variables
        lam2_chi = 16 *pi**2 *(fit_data['Fpi'] *fit_data['FK'])
        eps2_pi = fit_data['mpi']**2 / lam2_chi
        eps2_k = fit_data['mk']**2 / lam2_chi
        eps2_ss = fit_data['mss']**2 / lam2_chi
        eps2_ju = fit_data['mju']**2 / lam2_chi
        eps2_ru = fit_data['mru']**2 / lam2_chi
        eps2_sj = eps2_ru
        eps2_rs = fit_data['mrs']**2 / lam2_chi
        del2_pq = fit_data['a2DI'] / lam2_chi
        #del2_pq = p['a2DI'] / lam2_chi
        eps2_x = (4.0/3.0) *eps2_k - (1.0/3.0) *eps2_pi + del2_pq
        aw02 = fit_data['aw0']**2
        mpil = fit_data['MpiL']

        # Force output array to have the correct shape
        output = 0 *eps2_pi

        # Lattice artifact term
        output = output + aw02 *p['l_a2']
        #output = output + aw02 *np.sqrt(aw02) *p['l_vol']

        # Volume term
        output = output + np.exp(-mpil) / np.sqrt(mpil) *p['l_vol']

        # Order = 0
        output = (output + 1
                 + (del2_pq *(eps2_k - eps2_pi)) / (6.0 *(eps2_x - eps2_ss))
                 + (del2_pq)**2 / (24.0 *(eps2_x - eps2_pi))
                 - (del2_pq)**2 / (12.0 *(eps2_x - eps2_ss))
                 - del2_pq / 8.0
                 + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['l_5lam']
                 )
        if order == 0:
            return output

        # Order = 1
        output = (output +
                 + 1/2.0 *eps2_ju *p['l_ju']
                 + 1/8.0 *p['l_pi'] *(
                      eps2_pi
                    - del2_pq *(eps2_x + eps2_pi) / (eps2_x - eps2_pi)
                    + (del2_pq)**2 *eps2_x / (3.0 *(eps2_x - eps2_pi)**2)
                    - 4 *(del2_pq)**2 *eps2_pi / (3.0 *(eps2_x - eps2_pi) *(eps2_ss - eps2_pi))
                 )
                 - 1/2.0 *eps2_sj *p['l_sj']
                 + 1/4.0 *eps2_ru *p['l_ru']
                 - 1/4.0 *eps2_rs *p['l_rs']
                 + 1/4.0 *p['l_ss'] *(
                      eps2_ss
                    + del2_pq *(3 *(eps2_ss)**2 +2 *(eps2_k -eps2_pi) *eps2_x - 3 *eps2_ss *eps2_x) / (3.0 *(eps2_x - eps2_ss)**2)
                    - (del2_pq)**2 *(2 *(eps2_ss)**2 - eps2_x *(eps2_ss + eps2_pi)) / (3.0 *(eps2_x - eps2_ss)**2 *(eps2_ss - eps2_pi))
                 )
                 - 3/8.0 *eps2_x *p['l_x'] *(
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

    def buildprior(self, prior, mopt=None, extend=False):
        newprior = gv.BufferDict()
        order = self.order

        #newprior['a2DI'] = prior['a2DI']

        # Sort keys by expansion order
        keys_0 = ['l_5lam']
        keys_1 = ['l_ju', 'l_pi', 'l_rs', 'l_ru', 'l_sj', 'l_ss', 'l_x']
        keys_2 = []
        keys_lat = ['l_a2']#, 'c_3_a']
        keys_vol = ['l_vol']

        for key in keys_0:
            newprior[key] = prior[key]

        if order >= 1:
            for key in keys_1:
                newprior[key] = prior[key]

        if order >= 2:
            for key in keys_2:
                newprior[key] = prior[key]

        # Lattice artifacts
        for key in keys_lat:
            newprior[key] = prior[key]

        # Volume term
        for key in keys_vol:
           newprior[key] = prior[key]

        return newprior

    def builddata(self, data):
        return data['y']
