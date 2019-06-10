import lsqfit
import numpy as np
import gvar as gv


class fitter(object):

    def __init__(self, fit_data=None, prior=None, fit_type=None,
                 order_fit=None, order_latt_spacing=None, order_vol=None):

        # Renormalization momentum
        if order_fit is None:
            order_fit = 0
        if order_latt_spacing is None:
            order_latt_spacing = 2
        if order_vol is None:
            order_vol = 1

        if fit_type is None:
            fit_type = 'mix'

        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.order_fit = order_fit
        self.order_latt_spacing = order_latt_spacing
        self.order_vol = order_vol
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
                    order_fit=self.order_fit, order_latt_spacing=self.order_latt_spacing,
                    order_vol=self.order_vol, fit_type=self.fit_type))
        return models



class fk_fpi_model(lsqfit.MultiFitterModel):

    def __init__(self, datatag, fit_data, fit_type=None,
                 order_fit=None, order_latt_spacing=None, order_vol=None):
        super(fk_fpi_model, self).__init__(datatag)

        if order_fit is None:
            order_fit = 0
        if order_latt_spacing is None:
            order_latt_spacing = 2
        if order_vol is None:
            order_vol = 1

        if fit_type is None:
            fit_type = 'mix'

        # Fit data
        self.fit_data = fit_data
        self.order_fit = order_fit
        self.order_latt_spacing = order_latt_spacing
        self.order_vol = order_vol
        self.fit_type = fit_type

    def fitfcn(self, p, fit_data=None):
        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        # Lattice artifact terms
        output = (self.fitfcn_volume_corrections(p)
                 + self.fitfcn_latt_spacing_corections(p)
                 + self.fitfcn_mpia_corrections(p))

        # mixed-action/xpt fits
        if self.fit_type == 'mix':
            output = output + self.fitfcn_mixed_action(p)
        elif self.fit_type == 'xpt':
            output = output + self.fitfcn_xpt(p)

        return output

    def fitfcn_mpia_corrections(self, p):
        a = p['a']
        mpi = p['mpi']
        lam2_chi = 16 *np.pi**2 *(p['Fpi'] *p['FK'])
        output = (a *mpi /lam2_chi)**2 *p['c_mpia2']
        return output

    # Look at arxiv/1001.4692, section 2.4
    def fitfcn_volume_corrections(self, p):
        mpil = p['MpiL']
        output = np.exp(-mpil) / np.sqrt(mpil) *p['c_vol']
        return output

    def fitfcn_latt_spacing_corections(self, p):
        order = self.order_latt_spacing
        a = p['a']
        output = 0

        if order < 2:
            return output

        for j in range(order):
            output = output + a**(j+2) *p['c_a'][j]

        return output


    def fitfcn_mixed_action(self, p):

        # Constants
        pi = np.pi
        order = self.order_fit

        # Independent variables
        lam2_chi = 16 *pi**2 *(p['Fpi'] *p['FK'])
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_x = (4.0/3.0) *eps2_k - (1.0/3.0) *eps2_pi

        l_pi, l_k, l_x = [np.log(eps2) for eps2 in [eps2_pi, eps2_k, eps2_x]]

        output = (1.0
                 + 5.0/8.0 *(eps2_pi) *l_pi
                 - 1.0/4.0 *(eps2_k) *l_k
                 - 3.0/8.0 *(eps2_x) *l_x
                 + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['L_5_lam']
                 )

        return output

    def fitfcn_xpt(self, p):

        # Constants
        pi = np.pi
        order = self.order_fit

        # Independent variables
        lam2_chi = 16 *pi**2 *(p['Fpi'] *p['FK'])
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
                 + 4 *(eps2_k - eps2_pi) *(4*pi)**2 *p['L_5_lam']
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

    def buildprior(self, prior, mopt=None, extend=False):
        newprior = gv.BufferDict()
        order_fit = self.order_fit
        order_latt_spacing = self.order_latt_spacing
        order_vol = self.order_vol

        if True:
            fit_data = self.fit_data
            newprior['Fpi'] = fit_data['Fpi']
            newprior['FK'] = fit_data['FK']
            newprior['mpi'] = fit_data['mpi']
            newprior['mk'] = fit_data['mk']
            newprior['mss'] = fit_data['mss']
            newprior['mju'] = fit_data['mju']
            newprior['mru'] = fit_data['mru']
            newprior['mrs'] = fit_data['mrs']
            newprior['a2DI'] = fit_data['a2DI']

            newprior['a'] = fit_data['a']
            newprior['MpiL'] = fit_data['MpiL']

        # Sort keys by expansion order
        keys_0 = ['L_5_lam']
        keys_1 = []
        keys_2 = []
        keys_lat = ['c_a']
        keys_vol = ['c_vol']
        keys_other = ['c_mpia2']

        for key in keys_0:
            newprior[key] = prior[key]

        if order_fit >= 1:
            for key in keys_1:
                newprior[key] = prior[key]

        if order_fit >= 2:
            for key in keys_2:
                newprior[key] = prior[key]

        # Lattice artifacts
        newprior['c_a'] = np.empty(order_latt_spacing)
        for j in range(order_latt_spacing):
            newprior['c_a'][j] = prior['c_a'][j]

        # Volume term
        newprior['c_vol'] = np.empty(order_vol)
        for j in range(order_vol):
           newprior['c_vol'][j] = prior['c_vol'][j]

         # Other term
        for key in keys_other:
           newprior[key] = prior[key]

        return newprior

    def builddata(self, data):
        return data['y']
