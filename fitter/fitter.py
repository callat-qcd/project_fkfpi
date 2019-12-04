import lsqfit
import numpy as np
import gvar as gv
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import special_functions as sf

class fitter(object):

    def __init__(self, order, fit_type, F2, fit_data=None, prior=None, chain_fits=True):
        self.prior = prior
        self.fit_data = fit_data
        self.fit = None
        self.empbayes_fit = None
        self.order = order
        self.fit_type = fit_type
        self.chain_fits = chain_fits
        self.F2 = F2

    def _make_fitargs(self, z):
        y_data = self._make_y_data()
        prior = self._make_prior()

        #z['chiral'] = np.around(z['chiral'], 2)
        #z['spacing'] = np.around(z['spacing'], 2)
        z['chiral'] = np.abs(z['chiral'])
        z['spacing'] = np.abs(z['spacing'])

        for key in prior.keys():
            if key in ['A_p', 'A_k']:
                #prior[key] = prior[key] *z['chiral']
                prior[key] = gv.gvar(0, z['chiral'])
            if key in ['A_loga', 'A_a']:
                #prior[key] = prior[key] *z['spacing']
                prior[key] = gv.gvar(0, z['spacing'])

        fitfcn = self._make_models()[-1].fitfcn

        print(z)
        print(prior['A_k'], prior['A_a'])

        return dict(data=y_data, fcn=fitfcn, prior=prior)

    def _make_empbayes_fit(self):
        models = self._make_models()
        y_data = self._make_y_data()
        prior = self._make_prior()

        z0 = gv.BufferDict()
        z0['chiral'] = 1.0
        z0['spacing'] = 1.0

        fit, z = lsqfit.empbayes_fit(z0, self._make_fitargs)

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
                            F2=self.F2, order=self.order, fit_type=fit_type))
            return models


        order = self.order.copy()
        models = np.append(models, fk_fpi_model(datatag=self.fit_type+'_'+order['fit'],
                    F2=self.F2, order=order, fit_type=self.fit_type))

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


        newprior['L_1'] = prior['L_1']
        newprior['L_2'] = prior['L_2']
        newprior['L_3'] = prior['L_3']
        newprior['L_4'] = prior['L_4']
        newprior['L_5'] = prior['L_5']
        newprior['L_6'] = prior['L_6']
        newprior['L_7'] = prior['L_7']
        newprior['L_8'] = prior['L_8']


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

    def __init__(self, datatag, fit_type, order, F2):
        super(fk_fpi_model, self).__init__(datatag)

        self.fit_type = fit_type
        self.order = order
        self.F2 = F2


    def fitfcn(self, p, fit_data=None):
        if fit_data is not None:
            for key in fit_data.keys():
                p[key] = fit_data[key]

        for key in self.order['exclude']:
            p[key] = 0

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

        if self.order['fit'] in ['nnlo', 'nnnlo']:
            output = output + self.fitfcn_nnlo_cts(p)

        if self.order['fit'] in ['nnnlo']:
            output = output + self.fitfcn_nnnlo_cts(p)

        for key in self.order['exclude']:
            del(p[key])

        return output


    def fitfcn_nnlo_cts(self, p):
        lam2_chi = p['lam2_chi']
        eps2_a = (p['a/w0'])**2 / (4 *np.pi)
        eps2_pi = p['mpi']**2 / lam2_chi
        eps2_k = p['mk']**2 / lam2_chi
        eps2_eta = (4./3) *eps2_k  - (1./3) *eps2_pi

        fcn_l = lambda x : x *np.log(x)

        if self.F2 == 'FKFK':
            output = (3./2) *(
                + (3./8) * fcn_l(eps2_pi)
                + (3./4) * fcn_l(eps2_k)
                + (3./8) * fcn_l(eps2_eta)

                - 4 *(4 *np.pi)**2 *(
                    + eps2_pi *p['L_4']
                    + eps2_k *(2 *p['L_4'] + p['L_5'])
                )
            )
        elif self.F2 == 'FKFpi':
            output = 0
        elif self.F2 == 'FpiFpi':
            output = 0

        return output *(eps2_k - eps2_pi)

    def fitfcn_nnnlo_cts(self, p):
        return None

    def fitfcn_ma_ratio(self, p):
        return None

    def fitfcn_ma(self, p):
        return None

    def fitfcn_xpt_ratio(self, p):
        return None

    def fitfcn_xpt(self, p):

        # Constants
        #order_vol = self.order['vol']

        # Independent variables
        mpi = p['mpi']
        mk = p['mk']
        meta = np.sqrt((4./3) *(mk**2) - (1./3) *(mpi**2))

        lam2_chi = p['lam2_chi']
        eps2_pi = mpi**2 / lam2_chi
        eps2_k = mk**2 / lam2_chi
        eps2_eta = meta**2 / lam2_chi

        #L = p['L']
        #mu = np.sqrt(lam2_chi)
        #F2 = lam2_chi /(4*np.pi)**2

        fcn_l = lambda x : x *np.log(x)

        output = (
            + 1

            + (5./8) *fcn_l(eps2_pi)
            - (1./4) *fcn_l(eps2_k)
            - (3./8) *fcn_l(eps2_eta)

            + 4 *(eps2_k - eps2_pi) *(4 *np.pi)**2 *p['L_5']

            + eps2_k**2 *sf.fcn_FF(eps2_pi/eps2_k)

            + sf.fcn_Kr_j(1, mpi, mk, lam2_chi) *np.log(eps2_pi) *np.log(eps2_pi)
            + sf.fcn_Kr_j(2, mpi, mk, lam2_chi) *np.log(eps2_pi) *np.log(eps2_k)
            + sf.fcn_Kr_j(3, mpi, mk, lam2_chi) *np.log(eps2_pi) *np.log(eps2_eta)
            + sf.fcn_Kr_j(4, mpi, mk, lam2_chi) *np.log(eps2_k) *np.log(eps2_k)
            + sf.fcn_Kr_j(5, mpi, mk, lam2_chi) *np.log(eps2_k) *np.log(eps2_eta)
            + sf.fcn_Kr_j(6, mpi, mk, lam2_chi) *np.log(eps2_eta) *np.log(eps2_eta)

            + sf.fcn_Cr_j(1, mpi, mk, lam2_chi, p) *np.log(eps2_pi)
            + sf.fcn_Cr_j(2, mpi, mk, lam2_chi, p) *np.log(eps2_k)
            + sf.fcn_Cr_j(3, mpi, mk, lam2_chi, p) *np.log(eps2_eta)
        )

        return output

    def buildprior(self, prior, mopt=None, extend=False):
        return prior

    def builddata(self, data):
        return data
