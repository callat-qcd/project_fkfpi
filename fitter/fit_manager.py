import lsqfit
import numpy as np
import gvar as gv
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
import re
from mpl_toolkits.mplot3d.axes3d import Axes3D
import os

#from .fitter import fitter
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import fitter.fitter as fit
import fitter.special_functions as sf
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#import fitter
#import special_functions as sf

class fit_manager(object):

    def __init__(self, fit_data, phys_point_data, prior=None, model_info=None,
                abbrs=None, bias_correct=True, plot_bs_N=100, fast_sunset=False, **kwargs):

        # Default value for abbrs
        if abbrs is None:
            abbrs = fit_data.keys()


        # Add default values to model_info dict
        if model_info is None:
            model_info = {}

        model_info.setdefault('fit_type', 'xpt')
        model_info.setdefault('order', 'nnlo')
        model_info.setdefault('F2', 'FpiFpi')

        model_info.setdefault('include_FV', True)
        model_info.setdefault('exclude', [])
        model_info.setdefault('use_bijnens_central_value', True)

        model_info.setdefault('include_alpha_s', False)
        model_info.setdefault('include_log', False)
        model_info.setdefault('include_log2', False)
        model_info.setdefault('include_sunset', False)
        model_info.setdefault('include_latt_n3lo', False)
        model_info.setdefault('include_latt_n4lo', False)


        # Get lam_chi^2=renorm_scale^2 depending on choice of F^2
        if model_info['F2'] == 'FKFpi':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_kpi']
        elif model_info['F2'] == 'FpiFpi':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_pipi']
        elif model_info['F2'] == 'FKFK':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_kk']
        elif model_info['F2'] == 'F0F0':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_00']

        # Default prior
        if prior is None:
            print("Using default prior.")

            # From Bijnens
            gamma = np.array([3./32, 3./16, 0, 1./8, 3./8, 11./144, 0, 5./48])
            L_mu0 = np.array([3.72, 4.93, -30.7, 0.89, 3.77, 0.11, -3.4, 2.94])*10**-4

            # F2 is F^2 in this context; but F0, F1 are F_0, F_1
            if model_info['F2'] == 'FpiFpi':
                F1 = phys_point_data['Fpi']
            elif model_info['F2'] == 'FKFpi':
                F1 = np.sqrt(phys_point_data['FK'] *phys_point_data['Fpi'])
            elif model_info['F2'] == 'FKFK':
                F1 = phys_point_data['Fpi']

            F0 = 80 # MeV
            L_mu1 = gv.mean(L_mu0 - gamma/(4 *np.pi)**2 *np.log(F1/F0))
            if model_info['use_bijnens_central_value']:
                L_mu1 = gv.gvar(L_mu1, L_mu1/2)
            else:
                L_mu1 = gv.gvar(np.repeat(0, len(L_mu1)), L_mu1)

            # Have uncertainty of at least 10^-4
            for j in range(len(L_mu1)):
                if gv.sdev(L_mu1[j]) < 10**-4:
                    L_mu1[j] = gv.gvar(gv.mean(L_mu1[j]), 10**-4)

            prior = {}
            # Gasser-Leutwyler LECs
            for j, L_i in enumerate(L_mu1):
                prior['L_'+str(j+1)] = L_i

            # nlo Gasser-Leutwyler LECs
            prior['L_4'] = gv.gvar('0.0(0.005)')
            prior['L_5'] = gv.gvar('0.0(0.005)')

            # nlo polynomial
            prior['A_x'] = gv.gvar('0.0(1.0)')

            # Lattice spacing terms
            prior['A_a'] = gv.gvar('0.0(100.0)')#gv.gvar('0.0(5.0)')
            prior['A_loga'] = gv.gvar('0.0(100.0)')#gv.gvar('0.0(5.0)')
            prior['A_aa'] = gv.gvar('0.0(1000.0)')#gv.gvar('0.0(50.0)')
            prior['A_aaa'] = gv.gvar('0.0(10000.0)')#gv.gvar('0.0(50.0)')

            # nnlo terms
            prior['A_k'] = gv.gvar('0.0(10.0)')#gv.gvar('0.0(5.0)')
            prior['A_p'] = gv.gvar('0.0(10.0)')#gv.gvar('0.0(5.0)')

            # nnnlo terms
            prior['A_ak'] = gv.gvar('0.0(100.0)')
            prior['A_ap'] = gv.gvar('0.0(100.0)')
            prior['A_kk'] = gv.gvar('0.0(100.0)')
            prior['A_kp'] = gv.gvar('0.0(100.0)')
            prior['A_pp'] = gv.gvar('0.0(100.0)')

        # Convert bootstrapped data into gvar data
        bias_corrector = lambda arr : arr[1:] + (arr[0] - np.mean(arr[1:]))
        gv_data = {}
        plot_data = {}
        for abbr in abbrs:
            gv_data[abbr] = {}
            plot_data[abbr] = {}
            for key in ['FK', 'Fpi', 'mpi', 'mk', 'mss', 'mju', 'mjs', 'mru', 'mrs']:
                if bias_correct:
                    gv_data[abbr][key] = bias_corrector(fit_data[abbr][key])
                else:
                    gv_data[abbr][key] = fit_data[abbr][key]
                plot_data[abbr][key] = fit_data[abbr][key][:plot_bs_N]

            gv_data[abbr] = gv.dataset.avg_data(gv_data[abbr], bstrap=True)
            gv_data[abbr]['FK/Fpi'] = gv_data[abbr]['FK'] / gv_data[abbr]['Fpi']
            plot_data[abbr]['FK/Fpi'] = plot_data[abbr]['FK'] / plot_data[abbr]['Fpi']

            # Sunset term
            if fast_sunset and model_info['include_sunset']:
                gv_data[abbr]['sunset'] = sf.fcn_FF((gv_data[abbr]['mpi'] / gv_data[abbr]['mk'])**2)

            to_gvar = lambda arr : gv.gvar(arr[0], arr[1])
            for key in ['a2DI']:
                gv_data[abbr][key] = to_gvar(fit_data[abbr][key])
                plot_data[abbr][key] = to_gvar(fit_data[abbr][key])

            for key in ['aw0']:
                gv_data[abbr]['a/w0'] = to_gvar(fit_data[abbr]['aw0'])
                plot_data[abbr]['a/w0'] = to_gvar(fit_data[abbr]['aw0'])

            for key in ['L', 'alpha_s']:
                value = fit_data[abbr][key]
                gv_data[abbr][key] = gv.gvar(value, value/1000000.0)
                plot_data[abbr][key] = gv.gvar(value, value/1000000.0)

            if model_info['F2'] == 'FKFpi':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['FK'] *gv_data[abbr]['Fpi']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['FK'] *plot_data[abbr]['Fpi']
            elif model_info['F2'] == 'FKFK':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['FK'] *gv_data[abbr]['FK']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['FK'] *plot_data[abbr]['FK']
            elif model_info['F2'] == 'FpiFpi':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['Fpi'] *gv_data[abbr]['Fpi']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['Fpi'] *plot_data[abbr]['Fpi']
            elif model_info['F2'] == 'F0F0':
                latt_spacing = { # Taken from arxiv/1503.02769, Table VIII
                    'a15' : gv.gvar('0.1511(18)')*1.00,
                    'a12' : gv.gvar('0.1206(11)')*1.00,
                    'a09' : gv.gvar('0.0873(11)')*1.00,
                }
                F0 = 131.5 / np.sqrt(2)
                hbar_c = 197.32698045930
                gv_data[abbr]['lam2_chi'] = (4 * np.pi *latt_spacing[abbr[:3]] *F0 / hbar_c)**2
                plot_data[abbr]['lam2_chi'] = (4 * np.pi *latt_spacing[abbr[:3]] *F0 / hbar_c)**2


        # Set object values
        self.w0 = phys_point_data['w0']
        self.bs_N = 1
        self.plot_bs_N = plot_bs_N
        self.abbrs = sorted(abbrs)
        self.fit_data = gv_data
        self.plot_data = plot_data
        self.model_info = model_info

        self._input_prior = prior
        self._fit = None
        self._phys_point_data = phys_point_data
        #self._error_budget = None
        #self._posterior = None
        #self._prior = None

    @property
    def delta_su2(self):
        lam2_chi = (4 *np.pi *gv.gvar('80(20)'))**2 #Defn of lam2_chi for delta_su2 calc per FLAG 2019
        eps2_pi = (self._get_phys_point_data('mpi'))**2 / lam2_chi
        eps2_k = (self._get_phys_point_data('mk'))**2 / lam2_chi
        fkfpi = self.fk_fpi

        R = gv.gvar('35.7(2.6)') # From FLAG
        eps_su2 = np.sqrt(3)/(4.0 *R)

        delta = np.sqrt(3) *eps_su2 *(
            - (4.0/3.0) *(fkfpi - 1)
            + (2.0/6.0) *(eps2_k - eps2_pi - eps2_pi *np.log(eps2_k/eps2_pi))
        )
        return delta

    @property
    def error_budget(self):
        return self._get_error_budget()

    def _get_error_budget(self, verbose=False, **kwargs):
        output = {}

        fk_fpi = self.fk_fpi
        prior = self.prior
        posterior = self.posterior
        phys_point_data = self.phys_point_data

        inputs = {}
        inputs.update(prior)
        output['disc'] = fk_fpi.partialsdev(
            [prior[key] for key in ['A_a', 'A_aa', 'A_loga', 'A_aaa'] if key in prior]
        )
        output['chiral'] = fk_fpi.partialsdev(
            [prior[key] for key in (set(prior) - set(['A_a', 'A_aa', 'A_loga', 'A_aaa']))]
        )

        phys_point = {}
        phys_point['phys:mpi'] = phys_point_data['mpi']
        phys_point['phys:mk'] = phys_point_data['mk']
        phys_point['phys:lam2_chi'] = phys_point_data['lam2_chi']
        inputs.update(phys_point)
        output['pp_input'] = fk_fpi.partialsdev(phys_point)


        # Since the input data is correlated,
        # we only need to use a single variable as a proxy
        # for all of the variables; we use 'lam2_chi'
        input_data = {}
        input_data['stat'] = self._get_prior('lam2_chi')
        inputs.update(input_data)
        output['stat'] = fk_fpi.partialsdev(input_data)

        if verbose:
            # default kwargs
            if kwargs is None:
                kwargs = {}
            kwargs.setdefault('percent', False)
            kwargs.setdefault('ndecimal', 5)
            kwargs.setdefault('verify', True)

            return gv.fmt_errorbudget(outputs={'FK/Fpi' : fk_fpi}, inputs=inputs, **kwargs)

        else:
            return output


    @property
    def fit(self):
        if self._fit is None:
            print("Making fit...")

            start_time = time.time()
            temp_fit = self._make_fit()
            end_time = time.time()

            self._fit = temp_fit
            print("Time (s): ",  end_time - start_time)

        return self._fit

    @property
    def fit_info(self):
        fit_info = {
            'name' : self.model,
            'FK/Fpi' : self.fk_fpi,
            'delta_su2' : self.delta_su2,
            'logGBF' : self.fit.logGBF,
            'chi2/df' : self.fit.chi2 / self.fit.dof,
            'Q' : self.fit.Q,
            'vol' : self.model_info['include_FV'],
            'phys_point' : self.phys_point_data,
            'error_budget' : self.error_budget,
            'prior' : self.prior,
            'posterior' : self.posterior,
        }

        return fit_info

    # Returns names of LECs in prior/posterior
    @property
    def fit_keys(self):
        keys1 = list(self._input_prior.keys())
        keys2 = list(self.fit.p.keys())
        parameters = np.intersect1d(keys1, keys2)
        return parameters

    @property
    def fk_fpi(self):
        return self._extrapolate_to_phys_point()

    def _extrapolate_to_phys_point(self, include_su2_isospin_corrrection=False):
        output = self.fk_fpi_fit_fcn(fit_data=self.phys_point_data.copy())
        if include_su2_isospin_corrrection:
            output *= np.sqrt(1 + self.delta_su2) # include SU(2) isospin breaking correction

        return output

    @property
    def model(self):
        name = self.model_info['fit_type'] +'_'+ self.model_info['F2'] +'_'+ self.model_info['order']
        if self.model_info['include_log']:
            name = name + '_log'
        if self.model_info['include_log2']:
            name = name + '_logSq'
        if self.model_info['include_sunset']:
            name = name + '_sunset'
        if self.model_info['include_alpha_s']:
            name = name + '_alphaS'
        if self.model_info['include_latt_n3lo']:
            name = name + '_a4'
        if self.model_info['include_latt_n4lo']:
            name = name + '_a6'
        if self.model_info['include_FV']:
            name = name + '_FV'
        if self.model_info['use_bijnens_central_value']:
            name = name + '_bijnens'
        return name

    @property
    def phys_point_data(self):
        return self._get_phys_point_data()

    # need to convert to/from lattice units
    def _get_phys_point_data(self, parameter=None):
        if parameter is None:
            return self._phys_point_data
        elif parameter == 'FK/Fpi':
            # Physical point without su(2) isospin correction
            return self._phys_point_data['FK/Fpi_pm'] / np.sqrt(1 + self.delta_su2)
        else:
            return self._phys_point_data[parameter]

    @property
    def posterior(self):
        return self._get_posterior()

    # Returns dictionary with keys fit parameters, entries gvar results
    def _get_posterior(self, param=None):
        if param is None:
            #if self._posterior is None:
            #    self._posterior = {param : self.fit.p[param] for param in self.fit_keys}
            #return self._posterior
            return {param : self.fit.p[param] for param in self.fit_keys}
        elif param == 'all':
            return self.fit.p
        else:
            return self.fit.p[param]

    @property
    def prior(self):
        return self._get_prior()

    def _get_prior(self, param=None):
        if param is None:
            #if self._prior is None:
            #    self._prior = {param : self.fit.prior[param] for param in self.fit_keys}
            #return self._prior
            return {param : self.fit.prior[param] for param in self.fit_keys}
        elif param == 'all':
            return self.fit.prior
        else:
            return self.fit.prior[param]

    def __str__(self):
        output = "\nModel: %s\n" %(self.model)
        output = output + "\nFitted/[FLAG] values at physical point:"
        output = output + '\n\tF_K / F_pi = %s [%s]'%(
                            self.fk_fpi,
                            self._get_phys_point_data('FK/Fpi'))
        output = output + '\t(delta_su2 = %s)' %(self.delta_su2)
        output = output + "\n\n"


        output += 'Parameters:\n'
        my_str = self.fit.format(pstyle='m')
        for item in my_str.split('\n'):
            for key in self.fit_keys:
                re = key+' '
                if re in item:
                    output += item + '\n'

        output += '\n'
        output += self.fit.format(pstyle=None)

        sig_fig = lambda x : np.around(x, int(np.floor(-np.log10(x))+3)) # Round to 3 sig figs
        output += '\nError Budget (relative error):\n'
        for key in self.error_budget:
            output += "\t%s: %s" %(key, sig_fig(self.error_budget[key]/self.fk_fpi.mean))

        return output

    def _fmt_key_as_latex(self, key):
        convert = {
            # data parameters
            'a/w0' : r'$a$ (fm)',
            'L' : r'$L$ (fm)',

            'mpi' : r'$m_\pi$',
            'Fpi' : r'$F_\pi$',

            'mk' : r'$m_K$',
            'FK' : r'$F_K$',

            'FK/Fpi' : r'$F_K / F_\pi$',

            'V' : r'$e^{-m L} / \sqrt{m L}$',

            # lattice artifacts
            'c_2_a' : r'$c^{(2)}_a$',
            'c_3_a' : r'$c^{(3)}_a$'
        }

        if key in convert.keys():
            return convert[key]
        else:
            return key


    def _make_fit(self):
        prepped_data = self._make_fit_data()
        # Need to randomize prior in bayesian-bootstrap hybrid
        temp_prior = self._input_prior
        temp_fitter = fit.fitter(fit_data=prepped_data, prior=temp_prior, **self.model_info)
        return temp_fitter.get_fit()

    def _make_empbayes_fit(self):
        prepped_data = self._make_fit_data()
        temp_prior = self._input_prior
        temp_fitter = fit.fitter(fit_data=prepped_data, prior=temp_prior, **self.model_info)

        empbayes_fit = temp_fitter.get_empbayes_fit()
        self._fit = empbayes_fit
        return empbayes_fit

    def _make_fit_data(self):
        prepped_data = {}
        for parameter in ['mjs', 'mju', 'mk', 'mpi', 'mrs', 'mru', 'mss', 'lam2_chi']:
            prepped_data[parameter] = np.array([self.fit_data[abbr][parameter] for abbr in self.abbrs])
        for parameter in ['sunset']:
            if parameter in self.fit_data[self.abbrs[0]]:
                prepped_data['sunset'] = np.array([self.fit_data[abbr]['sunset'] for abbr in self.abbrs])
        for parameter in ['a/w0', 'a2DI', 'L', 'alpha_s']:
            prepped_data[parameter] = np.array([self.fit_data[abbr][parameter] for abbr in self.abbrs])

        prepped_data['y'] = [self.fit_data[abbr]['FK/Fpi'] for abbr in self.abbrs]

        return prepped_data

    def create_prior_from_fit(self):
        output = {}
        temp_prior = self._make_empbayes_fit().prior
        for key in self.fit_keys:
            if key in ['L_4', 'L_5']:
                output[key] = gv.gvar(0, 0.005)
            elif key in ['A_a', 'A_k', 'A_p', 'A_loga', 'A_aa', 'A_aaa']:
                mean = 0
                sdev = gv.sdev(temp_prior[key])
                output[key] = gv.gvar(mean, sdev)
            elif ['L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9']:
                if self.model_info['use_bijnens_central_value']:
                    output[key] = temp_prior[key]
                else:
                    output[key] = gv.gvar(0, gv.sdev(temp_prior[key]))
            else:
                print('Missing key when saving prior!', key)
        return output

    def extrapolate_to_ensemble(self, abbr):
        abbr_n = self.abbrs.index(abbr)

        try: # normal fit logic
            model_name = list(self.fit.data.keys())[-1]
            return self.fit.fcn(self._get_posterior('all'))[model_name][abbr_n]
        except AttributeError: # empbayes fit logic
            return self.fit.fcn(self._get_posterior('all'))[abbr_n]

    def fk_fpi_fit_fcn(self, fit_data=None, posterior=None, fit_type=None, debug=None):
        if fit_data is None:
            fit_data = self.phys_point_data.copy()
        if posterior is None:
            posterior = self.posterior.copy()
        if fit_type is None:
            fit_type = self.model_info['fit_type']

        if fit_type in ['ma', 'xpt']:
            fit_type = 'xpt'
        elif fit_type in ['ma-ratio', 'xpt-ratio']:
            fit_type = 'xpt-ratio'

        model_info = self.model_info.copy()
        model_info['fit_type'] = fit_type

        model = fit.fitter(**model_info)._make_models()[0]
        return model.fitfcn(p=posterior, fit_data=fit_data, debug=debug)


    def fmt_error_budget(self, **kwargs):
        return self._get_error_budget(verbose=True, **kwargs)


    def make_plots(self, show_error_ellipses=False):
        figs = [self.plot_fit_info()]
        figs.append(self.plot_fit_bar_graph())

        figs.append(self.plot_fit('mpi'))
        #squared = lambda x : x**2
        #figs.append(self.plot_parameters(xy_parameters=['mpi', 'FK/Fpi'],
        #            xlabel='$F_\pi$ (MeV)', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['FK', 'FK/Fpi'],
        #            xlabel='$F_K$ (MeV)', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['mpi', 'FK/Fpi'],
        #            xfcn=squared, xlabel='$m_\pi^2$ (MeV)$^2$', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['mk', 'FK/Fpi'],
        #            xfcn=squared, xlabel='$m_K^2$ (MeV)$^2$', color_parameter='a'))

        if show_error_ellipses:
            # Make error ellipses
            # Casts indixes of fit_keys as an upper-triangular matrix,
            # thereby allowing us to get all the 2-combinations
            # of the set of fit_keys
            fit_keys = list(self.posterior)
            rs,cs = np.triu_indices(len(fit_keys),1)
            for (i, j) in zip(rs, cs):
                figs.append(self.plot_error_ellipsis(x_key=fit_keys[i], y_key=fit_keys[j]))

        return figs

    # Takes keys from posterior (eg, 'L_5' and 'L_4')
    def plot_error_ellipsis(self, x_key, y_key):
        x = self._get_posterior(x_key)
        y = self._get_posterior(y_key)


        fig, ax = plt.subplots()

        corr = '{0:.3g}'.format(gv.evalcorr([x, y])[0,1])
        std_x = '{0:.3g}'.format(gv.sdev(x))
        std_y = '{0:.3g}'.format(gv.sdev(y))
        text = ('$R_{x, y}=$ %s\n $\sigma_x =$ %s\n $\sigma_y =$ %s'
                % (corr,std_x,std_y))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        C = gv.evalcov([x, y])
        eVe, eVa = np.linalg.eig(C)
        for e, v in zip(eVe, eVa.T):
            plt.plot([gv.mean(x)-1*np.sqrt(e)*v[0], 1*np.sqrt(e)*v[0] + gv.mean(x)],
                     [gv.mean(y)-1*np.sqrt(e)*v[1], 1*np.sqrt(e)*v[1] + gv.mean(y)],
                     'k-', lw=2)

        #plt.scatter(x-np.mean(x), y-np.mean(y), rasterized=True, marker=".", alpha=100.0/self.bs_N)
        #plt.scatter(x, y, rasterized=True, marker=".", alpha=100.0/self.bs_N)

        plt.grid()
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.xlabel(self._fmt_key_as_latex(x_key), fontsize = 24)
        plt.ylabel(self._fmt_key_as_latex(y_key), fontsize = 24)

        fig = plt.gcf()
        plt.close()
        return fig

    def plot_fit_info(self):
        plt.axis('off')
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        text = self.__str__()#.expandtabs()
        plt.text(0.0, 1.0, str(text),
                 fontsize=12, ha='left', va='top', family='monospace', bbox=props)

        plt.tight_layout()
        fig = plt.gcf()
        plt.close()

        return fig

    def plot_fit_bar_graph(self):

        y = 1
        labels = np.array([])

        # Physical point:
        for j in range(1):
            plt.axhline(y-2, ls ='-', color='C4')

            # FLAG
            data = self._get_phys_point_data('FK/Fpi')
            x = gv.mean(data)
            xerr = gv.sdev(data)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0, label='Latt/FLAG')
            labels = np.append(labels, str(""))
            y = y - 1

            # fit result
            fit_value = self.fk_fpi
            x = gv.mean(fit_value)
            xerr = gv.sdev(fit_value)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C2', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C3', elinewidth=10.0, label='Fit')
            y = y + 2



            labels = np.append(labels, str("Phys. point"))
            plt.axhline(y, ls='--')

            y = y + 1
            labels = np.append(labels, str(""))
            plt.axhline(y-1, ls ='-', color='C4')


        for abbr in self.abbrs:
            # fit result
            fit_value = self.extrapolate_to_ensemble(abbr)
            x = gv.mean(fit_value)
            xerr = gv.sdev(fit_value)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C2', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C3', elinewidth=10.0)
            y = y + 1

            # data
            data = self.fit_data[abbr]['FK/Fpi']
            x = gv.mean(data)
            xerr = gv.sdev(data)

            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)
            labels = np.append(labels, str(""))
            y = y + 1

            labels = np.append(labels, str(abbr))
            plt.axhline(y, ls='--')

            y = y + 1
            labels = np.append(labels, str(""))


        plt.legend()

        plt.yticks(1*list(range(len(labels))), labels, fontsize=15, rotation=45)
        plt.ylim(-1, y)
        plt.xlabel(self._fmt_key_as_latex('FK/Fpi'), fontsize = 24)

        fig = plt.gcf()
        plt.close()


        return fig


    # fit vs param (eg, param = 'mpi' or 'a')
    def plot_fit(self, param):

        # used to convert to phys units
        hbar_c = 197.32698045930
        if param == 'mpi':
            phys_params = ['mk']
            xlabel = '$\epsilon_\pi^2$'
        elif param == 'mk':
            phys_params = ['mpi']
            xlabel = '$\epsilon_K^2$'
        elif param == 'a':
            phys_params = ['mpi', 'mk']
            xlabel = '$\epsilon_a^2$'
        else:
            return None

        # Data for plot
        plot_data = {}

        # Convert to eps^2
        if param in ['mpi', 'mk']:
            plot_data['x'] = {abbr : ((self.plot_data[abbr][param] *hbar_c
                                       / (self.plot_data[abbr]['a/w0'] *self.w0))**2
                                       /self._get_phys_point_data('lam2_chi'))
                              for abbr in self.abbrs}
        elif param in ['a']:
            plot_data['x'] = {abbr : (self.plot_data[abbr]['a/w0'] / (4 *np.pi))**2 for abbr in self.abbrs}

        plot_data['y'] = self.shift_fk_fpi_for_phys_params(phys_params=phys_params)
        color_data = {abbr : self.plot_data[abbr]['a/w0'] *self.w0 for abbr in self.abbrs}

        # Get lattice spacings used in fit
        lattice_spacings = np.unique(self._make_fit_data()['a/w0']) *self.w0

        # Color lattice spacings
        cmap = matplotlib.cm.get_cmap('rainbow_r')
        colors = [cmap(j/len(lattice_spacings)) for j in range(len(lattice_spacings)+1)]

        # Plot fit
        for j, a in enumerate(sorted(np.append([gv.gvar('0(0)')], lattice_spacings), reverse=True)):

            # Get the range of m^2 (in phys units)
            if param in ['mpi', 'mk']:
                minimum = np.nanmin([np.nanmin(
                    np.sqrt([plot_data['x'][abbr] *self._get_phys_point_data('lam2_chi') for abbr in self.abbrs])
                ) for abbr in self.abbrs])
                maximum = np.nanmax([np.nanmax(
                    np.sqrt([plot_data['x'][abbr] *self._get_phys_point_data('lam2_chi')  for abbr in self.abbrs])
                ) for abbr in self.abbrs])
            # Get the range of a^2
            elif param in ['a']:
                minimum = 0
                maximum = np.nanmax([np.nanmax(
                    np.sqrt([plot_data['x'][abbr] *(4 *np.pi)**2 for abbr in self.abbrs])
                ) for abbr in self.abbrs])

            minimum = gv.mean(minimum)
            maximum = gv.mean(maximum)
            delta = maximum - minimum

            x = np.linspace(np.max((minimum - 0.05*delta, 0)), maximum + 0.05*delta)

            # Get phys point data, substituting x-data and current 'a' in loop
            prepped_data = self.phys_point_data.copy()
            prepped_data['a/w0'] = a/self.w0
            if param in ['mpi', 'mk']:
                prepped_data[param] = x
            elif param in ['a']:
                prepped_data['a/w0'] = x

            # Covert m -> eps^2
            if param in ['mpi', 'mk']:
                x = x**2 / self._get_phys_point_data('lam2_chi')
            elif param in ['a']:
                x = x**2 / (4 *np.pi)**2

            y = self.fk_fpi_fit_fcn(fit_data=prepped_data)

            pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
            plt.plot(pm(x, 0), pm(y, 0), '--', color=colors[j%len(colors)], label='$a=$%s (fm)'%(str(a)), rasterized=True)
            plt.fill_between(pm(x, 0), pm(y, -1), pm(y, 1), alpha=0.40, color=colors[j%len(colors)], rasterized=True)

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        minimum = 0
        maximum = np.max([gv.mean(color_data[abbr]) for abbr in self.abbrs])
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        length = self.plot_bs_N
        x = np.zeros(length * len(self.abbrs))
        y = np.zeros(length * len(self.abbrs))
        z = np.zeros(length * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*length:(j+1)*length] = gv.mean(plot_data['x'][abbr])
            y[j*length:(j+1)*length] = gv.mean(plot_data['y'][abbr])
            z[j*length:(j+1)*length] = gv.mean(color_data[abbr])

        # Plot data
        sc = plt.scatter(x, y, c=z, vmin=minimum, vmax=maximum,
                         cmap=cmap, rasterized=True, marker="o", alpha=10.0/self.plot_bs_N, edgecolor='black')

        # Plot FLAG result
        x_phys = self._get_phys_point_data(param)**2 / self._get_phys_point_data('lam2_chi')
        plt.axvline(gv.mean(x_phys), label='Phys point')
        y_phys = self._get_phys_point_data('FK/Fpi')

        plt.errorbar(x=gv.mean(x_phys), xerr=0,
                     y=gv.mean(y_phys), yerr=gv.sdev(y_phys), label='FLAG',
                    color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)

        # Plot labels
        plt.legend()
        plt.grid()
        plt.xlabel(xlabel, fontsize = 24)
        plt.ylabel('$F_K/F_\pi$', fontsize = 24)

        # Format colorbar
        color_bar = plt.colorbar(sc)
        color_bar.set_alpha(0.8)
        color_bar.draw_all()
        color_bar.set_label('$a$ (fm)', fontsize = 24)

        # Set xlim, ylim -- only works if xy_parameters[i] is a vector, not a scalar
        min_max = lambda x : [np.min(x), np.max(x)]
        try:
            xmin, xmax = min_max(np.concatenate([gv.mean(plot_data['x'][abbr]) for abbr in self.abbrs]))
            ymin, ymax = min_max(np.concatenate([gv.mean(plot_data['y'][abbr]) for abbr in self.abbrs]))
            xdelta = xmax - xmin
            ydelta = ymax - ymin
            plt.xlim(xmin-0.05*xdelta, xmax+0.05*xdelta) #xmin-0.05*xdelta
            #plt.ylim(ymin-0.05*ydelta, ymax+0.05*ydelta)
            #plt.ylim(1.04, 1.20)
        except ValueError:
            pass

        fig = plt.gcf()
        plt.close()

        return fig

    def plot_parameters(self, xy_parameters, color_parameter=None,
                        xfcn=None, xlabel=None, yfcn=None, ylabel=None):

        # used to convert to phys units
        hbar_c = 197.32698045930

        if xlabel is None:
            xlabel = self._fmt_key_as_latex(xy_parameters[0])
        if ylabel is None:
            ylabel = self._fmt_key_as_latex(xy_parameters[1])

        if xfcn is None:
            xfcn = lambda x : 1 * x
        if yfcn is None:
            yfcn = lambda y : 1 * y

        # Make plot data
        plot_data = {}
        myfcn = [xfcn, yfcn]
        for j, parameter in enumerate(xy_parameters):
            if parameter in ['FK/Fpi', 'FK / Fpi']:
                plot_data[j] = {abbr :  myfcn[j](self.plot_data[abbr]['FK'] /self.plot_data[abbr]['Fpi']) for abbr in self.abbrs}
            elif parameter in ['mpi', 'mju', 'mru', 'mk', 'mrs', 'mss', 'FK', 'Fpi']:
                # Convert to physical units
                plot_data[j] = {}
                for abbr in self.abbrs:
                    plot_data[j][abbr] = myfcn[j](self.plot_data[abbr][parameter] *hbar_c / (self.plot_data[abbr]['a/w0'] *self.w0))
            elif parameter in ['mpiL', 'mjuL', 'mruL', 'mkL', 'mrsL', 'mssL']:
                m_name = parameter[:-1]
                plot_data[j] = {}
                for abbr in self.abbrs:
                    plot_data[j][abbr] = myfcn[j](self.plot_data[abbr][m_name] *self.plot_data[abbr]['L'])
            else:
                plot_data[j] = {abbr :  myfcn[j](self.plot_data[abbr][parameter]) for abbr in self.abbrs}


        # Get data for color coding graph
        if color_parameter is None:
            color_parameter = 'a'

        if color_parameter in ['a']:
            color_data = {abbr : np.repeat(self.plot_data[abbr]['a/w0'] *self.w0, self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter in ['L']:
            color_data = {abbr : np.repeat(gv.mean(self.plot_data[abbr][color_parameter]), self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter == 'mpiL':
            color_data = {abbr : gv.mean(self.plot_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}
        else:
            color_data = {abbr : gv.mean(self.plot_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        min_max = lambda x : [np.min(x), np.max(x)]
        minimum, maximum = min_max(np.concatenate([gv.mean(color_data[abbr]) for abbr in self.abbrs]))
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        x = np.zeros(self.plot_bs_N * len(self.abbrs))
        y = np.zeros(self.plot_bs_N * len(self.abbrs))
        z = np.zeros(self.plot_bs_N * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*self.plot_bs_N:(j+1)*self.plot_bs_N] = gv.mean(plot_data[0][abbr])
            y[j*self.plot_bs_N:(j+1)*self.plot_bs_N] = gv.mean(plot_data[1][abbr])
            z[j*self.plot_bs_N:(j+1)*self.plot_bs_N] = gv.mean(color_data[abbr])

        # Plot data
        sc = plt.scatter(x, y, c=z, vmin=minimum, vmax=maximum,
                         cmap=cmap, rasterized=True, marker=".", alpha=100.0/self.bs_N)

        # Plot labels
        plt.grid()
        plt.xlabel(xlabel, fontsize = 24)
        plt.ylabel(ylabel, fontsize = 24)

        # Format colorbar
        color_bar = plt.colorbar(sc)
        color_bar.set_alpha(0.8)
        color_bar.draw_all()
        color_bar.set_label(self._fmt_key_as_latex(color_parameter), fontsize = 24)

        # Set xlim, ylim -- only works if xy_parameters[i] is a vector, not a scalar
        min_max = lambda x : [np.min(x), np.max(x)]
        try:
            xmin, xmax = min_max(np.concatenate([gv.mean(plot_data[0][abbr]) for abbr in self.abbrs]))
            ymin, ymax = min_max(np.concatenate([gv.mean(plot_data[1][abbr]) for abbr in self.abbrs]))
            xdelta = xmax - xmin
            ydelta = ymax - ymin
            plt.xlim(xmin-0.05*xdelta, xmax+0.05*xdelta) #xmin-0.05*xdelta
            #plt.ylim(ymin-0.05*ydelta, ymax+0.05*ydelta)
        except ValueError:
            pass

        fig = plt.gcf()
        plt.close()
        return fig

    # phys_params is a list: eg, ['mk'] or ['mk', 'mpi']
    def shift_fk_fpi_for_phys_params(self, phys_params, use_ratio=False):
        hbar_c = 197.32698045930 # MeV-fm
        lam2_chi = self._get_phys_point_data('lam2_chi')

        shifted_fkfpi = {}
        for abbr in self.abbrs:

            fkfpi_ens = self.plot_data[abbr]['FK'] / self.plot_data[abbr]['Fpi']
            fkfpi_fit_ens = self.extrapolate_to_ensemble(abbr)

            temp_data = {}
            for param in list(self.plot_data[abbr]):
                # Covert to phys units
                if param in ['FK', 'Fpi', 'mpi', 'mk']:
                    temp_data[param] = (self.plot_data[abbr][param] *hbar_c
                                        / (self.plot_data[abbr]['a/w0'] *self.w0))
                else:
                    temp_data[param] = self.plot_data[abbr][param]

            for param in phys_params:
                temp_data[param] = self._get_phys_point_data(param)

            temp_data['lam2_chi'] = self._get_phys_point_data('lam2_chi')

            fkfpi_fit_phys = self.fk_fpi_fit_fcn(fit_data=temp_data)

            if use_ratio:
                shifted_fkfpi[abbr] = fkfpi_ens *(fkfpi_fit_phys / fkfpi_fit_ens)
            else:
                shifted_fkfpi[abbr] = fkfpi_ens + fkfpi_fit_phys - fkfpi_fit_ens

        return shifted_fkfpi
