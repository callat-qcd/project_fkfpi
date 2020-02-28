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

    def __init__(self, fit_data, phys_point_data, prior=None, **kwargs):

        # Default values
        order = None
        fit_type = 'xpt'
        abbrs = fit_data.keys()
        include_su2_isospin_corrrection = False
        use_bijnens_central_value = True
        bias_correct = True
        plot_bs_N = 100
        fast_sunset = False
        F2 = 'FpiFpi'

        # Overwrite defaults using kwargs
        if 'fit_type' in kwargs:
            fit_type = kwargs['fit_type']

        if 'abbrs' in kwargs:
            abbrs = kwargs['abbrs']

        if 'include_su2_isospin_corrrection' in kwargs:
            include_su2_isospin_corrrection = kwargs['include_su2_isospin_corrrection']

        if 'use_bijnens_central_value' in kwargs:
            use_bijnens_central_value = kwargs['use_bijnens_central_value']

        if 'order' in kwargs:
            order = kwargs['order']

        if 'F2' in kwargs:
            F2 = kwargs['F2']
            # Set this now so get_phys_point_data('lam2_chi') works
            self.F2 = F2

        if 'bias_correct' in kwargs:
            bias_correct = kwargs['bias_correct']

        if 'plot_bs_N' in kwargs:
            plot_bs_N = kwargs['plot_bs_N']

        if 'fast_sunset' in kwargs:
            fast_sunset = kwargs['fast_sunset']

        # Get lam_chi^2=renorm_scale^2 depending on choice of F^2
        if F2 == 'FKFpi':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_kpi']
        elif F2 == 'FpiFpi':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_pipi']
        elif F2 == 'FKFK':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_kk']
        elif F2 == 'F0F0':
            phys_point_data['lam2_chi'] = phys_point_data['lam2_chi_00']

        # Add default values to order dict
        order_temp = {
            'fit' : 'nlo',
            'vol' : 1,
            'exclude' : [],

            # semi-nnlo corrections
            'include_alpha_s' : False,
            'include_log' : False,
            'include_log2' : False,
            'include_sunset' : False,

            # nnnlo corrections
            'include_latt_n3lo' : False,
        }
        if order is None:
            order = order_temp
        else:
            for key in order_temp.keys():
                if key not in order.keys():
                    order[key] = order_temp[key]

        # Default prior
        if prior is None:
            print("Using default prior.")

            # From Bijnens
            gamma = np.array([3./32, 3./16, 0, 1./8, 3./8, 11./144, 0, 5./48])
            L_mu0 = np.array([3.72, 4.93, -30.7, 0.89, 3.77, 0.11, -3.4, 2.94])*10**-4

            # F2 is F^2 in this context; but F0, F1 are F_0, F_1
            if F2 == 'FpiFpi':
                F1 = phys_point_data['Fpi']
            elif F2 == 'FKFpi':
                F1 = np.sqrt(phys_point_data['FK'] *phys_point_data['Fpi'])
            elif F2 == 'FKFK':
                F1 = phys_point_data['Fpi']

            F0 = 80 # MeV
            L_mu1 = gv.mean(L_mu0 - gamma/(4 *np.pi)**2 *np.log(F1/F0))
            if use_bijnens_central_value:
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
            prior['A_c'] = gv.gvar('0.0(0.05)')

            # Lattice spacing terms
            prior['A_a'] = gv.gvar('0.0(100.0)')#gv.gvar('0.0(5.0)')
            prior['A_loga'] = gv.gvar('0.0(100.0)')#gv.gvar('0.0(5.0)')
            prior['A_aa'] = gv.gvar('0.0(1000.0)')#gv.gvar('0.0(50.0)')

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
            if fast_sunset and order['include_sunset']:
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

            if F2 == 'FKFpi':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['FK'] *gv_data[abbr]['Fpi']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['FK'] *plot_data[abbr]['Fpi']
            elif F2 == 'FKFK':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['FK'] *gv_data[abbr]['FK']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['FK'] *plot_data[abbr]['FK']
            elif F2 == 'FpiFpi':
                gv_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *gv_data[abbr]['Fpi'] *gv_data[abbr]['Fpi']
                plot_data[abbr]['lam2_chi'] = (4 *np.pi)**2 *plot_data[abbr]['Fpi'] *plot_data[abbr]['Fpi']
            elif F2 == 'F0F0':
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
        self.phys_point_data = phys_point_data
        self.include_su2_isospin_corrrection = include_su2_isospin_corrrection
        self.use_bijnens_central_value = use_bijnens_central_value
        self.bs_N = 1
        self.F2 = F2
        self.plot_bs_N = plot_bs_N
        self.abbrs = sorted(abbrs)
        self.fit_data = gv_data
        self.plot_data = plot_data
        self.prior = prior
        self.posterior = None
        self.order = order
        self.fits = None
        self.fit_type = fit_type
        self.fast_sunset = fast_sunset

    @property
    def fit(self):
        if self.fits is None:
            self.bootstrap_fits()
        return self.fits[0]

    @property
    def fit_info(self):
        fit_info = {
            'name' : self.get_name(),
            'FK/Fpi' : self.extrapolate_to_phys_point(),
            'delta_su2' : self.get_delta_su2_correction(),
            'logGBF' : self.get_fit().logGBF,
            'chi2/df' : self.get_fit().chi2 / self.get_fit().dof,
            'Q' : self.get_fit().Q,
            'vol' : self.order['vol'],
            'phys_point' : self.get_phys_point_data(),
            'prior' : {},
            'posterior' : {},
            'error_budget' : self.get_error_budget()
        }

        for key in self.get_fit_keys():
            fit_info['prior'][key] = self.get_prior(key)
            fit_info['posterior'][key] = self.get_posterior(key)

        return fit_info


    def __str__(self):
        prior = self.prior
        output = "\nModel: %s" %(self.get_name())
        output = output + "\n\nFitting to %s \n" %(self.order['fit'])
        #output = output + " with lattice corrections O(%s) \n" %(self.order['latt_spacing'])
        output = output + " with volume corrections O(%s) \n" %(self.order['vol'])
        output = output + "Fitted/[FLAG] values at physical point (including SU(2) isospin corrections: %s):\n" %(self.include_su2_isospin_corrrection)
        output = output + '\nF_K / F_pi = %s [%s]'%(
                            self.extrapolate_to_phys_point(),
                            self.get_phys_point_data('FK/Fpi'))
        output = output + '   (delta_su2 = %s)' %(self.get_delta_su2_correction())
        output = output + "\n\n"

        fit_parameters = self.get_posterior()
        new_table = { key : [fit_parameters[key], prior[key]] for key in sorted(fit_parameters.keys())}

        output = output + gv.tabulate(new_table, ncol=2, headers=['Parameter', 'Result[0] / Prior[1]'])

        output = output + '\n---\nboot0 fit results:\n'
        output = output + self.get_fit().format(pstyle=None)

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

    def _make_fit(self, j):
        prepped_data = self._make_fit_data(j)
        # Need to randomize prior in bayesian-bootstrap hybrid
        temp_prior = self._randomize_prior(self.prior, j)
        temp_fitter = fit.fitter(fit_data=prepped_data, prior=temp_prior, F2=self.F2,
                        order=self.order, fit_type=self.fit_type, fast_sunset=self.fast_sunset)
        return temp_fitter.get_fit()

    def _make_empbayes_fit(self):
        prepped_data = self._make_fit_data(0)
        temp_prior = self._randomize_prior(self.prior, 0)
        temp_fitter = fit.fitter(fit_data=prepped_data, prior=temp_prior, F2=self.F2,
                        order=self.order, fit_type=self.fit_type, fast_sunset=self.fast_sunset)

        empbayes_fit = temp_fitter.get_empbayes_fit()
        self.fits = [empbayes_fit]
        return empbayes_fit

    def _make_fit_data(self, j):
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

    # This will consistently give the same results despite randomizing prior
    def _randomize_prior(self, prior, j):
        if j==0:
            return prior

        new_prior = {}
        for key in prior.keys():
            p_mean = gv.mean(prior[key])
            p_sdev = gv.sdev(prior[key])

            p_re = re.sub('[^A-Za-z0-9]+','',key)
            if len(p_re) > 6:
                p_re = p_re[0:6]
            s = int(p_re,36)

            np.random.seed(s+j)
            new_prior[key] = gv.gvar(np.random.normal(p_mean,scale=p_sdev), p_sdev)

        return new_prior

    def bootstrap_fits(self):
        print("Making fits...")
        start_time = time.time()
        self.fits = np.array([])

        for j in range(self.bs_N):
            temp_fit = self._make_fit(j)
            self.fits = np.append(self.fits, temp_fit)

            sys.stdout.write("\r{0}% complete".format(int((float(j+1)/self.bs_N)*100)))
            print("",)
            sys.stdout.flush()
        end_time = time.time()
        print("Time (s): ",  end_time - start_time)
        print("Compiling results...")

    def create_prior_from_fit(self):
        output = {}
        temp_prior = self._make_empbayes_fit().prior
        for key in self.get_fit_keys():
            if key in ['L_4', 'L_5']:
                output[key] = gv.gvar(0, 0.005)
            elif key in ['A_a', 'A_k', 'A_p', 'A_loga', 'A_aa']:
                mean = 0
                sdev = gv.sdev(temp_prior[key])
                output[key] = gv.gvar(mean, sdev)
            elif ['L_1', 'L_2', 'L_3', 'L_6', 'L_7', 'L_8', 'L_9']:
                if self.use_bijnens_central_value:
                    output[key] = temp_prior[key]
                else:
                    output[key] = gv.gvar(0, gv.sdev(temp_prior[key]))
            else:
                print('Missing key when saving prior!', key)
        return output

    def extrapolate_to_ensemble(self, abbr):
        abbr_n = self.abbrs.index(abbr)
        model_name = list(self.get_fit().data.keys())[-1] # pretty hacky way to get this

        if self.bs_N == 1:
            temp_fit = self.get_fit()
            return temp_fit.fcn(temp_fit.p)[model_name][abbr_n]

        else:
            results = []
            for j in self.bs_N:
                temp_fit = self.fits[j]
                results.append(gv.mean(temp_fit.fcn(temp_fit.p)[model_name][abbr_n]))
            return gv.gvar(np.mean(results), np.std(results))

    def extrapolate_to_phys_point(self, include_su2_isospin_corrrection=None):
        if include_su2_isospin_corrrection is None:
            include_su2_isospin_corrrection = self.include_su2_isospin_corrrection

        output = self.fk_fpi_fit_fcn(fit_data=self.get_phys_point_data().copy())
        if include_su2_isospin_corrrection:
            output *= np.sqrt(1 + self.get_delta_su2_correction()) # include SU(2) isospin breaking correction

        # Logic for frequentist vs bayesian fits
        try:
            return output[0]
        except TypeError:
            return output

    def fk_fpi_fit_fcn(self, fit_data=None, fit_parameters=None, fit_type=None, debug=None):
        if fit_data is None:
            fit_data = self.get_phys_point_data().copy()
        if fit_parameters is None:
            fit_parameters = self.get_posterior().copy()
        if fit_type is None:
            fit_type = self.fit_type

        if fit_type in ['ma', 'xpt']:
            fit_type = 'xpt'
        elif fit_type in ['ma-ratio', 'xpt-ratio']:
            fit_type = 'xpt-ratio'

        model = fit.fitter(order=self.order, fit_type=fit_type, F2=self.F2, fast_sunset=self.fast_sunset)._make_models()[0]
        return model.fitfcn(p=fit_parameters, fit_data=fit_data, debug=debug)

    def get_delta_su2_correction(self):
        lam2_chi = (4 *np.pi *gv.gvar('80(20)'))**2 #lam2_chi = self.get_phys_point_data('lam2_chi')
        eps2_pi = (self.get_phys_point_data('mpi'))**2 / lam2_chi
        eps2_k = (self.get_phys_point_data('mk'))**2 / lam2_chi
        fkfpi = self.extrapolate_to_phys_point()

        R = gv.gvar('35.7(2.6)') # From FLAG
        eps_su2 = np.sqrt(3)/(4.0 *R)

        delta = np.sqrt(3) *eps_su2 *(
            - (4.0/3.0) *(fkfpi - 1)
            + (2.0/6.0) *(eps2_k - eps2_pi - eps2_pi *np.log(eps2_k/eps2_pi))
        )
        return delta

    def get_error_budget(self, print_budget=False):
        output = {}

        fk_fpi = self.extrapolate_to_phys_point()
        prior = self.get_prior()
        posterior = self.get_posterior()
        phys_point_data = self.get_phys_point_data()

        inputs = {}
        inputs.update(prior)
        output['disc'] = fk_fpi.partialsdev(
            [prior[key] for key in ['A_a', 'A_aa', 'A_loga'] if key in prior]
        )
        output['chiral'] = fk_fpi.partialsdev(
            [prior[key] for key in (set(prior) - set(['A_a', 'A_aa', 'A_loga']))]
        )

        phys_point = {}
        phys_point['mpi'] = phys_point_data['mpi']
        phys_point['mk'] = phys_point_data['mk']
        phys_point['lam2_chi'] = phys_point_data['lam2_chi']
        inputs.update(phys_point)
        output['pp_input'] = fk_fpi.partialsdev(phys_point)


        # Since the input data is correlated,
        # we only need to use a single variable as a proxy
        # for all of the variables; we use 'lam2_chi'
        input_data = {}
        input_data['input_data'] = self.get_prior('lam2_chi')
        inputs.update(input_data)
        output['stat'] = fk_fpi.partialsdev(input_data)

        if print_budget:
            print('FK/Fpi =', fk_fpi, '\n')
            print(gv.fmt_errorbudget(outputs={'FK/Fpi' : fk_fpi}, inputs=inputs, percent=False, ndecimal=5, verify=True))

        return output

    def get_fit(self):
        if self.fits is None:
            self.bootstrap_fits()
        return self.fits[0]

    def get_fit_info(self):
        fit_info = {
            'name' : self.get_name(),
            'FK/Fpi' : self.extrapolate_to_phys_point(),
            'delta_su2' : self.get_delta_su2_correction(),
            'logGBF' : self.get_fit().logGBF,
            'chi2/df' : self.get_fit().chi2 / self.get_fit().dof,
            'Q' : self.get_fit().Q,
            'vol' : self.order['vol'],
            'phys_point' : self.get_phys_point_data(),
            'prior' : {},
            'posterior' : {},
            'error_budget' : self.get_error_budget()
        }

        for key in self.get_fit_keys():
            fit_info['prior'][key] = self.get_prior(key)
            fit_info['posterior'][key] = self.get_posterior(key)

        return fit_info


    # Returns keys of fit parameters
    def get_fit_keys(self):
        if self.fits is None:
            self.bootstrap_fits()

        keys1 = list(self.prior.keys())
        keys2 = list(self.get_fit().p.keys())
        parameters = np.intersect1d(keys1, keys2)
        return sorted(parameters)

    # Returns dictionary with keys fit parameters, entries gvar results
    def get_posterior(self, parameter=None):
        # If fitting only the central values
        if self.posterior is None:
            if self.fits is None:
                self.bootstrap_fits()

            keys1 = list(self.prior.keys())
            keys2 = list(self.get_fit().p.keys())
            parameters = np.intersect1d(keys1, keys2)

            self.posterior = {parameter : self.get_fit().p[parameter] for parameter in parameters}

        if parameter is not None:
            return self.posterior[parameter]
        else:
            return self.posterior

    def get_name(self):
        name = self.fit_type +'_'+ self.F2+'_'+self.order['fit']
        if self.order['include_log']:
            name = name + '_log'
        if self.order['include_log2']:
            name = name + '_logSq'
        if self.order['include_sunset']:
            name = name + '_sunset'
        if self.order['include_alpha_s']:
            name = name + '_alphaS'
        if self.order['include_latt_n3lo']:
            name = name + '_a4'
        if self.order['vol'] > 6:
            name = name + '_FV'
        if self.use_bijnens_central_value:
            name = name + '_bijnens'
        return name


    # need to convert to/from lattice units
    def get_phys_point_data(self, parameter=None):
        if parameter is None:
            return self.phys_point_data
        elif parameter == 'FK/Fpi':
            # Physical point without su(2) isospin correction
            return self.phys_point_data['FK/Fpi_pm'] / np.sqrt(1 + self.get_delta_su2_correction())
        else:
            return self.phys_point_data[parameter]

    # Returns dictionary with keys fit parameters, entries gvar results
    def get_posterior(self, parameter=None):
        # If fitting only the central values
        if self.posterior is None:
            if self.fits is None:
                self.bootstrap_fits()

            keys1 = list(self.prior.keys())
            keys2 = list(self.get_fit().p.keys())
            parameters = np.intersect1d(keys1, keys2)

            self.posterior = {parameter : self.get_fit().p[parameter] for parameter in parameters}

        if parameter is not None:
            return self.posterior[parameter]
        else:
            return self.posterior

    def get_prior(self, key=None):
        if key is None:
            output = {key : self.get_fit().prior[key] for key in self.get_fit_keys()}
            return output
        if key is not None:
            return self.get_fit().prior[key]


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
            fit_keys = list(self.get_posterior())
            rs,cs = np.triu_indices(len(fit_keys),1)
            for (i, j) in zip(rs, cs):
                figs.append(self.plot_error_ellipsis([fit_keys[i], fit_keys[j]]))

        return figs

    # Takes keys from posterior (eg, 'L_5' and 'L_4')
    def plot_error_ellipsis(self, x_key, y_key):
        x = self.get_posterior(x_key)
        y = self.get_posterior(y_key)


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
            data = self.get_phys_point_data('FK/Fpi')
            x = gv.mean(data)
            xerr = gv.sdev(data)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0, label='Latt/FLAG')
            labels = np.append(labels, str(""))
            y = y - 1

            # fit result
            fit_value = self.extrapolate_to_phys_point()
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
                                       /self.get_phys_point_data('lam2_chi'))
                              for abbr in self.abbrs}
        elif param in ['a']:
            plot_data['x'] = {abbr : (self.plot_data[abbr]['a/w0'] / (4 *np.pi))**2 for abbr in self.abbrs}

        plot_data['y'] = self.shift_fk_fpi_for_phys_params(phys_params=phys_params)
        color_data = {abbr : self.plot_data[abbr]['a/w0'] *self.w0 for abbr in self.abbrs}

        # Get lattice spacings used in fit
        lattice_spacings = np.unique(self._make_fit_data(0)['a/w0']) *self.w0

        # Color lattice spacings
        cmap = matplotlib.cm.get_cmap('rainbow_r')
        colors = [cmap(j/len(lattice_spacings)) for j in range(len(lattice_spacings)+1)]

        # Plot fit
        for j, a in enumerate(sorted(np.append([gv.gvar('0(0)')], lattice_spacings), reverse=True)):

            # Get the range of m^2 (in phys units)
            if param in ['mpi', 'mk']:
                minimum = np.nanmin([np.nanmin(
                    np.sqrt([plot_data['x'][abbr] *self.get_phys_point_data('lam2_chi') for abbr in self.abbrs])
                ) for abbr in self.abbrs])
                maximum = np.nanmax([np.nanmax(
                    np.sqrt([plot_data['x'][abbr] *self.get_phys_point_data('lam2_chi')  for abbr in self.abbrs])
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
            prepped_data = self.get_phys_point_data().copy()
            prepped_data['a/w0'] = a/self.w0
            if param in ['mpi', 'mk']:
                prepped_data[param] = x
            elif param in ['a']:
                prepped_data['a/w0'] = x

            # Covert m -> eps^2
            if param in ['mpi', 'mk']:
                x = x**2 / self.get_phys_point_data('lam2_chi')
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
        x_phys = self.get_phys_point_data(param)**2 / self.get_phys_point_data('lam2_chi')
        plt.axvline(gv.mean(x_phys), label='Phys point')
        y_phys = self.get_phys_point_data('FK/Fpi')

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
        lam2_chi = self.get_phys_point_data('lam2_chi')

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
                temp_data[param] = self.get_phys_point_data(param)

            temp_data['lam2_chi'] = self.get_phys_point_data('lam2_chi')

            fkfpi_fit_phys = self.fk_fpi_fit_fcn(fit_data=temp_data)

            if use_ratio:
                shifted_fkfpi[abbr] = fkfpi_ens *(fkfpi_fit_phys / fkfpi_fit_ens)
            else:
                shifted_fkfpi[abbr] = fkfpi_ens + fkfpi_fit_phys - fkfpi_fit_ens

        return shifted_fkfpi
