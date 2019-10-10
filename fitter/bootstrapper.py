import lsqfit
import numpy as np
import gvar as gv
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
import re
from mpl_toolkits.mplot3d.axes3d import Axes3D

from fitter import fitter

class bootstrapper(object):

    def __init__(self, fit_data, prior=None, abbrs=None, bs_N=None,
                 order=None, fit_type=None, F2=None, chain_fits=None, use_boot0=None):

        w0 = gv.gvar('0.175(10)') # Still needs to be determined, but something like this
        self.w0 = w0

        if fit_type is None:
            fit_type = 'ma-taylor'

        if F2 is None:
            F2 = 'FKFpi'
        self.F2 = F2

        if chain_fits is None:
            chain_fits = True

        if bs_N is None or bs_N==0:
            bs_N = len(fit_data[fit_data.keys()[0]]['mpi'])
        plot_bs_N = 100

        # Add default values to order dict
        order_temp = {
            'fit' : 'nlo',
            'vol' : 1,
            'exclude' : [],
            'include_log' : False,
            'include_log2' : True
        }
        if order is None:
            order = order_temp
        else:
            for key in order_temp.keys():
                if key not in order.keys():
                    order[key] = order_temp[key]

        if prior is None:
            print "Using default prior."
            prior = {
                # nlo terms
                'L_5' : '0.0(0.005)', #'0(0.001)',
                'L_4' : '0.0(0.005)',

                # nlo-log terms
                'A_loga' : '0.0(5.0)',

                # nnlo terms
                'A_a' : '0.0(5.0)', #'0(100)',
                'A_k' : '0.0(5.0)', #'0(1)',
                'A_p' : '0.0(5.0)', #'0(10)',

                # nnnlo terms
                'A_aa' : '0(5)', #'0(100000)',
                'A_ak' : '0(5)', #'0(1000)',
                'A_ap' : '0(5)',#'0(10000)',
                'A_kk' : '0(5)', #'0(10)',
                'A_kp' : '0(5)', #'0(100)',
                'A_pp' : '0(5)', #'0(1000)',
            }
            prior = gv.gvar(prior)

        if abbrs is None:
            abbrs = fit_data.keys()

        gv_data = {}
        plot_data = {}
        for abbr in abbrs:
            gv_data[abbr] = {}
            plot_data[abbr] = {}
            for key in ['FK', 'Fpi', 'mpi', 'mk', 'mss', 'mju', 'mjs', 'mru', 'mrs']:
                gv_data[abbr][key] = fit_data[abbr][key]
                plot_data[abbr][key] = fit_data[abbr][key][:plot_bs_N]

            gv_data[abbr] = gv.dataset.avg_data(gv_data[abbr], bstrap=True)
            gv_data[abbr]['FK/Fpi'] = gv_data[abbr]['FK'] / gv_data[abbr]['Fpi']
            plot_data[abbr]['FK/Fpi'] = plot_data[abbr]['FK'] / plot_data[abbr]['Fpi']

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


        self.bs_N = bs_N
        self.plot_bs_N = plot_bs_N
        self.abbrs = sorted(abbrs)
        #self.variable_names = sorted(fit_data[self.abbrs[0]].keys())
        self.fit_data = gv_data
        self.plot_data = plot_data
        self.prior = prior
        self.order = order
        self.fits = None
        self.chain_fits = chain_fits
        #self.bs_fit_parameters = None
        self.fit_type = fit_type

    def __str__(self):
        bs_fit_parameters = self.get_bootstrapped_fit_parameters()
        prior = self.prior
        output = "\n\nFit type: %s (F^2 = %s, bsN = %s)" %(self.fit_type, self.F2, self.bs_N)
        output = output + "\n\nFitting to %s \n" %(self.order['fit'])
        #output = output + " with lattice corrections O(%s) \n" %(self.order['latt_spacing'])
        output = output + " with volume corrections O(%s) \n" %(self.order['vol'])
        output = output + " chained: %s \n" %(self.chain_fits)
        output = output + "Fitted/[FLAG] values at physical point (including SU(2) isospin corrections):\n"
        output = output + '\nF_K / F_pi = %s [%s]'%(
                            self.extrapolate_to_phys_point(include_su2_isospin_corrrection=True),
                            self.get_phys_point_data('FK/Fpi_pm'))
        output = output + '   (delta_su2 = %s)' %(self.get_delta_su2_correction())
        output = output + "\n\n"

        fit_parameters = self.get_fit_parameters()
        new_table = { key : [fit_parameters[key], prior[key]] for key in sorted(fit_parameters.keys())}

        output = output + gv.tabulate(new_table, ncol=2, headers=['Parameter', 'Result[0] / Prior[1]'])

        output = output + '\n---\nboot0 fit results:\n'
        output = output + self.fits[0].format(pstyle=None)

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
        temp_fitter = fitter(fit_data=prepped_data, prior=temp_prior,
                        order=self.order, fit_type=self.fit_type, chain_fits=self.chain_fits)
        return temp_fitter.get_fit()

    def _make_empbayes_fit(self):
        prepped_data = self._make_fit_data(0)
        temp_prior = self._randomize_prior(self.prior, 0)
        temp_fitter = fitter(fit_data=prepped_data, prior=temp_prior,
                        order=self.order, fit_type=self.fit_type, chain_fits=self.chain_fits)
        return temp_fitter.get_empbayes_fit()

    def _make_fit_data(self, j):
        prepped_data = {}
        for parameter in ['mjs', 'mju', 'mk', 'mpi', 'mrs', 'mru', 'mss', 'lam2_chi']:
            prepped_data[parameter] = np.array([self.fit_data[abbr][parameter] for abbr in self.abbrs])
        #for parameter in ['w0']:
        #    prepped_data[parameter] = self.fit_data[abbr][parameter]
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
        print "Making fits..."
        start_time = time.time()
        self.fits = np.array([])

        for j in range(self.bs_N):
            temp_fit = self._make_fit(j)
            self.fits = np.append(self.fits, temp_fit)

            sys.stdout.write("\r{0}% complete".format(int((float(j+1)/self.bs_N)*100)))
            print "",
            sys.stdout.flush()
        end_time = time.time()
        print "Time (s): ",  end_time - start_time
        print "Compiling results..."

    def create_prior_from_fit(self):
        output = {}
        temp_parameters = self._make_empbayes_fit().prior
        for key in self.get_fit_keys():
            mean = 0
            if key in ['L_4', 'L_5']:
                #mean = gv.mean(self.get_fit_parameters(key))
                sdev = 0.001 #5 *gv.sdev(self.get_fit_parameters(key))
                output[key] = gv.gvar(mean, sdev)
            elif key in ['A_a','A_k','A_p']:
                #mean = gv.mean(self.get_fit_parameters(key))
                sdev = gv.sdev(temp_parameters[key])
                output[key] = gv.gvar(mean, sdev)
            elif key in ['A_aa', 'A_ak', 'A_ap', 'A_kk', 'A_kp', 'A_pp']:
                #mean = gv.mean(self.get_fit_parameters(key))
                sdev = gv.sdev(temp_parameters[key])
                output[key] = gv.gvar(mean, sdev)

        return output

    def extrapolate_to_ensemble(self, abbr):
        abbr_n = self.abbrs.index(abbr)
        model_name = self.fits[0].data.keys()[-1] # pretty hacky way to get this

        if self.bs_N == 1:
            temp_fit = self.fits[0]
            return temp_fit.fcn(temp_fit.p)[model_name][abbr_n]

        else:
            results = []
            for j in self.bs_N:
                temp_fit = self.fits[j]
                results.append(gv.mean(temp_fit.fcn(temp_fit.p)[model_name][abbr_n]))
            return gv.gvar(np.mean(results), np.std(results))

    def extrapolate_to_phys_point(self, include_su2_isospin_corrrection=False):
        output = self.fk_fpi_fit_fcn(self.get_phys_point_data())
        if include_su2_isospin_corrrection:
            output *= np.sqrt(1 + self.get_delta_su2_correction()) # include SU(2) isospin breaking correction

        # Logic for frequentist vs bayesian fits
        try:
            return output[0]
        except TypeError:
            return output

    def fk_fpi_fit_fcn(self, fit_data=None, fit_parameters=None, fit_type=None):
        if fit_data is None:
            fit_data = self.get_phys_point_data()
        if fit_parameters is None:
            fit_parameters = self.get_fit_parameters().copy()
        if fit_type is None:
            fit_type = self.fit_type

        for key in fit_data.keys():
            fit_parameters[key] = fit_data[key]
            temp_fit = self.fits[0]
            model_name = temp_fit.fcn(temp_fit.p).keys()[-1]

        return temp_fit.fcn(p=fit_parameters)[model_name]

    # Returns dictionary with keys fit parameters, entries bs results
    def get_bootstrapped_fit_parameters(self):
        # Make fits if they haven't been made yet
        if self.fits is None:
            self.bootstrap_fits()

        # Get all fit parameters
        keys1 = self.prior.keys()
        keys2 = self.fits[0].p.keys()
        parameters = np.intersect1d(keys1, keys2)

        # Make dictionary
        temp_fit = self.fits[0]
        bs_fit_parameters = {parameter : np.array([gv.mean(temp_fit.p[parameter])]).flatten() for parameter in parameters}

        # Populate it with bs fits
        for j in range(1, self.bs_N):
            temp_fit = self.fits[j]
            for parameter in parameters:
                bs_fit_parameters[parameter] = np.vstack((bs_fit_parameters[parameter],
                                                        np.array([gv.mean(temp_fit.p[parameter])]).flatten()))

        return bs_fit_parameters

    def get_delta_su2_correction(self):
        lam2_chi = self.get_phys_point_data('lam2_chi')
        eps2_pi = (self.get_phys_point_data('mpi'))**2 / lam2_chi
        eps2_k = (self.get_phys_point_data('mk'))**2 / lam2_chi
        fkfpi = self.extrapolate_to_phys_point()

        R = gv.gvar('35.7(2.6)') # From FLAG
        eps_su2 = np.sqrt(3)/(4.0 *R)

        delta = np.sqrt(3) *eps_su2 *(
            - (4.0/3.0) *(fkfpi - 1)
            + (4.0/3.0) *(eps2_k - eps2_pi - eps2_pi *np.log(eps2_k/eps2_pi))
        )
        return delta


    def get_fit_info(self):
        fit_type = self.fit_type+'_'+self.order['fit']
        fit_info = {
            'name' : self.fit_type+'_'+self.order['fit'],
            'fit' : self.extrapolate_to_phys_point(include_su2_isospin_corrrection=True),
            'delta_su2' : self.get_delta_su2_correction(),
            'logGBF' : self.fits[0].logGBF,
            'chi2/df' : self.fits[0].chi2 / self.fits[0].dof,
            'Q' : self.fits[0].Q,
            'vol corr' : self.order['vol'],
            #'latt corr' : self.order['latt_spacing']
        }
        for key in self.get_fit_keys():
            fit_info[key] = self.get_fit_parameters(key)

        return fit_info


    # Returns keys of fit parameters
    def get_fit_keys(self):
        if self.fits is None:
            self.bootstrap_fits()

        keys1 = self.prior.keys()
        keys2 = self.fits[0].p.keys()
        parameters = np.intersect1d(keys1, keys2)
        return sorted(parameters)

    # Returns dictionary with keys fit parameters, entries gvar results
    def get_fit_parameters(self, parameter=None):
        # If fitting only the central values
        if self.bs_N == 1:
            if self.fits is None:
                self.bootstrap_fits()

            keys1 = self.prior.keys()
            keys2 = self.fits[0].p.keys()
            parameters = np.intersect1d(keys1, keys2)

            fit_parameters = {parameter : self.fits[0].p[parameter] for parameter in parameters}
        else:
            fit_parameters = gv.dataset.avg_data(self.get_bootstrapped_fit_parameters(), bstrap=True)

        if parameter is not None:
            return fit_parameters[parameter]
        else:
            return fit_parameters

    # need to convert to/from lattice units
    def get_phys_point_data(self, parameter=None):
        phys_point_data = {
            'a/w0' : 0,
            'a' : 0,
            'L' : np.infty,
            'alpha_s' : 0, # Not sure, but not important since it comes with a^2

            'mpi' : gv.gvar('134.8(3)'), # '138.05638(37)'
            'mk' : gv.gvar('494.2(3)'), # '495.6479(92)'
            'mss' : gv.gvar('688.5(2.2)'), # Taken from arxiv/1303.1670

            'a2DI' : 0,
            'Fpi' : gv.gvar('91.9(3.5)'),
            'FK' : gv.gvar('110.38(64)'),
            'w0' : self.w0,

            'FK/Fpi_pm' : gv.gvar('1.1932(19)'), # FLAG, SU(2) isospin corrected value (arxiv/1902.08191, eqn 80)
        }
        # Or get mss, mrs with Gell-Mann-Oakes-Renner relations: arxiv/0505265 (3.45)
        mpi = phys_point_data['mpi']
        mk = phys_point_data['mk']
        phys_point_data['mss'] = np.sqrt(2 *(mk)**2 - (mpi)**2) *1.0000001 # prevents division by 0

        # ma pion
        phys_point_data['mju'] = phys_point_data['mpi']

        # ma kaon
        phys_point_data['mru'] = phys_point_data['mk']
        phys_point_data['mjs'] = phys_point_data['mk']

        # ma eta_s
        phys_point_data['mrs'] = phys_point_data['mss']

        FK = phys_point_data['FK']
        Fpi = phys_point_data['Fpi']
        if self.F2 == 'FKFpi':
            phys_point_data['lam2_chi'] = (4*np.pi)**2 *FK *Fpi
        elif self.F2 == 'FpiFpi':
            phys_point_data['lam2_chi'] = (4*np.pi)**2 *Fpi *Fpi
        elif self.F2 == 'FKFK':
            phys_point_data['lam2_chi'] = (4*np.pi)**2 *FK *FK

        if parameter is None:
            return phys_point_data
        else:
            return phys_point_data[parameter]

    def make_plots(self, show_error_ellipses=False, show_bootstrap_histograms=False):
        figs = [self.plot_fit_info()]
        figs.append(self.plot_fit_bar_graph())

        figs.append(self.plot_fit_vs_eps2pi())
        #squared = lambda x : x**2
        #figs.append(self.plot_parameters(xy_parameters=['mpi', 'FK/Fpi'],
        #            xlabel='$F_\pi$ (MeV)', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['FK', 'FK/Fpi'],
        #            xlabel='$F_K$ (MeV)', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['mpi', 'FK/Fpi'],
        #            xfcn=squared, xlabel='$m_\pi^2$ (MeV)$^2$', color_parameter='a'))
        #figs.append(self.plot_parameters(xy_parameters=['mk', 'FK/Fpi'],
        #            xfcn=squared, xlabel='$m_K^2$ (MeV)$^2$', color_parameter='a'))

        if show_bootstrap_histograms:
            # Make histograms
            fit_keys = self.get_fit_keys()
            for key in fit_keys:
                figs.append(self.plot_parameter_histogram(key))

        if show_error_ellipses:
            # Make error ellipses
            # Casts indixes of fit_keys as an upper-triangular matrix,
            # thereby allowing us to get all the 2-combinations
            # of the set of fit_keys
            fit_keys = self.get_bootstrapped_fit_parameters().keys()
            rs,cs = np.triu_indices(len(fit_keys),1)
            for (i, j) in zip(rs, cs):
                figs.append(self.plot_error_ellipsis([fit_keys[i], fit_keys[j]]))

        return figs

    # Takes an array (eg, ['l_ju', 'l_pi'])
    def plot_error_ellipsis(self, fit_keys):
        x = self.get_bootstrapped_fit_parameters()[fit_keys[0]]
        y = self.get_bootstrapped_fit_parameters()[fit_keys[1]]


        fig, ax = plt.subplots()

        corr = '{0:.3g}'.format(np.corrcoef(x, y)[0,1])
        std_x = '{0:.3g}'.format(np.std(x))
        std_y = '{0:.3g}'.format(np.std(y))
        text = ('$R_{x, y}=$ %s\n $\sigma_x =$ %s\n $\sigma_y =$ %s'
                % (corr,std_x,std_y))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        C = np.cov([x, y])
        eVe, eVa = np.linalg.eig(C)
        #print eVe, eVa
        for e, v in zip(eVe, eVa.T):
            plt.plot([np.mean(x)-1*np.sqrt(e)*v[0], 1*np.sqrt(e)*v[0] + np.mean(x)],
                     [np.mean(y)-1*np.sqrt(e)*v[1], 1*np.sqrt(e)*v[1] + np.mean(y)],
                     'k-', lw=2)

        #plt.scatter(x-np.mean(x), y-np.mean(y), rasterized=True, marker=".", alpha=100.0/self.bs_N)
        plt.scatter(x, y, rasterized=True, marker=".", alpha=100.0/self.bs_N)

        plt.grid()
        plt.gca().set_aspect('equal', adjustable='datalim')
        plt.xlabel(self._fmt_key_as_latex(fit_keys[0]), fontsize = 24)
        plt.ylabel(self._fmt_key_as_latex(fit_keys[1]), fontsize = 24)

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
        #print "----Plotting fit bar graph-----"

        y = 1
        labels = np.array([])

        # Physical point:
        for j in range(1):
            plt.axhline(y-2, ls ='-', color='C4')

            # FLAG
            data = self.get_phys_point_data('FK/Fpi_pm') / np.sqrt(1 + self.get_delta_su2_correction())
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
            data = self.fit_data[abbr]['FK'] / self.fit_data[abbr]['Fpi']
            x = gv.mean(np.mean(data))
            xerr = gv.sdev(data[0])

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

        plt.yticks(1*range(len(labels)), labels, fontsize=15, rotation=45)
        plt.ylim(-1, y)
        plt.xlabel(self._fmt_key_as_latex('FK/Fpi'), fontsize = 24)

        fig = plt.gcf()
        plt.close()

        #print "----Done plotting fit bar graph-----"

        return fig


    def plot_fit_vs_eps2pi(self):

        #print "-----Plotting fit----"
        # used to convert to phys units
        hbar_c = 197.327

        plot_data = {}
        plot_data[0] = {abbr : ((self.plot_data[abbr]['mpi'] *hbar_c / (self.plot_data[abbr]['a/w0'] *self.w0))**2 /self.get_phys_point_data('lam2_chi'))
                          for abbr in self.abbrs}
        plot_data[1] = {abbr : self.shift_fk_fpi_for_phys_mk(abbr, self.plot_data) for abbr in self.abbrs}

        color_data = {abbr : self.plot_data[abbr]['a/w0'] *self.w0 for abbr in self.abbrs}

        # Plot fit
        #colors = ['black', 'purple', 'green', 'red']
        colors = ['red', 'green', 'purple', 'black']
        lattice_spacings = np.unique(self._make_fit_data(0)['a/w0']) *self.w0
        for j, a in enumerate(sorted(np.append([gv.gvar('0(0)')], lattice_spacings), reverse=True)):

            # Get the range of pion masses (in phys units)
            minimum = np.nanmin([np.nanmin(
                np.sqrt([plot_data[0][abbr] *self.get_phys_point_data('lam2_chi') for abbr in self.abbrs])
            ) for abbr in self.abbrs])
            maximum = np.nanmax([np.nanmax(
                np.sqrt([plot_data[0][abbr] *self.get_phys_point_data('lam2_chi')  for abbr in self.abbrs])
            ) for abbr in self.abbrs])
            minimum = gv.mean(minimum)
            maximum = gv.mean(maximum)
            delta = maximum - minimum

            x = np.linspace(np.max((minimum - 0.05*delta, 0)), maximum + 0.05*delta)

            # Get phys point data, substituting x-data and current 'a' in loop
            prepped_data = self.get_phys_point_data()
            prepped_data['mpi'] = x
            prepped_data['a/w0'] = a/self.w0

            # Covert m_pi -> (eps_pi)^2
            x = x**2/self.get_phys_point_data('lam2_chi')
            y = self.fk_fpi_fit_fcn(fit_data=prepped_data)

            pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
            plt.plot(pm(x, 0), pm(y, 0), '--', color=colors[j], label='$a=$%s (fm)'%(str(a)), rasterized=True)
            plt.fill_between(pm(x, 0), pm(y, -1), pm(y, 1), alpha=0.40, color=colors[j], rasterized=True)

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        min_max = lambda x : [np.min(x), np.max(x)]
        minimum, maximum = min_max([gv.mean(color_data[abbr]) for abbr in self.abbrs])
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        length = self.plot_bs_N
        x = np.zeros(length * len(self.abbrs))
        y = np.zeros(length * len(self.abbrs))
        z = np.zeros(length * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*length:(j+1)*length] = gv.mean(plot_data[0][abbr])
            y[j*length:(j+1)*length] = gv.mean(plot_data[1][abbr])
            z[j*length:(j+1)*length] = gv.mean(color_data[abbr])

        # Plot data
        sc = plt.scatter(x, y, c=z, vmin=minimum, vmax=maximum,
                         cmap=cmap, rasterized=True, marker="o", alpha=10.0/self.plot_bs_N, edgecolor='black')

        # Plot FLAG result
        x_phys = self.get_phys_point_data('mpi')**2 / self.get_phys_point_data('lam2_chi')
        plt.axvline(gv.mean(x_phys), label='Phys point')
        y_phys = self.get_phys_point_data('FK/Fpi_pm') / np.sqrt(1 + self.get_delta_su2_correction())

        plt.errorbar(x=gv.mean(x_phys), xerr=0,
                     y=gv.mean(y_phys), yerr=gv.sdev(y_phys), label='FLAG',
                    color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)

        # Plot labels
        plt.legend()
        plt.grid()
        plt.xlabel('$\epsilon_\pi^2$', fontsize = 24)
        plt.ylabel('$F_K/F_\pi$', fontsize = 24)

        # Format colorbar
        color_bar = plt.colorbar(sc)
        color_bar.set_alpha(0.8)
        color_bar.draw_all()
        color_bar.set_label('$a$ (fm)', fontsize = 24)

        # Set xlim, ylim -- only works if xy_parameters[i] is a vector, not a scalar
        min_max = lambda x : [np.min(x), np.max(x)]
        try:
            xmin, xmax = min_max(np.concatenate([gv.mean(plot_data[0][abbr]) for abbr in self.abbrs]))
            ymin, ymax = min_max(np.concatenate([gv.mean(plot_data[1][abbr]) for abbr in self.abbrs]))
            xdelta = xmax - xmin
            ydelta = ymax - ymin
            plt.xlim(xmin-0.05*xdelta, xmax+0.05*xdelta) #xmin-0.05*xdelta
            #plt.ylim(ymin-0.05*ydelta, ymax+0.05*ydelta)
            plt.ylim(1.04, 1.20)
        except ValueError:
            pass

        fig = plt.gcf()
        plt.close()

        #print "----Done plotting fit----"
        return fig

    def plot_parameter_histogram(self, parameter):
        data = self.get_bootstrapped_fit_parameters()[parameter]
        fig, ax = plt.subplots()
        n, bins, patches = plt.hist(data, self.bs_N/10, normed=True, facecolor='green', alpha=0.75)

        mu = np.mean(data)
        sigma = np.std(data)

        # Overlay a gaussian
        y = matplotlib.mlab.normpdf(bins, mu, sigma)
        l = plt.plot(bins, y, 'r--', linewidth=1)

        plt.title("BS dist", fontsize=30)
        plt.xlabel("$p = $"+self._fmt_key_as_latex(parameter), fontsize=24)
        plt.ylabel('Frequency', fontsize=24)
        plt.grid(True)

        text = ('$\overline{p} (s_{\overline{p}}) = $ %s \n Prior: %s'
                % (gv.gvar(mu, sigma), self.prior[parameter]))

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)

        fig = plt.gcf()
        plt.close()

        return fig

    def plot_parameters(self, xy_parameters, color_parameter=None,
                        xfcn=None, xlabel=None, yfcn=None, ylabel=None):

        # used to convert to phys units
        hbar_c = 197.327

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


    def shift_fk_fpi_for_phys_mk_alt(self, abbr, data):
        pi = np.pi
        hbar_c = 197.327
        lam2_chi = self.get_phys_point_data('lam2_chi')

        fkfpi_ens = data[abbr]['FK'] / data[abbr]['Fpi']
        fkfpi_fit_ens = self.extrapolate_to_ensemble(abbr)

        temp_data = self.get_phys_point_data()
        temp_data['mpi'] = data[abbr]['mpi'] *hbar_c / (data[abbr]['a/w0'] *self.w0)
        temp_data['a/w0'] = data[abbr]['a/w0']
        fkfpi_fit_phys_vary_mpi = self.fk_fpi_fit_fcn(fit_data=temp_data)

        shifted_fkfpi = fkfpi_ens *(fkfpi_fit_phys_vary_mpi / fkfpi_fit_ens)

        return shifted_fkfpi

    def shift_fk_fpi_for_phys_mk(self, abbr, data):
        pi = np.pi
        hbar_c = 197.327
        lam2_chi = self.get_phys_point_data('lam2_chi')

        fkfpi_ens = data[abbr]['FK'] / data[abbr]['Fpi']
        fkfpi_fit_ens = self.extrapolate_to_ensemble(abbr)

        temp_data = self.get_phys_point_data()
        temp_data['mpi'] = data[abbr]['mpi'] *hbar_c / (data[abbr]['a/w0'] *self.w0)
        temp_data['a/w0'] = data[abbr]['a/w0']
        fkfpi_fit_phys_vary_mpi = self.fk_fpi_fit_fcn(fit_data=temp_data)

        shifted_fkfpi = fkfpi_ens + fkfpi_fit_phys_vary_mpi - fkfpi_fit_ens

        return shifted_fkfpi
