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
                 order=None, fit_type=None, F2=None):

        w0 = gv.gvar('0.17(10)') # Still needs to be determined, but something like this
        self.w0 = w0

        if fit_type is None:
            fit_type = 'ma-taylor'

        if F2 is None:
            F2 = 'FKFpi'
        self.F2 = F2

        if bs_N is None or bs_N==0:
            bs_N = len(fit_data[fit_data.keys()[0]]['mpi'])

        if order is None:
            order = {
                'fit' : 'nlo',
                'vol' : 1
            }

        if prior is None:
            print "Using default prior."
            prior = {
                # nlo terms
                'L_5' : '0.0002(1)' , #'0.000234 (100) ', #'0.00153 (12)', #'0.00(1)', #0.0002(1)
                'L_4' : '-0.0002(2) ', #'-0.000710 (24)', #'0.00(1)',

                # nnlo terms
                'A_a' : '0(1)',
                'A_x' : '0(1)',
                'A_k' : '0(1)',
                'A_p' : '0(1)',

                # nnnlo terms
                'A_aa' : '0(1)',
                'A_ax' : '0(1)',
                'A_ak' : '0(1)',
                'A_ap' : '0(1)',
                'A_xx' : '0(1)',
                'A_xk' : '0(1)',
                'A_xp' : '0(1)',
                'A_kk' : '0(1)',
                'A_kp' : '0(1)',
                'A_pp' : '0(1)',
            }
            prior = gv.gvar(prior)

        if abbrs is None:
            abbrs = fit_data.keys()

        # Resize array data to bs_N
        data = {}
        for abbr in abbrs:
            data[abbr] = {}
            for data_parameter in fit_data[abbr].keys():
                if data_parameter in ['Fpi', 'FK', 'mpi', 'mk', 'mss', 'mju', 'mjs', 'mru', 'mrs']:
                    means = fit_data[abbr][data_parameter][:bs_N]
                    unc = np.std(fit_data[abbr][data_parameter])
                    #print gv.gvar(means, np.repeat(unc, len(means)))
                    data[abbr][data_parameter] = gv.gvar(means, np.repeat(unc, len(means)))
                elif data_parameter in ['a2DI']:
                    to_gvar = lambda arr : gv.gvar(arr[0], arr[1])
                    data[abbr][data_parameter] = to_gvar(fit_data[abbr][data_parameter]) #np.repeat(to_gvar(fit_data[abbr][data_parameter]), bs_N)
                elif data_parameter in ['aw0']: # this is a/w0
                    to_gvar = lambda arr : gv.gvar(arr[0], arr[1])
                    #data[abbr]['w0'] = w0 #np.repeat(w0, bs_N)
                    #data[abbr]['a'] = to_gvar(fit_data[abbr]['aw0'])/w0 #np.repeat(to_gvar(fit_data[abbr][data_parameter])/w0, bs_N)
                    data[abbr]['a/w0'] = to_gvar(fit_data[abbr]['aw0'])
                elif data_parameter in ['MpiL']:
                    L = int(np.median(fit_data[abbr]['MpiL'] / fit_data[abbr]['mpi']))
                    data[abbr]['L'] = gv.gvar(L, L/100000.0) # Pretty certain about this value....


        for abbr in abbrs:
            hbar_c = 197.327
            a = data[abbr]['a/w0'] *w0
            data[abbr]['lam2_chi'] = self.get_phys_point_data('lam2_chi') *(a /hbar_c)**2

        self.bs_N = bs_N
        self.abbrs = sorted(abbrs)
        #self.variable_names = sorted(fit_data[self.abbrs[0]].keys())
        self.fit_data = data
        self.prior = prior
        self.order = order
        self.fits = None
        #self.bs_fit_parameters = None
        self.fit_type = fit_type

    def __str__(self):
        bs_fit_parameters = self.get_bootstrapped_fit_parameters()
        prior = self.prior
        output = "\n\nFit type: %s (F^2 = %s, bsN = %s)" %(self.fit_type, self.F2, self.bs_N)
        output = output + "\n\nFitting to %s \n" %(self.order['fit'])
        #output = output + " with lattice corrections O(%s) \n" %(self.order['latt_spacing'])
        output = output + " with volume corrections O(%s) \n" %(self.order['vol'])
        output = output + "Fitted/[FLAG] values at physical point:\n"
        output = output + '\nF_K / F_pi = %s [%s]\n'%(
                            self.extrapolate_to_phys_point(),
                            self.get_phys_point_data('FK/Fpi'))
        output = output + "\n"

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
                        order=self.order, fit_type=self.fit_type)
        return temp_fitter.get_fit()

    def _make_fit_data(self, j):
        prepped_data = {}
        for parameter in ['FK', 'Fpi', 'mjs', 'mju', 'mk', 'mpi', 'mrs', 'mru', 'mss']:
            prepped_data[parameter] = np.array([self.fit_data[abbr][parameter][j] for abbr in self.abbrs])
        #for parameter in ['w0']:
        #    prepped_data[parameter] = self.fit_data[abbr][parameter]
        for parameter in ['a/w0', 'a2DI', 'L', 'lam2_chi',]:
            prepped_data[parameter] =  np.array([self.fit_data[abbr][parameter] for abbr in self.abbrs])


        y = ([self.fit_data[abbr]['FK'][j]/self.fit_data[abbr]['Fpi'][j] for abbr in self.abbrs])
        prepped_data['y'] = y

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
        for key in self.get_fit_parameters().keys():
            mean = gv.mean(self.get_fit_parameters(key))
            sdev = gv.sdev(self.get_fit_parameters(key))
            unc = np.min((0.5, 5*sdev)) #0.25*gv.mean(self.get_fit_parameters()[key])

            output[key] = gv.gvar(mean, unc)

        return output

    def extrapolate_to_ensemble(self, abbr):
        to_gvar = lambda x : gv.gvar(np.median(gv.mean(x)), np.median(gv.sdev(x)))
        return to_gvar(self.fk_fpi_fit_fcn(self.fit_data[abbr]))

    def extrapolate_to_phys_point(self):
        output = self.fk_fpi_fit_fcn(self.get_phys_point_data())

        # Logic for frequentist vs bayesian fits
        try:
            return output[0]
        except TypeError:
            return output

    def fk_fpi_fit_fcn(self, fit_data=None, fit_parameters=None, fit_type=None):
        if fit_data is None:
            fit_data = self.get_phys_point_data()
        if fit_parameters is None:
            fit_parameters = self.get_fit_parameters()
        if fit_type is None:
            fit_type = self.fit_type

        model = fitter(order=self.order, fit_type=fit_type)._make_models()[-1]
        #print "model", model.datatag
        return model.fitfcn(p=fit_parameters, fit_data=fit_data)

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

    def get_fit_info(self):
        fit_type = self.fit_type+'_'+self.order['fit']
        fit_info = {
            'name' : self.fit_type+'_'+self.order['fit'],
            'fit' : self.extrapolate_to_phys_point(),
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

            'mpi' : gv.gvar('138.05638(37)'),
            'mju' : gv.gvar('138.05638(37)'),

            'mk' : gv.gvar('495.6479(92)'),
            'mru' : gv.gvar('495.6479(92)'),
            'mjs' : gv.gvar('495.6479(92)'),

            'mss' : gv.gvar('688.5(2.2)'), # arxiv/1303.1670
            'mrs' : gv.gvar('688.5(2.2)'),

            'a2DI' : 0, # Need to check this
            'Fpi' : gv.gvar('91.9(3.5)'),
            'FK' : gv.gvar('110.38(64)'),
            'w0' : self.w0,

            'FK/Fpi' : gv.gvar('1.1932(19)') # FLAG
        }

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
                         ecolor='C1', elinewidth=10.0, label='FLAG/Lat')
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

        return fig


    def plot_fit_vs_eps2pi(self):
        # used to convert to phys units
        hbar_c = 197.327

        plot_data = {
            0 :  {abbr : (self.fit_data[abbr]['mpi'] *hbar_c / (self.fit_data[abbr]['a/w0'] *self.w0)**2 /self.get_phys_point_data('lam2_chi'))
                              for abbr in self.abbrs},
            1 : {abbr : self.shift_fk_fpi_for_phys_mk(abbr) for abbr in self.abbrs}
        }

        color_data = {abbr : np.repeat(self.fit_data[abbr]['a/w0'] *self.w0, self.bs_N).ravel() for abbr in self.abbrs}

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        min_max = lambda x : [np.min(x), np.max(x)]
        minimum, maximum = min_max(np.concatenate([gv.mean(color_data[abbr]) for abbr in self.abbrs]))
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        x = np.zeros(self.bs_N * len(self.abbrs))
        y = np.zeros(self.bs_N * len(self.abbrs))
        z = np.zeros(self.bs_N * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[0][abbr])
            y[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[1][abbr])
            z[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(color_data[abbr])

        # Plot data
        sc = plt.scatter(x, y, c=z, vmin=minimum, vmax=maximum,
                         cmap=cmap, rasterized=True, marker=".", alpha=100.0/self.bs_N)


        # Plot FLAG result
        x_phys = self.get_phys_point_data('mpi')**2 / self.get_phys_point_data('lam2_chi')
        plt.axvline(gv.mean(x_phys), label='Phys point')
        y_phys = self.get_phys_point_data('FK/Fpi')

        plt.errorbar(x=gv.mean(x_phys), xerr=0,
                     y=gv.mean(y_phys), yerr=gv.sdev(y_phys), label='FLAG',
                    color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)

        # Plot fit
        colors = ['black', 'purple', 'green', 'red']
        lattice_spacings = np.unique(self._make_fit_data(0)['a/w0']) *self.w0
        for j, a in enumerate(sorted(np.append([gv.gvar('0(0)')], lattice_spacings))):

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
            prepped_data['a/w0'] = a *self.w0

            # Covert m_pi -> (eps_pi)^2
            x = x**2/self.get_phys_point_data('lam2_chi')
            y = self.fk_fpi_fit_fcn(fit_data=prepped_data)

            pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
            plt.plot(pm(x, 0), pm(y, 0), '--', color=colors[j], label='$a=$%s (fm)'%(str(a)), rasterized=True)
            plt.fill_between(pm(x, 0), pm(y, -1), pm(y, 1), alpha=0.20, color=colors[j], rasterized=True)

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
                plot_data[j] = {abbr :  myfcn[j](self.fit_data[abbr]['FK'] /self.fit_data[abbr]['Fpi']) for abbr in self.abbrs}
            elif parameter in ['mpi', 'mju', 'mru', 'mk', 'mrs', 'mss', 'FK', 'Fpi']:
                # Convert to physical units

                plot_data[j] = {}
                for abbr in self.abbrs:
                    plot_data[j][abbr] = myfcn[j](self.fit_data[abbr][parameter] *hbar_c / (self.fit_data[abbr]['a/w0'] *self.w0))
            else:
                plot_data[j] = {abbr :  myfcn[j](self.fit_data[abbr][parameter]) for abbr in self.abbrs}


        # Get data for color coding graph
        if color_parameter is None:
            color_parameter = 'a'

        if color_parameter in ['a']:
            color_data = {abbr : np.repeat(self.fit_data[abbr]['a/w0'] *self.w0, self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter in ['L']:
            color_data = {abbr : np.repeat(gv.mean(self.fit_data[abbr][color_parameter]), self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter == 'MpiL':
            color_data = {abbr : gv.mean(self.fit_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}
        else:
            color_data = {abbr : gv.mean(self.fit_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        min_max = lambda x : [np.min(x), np.max(x)]
        minimum, maximum = min_max(np.concatenate([gv.mean(color_data[abbr]) for abbr in self.abbrs]))
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        x = np.zeros(self.bs_N * len(self.abbrs))
        y = np.zeros(self.bs_N * len(self.abbrs))
        z = np.zeros(self.bs_N * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[0][abbr])
            y[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[1][abbr])
            z[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(color_data[abbr])

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


    def shift_fk_fpi_for_phys_mk_alt(self, abbr):
        hbar_c = 197.327
        lam2_chi = self.get_phys_point_data('lam2_chi')

        fkfpi_data = self.fit_data[abbr]['FK'] / self.fit_data[abbr]['Fpi']
        fkfpi_fit_mk_phys = self.extrapolate_to_phys_point()

        data = self.get_phys_point_data()
        data['mk'] = self.fit_data[abbr]['mk'] * hbar_c / (self.fit_data[abbr]['a/w0'] *self.w0)
        fkfpi_fit_mk_ens = self.fk_fpi_fit_fcn(data)

        return fkfpi_data *(fkfpi_fit_mk_phys / fkfpi_fit_mk_ens)

    def shift_fk_fpi_for_phys_mk(self, abbr):
        pi = np.pi
        hbar_c = 197.327
        lam2_chi = self.get_phys_point_data('lam2_chi')

        fkfpi_ens = self.fit_data[abbr]['FK'] / self.fit_data[abbr]['Fpi']

        to_phys = lambda m : m *hbar_c / (self.fit_data[abbr]['a/w0'] *self.w0)
        eps2_k = (to_phys(self.fit_data[abbr]['mk']))**2 / lam2_chi
        eps2_k_phys = (self.get_phys_point_data('mk'))**2 / lam2_chi
        eps2_pi = (to_phys(self.fit_data[abbr]['mpi']))**2 / lam2_chi
        eps2_pi_phys = (self.get_phys_point_data('mpi'))**2 / lam2_chi

        to_eta = lambda xk, xpi : (4.0/3.0) *xk - (1.0/3.0) *xpi
        eps2_x = to_eta(eps2_k, eps2_pi)
        eps2_x_phys = to_eta(eps2_k_phys, eps2_pi_phys)

        L_5 = self.get_fit_parameters('L_5')

        delta = lambda xphys, xens : (xphys *np.log(xphys) - xens *np.log(xens))
        shifted_fkfpi = (
            fkfpi_ens
            - (1.0/4.0) *delta(eps2_k_phys, eps2_k)
            - (3.0/8.0) *delta(eps2_x_phys, eps2_x)
            + 4 *(4*pi)**2 *(eps2_k_phys - eps2_k) *L_5
        )

        return shifted_fkfpi
