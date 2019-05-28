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

    def __init__(self, fit_data, prior, bs_N=None, order=None):

        if bs_N is None:
            bs_N = len(fit_data[fit_data.keys()[0]]['mpi'])

        if order is None:
            order = 0

        # Resize array data to bs_N
        for abbr in fit_data.keys():
            for data_parameter in fit_data[abbr].keys():
                try:
                    fit_data[abbr][data_parameter] = fit_data[abbr][data_parameter][:bs_N]
                except IndexError:
                    pass
                if data_parameter == 'a2DI':
                    temp = fit_data[abbr][data_parameter]
                    fit_data[abbr][data_parameter] = np.random.normal(temp[0], temp[1], bs_N)

        #abbrs = sorted(fit_data.keys())
        #to_gvar = lambda x : gv.gvar(x[0], x[1])
        #prior['a2DI'] = [to_gvar(fit_data[abbr]['a2DI']) for abbr in abbrs]

        self.bs_N = bs_N
        self.abbrs = sorted(fit_data.keys())
        self.variable_names = sorted(fit_data[self.abbrs[0]].keys())
        self.fit_data = fit_data
        self.prior = prior
        self.order = order

        self.w0 = 5.81743
        self.fits = None
        self.bs_fit_parameters = None
        self.fit_parameters = None

    def __str__(self):
        bs_fit_parameters = self.get_bootstrap_parameters()
        prior = self.prior

        output = "\n\nFitting to order %s \n" %(self.order)
        output = output + "Fitted/[Experimental] values at physical point:\n"
        output = output + '\n\t\tF_K / F_pi = %s [%s]\n'%(
                            self.fk_fpi_fit_fcn(self.get_phys_point_data()),
                            self.get_phys_point_data('FK')/self.get_phys_point_data('Fpi'))
        output = output + "\n"

        table = gv.dataset.avg_data(bs_fit_parameters, bstrap=True)
        new_table = { key : [table[key], prior[key]] for key in table.keys()}
        output = output + gv.tabulate(new_table, ncol=2, headers=['Parameter', 'Result[0] / Prior[1]'])

        return output

    # Returns dictionary with keys fit parameters, entries bs results
    def get_bootstrap_parameters(self):
        bs_fit_parameters = self.bs_fit_parameters
        if self.bs_fit_parameters is None:

            # Make fits if they haven't been made yet
            if self.fits is None:
                self.bootstrap_fits()

            # Get all fit parameters
            parameters = self.fits[0].p.keys()

            # Make dictionary if it doesn't exist yet
            if bs_fit_parameters is None:
                bs_fit_parameters = {parameter : np.array([]) for parameter in parameters}

            print "Getting bs fit parameters..."
            for j in range(self.bs_N):
                temp_fit = self.fits[j]
                for parameter in parameters:
                    bs_fit_parameters[parameter] = np.append(bs_fit_parameters[parameter], gv.mean(temp_fit.p[parameter]))

                sys.stdout.write("\r{0}% complete".format(int((float(j+1)/self.bs_N)*100)))
                print "",
                sys.stdout.flush()

            self.bs_fit_parameters = bs_fit_parameters
            return bs_fit_parameters
        else:
            return bs_fit_parameters

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


    def _make_fit(self, j):
        prepped_data = self._make_fit_data(j)
        # Need to randomize prior in bayesian-bootstrap hybrid
        temp_prior = self._randomize_prior(self.prior, j)
        temp_fitter = fitter(fit_data=prepped_data, prior=temp_prior, order=self.order)
        return temp_fitter.get_fit()

    def _make_fit_data(self, j):
        prepped_data = {}
        for parameter in ['mpi', 'mk', 'mju', 'mrs', 'mru', 'mss', 'a2DI', 'Fpi', 'FK']:
            prepped_data[parameter] = np.array([np.asscalar(self.fit_data[abbr][parameter][j]) for abbr in self.abbrs])

        y_mean = ([np.mean(self.fit_data[abbr]['FK']/self.fit_data[abbr]['Fpi'])
                   for abbr in self.abbrs])
        y_cov = np.diagflat([np.var(self.fit_data[abbr]['FK']/self.fit_data[abbr]['Fpi'])
                             for abbr in self.abbrs])
        prepped_data['y'] = gv.gvar(y_mean, y_cov)

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

    # need to convert to/from lattice units
    def get_phys_point_data(self, parameter=None):
        phys_point_data = {
            'a' : 0,
            'V' : 0,

            'mpi' : gv.gvar('138.05638(37)'),
            'mju' : gv.gvar('138.05638(37)'),

            'mk' : gv.gvar('495.6479(92)'),
            'mru' : gv.gvar('495.6479(92)'),

            'a2DI' : 0, # Need to check this
            'Fpi' : gv.gvar('91.9(3.5)'),
            'FK' : gv.gvar('110.38(64)')
        }
        phys_point_data['mss'] = np.sqrt((4*phys_point_data['mk']**2 - phys_point_data['mpi']**2)/3.0)
        phys_point_data['mss'] = phys_point_data['mss'] * 1.00000001
        phys_point_data['mrs'] = phys_point_data['mss']
        if parameter is None:
            return phys_point_data
        else:
            return phys_point_data[parameter]

    def fk_fpi_fit_fcn(self, fit_data=None, fit_parameters=None):
        if fit_data is None:
            fit_data = self.get_phys_point_data()
        if fit_parameters is None:
            fit_parameters = self.get_fit_parameters()

        model = fitter(order=self.order)._make_models(fit_data)[0]
        return model.fitfcn(p=fit_parameters, fit_data=fit_data)

    def extrapolate_to_phys_point(self):
        return self.fk_fpi_fit_fcn(self.get_phys_point_data())

    def extrapolate_to_ensemble(self, abbr):
        to_gvar = lambda x : gv.gvar(np.median(gv.mean(x)), np.median(gv.sdev(x)))
        return to_gvar(self.fk_fpi_fit_fcn(self.fit_data[abbr]))

    def _fmt_key_as_latex(self, key):
        convert = {
            # data parameters
            'a' : r'$a$ (fm)',
            'L' : r'$L$ (fm)',

            'mpi' : r'$m_\pi$ (MeV)',
            'Fpi' : r'$F_\pi$ (MeV)',

            'mk' : r'$m_K$ (MeV)',
            'FK' : r'$F_K$ (MeV)',

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

    def plot_fit_bar_graph(self):
        y = 0
        labels = np.array([])
        for abbr in self.abbrs:
            # data
            data = self.fit_data[abbr]['FK'] / self.fit_data[abbr]['Fpi']
            x = np.mean(data)
            xerr = np.std(data)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)

            labels = np.append(labels, str(""))
            y = y + 1

            # fit result
            fit_value = self.extrapolate_to_ensemble(abbr)
            x = gv.mean(fit_value)
            xerr = gv.sdev(fit_value)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C2', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C3', elinewidth=10.0)

            y = y + 1
            labels = np.append(labels, str(abbr))
            plt.axhline(y, ls='--')

            y = y + 1
            labels = np.append(labels, str(""))

        # Physical point:
        for j in range(1):
            plt.axhline(y-1, ls ='-', color='C4')
            data = self.get_phys_point_data('FK') / self.get_phys_point_data('FK')
            x = gv.mean(data)
            xerr = gv.sdev(data)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C0', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C1', elinewidth=10.0)

            labels = np.append(labels, str(""))
            y = y + 1

            # fit result
            fit_value = self.extrapolate_to_phys_point()
            x = gv.mean(fit_value)
            xerr = gv.sdev(fit_value)
            plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                         color='C2', marker='o', capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                         ecolor='C3', elinewidth=10.0)

            y = y + 1
            labels = np.append(labels, str("Phys. point"))
            plt.axhline(y, ls='--')

            y = y + 1
            labels = np.append(labels, str(""))
            plt.axhline(y-1, ls ='-', color='C4')


        plt.yticks(1*range(len(labels)), labels, fontsize=15, rotation=45)
        plt.ylim(-1, y)
        plt.xlabel(self._fmt_key_as_latex('FK/Fpi'), fontsize = 24)

        fig = plt.gcf()
        plt.close()

        return fig


    def plot_parameter_histogram(self, parameter):
        data = self.get_bootstrap_parameters()[parameter]
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


    # xy_parameters is an array (eg, ['mpi', 'Fpi'])
    def plot_parameters(self, xy_parameters, color_parameter=None,
                        xfcn=None, xlabel=None, yfcn=None, ylabel=None):
        if xfcn is None:
            xfcn = lambda x : 1 * x
            xlabel = self._fmt_key_as_latex(xy_parameters[0])
        else:
            # Must also specify xlabel
            xlabel = xlabel

        if yfcn is None:
            yfcn = lambda y : 1 * y
            ylabel = self._fmt_key_as_latex(xy_parameters[1])
        else:
            # Must also specify ylabel
            ylabel = ylabel

        plot_data = {}
        myfcn = [xfcn, yfcn]
        for j, parameter in enumerate(xy_parameters):
            if parameter in ['FK/Fpi', 'FK / Fpi']:
                plot_data[j] = {abbr :  myfcn[j](self.fit_data[abbr]['FK'] /self.fit_data[abbr]['Fpi']) for abbr in self.abbrs}
            else:
                plot_data[j] = {abbr :  myfcn[j](self.fit_data[abbr][parameter]) for abbr in self.abbrs}


        if color_parameter is None:
            color_parameter = 'a'

        if color_parameter in ['a']:
            color_data = {abbr : np.repeat(self.fit_data[abbr]['aw0'][0]/self.w0, self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter in ['L']:
            color_data = {abbr : np.repeat(gv.mean(self.fit_data[abbr][color_parameter]), self.bs_N).ravel() for abbr in self.abbrs}
        elif color_parameter == 'V':
            color_data = {abbr : gv.mean(self.fit_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}
        else:
            color_data = {abbr : gv.mean(self.fit_data[abbr][color_parameter]).ravel() for abbr in self.abbrs}

        # Color by lattice spacing/length
        cmap = matplotlib.cm.get_cmap('rainbow')
        min_max = lambda x : [np.min(x), np.max(x)]
        minimum, maximum = min_max(np.concatenate([color_data[abbr] for abbr in self.abbrs]))
        norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)

        # Get scatter plot & color data
        x = np.zeros(self.bs_N * len(self.abbrs))
        y = np.zeros(self.bs_N * len(self.abbrs))
        z = np.zeros(self.bs_N * len(self.abbrs))
        for j, abbr in enumerate(self.abbrs):
            x[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[0][abbr])
            y[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(plot_data[1][abbr])
            z[j*self.bs_N:(j+1)*self.bs_N] = color_data[abbr]

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
        '''
        min_max = lambda x : [np.min(x), np.max(x)]
        try:
            xmin, xmax = min_max(np.concatenate([gv.mean(plot_data[0][abbr]) for abbr in self.abbrs]))
            ymin, ymax = min_max(np.concatenate([gv.mean(plot_data[1][abbr]) for abbr in self.abbrs]))
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
        except ValueError:
            pass'''

        fig = plt.gcf()
        plt.close()
        return fig

    # Returns dictionary with keys fit parameters, entries gvar results
    def get_fit_parameters(self):
        if self.fit_parameters is None:
            return gv.dataset.avg_data(self.get_bootstrap_parameters(), bstrap=True)
        else:
            return self.fit_parameters

    # Returns keys of fit parameters
    def get_fit_keys(self):
        return sorted(self.get_bootstrap_parameters().keys())

    #################################################
    #### EVERYTHING BELOW INCOMPATIBLE WITH CLASS ###
    #################################################














    # Numerically performs a 1D interpolation
    #   f_pi = f_pi[m_pi]
    # at a given lattice spacing.
    # Used for plotting fits.
    def interpolate_data_by_lattice_spacing(self, a):
        m_pi = np.array([])
        f_pi = np.array([])
        for abbr in self.abbrs:
            # Use abbr name to get approx lattice spacing
            if abbr[1:3] == a:
                m_pi = np.append(m_pi, np.mean(self.fit_data[abbr]['m_pion']))
                f_pi = np.append(f_pi, np.mean(self.fit_data[abbr]['f_pion']))

        # Must sort x variable (ie, m_pi)
        # for interpolate to work
        index = np.argsort(m_pi)
        m_pi = m_pi[index]
        f_pi = f_pi[index]

        m_pi_out = np.arange(m_pi[0], m_pi[-1])
        output_dict = {
            'V' : 0,
            'm_pion' : m_pi_out,
            'f_pion' : np.interp(m_pi_out, m_pi, f_pi),
            'a' : int(a)/100.0
        }
        return output_dict





    # Takes an array (eg, ['m_omega_0', 'beta_6_omega'])
    def plot_error_ellipsis(self, fit_keys):
        x = self.get_bootstrap_parameters()[fit_keys[0]]
        y = self.get_bootstrap_parameters()[fit_keys[1]]


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

    # Plots m_omega vs m_pion with fits for each lattice spacing
    # xfcn acts on x variable in plot; eg, if squared = lambda x : x**2,
    # then squared is a lambda function that squares each pion mass
    def plot_fit(self, show_plot=True, xfcn=None, xlabel=None, title=None):

        if xfcn is None:
            xfcn = lambda x : 1 * x
            xlabel = self._fmt_key_as_latex('m_pion')
        else:
            # Must also specify xlabel
            xlabel = xlabel

        figs = np.array([])
        for m_hyperon in self.hyperon_mass_keys:
            # Plot data points
            # Color by lattice spacing
            cmap = matplotlib.cm.get_cmap('rainbow')
            min_max = lambda x : [np.min(x), np.max(x)]
            minimum, maximum = min_max([gv.mean(self.fit_data[abbr]['a']) for abbr in self.abbrs])
            norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)
            color_data = {abbr : np.repeat(gv.mean(self.fit_data[abbr]['a']), self.bs_N).ravel() for abbr in self.abbrs}

            # Get scatter plot & color data
            x = np.zeros(self.bs_N * len(self.abbrs))
            y = np.zeros(self.bs_N * len(self.abbrs))
            z = np.zeros(self.bs_N * len(self.abbrs))
            for j, abbr in enumerate(self.abbrs):
                x[j*self.bs_N:(j+1)*self.bs_N] = xfcn(gv.mean(self.fit_data[abbr]['m_pion']).ravel())
                y[j*self.bs_N:(j+1)*self.bs_N] = gv.mean(self.fit_data[abbr][m_hyperon]).ravel()
                z[j*self.bs_N:(j+1)*self.bs_N] = color_data[abbr]

            sc = plt.scatter(x, y, c=z, vmin=minimum, vmax=maximum,
                             cmap=cmap, rasterized=True, marker=".", alpha=100.0/self.bs_N)

            # Plot fits
            # Get different approx lattice spacings
            # These are taken from the abbr name -- there aren't the actual values!
            # If we were to use the actual values, we would need to round to anyway;
            # using the abbr name is a shortcut.
            lattice_spacings = np.unique([abbr[1:3] for abbr in self.abbrs])
            for a in lattice_spacings:
                interpolated_data = self.interpolate_data_by_lattice_spacing(a)
                interpolated_results = self.mass_fit_fcn(interpolated_data) # Generates m_omega, eg
                for m_hyperon in interpolated_results.keys():
                    interpolated_data[m_hyperon] = interpolated_results[m_hyperon]

                x = xfcn(interpolated_data['m_pion'])
                y = interpolated_data[m_hyperon]

                pm = lambda x, k : gv.mean(x) + k*gv.sdev(x)
                color = cmap(norm(int(a)/100.0)) # Convert string to int
                plt.plot(x, pm(y, 0), '--', color=color, label='$a=$'+a)
                plt.plot(x, pm(y, -1), x, pm(y, 1), color=color)
                plt.fill_between(x, pm(y, -1), pm(y, 1), facecolor=color, alpha = 0.10)

            # Draw line at physical value of l_exp
            plt.axvline(xfcn(gv.mean(self.get_phys_point_data('m_pion'))),
                        linestyle = '--', color='black', alpha=0.3, label="Physical point")

            x = xfcn(self.get_phys_point_data('m_pion'))
            plt.errorbar(x=gv.mean(x),
                         xerr=gv.sdev(x),
                         y=gv.mean(self.get_phys_point_data(m_hyperon)),
                         yerr=gv.sdev(self.get_phys_point_data(m_hyperon)),
                         fmt='o', capthick=2.0, elinewidth=5.0, label='Experiment')

            plt.errorbar(x=gv.mean(x),
                         xerr=gv.sdev(x),
                         y=gv.mean(self.mass_fit_fcn(self.get_phys_point_data())[m_hyperon]),
                         yerr=gv.sdev(self.mass_fit_fcn(self.get_phys_point_data()))[m_hyperon],
                         fmt='o', capthick=2.0, elinewidth=5.0, label='Calculated')

            # Set features for plot
            plt.legend()
            plt.grid()
            if title is None:
                if self.use_phys_f_pi:
                    plt.title("Best Fit (at $V=\infty$, physical $f_\pi =$ %s MeV)"
                               %(self.get_phys_point_data('f_pion')), fontsize=24)
                else:
                    plt.title("Best Fit (at $V=\infty$, interpolating $f_\pi = f_\pi(m_\pi)$)", fontsize=24)
            else:
                plt.title(title, fontsize=24)
            plt.xlabel(xlabel, fontsize=24)
            plt.ylabel(self._fmt_key_as_latex(m_hyperon), fontsize=24)

            fig = plt.gcf()
            if show_plot:
                plt.show()
            else:
                plt.close()
            figs = np.append(figs, fig)

        return figs


    def plot_contour(self, show_plot=True):

        figs = []
        for m_hyperon in self.hyperon_mass_keys:
            min_max = lambda x : np.array([np.min(x), np.max(x)])

            minimum, maximum = min_max([self.fit_data[abbr]['m_pion'] for abbr in self.abbrs])
            m_pi = np.linspace(minimum, maximum)

            minimum, maximum = min_max([self.fit_data[abbr]['f_pion'] for abbr in self.abbrs])
            f_pi = np.linspace(minimum, maximum)

            to_phy_data = lambda x, y : {'m_pion' : x, 'f_pion' : y, 'V' : 0,  'a' : 0}

            X, Y = np.meshgrid(m_pi, f_pi)

            # Generate omega mass for each point on meshgrid
            zs = np.array([self.mass_fit_fcn(to_phy_data(x, y))[m_hyperon] for x, y in zip(np.ravel(X), np.ravel(Y))])
            Z = gv.mean(zs.reshape(X.shape))
            Z_err = gv.sdev(zs.reshape(X.shape))

            cmap = matplotlib.cm.get_cmap('rainbow')
            min_max = lambda x : [np.min(x), np.max(x)]
            minimum, maximum = min_max(Z_err/Z)
            norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)
            color = cmap(norm(Z_err/Z))

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            ax.plot_surface(X, Y, Z, facecolors=color)

            m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
            m.set_array(color)
            color_bar = plt.colorbar(m, fraction=0.01, shrink=1.0)

            color_bar.set_alpha(0.8)
            color_bar.draw_all()
            color_bar.set_label('Relative Error', fontsize = 24)
            ax.set_xlabel(self._fmt_key_as_latex('m_pion'), fontsize=24)
            ax.set_ylabel(self._fmt_key_as_latex('f_pion'), fontsize=24)
            ax.set_zlabel(self._fmt_key_as_latex(m_hyperon), fontsize=24)

            plt.title('Fit contour at $V=\infty$, $a=0$', fontsize=24)

            fig = plt.gcf()
            if show_plot:
                plt.show()
            else:
                plt.close()
            figs = np.append(figs, fig)

        return figs

    def plot_best_fit_info(self):
        plt.axis('off')
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # place a text box in upper left in axes coords
        text = self.__str__()#.expandtabs()
        plt.text(1.0, 1.0, str(text),
                 fontsize=12, ha='right', va='top', family='monospace', bbox=props)

        plt.tight_layout()
        fig = plt.gcf()
        plt.close()

        return fig

    def make_plots(self):
        figs = self.plot_best_fit_info()

        if self.strangeness == 3:
            figs = np.append(figs,
                             self.plot_parameters(xy_parameters=['m_pion', 'm_omega'], color_parameter='a'))

        figs = np.append(figs, self.plot_fit(show_plot=False))

        # Make histograms
        fit_keys = self.get_fit_keys()
        for key in fit_keys:
            figs = np.append(figs,
                             self.plot_parameter_histogram(key))


        # Make error ellipses
        # Casts indixes of fit_keys as an upper-triangular matrix,
        # thereby allowing us to get all the 2-combinations
        # of the set of fit_keys
        rs,cs = np.triu_indices(len(fit_keys),1)
        for (i, j) in zip(rs, cs):
            figs = np.append(figs,
                             self.plot_error_ellipsis([fit_keys[i], fit_keys[j]]))

        return figs
