import lsqfit
import numpy as np
import gvar as gv
import time
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import scipy.stats as stats
from collections import OrderedDict

import fitter.fitter as fit
import fitter.data_loader as dl

class model_average(object):
    def __init__(self, fit_results):
        for model in fit_results.keys():
            fit_results[model]['FK/Fpi_pm'] = gv.gvar(fit_results[model]['FK/Fpi']) *np.sqrt(1 + gv.gvar(fit_results[model]['delta_su2']))
        self.fit_results = fit_results
        self.data_loader = dl.data_loader()

    def __str__(self):
        output  = 'FK/Fpi_pm: %s\n'%(self.average('FK/Fpi_pm'))
        output += '[FLAG:     %s]\n'%(self._get_phys_point_data()['FK/Fpi_pm'])

        output += '\n---\n'
        fk_fpi = self.average('FK/Fpi', split_unc=True)
        output += 'FK/Fpi:        %s \n' %(gv.gvar(fk_fpi[0], np.sqrt(fk_fpi[1]**2 + fk_fpi[2]**2)))
        output += 'delta_su(2):  %s \n' %(self.average('delta_su2'))

        sig_fig = lambda x : np.around(x, int(np.floor(-np.log10(x))+3)) if x>0 else x
        output += '\n---\n'
        output += 'Uncertainty: \n'
        output += '   Unexplained: % .5f \n' %(fk_fpi[1])
        output += '   Explained:   % .5f \n' %(fk_fpi[2])

        error_budget = self.error_budget()
        if error_budget['chiral'] is not None:
            output += '\n---\n'
            output += 'Error Budget: \n'
            output += '   Chiral:      % .5f \n' %(error_budget['chiral'])
            output += '   Phys Point:  % .5f \n' %(error_budget['pp_input'])
            output += '   Statistical: % .5f \n' %(error_budget['stat'])

        model_list = self.get_model_names(by_weight=True)
        weight = lambda model_k : np.exp(self.fit_results[model_k]['logGBF']) / np.sum([np.exp(self.fit_results[model_l]['logGBF']) for model_l in model_list])
        output += '\n---\n'
        output += 'Highest Weight: \n'
        for k in range(np.min([5, len(model_list)])):
            output += '  % .3f:  %s\n' %(weight(model_list[k]), model_list[k])

        return output

    def _get_model_info_from_name(self, name):
        return self.data_loader.get_model_info_from_name(name)

    def _get_fit_extrapolation(self, model):
        return gv.gvar(self.fit_results[model]['FK/Fpi'])

    def _get_fit_posterior(self, model):
        if 'posterior' in self.fit_results[model] and bool(self.fit_results[model]['posterior']):
            return self.fit_results[model]['posterior']
        else:
            return None

    def _get_fit_prior(self, model):
        if 'prior' in self.fit_results[model] and bool(self.fit_results[model]['prior']):
            return self.fit_results[model]['prior']
        else:
            return None

    def _get_phys_point_data(self, model=None):
        phys_point_data = {
            'a/w0' : 0,
            'a' : 0,
            'L' : np.infty,
            'alpha_s' : 0, # Not sure, but not important since it comes with a^2

            'mpi' : gv.gvar('134.8(3)'), # '138.05638(37)'
            'mk' : gv.gvar('494.2(3)'), # '495.6479(92)'

            'a2DI' : 0,
            'Fpi' : gv.gvar(130.2/np.sqrt(2), 0.8/np.sqrt(2)), # PDG
            'FK' : gv.gvar(155.5/np.sqrt(2), 0.7/np.sqrt(2)), # PDG
            'w0' : gv.gvar('0.175(10)'),

            'FK/Fpi_pm' : gv.gvar('1.1932(19)'), # FLAG, SU(2) isospin corrected value (arxiv/1902.08191, eqn 80)
        }

        if model is not None:
            model_info = self._get_model_info_from_name(model)

            FK = phys_point_data['FK']
            Fpi = phys_point_data['Fpi']
            if model_info['F2'] == 'FKFpi':
                phys_point_data['lam2_chi'] = (4*np.pi)**2 *FK *Fpi
            elif model_info['F2'] == 'FpiFpi':
                phys_point_data['lam2_chi'] = (4*np.pi)**2 *Fpi *Fpi
            elif model_info['F2'] == 'FKFK':
                phys_point_data['lam2_chi'] = (4*np.pi)**2 *FK *FK

        return phys_point_data

    def _param_keys_dict(self, param):
        if param == 'FK/Fpi':
            return '$F_K/F_\pi$'
        elif param == 'FK/Fpi_pm':
            return '$F^\pm_K/F^\pm_\pi$'
        elif param == 'delta_su2':
            return '$\delta_{SU(2)}$'

        elif param == 'L_5':
            return '$L_5$'
        elif param == 'L_4':
            return '$L_4$'

        elif param == 'A_k':
            return '$A_k$'
        elif param == 'A_p':
            return '$A_p$'
        elif param == 'A_a':
            return '$A_a$'
        elif param == 'A_loga':
            return '$A_{\log a}$'

        elif param == 'FKFK':
            return '$F_K^2$'
        elif param == 'FKFpi':
            return '$F_K F_\pi$'
        elif param == 'FpiFpi':
            return '$F_\pi^2$'

        else:
            return param

    def average(self, param=None, split_unc=False, include_unc=True):
        if param is None:
            param = 'FK/Fpi'

        y = {}
        for model in self.get_model_names():

            if param == 'FK/Fpi':
                y[model] = self._get_fit_extrapolation(model)

            elif self._get_fit_posterior(model) is not None and param in self._get_fit_posterior(model):
                y[model] =  self._get_fit_posterior(model)[param]

            # Error budget
            elif param.startswith('eb:') and 'error_budget' in self.fit_results[model]:
                y[model] = self.fit_results[model]['error_budget'][param.split(':')[-1]]

            elif param in self.fit_results[model]:
                y[model] = self.fit_results[model][param]

            else:
                y[model] = None


        # Only get results that aren't None
        nonempty_keys = []
        for model in self.get_model_names():
            if (y[model] is not np.nan) and (y[model]is not None):
                nonempty_keys.append(model)

        if nonempty_keys == []:
            return None

        # calculate P( M_k | D )
        prob_Mk_given_D = lambda model_k : (
             np.exp(self.fit_results[model_k]['logGBF']) / np.sum([np.exp(self.fit_results[model_l]['logGBF']) for model_l in nonempty_keys])
        )

        # Get central value
        expct_y = 0
        for model in nonempty_keys:
            expct_y += gv.mean(gv.gvar(y[model])) *prob_Mk_given_D(model)

        if not include_unc:
            return expct_y

        # Get variance
        if not split_unc:
            var_y = 0
            for model in nonempty_keys:
                var_y += gv.var(gv.gvar(y[model])) *prob_Mk_given_D(model)
            for model in nonempty_keys:
                var_y += (gv.mean(gv.gvar(y[model])))**2 *prob_Mk_given_D(model)

            var_y -= (expct_y)**2

            return gv.gvar(expct_y, np.sqrt(var_y))

        # Split statistics (unexplained var)/model selection (explained var)
        if split_unc:
            var_model = 0
            for model in nonempty_keys:
                var_model += gv.var(gv.gvar(y[model])) *prob_Mk_given_D(model)

            var_selection = 0
            for model in nonempty_keys:
                var_selection += (gv.mean(gv.gvar(y[model])))**2 *prob_Mk_given_D(model)
            var_selection -= (expct_y)**2

            return [expct_y, np.sqrt(var_model), np.sqrt(var_selection)]

    def error_budget(self):
        output = {}
        for key in ['chiral', 'pp_input', 'stat']:
            output[key] = self.average('eb:'+key, include_unc=False)

        return output

    def extrapolate_to_phys_point(self, model):
        data = self._get_phys_point_data(model)
        return self.fitfcn(model, data)

    def fitfcn(self, model, data, p=None):
        model_info = self._get_model_info_from_name(model).copy()

        if p is None:
            p = self._get_fit_posterior(model)

        if model_info['fit_type'] in ['xpt', 'ma']:
            model_info['fit_type'] = 'xpt'
            fitfcn = fit.fk_fpi_model(datatag='xpt', **model_info).fitfcn
        elif model_info['fit_type'] in ['xpt-ratio', 'ma-ratio']:
            model_info['fit_type'] = 'xpt-ratio'
            fitfcn = fit.fk_fpi_model(datatag='xpt-ratio', **model_info).fitfcn
        elif model_info=='poly':
            fitfcn = fit.fk_fpi_model(datatag='poly', **model_info).fitfcn

        #fitfcn = fit.fk_fpi_model(datatag='temp', **model_info).fitfcn

        return fitfcn(p=p, fit_data=data)

    def get_model_names(self, by_weight=False):
        if by_weight:
            models_list = self.get_model_names()
            temp = {model : self.fit_results[model]['logGBF'] for model in models_list}
            sorted_list = [model for model, logGBF
                           in sorted(temp.items(), key=lambda item: item[1], reverse=True)]
            return sorted_list

        else:
            return sorted(list(self.fit_results))

    def plot_comparison(self, param=None, other_results=None, title=None, xlabel=None,
                        show_model_avg=True):

        if param is None:
            param = 'FK/Fpi_pm'

        if title is None:
            title = ""
        if xlabel is None:
            xlabel = self._param_keys_dict(param)


        colors = ['salmon', 'darkorange', 'mediumaquamarine', 'orchid', 'navy']
        markers = ['^', 'o', 'v', '*']

        results_array = [self.fit_results]
        if other_results is not None:
            results_array.append(other_results)

        # These axes compare fits
        ax_fits = plt.axes([0.10,0.10,0.49,0.8])

        y=0
        labels = np.array([])
        for results in results_array:
            plt.axhline(y-0.5, ls='--')
            for name in sorted(results.keys()):

                param_value = None

                if param in results[name].keys():
                    param_value = gv.gvar(results[name][param])
                elif 'posterior' in results[name] and param in results[name]['posterior'].keys():
                    param_value = gv.gvar(results[name]['posterior'][param])
                else:
                    param_value = None

                if param_value is not None:

                    # Add band for FLAG
                    if name in ['FLAG']:
                        color = 'palevioletred'
                        x = gv.mean(param_value)
                        xerr = gv.sdev(param_value)
                        plt.axvspan(x-xerr, x+xerr, alpha=0.3, color=color, label='FLAG')

                        y = y + 1

                    else:
                        # Color by fit_type model
                        model_info = self._get_model_info_from_name(name)
                        if model_info['fit_type'] == 'ma':
                            color = colors[0]
                        elif model_info['fit_type'] == 'xpt':
                            color = colors[1]
                        elif model_info['fit_type'] == 'xpt-ratio':
                            color = colors[2]
                        elif model_info['fit_type'] == 'ma-ratio':
                            color = colors[3]
                        else:
                            color = colors[4]

                        # Marker by F^2
                        if model_info['F2'] == 'FKFK':
                            marker = markers[0]
                        elif model_info['F2'] == 'FKFpi':
                            marker = markers[1]
                        elif model_info['F2'] == 'FpiFpi':
                            marker = markers[2]
                        else:
                            marker = markers[3]

                        x = gv.mean(param_value)
                        xerr = gv.sdev(param_value)

                        plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                                     marker=marker, capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                                     color=color, elinewidth=5.0, label=name.split('_')[0])
                        y = y + 1
                        labels = np.append(labels, str(name))

        ymax = y

        # Show model average
        if show_model_avg:
            try:
                avg = self.average(param)
                pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
                plt.axvspan(pm(avg, -1), pm(avg, +1), alpha=0.3, color='cornflowerblue', label='model avg')
                plt.axvline(pm(avg, 0), ls='-.', color='m')
            except:
                pass

        # Get unique labels
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys())

        #plt.yticks(range(len(labels)), labels, fontsize=, rotation=65)
        plt.yticks([])
        plt.ylim(-1, y)

        #plt.title(particle+", "+abbr, fontsize=30)
        plt.xlabel(xlabel, fontsize=24)
        plt.title(title, fontsize=24)

        plt.grid(ls='--')

        # These axes compare the logGBFs
        ax_logGBF = plt.axes([0.60,0.10,0.09,0.8])

        # Get max logGBF
        logGBF_max = np.nanmax([gv.mean(gv.gvar(self.fit_results[model]['logGBF']))
                               for model in self.fit_results.keys()])

        y=0
        labels = np.array([])
        for results in [self.fit_results]:
            for name in sorted(results.keys()):
                model_info = self._get_model_info_from_name(name)
                if model_info['fit_type'] == 'ma':
                    color = colors[0]
                elif model_info['fit_type'] == 'xpt':
                    color = colors[1]
                elif model_info['fit_type'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['fit_type'] == 'ma-ratio':
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if model_info['F2'] == 'FKFK':
                    marker = markers[0]
                elif model_info['F2'] == 'FKFpi':
                    marker = markers[1]
                elif model_info['F2'] == 'FpiFpi':
                    marker = markers[2]
                else:
                    marker = markers[3]


                logGBF = gv.mean(gv.gvar(results[name]['logGBF']))
                x = np.exp(logGBF - logGBF_max)

                alpha = 1
                if x < 0.01:
                    alpha = 0

                #plt.axvline(x, ls='--', alpha=0.4)
                plt.scatter(x=x, y=y, color=color, marker=marker, alpha=alpha)
                y = y + 1
                labels = np.append(labels, str(name))

        for ti in np.arange(5)/4.0:
            plt.axvline(ti, ls='--', alpha=0.2)

        plt.yticks([])
        plt.ylim(-1, ymax)
        plt.xlabel("Norm\nGBF", fontsize=24)

        # These axes compare the Q-values
        ax_Q = plt.axes([0.70,0.10,0.09,0.8])

        y=0
        labels = np.array([])
        for results in [self.fit_results]:
            for name in sorted(results.keys()):
                model_info = self._get_model_info_from_name(name)
                if model_info['fit_type'] == 'ma':
                    color = colors[0]
                elif model_info['fit_type'] == 'xpt':
                    color = colors[1]
                elif model_info['fit_type'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['fit_type'] == 'ma-ratio':
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if model_info['F2'] == 'FKFK':
                    marker = markers[0]
                elif model_info['F2'] == 'FKFpi':
                    marker = markers[1]
                elif model_info['F2'] == 'FpiFpi':
                    marker = markers[2]
                else:
                    marker = markers[3]

                x = gv.mean(gv.gvar(results[name]['Q']))
                #plt.axvline(x, ls='--', alpha=0.4)
                plt.scatter(x=x, y=y, color=color, marker=marker)
                y = y + 1
                labels = np.append(labels, str(name))

        for ti in np.arange(5)/4.0:
            plt.axvline(ti, ls='--', alpha=0.2)

        plt.yticks([])
        plt.xlim(-0.05, 1.05)
        plt.ylim(-1, ymax)
        plt.xlabel("$Q$", fontsize=24)


        # These axes compare the reduced chi2
        ax_chi2 = plt.axes([0.80,0.10,0.1,0.8])

        y=0
        labels = np.array([])
        for results in [self.fit_results]:
            for name in sorted(results.keys()):
                model_info = self._get_model_info_from_name(name)
                if model_info['fit_type'] == 'ma':
                    color = colors[0]
                elif model_info['fit_type'] == 'xpt':
                    color = colors[1]
                elif model_info['fit_type'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['fit_type'] == 'ma-ratio':
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if model_info['F2'] == 'FKFK':
                    marker = markers[0]
                elif model_info['F2'] == 'FKFpi':
                    marker = markers[1]
                elif model_info['F2'] == 'FpiFpi':
                    marker = markers[2]
                else:
                    marker = markers[3]

                x = gv.mean(gv.gvar(results[name]['chi2/df']))
                #plt.axvline(x, ls='--', alpha=0.4)
                plt.scatter(x=x, y=y, color=color, marker=marker)
                y = y + 1
                labels = np.append(labels, str(name))

        for ti in np.arange(9)/4.0:
            plt.axvline(ti, ls='--', alpha=0.2)

        plt.yticks([])
        plt.xlim(0, 2)
        plt.ylim(-1, ymax)
        plt.xlabel(r"$\chi^2_\nu$", fontsize=24)

        fig = plt.gcf()
        plt.close()

        return fig

    # parameter = 'a', 'mpi', 'volume'
    def plot_fits(self, parameter):

        # Check that posterior contains covarianace matrix
        temp_model_name = list(self.get_model_names())[0]
        temp_po = self._get_fit_posterior(temp_model_name)
        temp_p = self._get_fit_prior(temp_model_name)

        if (temp_p is None) or (temp_po is None) or (gv.uncorrelated(temp_po['A_a'], temp_p['A_a'])):
            print('Correlations between gvar variables lost! Fit plots will be inaccurate.')
            return None

        colors = ['darkorange', 'mediumaquamarine', 'orchid', 'skyblue', 'silver']
        #colors = ['cyan', 'magenta', 'yellow', 'black', 'silver']

        data = self._get_phys_point_data()
        if parameter == 'a':
            xlabel = '$a$ (fm)'

            x = np.linspace(0, 0.8804, 50)
            data['a/w0'] = x
            x = x *gv.mean(self._get_phys_point_data()['w0'])

        elif parameter == 'mpi':
            xlabel = '$m_\pi$ (MeV)'
            plt.axvline(gv.mean(self._get_phys_point_data()['mpi']))

            x = np.linspace(10, 400, 50)
            data['mpi'] = x


        elif parameter == 'volume':
            xlabel = '$e^{-m_\pi L} / \sqrt{m_\pi L}$'

            x = np.linspace(3, 6, 50) /gv.mean(self._get_phys_point_data()['mpi'])
            data['L'] = x
            x = np.exp(-x *self._get_phys_point_data()['mpi']) / np.sqrt(x *self._get_phys_point_data()['mpi'])


        else:
            return None

        total_GBF = np.sum([np.exp(self.fit_results[model_l]['logGBF']) for model_l in self.fit_results.keys()])

        pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
        for name in self.fit_results.keys():

            model_info = self._get_model_info_from_name(name)
            data['lam2_chi'] = self._get_phys_point_data(name)['lam2_chi']
            if model_info['F2'] == 'FKFK':
                color = colors[0]
                F2 = self._param_keys_dict('FKFK')
            elif model_info['F2'] == 'FKFpi':
                color = colors[1]
                F2 = self._param_keys_dict('FKFpi')
            elif model_info['F2'] == 'FpiFpi':
                color = colors[2]
                F2 = self._param_keys_dict('FpiFpi')
            else:
                color = colors[3]

            y = self.fitfcn(name, data) *np.sqrt(1 + gv.gvar(self.fit_results[name]['delta_su2']))

            weight = np.exp(self.fit_results[name]['logGBF']) / total_GBF
            plt.fill_between(pm(x, 0), pm(y, -1), pm(y, 1), color=color, alpha=np.min([2*weight, 1]), rasterized=False)


        y = np.repeat(self._get_phys_point_data(name)['FK/Fpi_pm'],len(x))
        plt.fill_between(pm(x, 0), pm(y, -1), pm(y, 1), color=colors[3], alpha=0.5, rasterized=False)


        # Add legend
        kk_patch = mpatches.Patch(color=colors[0], alpha=0.5, label=self._param_keys_dict('FKFK'))
        kp_patch = mpatches.Patch(color=colors[1], alpha=0.5, label=self._param_keys_dict('FKFpi'))
        pp_patch = mpatches.Patch(color=colors[2], alpha=0.5, label=self._param_keys_dict('FpiFpi'))
        flag_patch = mpatches.Patch(color=colors[3], alpha=0.5, label='FLAG')
        plt.legend(handles=[flag_patch,kk_patch, kp_patch, pp_patch], fontsize=15)


        plt.xlim(np.min(gv.mean(x)), np.max(gv.mean(x)))
        plt.xlabel(xlabel, fontsize=24)
        plt.ylabel('$F^\pm_K/F^\pm_\pi$', fontsize=24)


        fk_fpi_avg = self.average()
        pm = lambda g, k : g.mean + k *g.sdev
        plt.ylim(pm(fk_fpi_avg, -5), pm(fk_fpi_avg, +5))

        fig = plt.gcf()
        plt.close()
        return fig

    # See self._get_model_info_from_name for possible values for 'vary_choice'
    def plot_histogram(self, param=None, title=None, xlabel=None, vary_choice='F2'):
        if param is None:
            param = 'FK/Fpi_pm'
        if xlabel is None:
            xlabel = self._param_keys_dict(param)+' (choose: '+vary_choice+')'
        if title is None:
            title = ""

        param_avg = self.average(param=param)
        pm = lambda g, k : g.mean + k *g.sdev
        x = np.linspace(pm(param_avg, -4), pm(param_avg, +4), 2000)

        # Determine ordering
        # Have larger contributions behind smaller contributions
        choices = np.unique([self._get_model_info_from_name(model)[vary_choice] for model in self.get_model_names()])
        temp_dict = {choice : 0 for choice in choices}
        for model in self.get_model_names():
            model_info = self._get_model_info_from_name(model)
            temp_dict[model_info[vary_choice]] += np.exp(self.fit_results[model]['logGBF'])
            choices = sorted(temp_dict, key=temp_dict.get, reverse=True)

        # Set colors
        #cmap = matplotlib.cm.get_cmap('gist_rainbow')
        #colors =  ['whitesmoke']
        #colors.extend([cmap(c) for c in np.linspace(0, 1, len(choices)+1)])
        #colors = ['whitesmoke', 'salmon', 'palegreen', 'lightskyblue', 'plum']
        colors = ['whitesmoke', 'crimson', 'springgreen', 'deepskyblue', 'magenta']

        for j, choice in enumerate(np.append(['All'], choices)):

            # read Bayes Factors
            logGBF_list = [self.fit_results[model]['logGBF'] for model in self.get_model_names()]

            # initiate a bunch of parameters
            y = 0
            y_list = []
            y_dict = dict()

            # weights
            w_lst = []
            wd = dict()

            # p. dist. fcn
            pdf = 0
            pdfdict = dict()

            # c. dist. fcn.
            cdf = 0
            cdfdict = dict()

            for model in self.get_model_names():
                model_info = self._get_model_info_from_name(model)

                r = np.nan
                if param == 'FK/Fpi':
                    r = self._get_fit_extrapolation(model)

                elif (self._get_fit_posterior(model) is not None) and (param in self._get_fit_posterior(model)):
                    r =  self._get_fit_posterior(model)[param]

                elif param in self.fit_results[model]:
                    r = gv.gvar(self.fit_results[model][param])

                else:
                    r = np.nan

                if r is not np.nan and r is not None and (str(model_info[vary_choice]) == choice or choice=='All'):
                    #r = gv.gvar(self.fit_results[model][param])
                    y_dict[model] = r

                    w = 1/sum(np.exp(np.array(logGBF_list)-self.fit_results[model]['logGBF']))
                    sqrtw = np.sqrt(w) # sqrt scales the std dev correctly
                    wd[model] = w
                    w_lst.append(w)

                    y += gv.gvar(w*r.mean,sqrtw*r.sdev)
                    y_list.append(r.mean)

                    p = stats.norm.pdf(x,r.mean,r.sdev)
                    pdf += w*p
                    pdfdict[model] = w*p


                    c = stats.norm.cdf(x,r.mean,r.sdev)
                    cdf += w*c
                    cdfdict[model] = w*c

            y_list = np.array(y_list)
            w_lst = np.array(w_lst)

            plot_params = {'x':x, 'pdf':pdf, 'pdfdict':pdfdict, 'cdf':cdf, 'cdfdict':cdfdict}

            gr = 1.618034333
            fs2_base = 3.50394
            lw = 0.5
            fs_l = 15
            fs_xy = 24
            ts = 15

            x = plot_params['x']
            ysum = plot_params['pdf']

            ydict = plot_params['pdfdict']
            cdf = plot_params['cdf']


            fig = plt.figure('result histogram')#,figsize=fig_size2)
            ax = plt.axes()

            for a in ydict.keys():
                ax.plot(x,ydict[a], color=colors[j], alpha=1.0, ls='dotted')


            ax.fill_between(x=x,y1=ysum,facecolor=colors[j], edgecolor='black',alpha=0.4,label=self._param_keys_dict(choice))
            ax.plot(x, ysum, color='k', alpha=1.0)
            if choice == 'All':
                # get 95% confidence
                lidx95 = abs(cdf-0.025).argmin()
                uidx95 = abs(cdf-0.975).argmin()
                ax.fill_between(x=x[lidx95:uidx95],y1=ysum[lidx95:uidx95],facecolor=colors[j],edgecolor='black',alpha=0.4)
                # get 68% confidence
                lidx68 = abs(cdf-0.158655254).argmin()
                uidx68 = abs(cdf-0.841344746).argmin()
                ax.fill_between(x=x[lidx68:uidx68],y1=ysum[lidx68:uidx68],facecolor=colors[j],edgecolor='black',alpha=0.4)
                # plot black curve over
                ax.errorbar(x=[x[lidx95],x[lidx95]],y=[0,ysum[lidx95]],color='black',lw=lw)
                ax.errorbar(x=[x[uidx95],x[uidx95]],y=[0,ysum[uidx95]],color='black',lw=lw)
                ax.errorbar(x=[x[lidx68],x[lidx68]],y=[0,ysum[lidx68]],color='black',lw=lw)
                ax.errorbar(x=[x[uidx68],x[uidx68]],y=[0,ysum[uidx68]],color='black',lw=lw)

                ax.errorbar(x=x,y=ysum,ls='-',color='black',lw=lw)



            leg = ax.legend(fontsize=fs_l, edgecolor='k',fancybox=False)
            ax.set_ylim(bottom=0)
            #ax.set_xlim([1.225,1.335])
            ax.set_xlabel(xlabel, fontsize=fs_xy)
            frame = plt.gca()
            frame.axes.get_yaxis().set_visible(False)
            ax.xaxis.set_tick_params(labelsize=ts,width=lw)

            # legend line width
            [ax.spines[key].set_linewidth(lw) for key in ax.spines]
            leg.get_frame().set_linewidth(lw)

        fig = plt.gcf()
        plt.close()
        return fig
