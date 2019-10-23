import lsqfit
import numpy as np
import gvar as gv
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
from collections import OrderedDict

from fitter import fitter

class model_average(object):
    def __init__(self, fit_results):
        self.fit_results = fit_results

    def _param_keys_dict(self, param):
        if param == 'fit':
            return '$F_K/F_\pi$'
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

        else:
            return param

    def average(self, fit_results=None, param=None):

        if param is None:
            param = 'fit'
        if fit_results is None:
            fit_results = self.fit_results

        # Only get results that aren't None
        nonempty_keys = []
        for model in fit_results.keys():
            if (fit_results[model][param] is not np.nan) and (fit_results[model][param] is not None):
                nonempty_keys.append(model)

        # calculate P( M_k | D )
        prob_Mk_given_D = lambda model_k : (
             np.exp(fit_results[model_k]['logGBF']) / np.sum([np.exp(fit_results[model_l]['logGBF']) for model_l in nonempty_keys])
        )

        # Get central value
        expct_y = 0
        for model in nonempty_keys:
            expct_y += gv.mean(gv.gvar(fit_results[model][param])) *prob_Mk_given_D(model)

        # Get variance
        var_y = 0
        for model in nonempty_keys:
            var_y += gv.var(gv.gvar(fit_results[model][param])) *prob_Mk_given_D(model)
        for model in nonempty_keys:
            var_y += (gv.mean(gv.gvar(fit_results[model][param])))**2 *prob_Mk_given_D(model)

        var_y -= (expct_y)**2

        return gv.gvar(expct_y, np.sqrt(var_y))



    def plot_comparison(self, fit_results=None, param=None, other_results=None, title=None, xlabel=None,
                        show_model_avg=True):

        if param is None:
            param = 'fit'
        if fit_results is None:
            fit_results = self.fit_results


        if title is None:
            title = ""
        if xlabel is None:
            xlabel = self._param_keys_dict(param)

        colors = ['salmon', 'darkorange', 'mediumaquamarine', 'orchid', 'silver']
        markers = ['^', 'o', 'v', '*']

        results_array = [fit_results]
        if other_results is not None:
            results_array.append(other_results)

        # These axes compare fits
        ax_fits = plt.axes([0.10,0.10,0.49,0.8])

        y=0
        labels = np.array([])
        for results in results_array:
            plt.axhline(y-0.5, ls='--')
            for name in results.keys():
                if param in results[name].keys():

                    # Add band for FLAG
                    if name in ['FLAG']:
                        color = 'palevioletred'
                        x = gv.mean(gv.gvar(results[name][param]))
                        xerr = gv.sdev(gv.gvar(results[name][param]))
                        plt.axvspan(x-xerr, x+xerr, alpha=0.3, color=color, label='FLAG')

                        y = y + 1

                    else:
                        # Color by base model
                        if 'ma_' in name:
                            color = colors[0]
                        elif 'xpt_' in name:
                            color = colors[1]
                        elif 'xpt-ratio_' in name:
                            color = colors[2]
                        elif 'ma-ratio_' in name:
                            color = colors[3]
                        else:
                            color = colors[4]

                        # Marker by F^2
                        if '_FKFK_' in name:
                            marker = markers[0]
                        elif '_FKFpi_' in name:
                            marker = markers[1]
                        elif '_FpiFpi_' in name:
                            marker = markers[2]
                        else:
                            marker = markers[3]

                        x = gv.mean(gv.gvar(results[name][param]))
                        xerr = gv.sdev(gv.gvar(results[name][param]))
                        plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                                     marker=marker, capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                                     color=color, elinewidth=5.0, label=name.split('_')[0])
                        y = y + 1
                        labels = np.append(labels, str(name))
                        
        ymax = y

        # Show model average
        if show_model_avg:
            avg = self.average(fit_results, param)
            pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
            plt.axvspan(pm(avg, -1), pm(avg, +1), alpha=0.3, color='cornflowerblue', label='model avg')
            plt.axvline(pm(avg, 0), ls='-.', color='m')

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
        logGBF_max = np.nanmax([gv.mean(gv.gvar(fit_results[model]['logGBF']))
                               for model in fit_results.keys()])

        y=0
        labels = np.array([])
        for results in [fit_results]:
            for name in results.keys():
                if 'ma_' in name:
                    color = colors[0]
                elif 'xpt_' in name:
                    color = colors[1]
                elif 'xpt-ratio_' in name:
                    color = colors[2]
                elif 'ma-ratio_' in name:
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if '_FKFK_' in name:
                    marker = markers[0]
                elif '_FKFpi_' in name:
                    marker = markers[1]
                elif '_FpiFpi_' in name:
                    marker = markers[2]
                else:
                    marker = markers[3]


                logGBF = gv.mean(gv.gvar(results[name]['logGBF']))
                x = np.exp(logGBF - logGBF_max)

                #plt.axvline(x, ls='--', alpha=0.4)
                plt.scatter(x=x, y=y, color=color, marker=marker)
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
        for results in [fit_results]:
            for name in results.keys():
                if 'ma_' in name:
                    color = colors[0]
                elif 'xpt_' in name:
                    color = colors[1]
                elif 'xpt-ratio_' in name:
                    color = colors[2]
                elif 'ma-ratio_' in name:
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if '_FKFK_' in name:
                    marker = markers[0]
                elif '_FKFpi_' in name:
                    marker = markers[1]
                elif '_FpiFpi_' in name:
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
        for results in [fit_results]:
            for name in results.keys():
                if 'ma_' in name:
                    color = colors[0]
                elif 'xpt_' in name:
                    color = colors[1]
                elif 'xpt-ratio_' in name:
                    color = colors[2]
                elif 'ma-ratio_' in name:
                    color = colors[3]
                else:
                    color = colors[4]

                # Marker by F^2
                if '_FKFK_' in name:
                    marker = markers[0]
                elif '_FKFpi_' in name:
                    marker = markers[1]
                elif '_FpiFpi_' in name:
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
