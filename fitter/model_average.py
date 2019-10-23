import lsqfit
import numpy as np
import gvar as gv
import time
import matplotlib
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
from collections import OrderedDict

from fitter import fitter

class model_average(object):
    def __init__(self, fit_results):
        self.fit_results = fit_results

    def _get_model_info_from_name(self, name):
        output_dict = {}
        output_dict['base'] = name.split('_')[0]
        output_dict['F2'] = name.split('_')[1]
        output_dict['include_alphaS'] = False
        output_dict['include_FV'] = False
        output_dict['include_logSq'] = False

        if '_FV' in name:
            output_dict['include_FV'] = True
        if '_alphaS' in name:
            output_dict['include_alphaS'] = True
        if '_logSq' in name:
            output_dict['include_logSq'] = True

        return output_dict

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

        elif param == 'FKFK':
            return '$F_K^2$'
        elif param == 'FKFpi':
            return '$F_K F_\pi$'
        elif param == 'FpiFpi':
            return '$F_\pi^2$'

        else:
            return param

    def average(self, param=None):
        if param is None:
            param = 'fit'


        # Only get results that aren't None
        nonempty_keys = []
        for model in self.fit_results.keys():
            if (self.fit_results[model][param] is not np.nan) and (self.fit_results[model][param] is not None):
                nonempty_keys.append(model)

        # calculate P( M_k | D )
        prob_Mk_given_D = lambda model_k : (
             np.exp(self.fit_results[model_k]['logGBF']) / np.sum([np.exp(self.fit_results[model_l]['logGBF']) for model_l in nonempty_keys])
        )

        # Get central value
        expct_y = 0
        for model in nonempty_keys:
            expct_y += gv.mean(gv.gvar(self.fit_results[model][param])) *prob_Mk_given_D(model)

        # Get variance
        var_y = 0
        for model in nonempty_keys:
            var_y += gv.var(gv.gvar(self.fit_results[model][param])) *prob_Mk_given_D(model)
        for model in nonempty_keys:
            var_y += (gv.mean(gv.gvar(self.fit_results[model][param])))**2 *prob_Mk_given_D(model)

        var_y -= (expct_y)**2

        return gv.gvar(expct_y, np.sqrt(var_y))



    def plot_comparison(self, param=None, other_results=None, title=None, xlabel=None,
                        show_model_avg=True):

        if param is None:
            param = 'fit'


        if title is None:
            title = ""
        if xlabel is None:
            xlabel = self._param_keys_dict(param)

        colors = ['salmon', 'darkorange', 'mediumaquamarine', 'orchid', 'silver']
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
                        model_info = self._get_model_info_from_name(name)
                        if model_info['base'] == 'ma':
                            color = colors[0]
                        elif model_info['base'] == 'xpt':
                            color = colors[1]
                        elif model_info['base'] == 'xpt-ratio':
                            color = colors[2]
                        elif model_info['base'] == 'ma-ratio':
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
            avg = self.average(param)
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
        logGBF_max = np.nanmax([gv.mean(gv.gvar(self.fit_results[model]['logGBF']))
                               for model in self.fit_results.keys()])

        y=0
        labels = np.array([])
        for results in [self.fit_results]:
            for name in results.keys():
                model_info = self._get_model_info_from_name(name)
                if model_info['base'] == 'ma':
                    color = colors[0]
                elif model_info['base'] == 'xpt':
                    color = colors[1]
                elif model_info['base'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['base'] == 'ma-ratio':
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
        for results in [self.fit_results]:
            for name in results.keys():
                model_info = self._get_model_info_from_name(name)
                if model_info['base'] == 'ma':
                    color = colors[0]
                elif model_info['base'] == 'xpt':
                    color = colors[1]
                elif model_info['base'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['base'] == 'ma-ratio':
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
            for name in results.keys():
                model_info = self._get_model_info_from_name(name)
                if model_info['base'] == 'ma':
                    color = colors[0]
                elif model_info['base'] == 'xpt':
                    color = colors[1]
                elif model_info['base'] == 'xpt-ratio':
                    color = colors[2]
                elif model_info['base'] == 'ma-ratio':
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


    def plot_histogram(self, param, title=None, xlabel=None):

        if title is None:
            title = ""
        if xlabel is None:
            xlabel = self._param_keys_dict(param)

        colors = ['salmon', 'darkorange', 'mediumaquamarine', 'orchid', 'silver']
        for j, F2 in enumerate(['All', 'FKFK', 'FKFpi', 'FpiFpi']):

            # read Bayes Factors
            logGBF_list = [self.fit_results[model]['logGBF'] for model in self.fit_results.keys()]

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


            # for plotting
            temp_array = np.array([self.fit_results[model][param] for model in self.fit_results.keys()])
            temp_array = gv.gvar(temp_array[temp_array != 'nan'])
            pm_arr = lambda arr, k : gv.mean(temp_array) + k *gv.sdev(temp_array)
            min_x = np.min(pm_arr(temp_array, -2))
            max_x = np.max(pm_arr(temp_array, 2))
            x = np.linspace(min_x, max_x, 2000)


            for model in self.fit_results.keys():
                model_info = self._get_model_info_from_name(model)
                if self.fit_results[model][param] is not np.nan and (model_info['F2'] == F2 or F2=='All'):
                    r = gv.gvar(self.fit_results[model][param])
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

            # '-','--','-.',':'
            # #ec5d57 #70bf41 #51a7f9

            fig = plt.figure('result histogram')#,figsize=fig_size2)
            ax = plt.axes()
            ax.fill_between(x=x,y1=ysum,facecolor=colors[j], edgecolor='black',alpha=0.4,label=self._param_keys_dict(F2))
            if F2 == 'All':
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

            for a in ydict.keys():
                ax.plot(x,ydict[a], color=colors[j])

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
