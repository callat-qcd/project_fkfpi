import numpy as np
import gvar as gv
import matplotlib
import matplotlib.pyplot as plt
from collections import OrderedDict

def model_average(results, param):

    # Only get results that aren't None
    nonempty_keys = []
    for model in results.keys():
        if (results[model][param] is not np.nan) and (results[model][param] is not None):
            nonempty_keys.append(model)

    # calculate P( M_k | D )
    prob_Mk_given_D = lambda model_k : (
         np.exp(results[model_k]['logGBF']) / np.sum([np.exp(results[model_l]['logGBF']) for model_l in nonempty_keys])
    )

    # Get central value
    expct_y = 0
    for model in nonempty_keys:
        expct_y += gv.mean(gv.gvar(results[model][param])) *prob_Mk_given_D(model)

    # Get variance
    var_y = 0
    for model in nonempty_keys:
        var_y += gv.var(gv.gvar(results[model][param])) *prob_Mk_given_D(model)
    for model in nonempty_keys:
        var_y += (gv.mean(gv.gvar(results[model][param])))**2 *prob_Mk_given_D(model)

    var_y -= (expct_y)**2

    return gv.gvar(expct_y, np.sqrt(var_y))


def plot_comparison(param, fit_results, other_results=None, title=None, xlabel=None,
                    show_model_avg=False):

    if title is None:
        title = ""
    if xlabel is None:
        xlabel = ""

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

            # Add band for FLAG
            if name in ['FLAG']:
                x = gv.mean(gv.gvar(results[name][param]))
                xerr = gv.sdev(gv.gvar(results[name][param]))
                plt.axvspan(x-xerr, x+xerr, alpha=0.3, color='g', label='FLAG')

                y = y + 1

            else:
                if 'ma_' in name:
                    fmt = 'ro'
                elif 'xpt_' in name:
                    fmt = 'co'
                elif 'xpt-ratio_' in name:
                    fmt = 'go'
                elif 'ma-ratio_' in name:
                    fmt = 'mo'
                else:
                    fmt = 'bo'

                x = gv.mean(gv.gvar(results[name][param]))
                xerr = gv.sdev(gv.gvar(results[name][param]))
                plt.errorbar(x=x, y=y, xerr=xerr, yerr=0.0,
                             fmt = fmt, capsize=0.0, mec='white', ms=10.0, alpha=0.6,
                             ecolor='b', elinewidth=5.0, label=name.split('_')[0])
                y = y + 1
                labels = np.append(labels, str(name))




    ymax = y

    # Show model average
    if show_model_avg:
        avg = model_average(fit_results, param)
        pm = lambda g, k : gv.mean(g) + k*gv.sdev(g)
        plt.axvspan(pm(avg, -1), pm(avg, +1), alpha=0.3, color='m', label='model avg')
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
                color = 'r'
            elif 'xpt_' in name:
                color = 'c'
            elif 'xpt-ratio_' in name:
                color = 'g'
            elif 'ma-ratio_' in name:
                color = 'm'
            else:
                color = 'b'

            logGBF = gv.mean(gv.gvar(results[name]['logGBF']))
            x = np.exp(logGBF - logGBF_max)

            #plt.axvline(x, ls='--', alpha=0.4)
            plt.scatter(x=x, y=y, color=color)
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
                color = 'r'
            elif 'xpt_' in name:
                color = 'c'
            elif 'xpt-ratio_' in name:
                color = 'g'
            elif 'ma-ratio_' in name:
                color = 'm'
            else:
                color = 'b'

            x = gv.mean(gv.gvar(results[name]['Q']))
            #plt.axvline(x, ls='--', alpha=0.4)
            plt.scatter(x=x, y=y, color=color)
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
                color = 'r'
            elif 'xpt_' in name:
                color = 'c'
            elif 'xpt-ratio_' in name:
                color = 'g'
            elif 'ma-ratio_' in name:
                color = 'm'
            else:
                color = 'b'

            x = gv.mean(gv.gvar(results[name]['chi2/df']))
            #plt.axvline(x, ls='--', alpha=0.4)
            plt.scatter(x=x, y=y, color=color)
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
