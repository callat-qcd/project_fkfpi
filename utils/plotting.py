import matplotlib.pyplot as plt
import numpy as np
import os, copy
import gvar as gv
# import chipt lib for fit functions
import chipt


class ExtrapolationPlots:

    def __init__(self, model, model_list, fitEnv, fit_result, switches):
        if not os.path.exists('figures'):
            os.makedirs('figures')
        # Figure formatting
        self.fig_width = 6.75 # in inches, 2x as wide as APS column
        self.gr        = 1.618034333 # golden ratio
        self.fig_size  = (self.fig_width, self.fig_width / self.gr)
        self.fig_size2 = (self.fig_width, self.fig_width * 1.6)
        self.plt_axes  = [0.14,0.14,0.855,0.855]
        self.fs_text   = 20 # font size of text
        self.fs_leg    = 16 # legend font size
        self.mrk_size  = '5' # marker size
        self.tick_size = 16 # tick size
        self.lw        = 1 # line width

        self.colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9', 'a06':'#00FFFF'}
        self.shapes = {'m400':'h', 'm350':'p', 'm310':'s', 'm220':'^', 'm130':'o', 'm135':'*'}
        self.labels = {
            'a15m400':'', 'a15m350':'', 'a15m310':'a15', 'a15m220':'','a15m135XL':'',
            'a12m400':'', 'a12m350':'', 'a12m310':'a12', 'a12m220':'', 'a12m130':'',
            'a12m220L':'', 'a12m220S':'',
            'a09m400':'', 'a09m350':'', 'a09m310':'a09', 'a09m220':'', 'a09m135':'',
            'a06m310L':'a06',
            }
        self.dx_cont = {
            'a15m400'  :0.0016, 'a12m400':0.0016, 'a09m400':0.0016,
            'a15m350'  :0.0008, 'a12m350':0.0008, 'a09m350':0.0008,
            'a15m310'  :0.,     'a12m310':0.,     'a09m310':0.,     'a06m310L':0.,
            'a15m220'  :-0.0008,'a12m220':-0.0008,'a09m220':-0.0008,
            'a15m135XL':-.0016, 'a12m130':-0.0016,'a09m135':-0.0016,
            'a12m220L' :-0.0012,'a12m220S':-0.0004,
        }

        self.model      = model
        self.FF         = model.split('_')[-1]
        self.fv         = 'FV' in model
        self.model_list = model_list
        self.fitEnv     = fitEnv
        self.fit_result = fit_result
        #for k, v in self.fit_results.items():
        #    print(k,v)
        self.switches   = switches

        # create fit functions for original and shifted points
        self.og_fit     = chipt.FitModel(self.model_list, _fv=self.fv, _FF=self.FF)
        self.shift_list = list(model_list)
        if 'ma' in model_list[0]:
            self.shift_list[0] = self.shift_list[0].replace('ma','xpt')
        self.shift_fit  = chipt.FitModel(self.shift_list, _fv=False, _FF=self.FF)


    def plot_vs_eps_asq(self,shift_points):
        self.shift_xp = copy.deepcopy(shift_points)
        for k in self.fit_result.p:
            if isinstance(k,str):
                self.shift_xp['p'][k] = self.fit_result.p[k]
        y_plot = []
        x_plot = []
        a_range = np.sqrt(np.arange(0, .16**2, .16**2 / 50))
        for a_fm in a_range:
            self.shift_xp['p']['aw0'] = a_fm / self.shift_xp['p']['w0']
            x_plot.append(self.shift_xp['p']['aw0']**2 / 4 / np.pi)
            y_plot.append(self.fitEnv._fit_function(self.shift_fit, self.shift_xp['x'], self.shift_xp['p']))
        x  = np.array([k.mean for k in x_plot])
        y  = np.array([k.mean for k in y_plot])
        dy = np.array([k.sdev for k in y_plot])

        if self.switches['milc_compare']:
            figsize = self.fig_size2
        else:
            figsize = self.fig_size

        self.fig_cont = plt.figure('FKFpi_vs_ea_'+self.model,figsize=figsize)
        self.ax_cont  = plt.axes(self.plt_axes)
        self.ax_cont.fill_between(x, y-dy, y+dy, color='#b36ae2', alpha=0.4)

        self.plot_data(p_type='ea')
        handles, labels = self.ax_cont.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        self.ax_cont.legend(handles, labels, ncol=4, fontsize=self.fs_leg)

        self.ax_cont.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=self.fs_text)
        self.ax_cont.set_ylabel(r'$F_K / F_\pi$',fontsize=self.fs_text)
        self.ax_cont.set_xlim(0,.065)
        self.ax_cont.set_ylim(1.135, 1.218)

        if self.switches['save_figs']:
            plt.savefig('figures/'+'FKFpi_vs_ea_'+self.model+'.pdf',transparent=True)

    def plot_vs_eps_pi(self,shift_points):
        self.shift_xp = copy.deepcopy(shift_points)
        eps_pisq_phys = (gv.gvar(self.shift_xp['p']['mpi'] / self.shift_xp['p']['Lchi_'+self.FF]))**2
        for k in self.fit_result.p:
            if isinstance(k,str):
                self.shift_xp['p'][k] = self.fit_result.p[k]
        y_plot = dict()
        y_plot['a15'] = []
        y_plot['a12'] = []
        y_plot['a09'] = []
        y_plot['a06'] = []
        y_plot['a00'] = []
        x_plot = []
        mpi_range = np.sqrt(np.arange(100, 411**2, 411**2/200))
        for a_mpi in mpi_range:
            x_plot.append(a_mpi**2 / (self.shift_xp['p']['Lchi_'+self.FF])**2)
            self.shift_xp['p']['mpi'] = a_mpi
            self.shift_xp['p']['aw0'] = 0
            y_plot['a00'].append(self.fitEnv._fit_function(self.shift_fit, self.shift_xp['x'], self.shift_xp['p']))
            for aa in ['a15','a12','a09','a06']:
                if aa == 'a06':
                    self.shift_xp['p']['aw0'] = self.fit_result.p[(aa+'m310L','aw0')]
                else:
                    self.shift_xp['p']['aw0'] = self.fit_result.p[(aa+'m310','aw0')]
                y_plot[aa].append(self.fitEnv._fit_function(self.shift_fit, self.shift_xp['x'], self.shift_xp['p']))
        x = np.array([k.mean for k in x_plot])
        y  = dict()
        dy = dict()
        y['a00']  = np.array([k.mean for k in y_plot['a00']])
        dy['a00'] = np.array([k.sdev for k in y_plot['a00']])
        for aa in ['a15','a12','a09','a06']:
            y[aa]  = np.array([k.mean for k in y_plot[aa]])
            dy[aa] = np.array([k.sdev for k in y_plot[aa]])

        self.fig_x = plt.figure('FKFpi_vs_epi_'+self.model, figsize=self.fig_size)
        self.ax_x  = plt.axes(self.plt_axes)
        self.ax_x.fill_between(x, y['a00']-dy['a00'], y['a00']+dy['a00'], color='#b36ae2',alpha=0.4)
        for aa in ['a15','a12','a09','a06']:
            self.ax_x.plot(x, y[aa], color=self.colors[aa])

        # plot physical eps_pi**2
        self.ax_x.axvline(eps_pisq_phys.mean,linestyle='--',color='#a6aaa9')
        self.ax_x.axvspan(eps_pisq_phys.mean -eps_pisq_phys.sdev, eps_pisq_phys.mean +eps_pisq_phys.sdev,
            alpha=0.4, color='#a6aaa9')
        # plot data
        self.plot_data(p_type='epi')
        handles, labels = self.ax_x.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        self.ax_x.legend(handles, labels, ncol=1, fontsize=self.fs_leg)

        # labels
        if self.FF == 'PP':
            self.ax_x.set_xlabel(r'$\epsilon_\pi^2 = (m_\pi / 4\pi F_\pi)^2$',fontsize=self.fs_text)
            self.ax_x.set_xlim(0,.094)
        elif self.FF == 'PK':
            self.ax_x.set_xlabel(r'$\epsilon_\pi^2 = (m_\pi / 4\pi)^2 /(F_\pi F_K)$',fontsize=self.fs_text)
            self.ax_x.set_xlim(0,.088)
        elif self.FF == 'KK':
            self.ax_x.set_xlabel(r'$\epsilon_\pi^2 = (m_\pi / 4\pi F_K)^2$',fontsize=self.fs_text)
            self.ax_x.set_xlim(0,.084)
        self.ax_x.set_ylabel(r'$F_K / F_\pi$',fontsize=self.fs_text)
        self.ax_x.set_ylim(1.02, 1.218)

        if self.switches['save_figs']:
            plt.savefig('figures/'+'FKFpi_vs_epi_'+self.model+'.pdf',transparent=True)


    def plot_data(self, p_type, offset=False, raw=False):
        y_shift = self.shift_data(p_type=p_type)
        for a_ens in self.switches['ensembles']:
            if a_ens in self.switches['ensembles_fit']:
                c = self.colors[a_ens.split('m')[0]]
            else:
                c = 'k'
            s = self.shapes['m'+a_ens.split('m')[1][0:3]]
            if p_type == 'ea':
                x  = self.fit_result.p[(a_ens, 'aw0')]**2 / 4 / np.pi
                dx = self.dx_cont[a_ens]
            elif p_type == 'epi':
                x  = (self.fit_result.p[(a_ens, 'mpi')] / self.fit_result.p[(a_ens, 'Lchi_'+self.FF)])**2
                dx = 0
            label = self.labels[a_ens]
            if p_type == 'ea':
                mfc = c
            elif p_type == 'epi':
                mfc = 'None'
            y = self.fit_result.y[a_ens] + y_shift[a_ens]
            if p_type == 'ea':
                self.ax_cont.errorbar(x=x.mean+dx, y=y.mean,xerr=x.sdev, yerr=y.sdev,
                    marker=s, color=c, mfc=mfc, linestyle='None', label=label)
            elif p_type == 'epi':
                self.ax_x.errorbar(x=x.mean+dx, y=y.mean,xerr=x.sdev, yerr=y.sdev,
                    marker=s, color=c, mfc=mfc, linestyle='None', label=label)

    def shift_data(self, p_type):
        y_shift = dict()
        if self.switches['debug_shift']:
            print('%9s   %11s   y_shift[ens]' %('ensemble', 'y[ens]'))
            print('---------------------------------------------------------------')
        for a_ens in self.switches['ensembles']:
            og_priors = dict()
            for k, v in self.fitEnv.p.items():
                if type(k) == tuple and k[0] == a_ens:
                    #print(a_ens,k,v)
                    if k in self.fit_result.p:
                    #if a_ens in self.switches['ensembles_fit']:
                        og_priors[k[1]] = self.fit_result.p[k] # the x-params which are priors
                    else:
                        # if it is from an excluded ensemble - get from original data
                        og_priors[k[1]] = v
            for k in self.fit_result.p:
                if isinstance(k,str):# grab the LECs from the fit results
                    og_priors[k] = self.fit_result.p[k]    # the LECs of the fit
            self.shift_xp['p']['aw0'] = self.fit_result.p[(a_ens,'aw0')]
            if p_type == 'epi':
                self.shift_xp['p']['aw0'] = self.fit_result.p[(a_ens,'aw0')]
                self.shift_xp['p']['mpi'] = self.fit_result.p[(a_ens,'mpi')] / self.fit_result.p[(a_ens, 'Lchi_'+self.FF)]
                self.shift_xp['p']['mk']  = self.shift_xp['p']['mk'] / self.shift_xp['p']['Lchi_'+self.FF]
                self.shift_xp['p']['Lchi_'+self.FF] = 1
                self.shift_xp['x']['alphaS'] = self.fit_result.x[a_ens]['alphaS']
            og_y    = self.fitEnv._fit_function(self.og_fit,    self.fitEnv.x[a_ens], og_priors)
            shift_y = self.fitEnv._fit_function(self.shift_fit, self.shift_xp['x'],   self.shift_xp['p'])
            y_shift[a_ens] = shift_y - og_y
            if self.switches['debug_shift']:
                print('%9s   %11s   %s' %(a_ens, og_y, shift_y))
        return y_shift



def plot_vs_eps_asq(model, fitEnv, fit_result, shift_point, model_list, FF, fv, switches):
    shift_xp   = dict(shift_point)
    shift_list = list(model_list)
    if 'ma' in shift_list[0]:
        shift_list[0] = shift_list[0].replace('ma','xpt')
    shift_fit  = chipt.FitModel(shift_list, _fv=False, _FF=FF)
    y_plot = []
    x_plot = []
    # make equally space in a**2 points
    a_range = np.sqrt(np.arange(0, .16**2, .16**2 / 50))
    for a_fm in a_range:
        shift_xp['p']['aw0'] = a_fm / shift_xp['p']['w0']
        x_plot.append(shift_xp['p']['aw0']**2 / 4 / np.pi)
        y_plot.append(fitEnv._fit_function(shift_fit, shift_xp['x'],   shift_xp['p']))
    x  = np.array([k.mean for k in x_plot])
    y  = np.array([k.mean for k in y_plot])
    dy = np.array([k.sdev for k in y_plot])

    fig = plt.figure('FKFpi_vs_ea_'+model,figsize=fig_size)
    ax  = plt.axes(plt_axes)
    ax.fill_between(x, y-dy, y+dy, color='#b36ae2', alpha=0.4)

    ax = plot_data_phys_mass(ax, fitEnv, fit_result, shift_point, model_list, FF, fv, switches)

    ax.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=fs_text)
    ax.set_ylabel(r'$F_K / F_\pi$',fontsize=fs_text)
    ax.set_xlim(0,.065)
    ax.set_ylim(1.14, 1.228)


def plot_data_phys_mass(ax, fitEnv, fit_result, shift_point, model_list, FF, fv, switches):
    og_fit     = chipt.FitModel(model_list, _fv=fv, _FF=FF)
    shift_xp   = dict(shift_point)
    shift_list = list(model_list)
    if 'ma' in shift_list[0]:
        shift_list[0] = shift_list[0].replace('ma','xpt')
    shift_fit  = chipt.FitModel(shift_list, _fv=False, _FF=FF)
    y_shift = dict()
    if switches['debug_shift']:
        print('%9s   %11s   y_shift[ens]' %('ensemble', 'y[ens]'))
        print('---------------------------------------------------------------')
    for a_ens in fit_result.x:
        # populate original priors
        og_priors    = dict()
        for k, v in fit_result.p.items():
            if type(k) == tuple and k[0] == a_ens:
                og_priors[k[1]] = v # the x-params which are priors
            else:
                og_priors[k] = v    # the LECs of the fit
                shift_xp['p'][k] = v
        # set shifted a/w0 for a given ensemble
        shift_xp['p']['aw0'] = fit_result.p[(a_ens,'aw0')]
        og_y    = fitEnv._fit_function(og_fit,    fit_result.x[a_ens], og_priors)
        shift_y = fitEnv._fit_function(shift_fit, shift_xp['x'],   shift_xp['p'])
        y_shift[a_ens] = shift_y.mean - og_y.mean
        if switches['debug_shift']:
            print('%9s   %11s   %s' %(a_ens, og_y, shift_y))
        # plot data
        s = shapes['m'+a_ens.split('m')[1][0:3]]
        x = shift_xp['p']['aw0']**2 / 4 / np.pi
        label = labels[a_ens]
        if a_ens in switches['ensembles_fit']:
            c = colors[a_ens.split('m')[0]]
        else:
            c = 'k'
        mfc = c
        dx = dx_cont[a_ens]

        ax.errorbar(x=x.mean+dx, y=fit_result.y[a_ens].mean+y_shift[a_ens], xerr=x.sdev, yerr=fit_result.y[a_ens].sdev,
                    marker=s, color=c, mfc=mfc, linestyle='None', label=label)

    return ax
