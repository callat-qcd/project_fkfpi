from __future__ import print_function
import os, sys
import numpy as np
import tables as h5
import lsqfit
import gvar as gv
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_style("ticks")
import chipt_awl as xpt

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

ens_long = {
    'a15m400'  :'l1648f211b580m0217m065m838',
    'a15m350'  :'l1648f211b580m0166m065m838',
    'a15m310'  :'l1648f211b580m013m065m838',
    'a15m220'  :'l2448f211b580m0064m0640m828',
    'a15m130'  :'l3248f211b580m00235m0647m831',
    'a15m135XL':'l4864f211b580m002426m06730m8447',
    'a12m400'  :'l2464f211b600m0170m0509m635',
    'a12m350'  :'l2464f211b600m0130m0509m635',
    'a12m310'  :'l2464f211b600m0102m0509m635',
    'a12m220'  :'l3264f211b600m00507m0507m628',
    'a12m220S' :'l2464f211b600m00507m0507m628',
    'a12m220L' :'l4064f211b600m00507m0507m628',
    'a12m130'  :'l4864f211b600m00184m0507m628',
    'a09m400'  :'l3264f211b630m0124m037m440',
    'a09m350'  :'l3264f211b630m00945m037m440',
    'a09m310'  :'l3296f211b630m0074m037m440',
    'a09m220'  :'l4896f211b630m00363m0363m430',
}
r_a = {
    'a15m400'  :gv.gvar(2.059,0.023),
    'a15m350'  :gv.gvar(2.059,0.023),
    'a15m310'  :gv.gvar(2.059,0.023),
    'a15m220'  :gv.gvar(2.073,0.013),
    'a15m130'  :gv.gvar(2.089,0.008),
    'a15m135XL':gv.gvar(2.089,0.008),
    'a12m400'  :gv.gvar(2.575,0.017),
    'a12m350'  :gv.gvar(2.575,0.017),
    'a12m310'  :gv.gvar(2.575,0.017),
    'a12m220S' :gv.gvar(2.585,0.019),
    'a12m220'  :gv.gvar(2.626,0.013),
    'a12m220L' :gv.gvar(2.614,0.009),
    'a12m130'  :gv.gvar(2.608,0.008),
    'a09m400'  :gv.gvar(3.499,0.024),
    'a09m350'  :gv.gvar(3.499,0.024),
    'a09m310'  :gv.gvar(3.499,0.024),
    'a09m220'  :gv.gvar(3.566,0.014)
    }

def format_h5_data(switches,data):
    x = dict()
    y = dict()
    p = dict()
    if switches['bs_bias']:
        print('Shifting BS data to boot0')

    for ens in switches['ensembles']:
        #print(ens)
        x[ens] = dict()
        data_dict = dict()
        for m in ['mpi','mk','mss','mju','mjs','mrs','mru']:
            data_dict[m] = data.get_node('/'+ens+'/'+m).read()
            #print('  %s = %.5f +- %.5f' %(m,data_dict[m].mean(),data_dict[m].std()))
        for f in ['FK','Fpi']:
            data_dict[f] = data.get_node('/'+ens+'/'+f).read()
        if switches['debug']:
            fkfpi = data_dict['FK']/data_dict['Fpi']
            print('  FK/Fpi = %.5f +- %.5f' %(fkfpi.mean(),fkfpi.std()))
        if switches['bs_bias']:
            data_bs = dict()
            for d in data_dict:
                data_bs[d] = data_dict[d][1:]
                if switches['debug']:
                    print(d,data_bs[d].shape)
            gvdata = gv.dataset.avg_data(data_bs,bstrap=True)
            if switches['debug']:
                gvdata_copy = dict(gvdata)
            for d in data_bs:
                gvdata[d] = gvdata[d] + (data_dict[d][0] - gvdata[d].mean)
                if switches['debug']:
                    print(d,gvdata[d].mean,gvdata_copy[d].mean)
        else:
            gvdata = gv.dataset.avg_data(data_dict,bstrap=True)
        if switches['debug']:
            print(gvdata)
        if switches['debug']:
            print('data_dict')
            print(data_dict)

        y[ens] = gvdata['FK']/gvdata['Fpi']
        print(ens,y[ens])

        x[ens]['mpiL'] = gvdata['mpi'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mkL']  = gvdata['mk'].mean  * data.get_node('/'+ens+'/L').read()
        x[ens]['mssL'] = gvdata['mss'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mjuL'] = gvdata['mju'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mjsL'] = gvdata['mjs'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mruL'] = gvdata['mru'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mrsL'] = gvdata['mrs'].mean * data.get_node('/'+ens+'/L').read()

        # MASSES
        p[(ens,'mpi')] = gvdata['mpi']
        p[(ens,'mk')]  = gvdata['mk']
        p[(ens,'mss')] = gvdata['mss']
        p[(ens,'mju')] = gvdata['mju']
        p[(ens,'mjs')] = gvdata['mjs']
        p[(ens,'mru')] = gvdata['mru']
        p[(ens,'mrs')] = gvdata['mrs']
        p[(ens,'Lchi_PP')] = 4 * np.pi * gvdata['Fpi']
        p[(ens,'Lchi_PK')] = 4 * np.pi * np.sqrt(gvdata['FK'] * gvdata['Fpi'])
        p[(ens,'Lchi_KK')] = 4 * np.pi * gvdata['FK']
        # HISQ params
        aw0 = data.get_node('/'+ens+'/aw0').read()
        p[(ens,'aw0')] = gv.gvar(aw0[0],aw0[1])
        a2di = data.get_node('/'+ens+'/a2DI').read()
        p[(ens,'a2DI')] = gv.gvar(a2di[0],a2di[1])
        x[ens]['alphaS'] = data.get_node('/'+ens+'/alpha_s').read()
    return {'x':x, 'y':y, 'p':p}

def fkfpi_phys(x_phys,fit):
    # use metaSq = 4/3 mK**2 - 1/3 mpi**2
    print('prediction from LQCD')
    #meta = np.sqrt(4./3 * x_phys['mka']**2 -1./3 * x_phys['mpi']**2)
    #fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
    #    meta=meta,Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    Fitc = xpt.Fit(switches={'ansatz':{'type':'xpt'}})
    for k in fit.p:
        if k in x_phys:
            pass
        else:
            x_phys[k] = fit.p[k]
    fkp = Fitc.fit_function({'Lchi':x_phys['Lchi']},x_phys)
    su2 = dsu2(FKpi=fkp,mpi=x_phys['mpi'],mk=x_phys['mka'],F0=x_phys['F0'])
    print('FK / Fpi = ',fkp)
    print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))
    # use meta = meta_pdg
    #print('prediction from LQCD + meta_PDG')
    #fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
    #    meta=x_phys['meta'],Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    #print('FK / Fpi = ',fkp)
    #print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))
    return fkp*np.sqrt(1+su2)

def dsu2(FKpi,mpi,mk,F0):
    R = gv.gvar(35.7,np.sqrt(1.9**2 + 1.8**2))
    eps_su2 = np.sqrt(3) / 4 / R
    d = np.sqrt(3) * eps_su2 * (\
        -4./3 * (FKpi-1) \
        + 2./(6*(4*np.pi)**2*F0**2)*(mk**2-mpi**2-mpi**2*np.log(mk**2/mpi**2))\
        )
    return d

def plot_data(ax,e,x,y):
    colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9'}
    shapes = {'m400':'h', 'm350':'p', 'm310':'s', 'm220':'^', 'm130':'o', 'm135':'*'}
    labels = {
        'a15m400':'', 'a15m350':'', 'a15m310':'', 'a15m220':'','a15m135XL':'a15',
        'a12m400':'', 'a12m350':'', 'a12m310':'', 'a12m220':'', 'a12m130':'a12',
        'a12m220L':'', 'a12m220S':'',
        'a09m400':'', 'a09m350':'', 'a09m310':'', 'a09m220':'a09',
    }
    c = colors[e.split('m')[0]]
    s = shapes['m'+e.split('m')[1][0:3]]
    ax.errorbar(x=x.mean,y=y.mean,xerr=x.sdev,yerr=y.sdev,
        marker=s,color=c,linestyle='None',label=labels[e])

    return ax

def perform_analysis(switches,priors):
    fit_results = dict()
    for base_model in switches['ansatz']['models']:
        for FPK in ['PP','PK','KK']:
            model = base_model +'_'+FPK
            switches['scale'] = FPK
            switches['ansatz']['model'] = model
            print('EFT: ',model)
            x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
            y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
            p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
            for key in priors:
                p_e[key] = priors[key]
            if not switches['default_priors']:
                try:
                    for key in ['p_4','k_4','s_4','saS_4']:
                        p_e[key] = gv.gvar(0,ip.nnlo_width[base_model][FPK][key])
                except Exception as e:
                    print('WARNING: non optimized NNLO prior widths')
                    print(str(e))
            else:
                print('using default prior widths')
            d_e = dict()
            d_e['x'] = x_e
            d_e['y'] = y_e
            d_e['p'] = p_e
            fit_e = xpt.Fit(switches,xyp_init=d_e)
            fit_e.fit_data()
            fit_results[model] = fit_e
            if switches['print_fit']:
                print(fit_e.fit.format(maxline=True))

            fit_e.report_phys_point(phys_point)

            if switches['make_plots']:
                epi_range = dict()
                epi_range['Lchi'] = phys_point['Lchi']
                epi_range['mk']   = phys_point['mk']
                epi_range['mpi']  = np.sqrt(np.arange(1,400**2,400**2 // 100))
                epi_range['mpi_phys'] = phys_point['mpi']

                fig_vs_epi = plt.figure('FKFpi_vs_epi_'+model)
                ax = plt.axes([.12, .12, .85, .85])
                ax = fit_e.vs_epi(epi_range, ax)
                y_shift = fit_e.shift_phys_epsK(phys_point)
                for e in switches['ensembles_fit']:
                    x = fit_e.fit.p[(e,'p2')]
                    y = y_e[e] + y_shift[e]
                    ax = plot_data(ax, e, x, y)
                ax.legend()
                ax.set_xlabel(r'$\epsilon_\pi^2$',fontsize=16)
                ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
                ax.set_xlim(0,.1)
                ax.set_ylim(1.05, 1.25)
                plt.savefig('figures/vs_episq_'+model+'.pdf',transparent=True)

    for base_model in switches['ansatz']['models']:
        for FPK in ['PP','PK','KK']:
            model = base_model +'_'+FPK
            L5_rho  = fit_results[model].fit.p['L5']
            L5_rho += 3./8 * 1/(4*np.pi)**2 * np.log(phys_point['Lchi']/770.)
            print('%25s: L5(m_rho) = %s' %(model,L5_rho))


def nnlo_prior_scan(switches,priors):
    logGBF_array = []
    model = switches['nnlo_priors_model']
    switches['ansatz']['model'] = model
    print(model,'Prior width study')
    x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
    y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
    p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
    for key in priors:
        p_e[key] = priors[key]
    d_e = dict()
    d_e['x'] = x_e
    d_e['y'] = y_e

    if switches['prior_group']:
        p_range = priors['p_range']
        a_range = priors['a_range']
        z = np.zeros([len(a_range),len(p_range)])
        tot = len(a_range) * len(p_range)
        i_t = 0
        for i_s,s4 in enumerate(a_range):
            p_e['s_4']   = gv.gvar(0,s4)
            if 'alphaS' in model:
                p_e['saS_4'] = gv.gvar(0,s4)
            for i_p,p4 in enumerate(p_range):
                p_e['p_4'] = gv.gvar(0,p4)
                p_e['k_4'] = gv.gvar(0,p4)
                d_e['p'] = p_e
                fit_e = xpt.Fit(switches,xyp_init=d_e)
                fit_e.fit_data()
                tmp = [ p_e[k].sdev for k in ['p_4','s_4'] ]
                tmp.append(fit_e.fit.logGBF)
                logGBF_array.append(tmp)
                sys.stdout.write('%4d out of %d\r' %(i_t,tot))
                sys.stdout.flush()
                i_t += 1
                z[i_s,i_p] = fit_e.fit.logGBF
    else:
        print('individual prior width study not supported [yet]')

    logGBF_array = np.array(logGBF_array)
    logGBF_max = np.argmax(logGBF_array[:,-1])
    print('optimal prior widths for %s' %model)
    logGBF_optimal = logGBF_array[logGBF_max]
    tmp = ''
    for i_k,k in enumerate(['p_4','s_4']):
        tmp += '%s = %.2f ' %(k,logGBF_optimal[i_k])
    print(tmp,'logGBF = ',logGBF_optimal[-1])

    lgbf = logGBF_array[:,-1]
    w = np.exp(lgbf - logGBF_optimal[-1])
    #w = w / w.sum()
    logGBF_w = np.copy(logGBF_array)
    logGBF_w[:,-1] = w

    if switches['prior_group']:
        z = np.exp(z - logGBF_optimal[-1])
        cmap_old = sns.cubehelix_palette(8, as_cmap=True)
        levels = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(sns.color_palette("BrBG_r", len(levels)).as_hex())
        fig, ax = plt.subplots(constrained_layout=True)
        CS = ax.contour(p_range,a_range,z,levels=levels, cmap=cmap)
        ax.clabel(CS, CS.levels[0::2], fmt='%2.1f', colors='k', fontsize=12)
        ax.plot(logGBF_optimal[0], logGBF_optimal[1],marker='X',color='r',markersize=20)
        ax.set_xlabel(r'$\sigma_{nnlo,\chi{\rm PT}}$',fontsize=16)
        ax.set_ylabel(r'$\sigma_{nnlo,a^2}$',fontsize=16)
        fig_name = 'prior_width_'+switches['ansatz']['model']+'_'+switches['scale']
        if not os.path.exists('figures'):
            os.makedirs('figures')
        plt.savefig('figures/'+fig_name+'.pdf',transparent=True)

def nlo_report(switches,priors,FV=True):
    latex = {
        'ma_nlo':'ma nlo',   'ma_nlo_FV':'ma nlo w/FV',
        'xpt_nlo':'xpt nlo', 'xpt_nlo_FV':'xpt nlo w/FV'}
    marker = {
        'ma_nlo':'s',  'ma_nlo_FV':'o',
        'xpt_nlo':'*', 'xpt_nlo_FV':'d'}
    color = {
        'ma_nlo':'r',  'ma_nlo_FV':'b',
        'xpt_nlo':'m', 'xpt_nlo_FV':'g'}
    if FV:
        models = ['ma_nlo','ma_nlo_FV','xpt_nlo','xpt_nlo_FV']
    else:
        models = ['ma_nlo','xpt_nlo']

    fit_results = dict()
    for model in models:
        switches['ansatz']['model'] = model
        print('EFT: ',model)
        model_result = dict()
        for e in switches['ensembles']:
            switches['ensembles_fit'] = [e]
            x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
            y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
            p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
            p_e['L5'] = priors['L5']
            d_e = dict()
            d_e['x'] = x_e
            d_e['y'] = y_e
            d_e['p'] = p_e
            fit_e = xpt.Fit(switches,xyp_init=d_e)
            fit_e.fit_data()
            model_result[e] = fit_e
        fit_results[model] = model_result
    fig_name = 'nlo_report'
    if FV:
        fig_name += '_FV'

    fig = plt.figure(fig_name)
    ax  = plt.axes([.13, .13, .84, .84])

    if len(models) == 2:
        print("%8s & %13s & %13s & \\\\" \
            %('ensemble',latex[models[0]], latex[models[1]]))
    else:
        print("%8s & %13s & %13s & %13s & %13s\\\\" \
            %('ensemble',latex[models[0]], latex[models[1]],latex[models[2]], latex[models[3]]))
    print("\\hline")
    for i_e,e in enumerate(switches['ensembles']):
        s = "%8s" %e
        for model in models:
            s += " & %13s" %str(fit_results[model][e].fit.p['L5'])
        s += "\\\\"
        print(s)

        for i_m,model in enumerate(models):
            if i_e == 0:
                label=latex[models[i_m]]
            else:
                label=''
            y = len(switches['ensembles']) - i_e + 0.1*i_m
            ax.errorbar(x=fit_results[model][e].fit.p['L5'].mean,y=y,
                xerr=fit_results[model][e].fit.p['L5'].sdev,linestyle='None',
                marker=marker[model],mfc='None',color=color[model],label=label)
    plt.yticks(np.arange(len(switches['ensembles']),0,-1),tuple(switches['ensembles']))
    ax.set_xlabel(r'$L_5$',fontsize=16)
    ax.legend(loc=1,fontsize=16)
    ax.set_xlim(-0.0002,0.001)
    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.savefig('figures/'+fig_name+'.pdf',transparent=True)

def fit_checker(switches,priors):
    fit_results = dict()
    for model in switches['ansatz']['models']:
        switches['ansatz']['model'] = model
        print('EFT: ',model)
        x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
        y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
        p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
        for key in priors:
            p_e[key] = priors[key]
        if not switches['default_priors']:
            try:
                for key in ['p_4','k_4','s_4','saS_4']:
                    p_e[key] = gv.gvar(0,ip.nnlo_width[model][switches['scale']][key])
            except Exception as e:
                print('WARNING: non optimized NNLO prior widths')
                print(str(e))
        else:
            print('using default prior widths')
        d_e = dict()
        d_e['x'] = x_e
        d_e['y'] = y_e
        d_e['p'] = p_e
        fit_e = xpt.Fit(switches,xyp_init=d_e)
        fit_e.fit_data()
        fit_results[model] = fit_e
        if switches['print_fit']:
            print(fit_e.fit.format(maxline=True))

    print('model & PK& chisq/dof& logGBF& $L_5$& $k_4$& $p_4$& $s_4$& $F_K/F_\pi$\\\\')
    for model in switches['ansatz']['models']:
        chi2dof,dof,logGBF,L5,k4,p4,s4,FKFpi = fit_results[model].check_fit(phys_point)
        print('%34s& %s& %.3f [%2d]& %.2f& %14s& %8s& %8s& %8s& %10s\\\\'\
            %(model.replace('_','\_'),switches['scale'],chi2dof,dof,logGBF,L5,k4,p4,s4,FKFpi))


if __name__ == "__main__":
    import input_params as ip
    print("python version    :", sys.version)
    #print("pandas version:", pd.__version__)
    print("numpy  version    :", np.__version__)
    print("matplotlib version:", matplotlib.__version__)
    print("gvar   version    :", gv.__version__)
    print("lsqfit version    :", lsqfit.__version__)
    print('')

    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.ion()

    # Load input params
    switches   = ip.switches
    priors     = ip.priors
    phys_point = ip.phys_point

    # Load data
    data  = h5.open_file('FK_Fpi_data.h5','r')
    gv_data = format_h5_data(switches,data)
    data.close()

    if switches['do_analysis']:
        perform_analysis(switches,priors)

    if switches['check_fit']:
        fit_checker(switches,priors)

    if switches['nnlo_priors']:
        nnlo_prior_scan(switches,priors)

    if switches['nlo_fv_report']:
        nlo_report(switches,priors,FV=True)
    if switches['nlo_report']:
        nlo_report(switches,priors,FV=False)

    plt.ioff()
    if run_from_ipython():
        plt.show(block=False)
    else:
        plt.show()
