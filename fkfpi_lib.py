from __future__ import print_function
import os, sys, shutil
import numpy as np
import scipy.stats as stats
import tables as h5
import pandas as pd
import lsqfit
import gvar as gv
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_style("ticks")
import yaml
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

a_ens = {
    'a15':gv.gvar(0.1509,0.0013),
    'a12':gv.gvar(0.1206,0.0010),
    'a09':gv.gvar(0.0875,0.0008),
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
        print("%9s %s" %(ens,y[ens]))

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
        p[(ens,'Lchi_F0')] = 4 * np.pi * gv.gvar(80,1) / 197.3 * a_ens[ens[0:3]]
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
        + 2. * (mk**2 -mpi**2 -mpi**2*np.log(mk**2/mpi**2)) / (6*(4*np.pi)**2*F0**2)
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

def perform_analysis(switches,priors,phys_point):
    fit_results = dict()
    for base_model in switches['ansatz']['models']:
        if switches['optimized_priors']:
            if not os.path.exists('data/saved_prior_search.yaml'):
                print('ERROR: you asked for optimized priors but')
                print('data/saved_prior_search.yaml')
                print('does not exist - exiting')
                sys.exit()
            with open('data/saved_prior_search.yaml','r') as fin:
                prior_grid = yaml.safe_load(fin.read())
        for FPK in switches['scales']:
            model = base_model +'_'+FPK
            switches['scale'] = FPK
            switches['ansatz']['model'] = model
            print('EFT: ',model)
            x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
            y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
            p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
            for key in priors:
                p_e[key] = priors[key]
            if switches['optimized_priors']:
                df = pd.DataFrame(prior_grid[model])
                sp_max = df.stack().idxmax()
                print('      setting prior widths:')
                print('      (s_4, saS_4) = %s; (p_4, k_4) = %s' %(sp_max[0],sp_max[1]))
                for key in ['s_4','saS_4']:
                    if key in p_e:
                        p_e[key] = gv.gvar(0,float(sp_max[0]))
                for key in ['p_4','k_4']:
                    if key in p_e:
                        p_e[key] = gv.gvar(0,float(sp_max[1]))
            else:
                print('      using default prior widths')
            d_e = dict()
            d_e['x'] = x_e
            d_e['y'] = y_e
            d_e['p'] = p_e
            fit_e = xpt.Fit(switches,xyp_init=d_e)
            fit_e.fit_data()
            fit_results[model] = fit_e
            if switches['print_fit']:
                print(fit_e.fit.format(maxline=True))
            if switches['debug_phys']:
                fit_e.report_phys_point(phys_point)['phys']

            print('DEBUG: error budget')
            tmp = fit_e.report_phys_point(phys_point)['phys']
            Lchi_phys = phys_point['Lchi_'+FPK]
            p_phys = dict()
            p_phys[('phys','p2')] = phys_point['mpi']**2 / Lchi_phys**2
            p_phys[('phys','k2')] = phys_point['mk']**2 / Lchi_phys**2
            p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
            print(tmp.partialsdev(fit_e.fit.y))
            print(tmp.partialsdev(fit_e.fit.y,p_phys[('phys','p2')]))
            print(tmp.partialsdev(fit_e.fit.y,p_phys[('phys','p2')],p_phys[('phys','k2')]))
            print(tmp.partialsdev(fit_e.fit.y,p_phys[('phys','p2')],p_phys[('phys','k2')],p_phys[('phys','e2')]))

            if switches['make_plots']:
                if base_model in ['xpt_nnlo_FV','xpt_nnlo_FV_a4']:#, 'xpt_nnlo_FV_logSq', 'xpt-ratio_nnlo_FV', 'xpt-ratio_nnlo_FV_logSq']:
                    ea_range = dict()
                    ea_range['Lchi'] = phys_point['Lchi_'+FPK]
                    ea_range['mpi']  = phys_point['mpi']
                    ea_range['mk']   = phys_point['mk']
                    ea_range['a']    = np.sqrt(np.arange(0,.16**2,.16**2/100))
                    ea_range['w0']   = phys_point['w0']

                    fig_vs_ea = plt.figure('FKFpi_vs_ea_'+model)
                    ax = plt.axes([.12, .12, .85, .85])
                    #print(fit_e.fit.format(maxline=True))
                    ax = fit_e.vs_ea(ea_range, ax)
                    # add data
                    y_shift = fit_e.shift_phys_mass(phys_point)
                    for e in switches['ensembles_fit']:
                        x = fit_e.fit.p[(e,'a2')]
                        y = y_e[e] + y_shift[e]
                        #print(e,y_e[e],y_shift[e])
                        ax = plot_data(ax, e, x, y)
                    ax.legend(ncol=3)
                    ax.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=16)
                    ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
                    ax.set_xlim(0,.065)
                    ax.set_ylim(1.135, 1.205)
                    plt.savefig('figures/vs_epasq_'+model+'.pdf',transparent=True)
                if False:
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
    if not switches['debug']:
        model_avg,model_var = bayes_model_avg(switches,fit_results,phys_point)
        #dsu2(FKpi,mpi,mk,F0)
        print('SU(2) correction',dsu2(model_avg,phys_point['mpi'],phys_point['mk'],phys_point['F0']))

    '''
    for base_model in switches['ansatz']['models']:
        for FPK in switches['scales']:
            model = base_model +'_'+FPK
            if fit_results[model].
            L5_rho  = fit_results[model].fit.p['L5']
            L5_rho += 3./8 * 1/(4*np.pi)**2 * np.log(phys_point['Lchi']/770.)
            print('%25s: L5(m_rho) = %s' %(model,L5_rho))
    '''
    return fit_results

def bayes_model_avg(switches,results,phys_point):
    logGBF_list = []
    r_list      = []
    w_list      = []
    models      = []
    FKFpi       = 0.
    pdf         = 0
    cdf         = 0
    pdf_pp      = 0
    pdf_pk      = 0
    pdf_kk      = 0
    x = np.arange(1.15,1.2101,.0001)
    for model in results:
        logGBF_list.append(results[model].fit.logGBF)
        FPK=model.split('_')[-1]
        results[model].switches['scale'] = FPK
        r_list.append(results[model].report_phys_point(phys_point)['phys'])
        if switches['debug']:
            print('\nDEBUG:',model,results[model].switches['scale'])
    print('\nAll Models')
    print("%33s %5s %8s %s" %('model','Q','w','FK/Fpi'))
    print('-------------------------------------------------------------------')
    for i_m,model in enumerate(results):
        models.append(model)
        r = r_list[i_m]
        w = np.exp(logGBF_list[i_m])
        w = w / np.sum(np.exp(np.array(logGBF_list)))
        w_list.append(w)
        print("%33s %.3f %.2e %s" %(model,results[model].fit.Q,w,r))
        FKFpi += gv.gvar(w*r.mean, np.sqrt(w)*r.sdev)
        p = stats.norm.pdf(x,r.mean,r.sdev)
        pdf += w * p
        if '_PP' in model:
            pdf_pp += w * p
        elif '_PK' in model:
            pdf_pk += w * p
        elif '_KK' in model:
            pdf_kk += w * p
        c = stats.norm.cdf(x,r.mean,r.sdev)
        cdf += w * c
    w_list = np.array(w_list)
    r_list = np.array(r_list)
    model_var = np.sum(w_list * np.array([ r.mean**2 for r in r_list])) - FKFpi.mean**2

    print('\nModels with Q > 0.05')
    print("%33s %4s %5s %s" %('model','Q','w','FK/Fpi'))
    print('-------------------------------------------------------------------')
    for i_m,model in enumerate(models):
        if results[model].fit.Q > 0.05:
            print("%33s %.2f %.3f %s" %(model,results[model].fit.Q,w_list[i_m],r_list[i_m]))
        else:
            if w_list[i_m] > 0.05:
                print("%33s %.2f %.3f %s" %(model,results[model].fit.Q,w_list[i_m],'NEGLECTED'))

    print('\nFull average with %d Models' %len(results))
    print("%s +- %.4f" %(FKFpi,np.sqrt(model_var)))

    fig = plt.figure('hist')
    ax = plt.axes([0.11,0.12,0.87,0.87])
    x0 = FKFpi.mean
    dx = np.sqrt(model_var + FKFpi.sdev**2)
    #ax.axvspan(x0-dx,x0+dx,color='k',alpha=0.05)
    #ax.axvline(x0-dx,color='k',linestyle='-.')
    #ax.axvline(x0+dx,color='k',linestyle='-.')
    #dx = FKFpi.sdev
    #ax.axvspan(x0-dx,x0+dx,color='k',alpha=0.2)
    #ax.axvline(x0,color='k')
    #ax.axvline(x0-dx,color='k',linestyle='--')
    #ax.axvline(x0+dx,color='k',linestyle='--')
    ax.plot(x,pdf,color='k')
    ax.fill_between(x=x,y1=pdf,color='k',alpha=0.1)
    # 95%
    lidx95 = abs(cdf-0.025).argmin()
    uidx95 = abs(cdf-0.975).argmin()
    ax.fill_between(x=x[lidx95:uidx95],y1=pdf[lidx95:uidx95],facecolor='k',edgecolor='k',alpha=0.1)
    #68%
    lidx68 = abs(cdf-0.158655254).argmin()
    uidx68 = abs(cdf-0.841344746).argmin()
    ax.fill_between(x=x[lidx68:uidx68],y1=pdf[lidx68:uidx68],facecolor='k',edgecolor='k',alpha=0.1)
    # black lines
    ax.errorbar(x=[x[lidx95],x[lidx95]],y=[0,pdf[lidx95]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[uidx95],x[uidx95]],y=[0,pdf[uidx95]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[lidx68],x[lidx68]],y=[0,pdf[lidx68]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[uidx68],x[uidx68]],y=[0,pdf[uidx68]],color='k')#,lw=0.5)

    # split Fpi, FK and Fpi FK results
    ax.fill_between(x=x,y1=pdf_kk,color='b',alpha=0.6,label=r'$F^2 \rightarrow F_K^2$')
    ax.fill_between(x=x,y1=pdf_pk,color='g',alpha=0.6,label=r'$F^2 \rightarrow F_\pi F_K$')
    ax.fill_between(x=x,y1=pdf_pp,color='r',alpha=0.6,label=r'$F^2 \rightarrow F_\pi^2$')
    ax.set_xlim([1.1575,1.2075])
    ax.set_ylim(ymin=0)
    ax.set_xlabel(r'$F_K / F_\pi$',fontsize=16)
    ax.legend()
    ax.set_ylabel(r'Bayes Model Avg PDF',fontsize=16)
    if switches['optimized_priors']:
        ax.text(0.05,0.95, 'optimized prior widths',verticalalignment='center',horizontalalignment='left',transform=ax.transAxes,fontsize=16)
        plt.savefig('figures/model_avg_hist_logGBF_optimal_priors.pdf',transpareent=True)
    else:
        ax.text(0.05,0.95, 'default prior widths',verticalalignment='center',horizontalalignment='left',transform=ax.transAxes,fontsize=16)
        plt.savefig('figures/model_avg_hist_default_priors.pdf',transpareent=True)

    return FKFpi, model_var

def bma(switches,result,isospin):
    # read Bayes Factors
    logGBF_list = []
    for a in switches['ansatz']['type']:
        logGBF_list.append(result[a]['fit'].logGBF)
    # initiate a bunch of parameters
    # gA
    gA = 0
    gA_lst = []
    gA_dict = dict()
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
    x = np.linspace(1.222,1.352,13000)
    # error breakdown for each model
    model_error = error_budget(switches,result)
    model_budget = {k:0 for k in model_error[list(model_error.keys())[0]]['std'].keys()}
    for a in switches['ansatz']['type']:
        r = result[a]['phys']['result']
        gA_dict[a] = r
        w = 1/sum(np.exp(np.array(logGBF_list)-result[a]['fit'].logGBF))
        sqrtw = np.sqrt(w) # sqrt scales the std dev correctly
        wd[a] = w
        w_lst.append(w)
        gA += gv.gvar(w*r.mean,sqrtw*r.sdev)
        gA_lst.append(r.mean)
        p = stats.norm.pdf(x,r.mean,r.sdev)
        pdf += w*p
        pdfdict[a] = w*p
        c = stats.norm.cdf(x,r.mean,r.sdev)
        cdf += w*c
        cdfdict[a] = w*c
        # error breakdown
        model_std = model_error[a]['std']
        model_budget = {k:model_budget[k]+w*model_std[k]**2 for k in model_std} # variance breakdown of model average
    gA_lst = np.array(gA_lst)
    w_lst = np.array(w_lst)
    model_var = np.sum(w_lst*gA_lst**2) - gA.mean**2
    final_error = np.sqrt(gA.sdev**2 + isospin**2)
    model_budget['isospin'] = isospin**2
    model_budget['model'] = model_var
    model_budget['total'] = model_budget['total']+model_budget['isospin']+model_budget['model'] # add in quadrature isospin and model variance
    pct_budget = {k:[np.sqrt(model_budget[k])/gA.mean*100] for k in model_budget} # percent breakdown of model average
    error = {'E(gA)': gA.mean, 's(gA)': final_error, 's(Mk)': np.sqrt(model_var), 'weights': wd, 'error_budget': model_budget, 'pct_budget': pct_budget, 'gA_dict':gA_dict}
    plot_params = {'x':x, 'pdf':pdf, 'pdfdict':pdfdict, 'cdf':cdf, 'cdfdict':cdfdict}
    return error, plot_params

def nnlo_prior_scan(switches,priors):
    for base_model in switches['ansatz']['models']:
        for FPK in switches['scales']:
            if os.path.exists('data/saved_prior_search.yaml'):
                with open('data/saved_prior_search.yaml','r') as fin:
                    prior_grid = yaml.safe_load(fin.read())
                shutil.copyfile('data/saved_prior_search.yaml','data/saved_prior_search.yaml.bak')
            else:
                prior_grid = dict()
            model = base_model +'_'+FPK
            switches['scale'] = FPK
            switches['ansatz']['model'] = model
            print('Prior width study: ',model)

            logGBF_array = []

            if model not in prior_grid:
                prior_grid[model] = dict()
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
                if switches['refine_prior'] and len(prior_grid[model]) > 0:
                    df = pd.DataFrame(prior_grid[model])
                    print('    maximum logGBF = %.4f' %(df.stack().max()))
                    sp_max = df.stack().idxmax()
                    print('    (s_4, p_4) = ', sp_max[0],sp_max[1])
                    new_s = []
                    sf = float("%.1f" %float(sp_max[0]))
                    for i in np.arange(sf-0.5,sf+0.5,.1):
                        ii = float("%.1f" %i)
                        if ii not in a_range and ii > 0.:
                            new_s.append(ii)
                    new_s = np.array(new_s)
                    a_range = np.concatenate((a_range,new_s))
                    new_p = []
                    pf = float(sp_max[1])
                    for i in np.arange(pf-0.5,pf+0.5,.1):
                        ii = float("%.1f" %i)
                        if ii not in p_range and ii > 0.:
                            new_p.append(ii)
                    new_p = np.array(new_p)
                    p_range = np.concatenate((p_range,new_p))
                    print(a_range)
                    print(p_range)
                    #sys.exit()
                z = np.zeros([len(a_range),len(p_range)])
                tot = len(a_range) * len(p_range)
                i_t = 0
                for i_s,s4 in enumerate(a_range):
                    s4_s = str(s4)
                    if s4_s not in prior_grid[model]:
                        prior_grid[model][s4_s] = dict()
                    p_e['s_4']   = gv.gvar(0,s4)
                    if 'alphaS' in model:
                        p_e['saS_4'] = gv.gvar(0,s4)
                    for i_p,p4 in enumerate(p_range):
                        p4_s = str(p4)
                        sys.stdout.write('%4d out of %d, s_4 = %s p_4 = %s\r' %(i_t,tot,s4_s,p4_s))
                        sys.stdout.flush()
                        if p4_s not in prior_grid[model][s4_s]:
                            p_e['p_4'] = gv.gvar(0,p4)
                            p_e['k_4'] = gv.gvar(0,p4)
                            d_e['p'] = p_e
                            fit_e = xpt.Fit(switches,xyp_init=d_e)
                            fit_e.fit_data()
                            tmp = [ p_e[k].sdev for k in ['p_4','s_4'] ]
                            tmp.append(fit_e.fit.logGBF)
                            logGBF_array.append(tmp)
                            prior_grid[model][s4_s][p4_s] = float(fit_e.fit.logGBF)
                            i_t += 1
                            z[i_s,i_p] = fit_e.fit.logGBF
                        else:
                            tmp = [ p4, s4 ]
                            tmp.append(prior_grid[model][s4_s][p4_s])
                            logGBF_array.append(tmp)
                            i_t += 1
                            z[i_s,i_p] = prior_grid[model][s4_s][p4_s]
            else:
                print('individual prior width study not supported [yet]')

            df = pd.DataFrame(prior_grid[model])
            print('    maximum logGBF', df.stack().max())
            sp_max = df.stack().idxmax()
            print('    (s_4, p_4) = ', sp_max[0],sp_max[1])
            logGBF_array = np.array(logGBF_array)
            logGBF_max = np.argmax(logGBF_array[:,-1])
            print('optimal prior widths for %s' %model)
            logGBF_optimal = logGBF_array[logGBF_max]
            tmp = ''
            for i_k,k in enumerate(['p_4','s_4']):
                tmp += '%s = %.2f ' %(k,logGBF_optimal[i_k])
            print(tmp,'logGBF = %.4f\n' %(logGBF_optimal[-1]))

            lgbf = logGBF_array[:,-1]
            w = np.exp(lgbf - logGBF_optimal[-1])
            #w = w / w.sum()
            logGBF_w = np.copy(logGBF_array)
            logGBF_w[:,-1] = w

            prior_file = open('data/saved_prior_search.yaml', 'w')
            yaml.dump(prior_grid, prior_file)
            prior_file.close()


    sys.exit()

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
        if switches['optimized_priors']:
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

def check_fit_function(switches,check_point):
    print('Performing check of fit function')
    for k in check_point:
        print('%s = %f' %(k,check_point[k]))
    switches['ensembles_fit'] = ['check']

    x_check = dict()
    x_check['check'] = {k:np.inf for k in ['mpiL','mkL']}
    x_check['check']['alphaS'] = 0.

    y_check = dict()
    y_check['check'] = 1.19897

    p_check = dict()
    for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
        p_check[k] = 0.
    for key in ['L1','L2','L3','L4','L5','L6','L7','L8','k_4','p_4']:
        p_check[key] = check_point[key]

    d_e = dict()
    d_e['x'] = x_check
    d_e['y'] = y_check

    for FPK in switches['scales']:
        model = 'xpt_nnlo_'+FPK
        switches['scale'] = FPK
        switches['ansatz']['model'] = model
        if FPK == 'PP':
            Lchi_check = 4*np.pi * check_point['Fpi']
        elif FPK == 'PK':
            Lchi_check = 4*np.pi * np.sqrt(check_point['Fpi'] * check_point['FK'])
        elif FPK == 'KK':
            Lchi_check = 4*np.pi * check_point['FK']
        p_check[('check','Lchi_'+FPK)] = Lchi_check
        d_e['p'] = p_check
        #print('debug p',d_e['p'])

        fit_check = xpt.Fit(switches,xyp_init=d_e)
        fit_check.check_fit_function(check_point)

        # ratio
        model = 'xpt-ratio_nnlo_'+FPK
        switches['ansatz']['model'] = model

        #Lchi_check = 4*np.pi * check_point['Fpi']
        #p_check[('check','Lchi_'+FPK)] = Lchi_check
        #d_e['p'] = p_check

        fit_check = xpt.Fit(switches,xyp_init=d_e)
        fit_check.check_fit_function(check_point)


if __name__ == "__main__":
    import input_params as ip
    print("python     version:", sys.version)
    #print("pandas version:", pd.__version__)
    print("numpy      version:", np.__version__)
    print("pandas     version:", pd.__version__)
    print("matplotlib version:", matplotlib.__version__)
    print("gvar       version:", gv.__version__)
    print("lsqfit     version:", lsqfit.__version__)
    print("yaml       version:", yaml.__version__)
    print('')

    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.ion()

    # Load input params
    switches   = ip.switches
    priors     = ip.priors
    phys_point = ip.phys_point
    check_fit  = ip.check_fit

    # check fit
    if switches['check_fit']:
        check_fit_function(switches,check_fit)
        sys.exit()

    # Load data
    data    = h5.open_file('FK_Fpi_data.h5','r')
    gv_data = format_h5_data(switches,data)
    data.close()

    if switches['do_analysis']:
        fit_results = perform_analysis(switches,priors,phys_point)

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
