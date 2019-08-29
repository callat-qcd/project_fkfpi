from __future__ import print_function
import sys
import numpy as np
import tables as h5
import lsqfit
import gvar as gv
import matplotlib.pyplot as plt
import chipt_awl as xpt

ens_long = {
    'a15m400' :'l1648f211b580m0217m065m838',
    'a15m350' :'l1648f211b580m0166m065m838',
    'a15m310' :'l1648f211b580m013m065m838',
    'a15m220' :'l2448f211b580m0064m0640m828',
    'a15m130' :'l3248f211b580m00235m0647m831',
    'a12m400' :'l2464f211b600m0170m0509m635',
    'a12m350' :'l2464f211b600m0130m0509m635',
    'a12m310' :'l2464f211b600m0102m0509m635',
    'a12m220' :'l3264f211b600m00507m0507m628',
    'a12m220S':'l2464f211b600m00507m0507m628',
    'a12m220L':'l4064f211b600m00507m0507m628',
    'a12m130' :'l4864f211b600m00184m0507m628',
    'a09m400' :'l3264f211b630m0124m037m440',
    'a09m350' :'l3264f211b630m00945m037m440',
    'a09m310' :'l3296f211b630m0074m037m440',
    'a09m220' :'l4896f211b630m00363m0363m430',
}
r_a = {
    'a15m400' :gv.gvar(2.059,0.023),
    'a15m350' :gv.gvar(2.059,0.023),
    'a15m310' :gv.gvar(2.059,0.023),
    'a15m220' :gv.gvar(2.073,0.013),
    'a15m130' :gv.gvar(2.089,0.008),
    'a12m400' :gv.gvar(2.575,0.017),
    'a12m350' :gv.gvar(2.575,0.017),
    'a12m310' :gv.gvar(2.575,0.017),
    'a12m220S':gv.gvar(2.585,0.019),
    'a12m220' :gv.gvar(2.626,0.013),
    'a12m220L':gv.gvar(2.614,0.009),
    'a12m130' :gv.gvar(2.608,0.008),
    'a09m400' :gv.gvar(3.499,0.024),
    'a09m350' :gv.gvar(3.499,0.024),
    'a09m310' :gv.gvar(3.499,0.024),
    'a09m220' :gv.gvar(3.566,0.014)
    }
L_ens = {
    'a15m400':16,'a15m350':16,'a15m310':16,'a15m220':24,'a15m130':32,
    'a12m400':24,'a12m350':24,'a12m310':24,'a12m220':32,'a12m220S':24,'a12m220L':40,'a12m130':48,
    'a09m400':32,'a09m350':32,'a09m310':32,'a09m220':48,}

def format_h5_data(switches,data):
    x = dict()
    y = dict()
    p = dict()

    for ens in switches['ensembles']:
        #print(ens)
        x[ens] = dict()
        data_dict = dict()
        for m in ['mpi','mk','mss','mju','mjs','mrs','mru']:
            data_dict[m] = data.get_node('/'+ens+'/'+m).read()
            #print('  %s = %.5f +- %.5f' %(m,data_dict[m].mean(),data_dict[m].std()))
        for f in ['FK','Fpi','Fss']:
            data_dict[f] = data.get_node('/'+ens+'/'+f).read()
        if switches['debug']:
            fkfpi = data_dict['FK']/data_dict['Fpi']
            print('  FK/Fpi = %.5f +- %.5f' %(fkfpi.mean(),fkfpi.std()))
        gvdata = gv.dataset.avg_data(data_dict,bstrap=True)
        if switches['debug']:
            print(gvdata)
        if switches['debug']:
            print('data_dict')
            print(data_dict)

        y[ens] = gvdata['FK']/gvdata['Fpi']

        x[ens]['mpiL'] = gvdata['mpi'].mean * L_ens[ens]
        x[ens]['mkL']  = gvdata['mk'].mean  * L_ens[ens]
        x[ens]['mssL'] = gvdata['mss'].mean * L_ens[ens]
        x[ens]['mjuL'] = gvdata['mju'].mean * L_ens[ens]
        x[ens]['mjsL'] = gvdata['mjs'].mean * L_ens[ens]
        x[ens]['mruL'] = gvdata['mru'].mean * L_ens[ens]
        x[ens]['mrsL'] = gvdata['mrs'].mean * L_ens[ens]

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


if __name__ == "__main__":
    import input_params as ip
    print("python version:", sys.version)
    #print("pandas version:", pd.__version__)
    print("numpy  version:", np.__version__)
    print("gvar   version:", gv.__version__)
    print("lsqfit version:", lsqfit.__version__)

    # Load input params
    switches   = ip.switches
    priors     = ip.priors
    phys_p     = ip.phys_p
    flag_FKFpi = ip.flag_FKFpi

    # Load data
    data  = h5.open_file('FK_Fpi_data.h5','r')
    gv_data = format_h5_data(switches,data)
    data.close()

    # do analysis
    models = [
        'xpt_nnnlo_FV','xpt_nnnlo_FV_alphaS',
        'xpt_nnlo_FV','xpt_nnlo_FV_alphaS',
        'ma_nnlo_FV_alphaS',
        'xpt-ratio_nnlo_FV_alphaS']
    models = ['xpt_nnnlo_FV','xpt_nnlo_FV_alphaS']
    fit_results = dict()
    for model in models:
        switches['ansatz']['model'] = model
        print('EFT: ',model)
        x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
        y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
        p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
        for key in priors:
            p_e[key] = priors[key]
        d_e = dict()
        d_e['x'] = x_e
        d_e['y'] = y_e
        d_e['p'] = p_e
        fit_e = xpt.Fit(switches,xyp_init=d_e)
        fit_e.fit_data()

        x_phys = dict()
        x_phys['phys'] = {k:np.inf for k in ['mpiL','mkL']}
        x_phys['phys']['alphaS'] = 0.
        y_phys = dict()
        y_phys['phys'] = phys_p['FK'] / phys_p['Fpi']
        p_phys = dict()
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p_phys[k] = 0.
        Lchi_phys = phys_p['Lchi']
        p_phys[('phys','p2')] = phys_p['mpi']**2 / Lchi_phys**2
        p_phys[('phys','k2')] = phys_p['mk']**2  / Lchi_phys**2
        p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
        p_phys[('phys','a2')] = 0.
        for k in fit_e.fit.p:
            if isinstance(k,str):
                print(k,fit_e.fit.p[k])
        for k in ['L5','L4','k_4','p_4','kp_6','k_6','p_6']:
            if k in fit_e.fit.p:
                p_phys[k] = fit_e.fit.p[k]
        fit_e.fv = False
        if 'ratio' in model:
            eft = 'xpt-ratio'
        else:
            eft = 'xpt'
        order = model.split('_')[1]
        fit_e.eft   = eft
        fit_e.order = order
        print('chi2/dof [dof] = %.2f [%d]    Q = %.2e    logGBF = %.3f' \
            %(fit_e.fit.chi2/fit_e.fit.dof,fit_e.fit.dof,fit_e.fit.Q,fit_e.fit.logGBF))
        print(fit_e.fit_function(x_phys,p_phys),'\n')


    if switches['nlo_fv_report']:
        models = ['ma_nlo','ma_nlo_FV','xpt_nlo','xpt_nlo_FV']
        mod = {
            'ma_nlo':'ma_nlo',  'ma_nlo_FV':'ma_nlo',
            'xpt_nlo':'xpt_nlo','xpt_nlo_FV':'xpt_nlo'}
        latex = {
            'ma_nlo':'ma nlo',   'ma_nlo_FV':'ma nlo w/FV',
            'xpt_nlo':'xpt nlo', 'xpt_nlo_FV':'xpt nlo w/FV'}
        marker = {
            'ma_nlo':'s',  'ma_nlo_FV':'o',
            'xpt_nlo':'*', 'xpt_nlo_FV':'d'}
        color = {
            'ma_nlo':'r',  'ma_nlo_FV':'b',
            'xpt_nlo':'m', 'xpt_nlo_FV':'g'}
        fit_results = dict()
        for model in models:
            switches['ansatz']['model'] = mod[model]
            print('EFT: ',mod[model])
            model_result = dict()
            for e in switches['ensembles']:
                switches['ensembles_fit'] = [e]
                if 'FV' in model:
                    switches['ansatz']['FV'] = True
                else:
                    switches['ansatz']['FV'] = False
                x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
                y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
                p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
                p_e['L4'] = gv.gvar(0,5.e-3)#priors['L4']
                p_e['L5'] = priors['L5']
                d_e = dict()
                d_e['x'] = x_e
                d_e['y'] = y_e
                d_e['p'] = p_e
                fit_e = fit_data(switches,d_e)
                model_result[e] = fit_e
            fit_results[model] = model_result
        plt.ion()
        fig = plt.figure('nlo_report_FV')
        ax  = plt.axes([.12, .12, .85, .85])

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
                s += " & %13s" %str(fit_results[model][e].p['L5'])
            s += "\\\\"
            print(s)

            for i_m,model in enumerate(models):
                if i_e == 0:
                    label=latex[models[i_m]]
                else:
                    label=''
                y = len(switches['ensembles']) - i_e + 0.1*i_m
                ax.errorbar(x=fit_results[model][e].p['L5'].mean,y=y,
                    xerr=fit_results[model][e].p['L5'].sdev,linestyle='None',
                    marker=marker[model],mfc='None',color=color[model],label=label)
        plt.yticks(np.arange(len(switches['ensembles']),0,-1),tuple(switches['ensembles']))
        ax.set_xlabel(r'$L_5$',fontsize=16)
        ax.legend(loc=1,fontsize=16)
        ax.set_xlim(-0.0001,0.001)
        plt.savefig('nlo_report_FV.pdf',transparent=True)
        plt.ioff()
        plt.show()


    if switches['nlo_report']:
        models = ['ma_nlo','ma-ratio_nlo','xpt_nlo','xpt-ratio_nlo']
        models = ['ma_nlo','xpt_nlo']
        latex = {'ma_nlo':'ma nlo', 'ma-ratio_nlo':'ma-r nlo',
            'xpt_nlo':'xpt nlo', 'xpt-ratio_nlo':'xpt-r nlo'}
        marker = {'ma_nlo':'s', 'ma-ratio_nlo':'o',
            'xpt_nlo':'d', 'xpt-ratio_nlo':'*'}
        color = {'ma_nlo':'r', 'ma-ratio_nlo':'g',
            'xpt_nlo':'b', 'xpt-ratio_nlo':'magenta'}
        fit_results = dict()
        #for model in ['ma_nlo','ma-Kfunc_nlo']:
        #for model in ['xpt_nlo','ma_nlo']:
        for model in models:
            switches['ansatz']['model'] = model
            print('EFT: ',model)
            model_result = dict()
            for e in switches['ensembles']:
                switches['ensembles_fit'] = [e]
                x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
                y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
                p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
                p_e['L4'] = gv.gvar(0,5.e-3)#priors['L4']
                p_e['L5'] = priors['L5']
                d_e = dict()
                d_e['x'] = x_e
                d_e['y'] = y_e
                d_e['p'] = p_e
                #print('fit_e')
                #print(d_e['x'])
                fit_e = fit_data(switches,d_e)
                #print(fit_e.format(maxline=True))
                model_result[e] = fit_e
            fit_results[model] = model_result

        plt.ion()
        fig = plt.figure('nlo_report')
        ax  = plt.axes([.12, .12, .85, .85])

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
                s += " & %13s" %str(fit_results[model][e].p['L5'])
            s += "\\\\"
            print(s)

            for i_m,model in enumerate(models):
                if i_e == 0:
                    label=latex[models[i_m]]
                else:
                    label=''
                y = len(fit_results['xpt_nlo']) - i_e + 0.1*i_m
                ax.errorbar(x=fit_results[model][e].p['L5'].mean,y=y,
                    xerr=fit_results[model][e].p['L5'].sdev,linestyle='None',
                    marker=marker[model],mfc='None',color=color[model],label=label)
        plt.yticks(np.arange(len(switches['ensembles']),0,-1),tuple(switches['ensembles']))
        ax.set_xlabel(r'$L_5$',fontsize=16)
        ax.legend(loc=1,fontsize=16)
        if len(models) == 2:
            ax.set_xlim(-0.0001,0.0008)
        else:
            ax.set_xlim(-0.0015,0.004)
        plt.savefig('nlo_report.pdf',transparent=True)
        plt.ioff()
        plt.show()
