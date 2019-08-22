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
    priors = dict()

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

        #y[ens] = np.array(gvdata['FK']/gvdata['Fpi'])
        y[ens] = gvdata['FK']/gvdata['Fpi']

        x[ens]['mpiL'] = gvdata['mpi'].mean * L_ens[ens]
        x[ens]['mkL']  = gvdata['mk'].mean  * L_ens[ens]
        x[ens]['mssL'] = gvdata['mss'].mean * L_ens[ens]
        x[ens]['mjuL'] = gvdata['mju'].mean * L_ens[ens]
        x[ens]['mjsL'] = gvdata['mjs'].mean * L_ens[ens]
        x[ens]['mruL'] = gvdata['mru'].mean * L_ens[ens]
        x[ens]['mrsL'] = gvdata['mrs'].mean * L_ens[ens]
        x[ens]['Lchi_PK'] = 4*np.pi*np.sqrt(gvdata['FK'].mean * gvdata['Fpi'].mean)
        x[ens]['Lchi_PP'] = 4 *np.pi *gvdata['Fpi'].mean
        x[ens]['Lchi_KK'] = 4 *np.pi *gvdata['FK'].mean

        # MASSES
        priors[(ens,'mpi')] = gvdata['mpi']
        priors[(ens,'mk')]  = gvdata['mk']
        priors[(ens,'mss')] = gvdata['mss']
        priors[(ens,'mju')] = gvdata['mju']
        priors[(ens,'mjs')] = gvdata['mjs']
        priors[(ens,'mru')] = gvdata['mru']
        priors[(ens,'mrs')] = gvdata['mrs']
        # HISQ params
        aw0 = data.get_node('/'+ens+'/aw0').read()
        priors[(ens,'aw0')] = gv.gvar(aw0[0],aw0[1])
        a2di = data.get_node('/'+ens+'/a2DI').read()
        priors[(ens,'a2DI')] = gv.gvar(a2di[0],a2di[1])
        x[ens]['alpha_S'] = data.get_node('/'+ens+'/alpha_s').read()
    return {'x':x, 'y':y, 'p':priors}

    """
    ''' MOVE FV CORRECTION TO FIT CLASS in XPT LIB '''
    if switches['ansatz']['FV']:
        fv_params = dict()
        fv_params['mjuL'] = mjuL
        fv_params['mpiL'] = mpiL
        fv_params['Lchi'] = Lchi
        fv_mns_inv = xpt.fv_correction(switches,fv_params,priors)
        #print('DEBUG FV:',fv_mns_inv)
        y = y - fv_mns_inv
    """

def fit_data(switches,data):
    x = data['x']
    y = data['y']
    p = data['p']
    Fitc = xpt.Fit(switches)
    #tmp = dict()
    #for p in prior:
    #    tmp[p] = prior[p].mean
        #print(p,prior[p].mean)
    #e_tmp = 'a15m310'
    #tmp[(e_tmp,'a2DI')] = 1.
    #tmp['L5'] = 0.001
    #for p in tmp:
    #    print(p,tmp[p])
    #print(x[e_tmp]['Lchi_PK'])
    #print('FK_Fpi_fit_function')
    #print(Fitc.fit_function(x,tmp))
    fit = lsqfit.nonlinear_fit(data=(x,y),prior=p,fcn=Fitc.fit_function)
    return fit

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
    #for e in gv_data['x']:
    #    print("%8s" %e)
    #    print(gv_data['x'][e])
    #print(gv_data['x'])
    #print(gv_data['p'][('a12m220','a2DI')])

    '''
    1 - build prior filter function
    2 - add FV correction to Fit class
    3 - fix counter term for ratio
    4 - make FK and Fpi fit functions
        - xpt
        - MA
    '''

    if switches['nlo_report']:
        fit_results = dict()
        #for model in ['ma_nlo','ma-Kfunc_nlo']:
        for model in ['xpt-taylor_nlo','xpt_nlo']:
            switches['ansatz']['model'] = model
            print(switches['ansatz']['model'])
            model_result = dict()
            for e in switches['ensembles']:
                e_lst = [e]
                x_e = {k:gv_data['x'][k] for k in e_lst}
                y_e = {k:gv_data['y'][k] for k in e_lst}
                p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in e_lst}
                p_e['L4'] = gv.gvar(0,1e-10)#priors['L4']
                p_e['L5'] = priors['L5']
                d_e = dict()
                d_e['x'] = x_e
                d_e['y'] = y_e
                d_e['p'] = p_e
                print('fit_e')
                fit_e = fit_data(switches,d_e)
                print(fit_e.format(maxline=True))
                model_result[e] = fit_e
            fit_results[model] = model_result
        print("%8s & %13s & %13s\\\\" %('ensemble','xpt nlo','ma nlo'))
        print("\\hline")
        for e in fit_results['xpt_nlo']:
            s = "%8s" %e
            for model in ['xpt_nlo','ma_nlo']:
                s += " & %13s" %str(fit_results[model][e].p['L5'])
            s += "\\\\"
            print(s)
