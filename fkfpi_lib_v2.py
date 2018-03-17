from __future__ import print_function
import numpy as np
import pandas as pd
import lsqfit
import gvar as gv
import matplotlib.pyplot as plt
import chipt_awl as xpt

ens_long = {
    'a15m310' :'l1648f211b580m013m065m838',
    'a15m220' :'l2448f211b580m0064m0640m828',
    'a15m130' :'l3248f211b580m00235m0647m831',
    'a12m310' :'l2464f211b600m0102m0509m635',
    'a12m220' :'l3264f211b600m00507m0507m628',
    'a12m220S':'l2464f211b600m00507m0507m628',
    'a12m220L':'l4064f211b600m00507m0507m628',
    'a12m130' :'l4864f211b600m00184m0507m628',
    'a09m310' :'l3296f211b630m0074m037m440',
    'a09m220' :'l4896f211b630m00363m0363m430',
}
r_a = {
    'a15m310' :gv.gvar(2.059,0.023),
    'a15m220' :gv.gvar(2.073,0.013),
    'a15m130' :gv.gvar(2.089,0.008),
    'a12m310' :gv.gvar(2.575,0.017),
    'a12m220S':gv.gvar(2.585,0.019),
    'a12m220' :gv.gvar(2.626,0.013),
    'a12m220L':gv.gvar(2.614,0.009),
    'a12m130' :gv.gvar(2.608,0.008),
    'a09m310' :gv.gvar(3.499,0.024),
    'a09m220' :gv.gvar(3.566,0.014)
    }
L_ens = {'a15m310':16,'a15m220':24,'a15m130':32,
    'a12m310':24,'a12m220':32,'a12m220S':24,'a12m220L':40,'a12m130':48,
    'a09m310':32,'a09m220':48,}

def decay_constant(s,mval,bsdata):
    bsdata = {k:np.array(bsdata[k]) for k in bsdata}
    Fpi = bsdata['z0p_pion']*2.*(mval['mq1']+bsdata['mresl'])/bsdata['e0_pion']**(3./2.)
    Fka = bsdata['z0p_kaon']*(mval['mq1']+mval['mq2']+bsdata['mresl']+bsdata['mress'])/bsdata['e0_kaon']**(3./2.)
    Fss = bsdata['z0p_etas']*2.*(mval['mq2']+bsdata['mress'])/bsdata['e0_etas']**(3./2.)
    if s['scale'] == 'PP':
        Lchi = 4 * np.pi * Fpi
    elif s['scale'] == 'PK':
        Lchi = 4 * np.pi * np.sqrt(Fpi*Fka)
    elif s['scale'] == 'KK':
        Lchi = Lchi = 4 * np.pi * Fka
    r = dict()
    r['Lchi'] = np.average(Lchi)
    r['Fka/Fpi'] = Fka/Fpi
    r['mpi'] = bsdata['e0_pion']
    r['mka'] = bsdata['e0_kaon']
    r['mss'] = bsdata['e0_etas']
    r['mju'] = bsdata['phi_ju']
    return r

def format_data(switches,data,mixed_data,hisq_params,priors):
    x    = list()
    y    = list()
    elist= list()
    mpi  = list()
    mpiL = list()
    mka  = list()
    mss  = list()
    aw0  = list()
    a2di = list()
    mju  = list()
    mjuL = list()
    mjs  = list()
    mru  = list()
    mrs  = list()
    Lchi = list()
    for ens in switches['ensemble']:
        e = ens_long[ens]
        elist.append(ens)
        #print(ens)
        # get from postgre csv dump
        # get mres, E0 and Z0p
        data2pt = data.sort_values(by='nbs').query("ensemble=='%s'" %e)[['e0_pion','z0p_pion','e0_kaon','z0p_kaon','e0_etas','z0p_etas','mresl','mress']].to_dict(orient='list')
        if ens in ['a12m220S','a12m220L']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l3264f211b600m00507m0507m628' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        else:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='%s' and tag=='%s'" %(e,t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        bsdata = dict(data2pt,**datamix)
        # get input light quark masses
        mval = data.query("ensemble=='%s'" %e)[['mq1','mq2']].iloc[0].to_dict()
        # get hisq params
        # milc file
        hp = hisq_params.query("ensemble=='%s'" %e)
        r = decay_constant(switches,mval,bsdata)
        if switches['ansatz']['FV']:
            r['mjuL'] = np.array(r['mju']*L_ens[ens])
            r['mpiL'] = np.array(r['mpi']*L_ens[ens])
            a2di = gv.gvar(hp['a2DI_mean'].iloc[0],hp['a2DI_sdev'].iloc[0])/r_a[ens]**2
            r['a2di'] = np.random.normal(a2di.mean,a2di.sdev,len(r['Fka/Fpi']))
            fv_mns_inv = xpt.fv_correction(switches,r)
        #print('DEBUG FV:',fv_mns_inv)
        y = r['Fka/Fpi'] - fv_mns_inv
        print(np.std(fv_mns_inv))
        gvdata = gv.dataset.avg_data({'y_iv':y,'dfv':fv_mns_inv,'y_fv':r['Fka/Fpi']},median=False,spread=True)
        print(gvdata)
        print(gv.evalcorr(gvdata))
        #print(np.average(r['Fka/Fpi']),np.std(r['Fka/Fpi']))
        raise SystemExit
        data_dict = decay_constant(switches,mval,gvdata)
        Lchi.append(data_dict['x']['Lchi'])
        x.append(data_dict['x']['Lchi'])
        y.append(data_dict['y']['Fka/Fpi'])
        mpi.append(data_dict['p']['mpi'])
        mpiL.append(data_dict['p']['mpi'].mean * L_ens[ens])
        mjuL.append(gvdata['phi_ju'].mean * L_ens[ens])
        mka.append(data_dict['p']['mka'])
        mss.append(data_dict['p']['mss'])
        mju.append(gvdata['phi_ju'])
        mjs.append(gvdata['phi_js'])
        mru.append(gvdata['phi_ru'])
        mrs.append(gvdata['phi_rs'])
        # milc file
        hp = hisq_params.query("ensemble=='%s'" %e)
        aw0_ens = gv.gvar(hp['aw0_mean'].iloc[0],hp['aw0_sdev'].iloc[0])
        aw0.append(aw0_ens)
        a2di.append(gv.gvar(hp['a2DI_mean'].iloc[0],hp['a2DI_sdev'].iloc[0])/r_a[ens]**2)
    priors['mpi'] = np.array(mpi)
    priors['mka'] = np.array(mka)
    priors['mss'] = np.array(mss)
    priors['aw0'] = np.array(aw0)
    priors['mju'] = np.array(mju)
    priors['mjs'] = np.array(mjs)
    priors['mru'] = np.array(mru)
    priors['mrs'] = np.array(mrs)
    if switches['ansatz']['a2dm'] == 'avg':
        a2dm_ju = priors['mju']**2 - priors['mpi']**2
        a2dm_sj = priors['mjs']**2 - priors['mka']**2
        a2dm_ru = priors['mru']**2 - priors['mka']**2
        a2dm_rs = priors['mrs']**2 - priors['mss']**2
        priors['a2dm'] = 1./4*(a2dm_ju +a2dm_sj +a2dm_ru +a2dm_rs)
    elif switches['ansatz']['a2dm'] == 'individual':
        priors['a2dm'] = priors['mju']**2-priors['mpi']**2
    y = np.array(y)
    priors['a2di'] = np.array(a2di)
    #priors['a2dm'] = np.array(a2dm)
    mpiL = np.array(mpiL)
    Lchi = np.array(Lchi)
    mjuL = np.array(mjuL)
    if switches['ansatz']['FV']:
        fv_params = dict()
        fv_params['mjuL'] = mjuL
        fv_params['mpiL'] = mpiL
        fv_params['Lchi'] = Lchi
        fv_mns_inv = xpt.fv_correction(switches,fv_params,priors)
        #print('DEBUG FV:',fv_mns_inv)
        y = y - fv_mns_inv
    return {'x':{'Lchi':np.array(x),'mpiL':np.array(mpiL),'elist':np.array(elist)}, 'y': y, 'p': priors}

def fit_data(switches,data,phys_params):
    x = data['x']
    y = data['y']
    Fitc = xpt.Fit(switches)
    prior = data['p']
    fit = lsqfit.nonlinear_fit(data=(x,y),prior=prior,fcn=Fitc.fit_function)
    return fit

def fkfpi_nlo(mpi,mka,meta,Lchi,L5):
    r =  1.
    r += 5./8 * mpi**2 / Lchi**2 * np.log(mpi**2/Lchi**2)
    r -= 1./4 * mka**2 / Lchi**2 * np.log(mka**2/Lchi**2)
    r -= 3./8 * meta**2 / Lchi**2 * np.log(meta**2 / Lchi**2)
    r += 4. * (mka**2 - mpi**2) / Lchi**2 * (4*np.pi)**2 * L5
    return r

def fkfpi_phys(x_phys,fit):
    # use metaSq = 4/3 mK**2 - 1/3 mpi**2
    print('prediction from LQCD')
    meta = np.sqrt(4./3 * x_phys['mka']**2 -1./3 * x_phys['mpi']**2)
    fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
        meta=meta,Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    su2 = dsu2(FKpi=fkp,mpi=x_phys['mpi'],mk=x_phys['mka'],F0=x_phys['F0'])
    print('FK / Fpi = ',fkp)
    print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))
    # use meta = meta_pdg
    print('prediction from LQCD + meta_PDG')
    fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
        meta=x_phys['meta'],Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    print('FK / Fpi = ',fkp)
    print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))

def dsu2(FKpi,mpi,mk,F0):
    R = gv.gvar(35.7,np.sqrt(1.9**2 + 1.8**2))
    d = np.sqrt(3) * np.sqrt(3) / 4 / R * (\
        -4./3 * (FKpi-1) \
        + 2./(6*(4*np.pi)**2*F0**2)*(mk**2-mpi**2-mpi**2*np.log(mk**2/mpi**2))\
        )
    return d
