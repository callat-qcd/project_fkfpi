from __future__ import print_function
import numpy as np
import pandas as pd
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
L_ens = {'a15m400':16,'a15m350':16,'a15m310':16,'a15m220':24,'a15m130':32,
    'a12m400':24,'a12m350':24,'a12m310':24,'a12m220':32,'a12m220S':24,'a12m220L':40,'a12m130':48,
    'a09m400':32,'a09m350':32,'a09m310':32,'a09m220':48,}

def decay_constant(s,mval,gvdata):
    Fpi = gvdata['z0p_pion']*2*(mval['mq1']+gvdata['mresl'])/gvdata['e0_pion']**(3./2.)
    Fka = gvdata['z0p_kaon']*(mval['mq1']+mval['mq2']+gvdata['mresl']+gvdata['mress'])/gvdata['e0_kaon']**(3./2.)
    if s['scale'] == 'PP':
        Lchi = 4 * np.pi * Fpi
    elif s['scale'] == 'PK':
        Lchi = 4 * np.pi * np.sqrt(Fpi*Fka)
    elif s['scale'] == 'KK':
        Lchi = Lchi = 4 * np.pi * Fka
    x = dict()
    x['Lchi'] = Lchi.mean
    y = dict()
    y['Fka/Fpi'] = Fka/Fpi
    p = dict()
    p['mpi'] = gvdata['e0_pion']
    p['mka'] = gvdata['e0_kaon']
    return {'x': x, 'y': y, 'p': p}

def format_data(switches,data,mixed_data,hisq_params,priors):
    x    = list()
    y    = list()
    mpi  = list()
    mpiL = list()
    mka  = list()
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
        #print(ens)
        # get from postgre csv dump
        # get mres, E0 and Z0p
        data2pt = data.sort_values(by='nbs').query("ensemble=='%s'" %e)[['e0_pion','z0p_pion','e0_kaon','z0p_kaon','mresl','mress']].to_dict(orient='list')
        if ens in ['a12m220S','a12m220L']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l3264f211b600m00507m0507m628' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        elif ens in ['a15m400','a15m350']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l1648f211b580m013m065m838' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        elif ens in ['a12m400','a12m350']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l2464f211b600m0102m0509m635' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        elif ens in ['a09m400','a09m350']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l3296f211b630m0074m037m440' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        else:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='%s' and tag=='%s'" %(e,t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        datamerge = dict(data2pt,**datamix)
        gvdata = gv.dataset.avg_data(datamerge,bstrap=True)
        mval = data.query("ensemble=='%s'" %e)[['mq1','mq2']].iloc[0].to_dict()
        data_dict = decay_constant(switches,mval,gvdata)
        Lchi.append(data_dict['x']['Lchi'])
        x.append(data_dict['x']['Lchi'])
        y.append(data_dict['y']['Fka/Fpi'])
        mpi.append(data_dict['p']['mpi'])
        mpiL.append(data_dict['p']['mpi'].mean * L_ens[ens])
        mjuL.append(gvdata['phi_ju'].mean * L_ens[ens])
        mka.append(data_dict['p']['mka'])
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
    mpiL = np.array(mpiL)
    Lchi = np.array(Lchi)
    priors['a2di'] = np.array(a2di)
    #priors['a2dm'] = np.array(a2dm)
    mjuL = np.array(mjuL)
    if switches['ansatz']['FV']:
        fv_params = dict()
        fv_params['mjuL'] = mjuL
        fv_params['mpiL'] = mpiL
        fv_params['Lchi'] = Lchi
        fv_mns_inv = xpt.fv_correction(switches,fv_params,priors)
        #print('DEBUG FV:',fv_mns_inv)
        y = y - fv_mns_inv
    return {'x':{'Lchi':np.array(x),'mpiL':np.array(mpiL)}, 'y': y, 'p': priors}

def fit_data(switches,data,phys_params):
    x = data['x']
    y = data['y']
    Fitc = xpt.Fit(switches)
    prior = data['p']
    fit = lsqfit.nonlinear_fit(data=(x,y),prior=prior,fcn=Fitc.fit_function)
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
    d = np.sqrt(3) * np.sqrt(3) / 4 / R * (\
        -4./3 * (FKpi-1) \
        + 2./(6*(4*np.pi)**2*F0**2)*(mk**2-mpi**2-mpi**2*np.log(mk**2/mpi**2))\
        )
    return d
