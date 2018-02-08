from __future__ import print_function
import numpy as np
import scipy.special as spsp
import pandas as pd
import lsqfit
import gvar as gv
import matplotlib.pyplot as plt

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

def decay_constant(s,mval,gvdata):
    Fpi = gvdata['z0p_pion']*2*(mval['mq1']+gvdata['mresl'])/gvdata['e0_pion']**(3./2.)
    Fka = gvdata['z0p_kaon']*(mval['mq1']+mval['mq2']+gvdata['mresl']+gvdata['mress'])/gvdata['e0_kaon']**(3./2.)
    Fss = gvdata['z0p_etas']*2*(mval['mq2']+gvdata['mress'])/gvdata['e0_etas']**(3./2.)
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
    p['mss'] = gvdata['e0_etas']
    return {'x': x, 'y': y, 'p': p}

def format_data(switches,data,mixed_data,hisq_params,priors):
    x    = list()
    y    = list()
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
        #print(ens)
        # get from postgre csv dump
        # get mres, E0 and Z0p
        data2pt = data.sort_values(by='nbs').query("ensemble=='%s'" %e)[['e0_pion','z0p_pion','e0_kaon','z0p_kaon','e0_etas','z0p_etas','mresl','mress']].to_dict(orient='list')
        if ens in ['a12m220S','a12m220L']:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='l3264f211b600m00507m0507m628' and tag=='%s'" %(t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        else:
            datamix = {t:np.squeeze(mixed_data.sort_values(by='nbs').query("ensemble=='%s' and tag=='%s'" %(e,t))[['E0']].as_matrix()) for t in ['phi_ju','phi_js','phi_ru','phi_rs']}
        datamerge = dict(data2pt,**datamix)
        gvdata = gv.dataset.avg_data(datamerge,bstrap=True)
        #print(gvdata)
        mval = data.query("ensemble=='%s'" %e)[['mq1','mq2']].iloc[0].to_dict()
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
        #a2dm.append(gv.gvar(hp['a2dm_mean'].iloc[0],hp['a2dm_sdev'].iloc[0])*aw0_ens**2)
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
    mpiL = np.array(mpiL)
    Lchi = np.array(Lchi)
    priors['a2di'] = np.array(a2di)
    #priors['a2dm'] = np.array(a2dm)
    mjuL = np.array(mjuL)
    if switches['ansatz']['FV']:
        cn = np.array([6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24])
        n_mag = np.sqrt(np.arange(1,len(cn)+1,1))
        k1ju = np.zeros_like(mpiL)
        k1pi = np.zeros_like(mpiL)
        k0pi = np.zeros_like(mpiL)
        k2pi = np.zeros_like(mpiL)
        for i,mL in enumerate(mpiL):
            k1ju[i] = np.sum(cn / (n_mag * mjuL[i]) * spsp.kn(1,n_mag * mjuL[i]))
            k1pi[i] = np.sum(cn / (n_mag * mL) * spsp.kn(1,n_mag * mL))
            k0pi[i] = np.sum(cn * spsp.kn(0,n_mag * mL))
            k2pi[i] = np.sum(cn * spsp.kn(2,n_mag * mL))
        p2 = (priors['mpi']/ Lchi)**2
        if switches['ansatz']['type'] == 'xpt':
            fv_mns_inv = 5./2 * p2 * k1pi
        elif switches['ansatz']['type'] == 'MA':
            ju  = priors['mju']**2 / Lchi**2
            k2  = priors['mka']**2 / Lchi**2
            x2  = 4./3 * k2 - p2/3 + priors['a2di'] / Lchi**2
            s2  = priors['mss']**2 / Lchi**2
            dju = priors['a2di'] / Lchi**2
            drs = priors['a2di'] / Lchi**2
            fv_mns_inv  = 2 * ju * k1ju + 0.5 * p2 * k1pi
            fv_mns_inv -= dju/8 * (2*k1pi - k0pi - k2pi)
            fv_mns_inv -= dju * p2 / (x2-p2) * k1pi
            fv_mns_inv += dju**2 / 24 * (4*p2*k1pi/(x2-p2)**2 \
                +(2*k1pi - k0pi -k2pi)/(x2-p2))
            fv_mns_inv += 2./3 * dju * drs * p2 * k1pi/(x2-p2)/(s2-p2)
        #print('DEBUG FV:',fv_mns_inv)
        y = y - fv_mns_inv
    return {'x':{'Lchi':np.array(x),'mpiL':np.array(mpiL)}, 'y': y, 'p': priors}

def fit_data(switches,data,phys_params):
    x = data['x']
    y = data['y']
    Fitc = Fit(switches)
    prior = data['p']
    fit = lsqfit.nonlinear_fit(data=(x,y),prior=prior,fcn=Fitc.fit_function)
    return fit

class Fit(object):
    def __init__(self,switches):
        self.switches = switches
        #self.n = switches['ansatz']['truncation']
    def fit_function(self,x,p):
        if self.switches['ansatz']['type'] == 'xpt':
            metaSq = 4./3 * p['mka']**2 - 1./3 * p['mpi']**2
            r =  1.
            r += 5./8 * p['mpi']**2 / x['Lchi']**2 * np.log(p['mpi']**2/x['Lchi']**2)
            r -= 1./4 * p['mka']**2 / x['Lchi']**2 * np.log(p['mka']**2/x['Lchi']**2)
            r -= 3./8 * metaSq / x['Lchi']**2 * np.log(metaSq / x['Lchi']**2)
            r += 4. * (p['mka']**2 - p['mpi']**2) / x['Lchi']**2 * (4*np.pi)**2 * p['L5']
            r += p['aw0']**2 * (p['mka']**2 - p['mpi']**2) / x['Lchi']**2 * p['s2'] * (4*np.pi)
            return r
        elif self.switches['ansatz']['type'] == 'MA':
            p2  = p['mpi']**2 / x['Lchi']**2
            k2  = p['mka']**2  / x['Lchi']**2
            s2  = p['mss']**2 / x['Lchi']**2
            if self.switches['ansatz']['a2dm'] == 'avg':
                ju = p2 + p['a2dm'] / x['Lchi']**2
                sj = k2 + p['a2dm'] / x['Lchi']**2
                ru = k2 + p['a2dm'] / x['Lchi']**2
                rs = s2 + p['a2dm'] / x['Lchi']**2
            elif self.switches['ansatz']['a2dm'] == 'individual':
                ju = (p['mju'] / x['Lchi'])**2
                sj = (p['mjs'] / x['Lchi'])**2
                ru = (p['mru'] / x['Lchi'])**2
                rs = (p['mrs'] / x['Lchi'])**2
            #print('ju:',ju)
            x2  = 4./3 * k2 - p2/3 + p['a2di'] / x['Lchi']**2
            dju = p['a2di'] / x['Lchi']**2
            drs = p['a2di'] / x['Lchi']**2
            r  = 1.
            # (-) pion log terms
            r += ju * np.log(ju)
            r += 0.5 * ru * np.log(ru)
            # (+) kaon log terms
            r += -0.5 * ju * np.log(ju)
            r += -1./4 * ru * np.log(ru)
            r += -0.5*sj * np.log(sj)
            r += -1./4 * rs * np.log(rs)
            r += -dju / 8
            r +=  dju**2 / 24 / (x2 - p2)
            r +=  drs * (k2-p2) / 6 / (x2-s2)
            r += -dju * drs / 12 / (x2 - s2)
            r += np.log(p2)/24 * (3*p2 \
                - 3*dju*(x2+p2)/(x2-p2) \
                + dju**2 * x2/(x2-p2)**2\
                -4*dju*drs*p2/(x2-p2)/(s2-p2)\
                )
            r += -x2/24 * np.log(x2)*(9 \
                -6*dju/(x2-p2) \
                + dju**2/(x2-p2)**2\
                +drs*(4*(k2-p2)+6*(s2-x2))/(x2-s2)**2 \
                -2*dju*drs*(2*s2-p2-x2)/(x2-s2)**2/(x2-p2)\
                )
            r += np.log(s2)/12 * (3*s2 \
                +drs*(3*s2**2 + 2*(k2-p2)*x2 -3*s2*x2)/(x2-s2)**2\
                -dju*drs*(2*s2**2 - x2*(s2+p2))/(x2-s2)**2 / (s2-p2)\
                )
            # counter term
            r += 4 * (k2 - p2) * (4*np.pi)**2 * p['L5']
            #r += p['aw0']**2 * (p['mka']**2 - p['mpi']**2) / x['Lchi']**2 * p['s2'] * (4*np.pi)
            return r

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
