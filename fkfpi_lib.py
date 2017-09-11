import numpy as np
import pandas as pd
import lsqfit
import gvar as gv
import matplotlib.pyplot as plt

def decay_constant(mval,gvdata):
    fpi = gvdata['z0p_pion']*np.sqrt(2.)*(2.*mval['mq1']+2.*gvdata['mresl'])/gvdata['e0_pion']**(3./2.)
    fka = gvdata['z0p_kaon']*np.sqrt(2.)*(mval['mq1']+mval['mq2']+gvdata['mresl']+gvdata['mress'])/gvdata['e0_kaon']**(3./2.)
    fss = gvdata['z0p_etas']*np.sqrt(2.)*(2.*mval['mq2']+2.*gvdata['mress'])/gvdata['e0_etas']**(3./2.)
    L = 4.*np.pi*np.sqrt(fka*fpi)
    x = dict()
    x['L'] = L.mean
    y = dict()
    y['fka/fpi'] = fka/fpi
    p = dict()
    p['mpi'] = gvdata['e0_pion']
    p['mka'] = gvdata['e0_kaon']
    p['mss'] = gvdata['e0_etas']
    return {'x': x, 'y': y, 'p': p}

def format_data(switches,data,hisq_params,priors):
    x = list()
    y = list()
    mpi = list()
    mka = list()
    mss = list()
    aw0 = list()
    a2di = list()
    a2dm = list()
    for e in switches['ensemble']:
        # get from postgre csv dump
        gvdata = gv.dataset.avg_data(data.query("ensemble=='%s'" %e)[['e0_pion','z0p_pion','e0_kaon','z0p_kaon','e0_etas','z0p_etas','mresl','mress']].to_dict(orient='list'),bstrap=True)
        mval = data.query("ensemble=='%s'" %e)[['mq1','mq2']].iloc[0].to_dict()
        data_dict = decay_constant(mval,gvdata)
        x.append(data_dict['x']['L'])
        y.append(data_dict['y']['fka/fpi'])
        mpi.append(data_dict['p']['mpi'])
        mka.append(data_dict['p']['mka'])
        mss.append(data_dict['p']['mss'])
        # milc file
        hp = hisq_params.query("ensemble=='%s'" %e)
        aw0.append(gv.gvar(hp['aw0_mean'].iloc[0],hp['aw0_sdev'].iloc[0]))
        a2di.append(gv.gvar(hp['a2DI_mean'].iloc[0],hp['a2DI_sdev'].iloc[0]))
        a2dm.append(gv.gvar(hp['a2dm_mean'].iloc[0],hp['a2dm_sdev'].iloc[0]))
    priors['mpi'] = np.array(mpi)
    priors['mka'] = np.array(mka)
    priors['mss'] = np.array(mss)
    priors['aw0'] = np.array(aw0)
    if switches['ansatz']['type'] == 'MA':
        priors['a2di'] = np.array(a2di)
        priors['a2dm'] = np.array(a2dm)
    return {'x': np.array(x), 'y': np.array(y), 'p': priors}

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
        self.n = switches['ansatz']['truncation']
    def get_priors(self,p):
        prior = dict()
        prior['L5'] = p['L5']
        prior['mpi'] = p['mpi']
        prior['mka'] = p['mka']
        prior['mss'] = p['mss']
        prior['aw0'] = p['aw0']
        for k in p:
            if int(k[1:]) <= self.n:
                prior[k] = p[k]
            else: pass
        return prior
    def frac(self,num,den):
        return num/den
    def xlog(self,num,den):
        return np.log(self.frac(num,den))
    def fit_function(self,x,p):
        if self.switches['ansatz']['type'] == 'xpt':
            r = 1
            r += 5./8.*self.frac(p['mpi']**2,x**2)*self.xlog(p['mpi']**2,x**2)
            r += -1./4.*self.frac(p['mka']**2,x**2)*self.xlog(p['mka']**2,x**2)
            r += -3./8.*self.frac(p['mss']**2,x**2)*self.xlog(p['mss']**2,x**2)
            r += (4.*np.pi)**2*self.frac(4.*(p['mka']**2-p['mpi']**2),x**2)*(p['L5']+p['aw0']**2*p['s2'])
            return r
        elif self.switches['ansatz']['type'] == 'MA':
            C4cls = C4(x,p)
            logs = C4cls()
            print("logs:",logs)
            r = 1
            r += logs
            r += (4.*np.pi)**2*self.frac(4.*(p['mka']**2-p['mpi']**2),x**2)*p['L5'] #+p['aw0']**2*p['s2'])
            return r 

class C4(object):
    # Eq. C4
    def __init__(self,x,p):
        # m_etas^2 = 4/3 * mka^2 - 1/3 mpi^2 # gellmann-okubo
        # mu = 4*pi * sqrt(Fpi Fka)
        self.mu2 = x**2
        self.eka2 = p['mka']**2/self.mu2
        self.epi2 = p['mpi']**2/self.mu2
        self.ess2 = p['mss']**2/self.mu2
        self.eju2 = (p['mpi']**2+p['a2dm'])/self.mu2
        self.eru2 = (p['mka']**2+p['a2dm'])/self.mu2
        self.esj2 = self.eru2
        self.ers2 = (p['mss']**2+p['a2dm'])/self.mu2
        self.dpq2 = p['a2di']/self.mu2
        # construct dependencies
        self.eX2 = 4.*self.eka2/3.-self.epi2/3.+self.dpq2 # Gellmann-Okubo relation
        # initialize contributions
        self.malog = 0
        self.pilog = 0
        self.sslog = 0
        self.etlog = 0
        self.pqcnt = 0
    def __call__(self):
        self.maLog()
        self.piLog()
        self.ssLog()
        self.etLog()
        self.pqCnt()
        result = self.malog + self.pilog + self.sslog + self.etlog + self.pqcnt
        return result
    def frac(self,num,den):
        return num/den
    def maLog(self):
        # mixed action logs
        cju = 1./2. * self.eju2*np.log(self.eju2)
        csj = 1./2. * self.esj2*np.log(self.esj2)
        cru = 1./4. * self.eru2*np.log(self.eru2)
        crs = 1./4. * self.ers2*np.log(self.ers2)
        self.malog = cju-csj+cru-crs
    def piLog(self):
        # pion logs
        c = 1./8. * np.log(self.epi2) 
        c1 = self.epi2
        c2 = self.frac(self.dpq2*(self.eX2+self.epi2),self.eX2-self.epi2)
        c3 = self.frac(self.dpq2**2*self.eX2,3.*(self.eX2-self.epi2)**2)
        c4 = self.frac(4.*self.dpq2**2*self.epi2,3.*(self.eX2-self.epi2)*(self.ess2-self.epi2))
        self.pilog = c*(c1-c2+c3-c4)
    def ssLog(self):
        # phi_ss log
        c = 1./4. * np.log(self.ess2)
        c1 = self.ess2
        c2 = self.frac(self.dpq2*(3.*self.ess2*2+2.*(self.eka2-self.epi2)*self.eX2-3.*self.ess2*self.eX2),3.*(self.eX2-self.ess2)**2)
        c3 = self.frac(self.dpq2**2*(self.ess2**2-self.eX2*(self.ess2+self.epi2)),3.*(self.eX2-self.ess2)**2*(self.ess2-self.epi2))
        self.sslog = c*(c1+c2-c3)
    def etLog(self):
        # eta logs
        print(self.eX2)
        c = 3./8.*self.eX2*np.log(self.eX2)
        c1 = 1.
        c2 = self.frac(2.*self.dpq2,3.*(self.eX2-self.epi2))
        c3 = self.frac(self.dpq2*(4.*(self.eka2-self.epi2)+6.*(self.ess2-self.eX2)),9.*(self.eX2-self.ess2)**2)
        c4 = self.frac(self.dpq2**2.,9.*(self.eX2-self.ess2)**2)
        c5 = self.frac(2.*self.dpq2**2*(2.*self.ess2-self.epi2-self.eX2),9.*(self.eX2-self.ess2)**2*(self.eX2-self.epi2))
        self.etlog = -c*(c1-c2+c3+c4-c5)
    def pqCnt(self):
        c1 = self.frac(self.dpq2*(self.eka2-self.epi2),6.*(self.eX2-self.ess2))
        c2 = self.frac(self.dpq2**2,24.*self.eX2-self.epi2)
        c3 = self.frac(self.dpq2**2,12.*self.eX2-self.ess2)
        c4 = self.frac(self.dpq2,8.)
        self.pqcnt = c1+c2-c3-c4
