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
    for e in switches['ensemble']:
        gvdata = gv.dataset.avg_data(data.query("ensemble=='%s'" %e)[['e0_pion','z0p_pion','e0_kaon','z0p_kaon','e0_etas','z0p_etas','mresl','mress']].to_dict(orient='list'),bstrap=True)
        mval = data.query("ensemble=='%s'" %e)[['mq1','mq2']].iloc[0].to_dict()
        data_dict = decay_constant(mval,gvdata)
        x.append(data_dict['x']['L'])
        y.append(data_dict['y']['fka/fpi'])
        mpi.append(data_dict['p']['mpi'])
        mka.append(data_dict['p']['mka'])
        mss.append(data_dict['p']['mss'])
        hp = hisq_params.query("ensemble=='%s'" %e)
        aw0.append(gv.gvar(hp['aw0_mean'].iloc[0],hp['aw0_sdev'].iloc[0]))
    priors['mpi'] = np.array(mpi)
    priors['mka'] = np.array(mka)
    priors['mss'] = np.array(mss)
    priors['aw0'] = np.array(aw0)
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
            r += p['
            return r
        elif self.switches['ansatz']['type'] == 'MA':
            pass

class C4(object):
    # Eq. C4
    def __init__(self,Fka,Fpi,mka,mpi,mss,a2dm,a2di):
        # m_etas^2 = 4/3 * mka^2 - 1/3 mpi^2 # gellmann-okubo
        # mu = 4*pi * sqrt(Fpi Fka)
        self.Fka = Fka
        self.Fpi = Fpi
        self.mka = mka
        self.mpi = mpi
        self.mss = mss
        self.a2dm = a2dm
        self.a2di = a2di
        # construct dependencies
        self.met = np.sqrt(4.*self.mka**2/3. - self.mpi**2/3.) # Gellmann-Okubo relation
        self.mu = 4.*np.pi*np.sqrt(self.Fpi*self.Fka) # scale
        self.fpi = np.sqrt(2*self.Fpi*self.Fka)
        # initialize contributions
        self.ma1 = 0
        self.ma2 = 0
        self.ma3 = 0
        self.ma4 = 0
        self.ma7 = 0
    def __call__(self):
        self.l1()
        self.l2()
        self.l3()
        self.l4()
        self.l7()
        result = self.ma1 + self.ma2 + self.ma3 + self.ma4 + self.ma7
        return result
    def frac(self,num,den):
        return num/den
    def clog(self,num,den):
        return np.log(self.frac(num,den))
    def l1(self):
        # line one of C4
        c1 = self.frac(self.mpi**2+self.a2dm,(4.*np.pi*self.fpi)**2) * self.clog(self.mpi**2+self.a2dm,self.mu**2)
        c2 = self.frac(self.mpi**2,(4.*np.pi*self.fpi)**2) * self.clog(self.mpi**2,self.mu**2)
        c3 = self.frac(self.mka**2+self.a2dm,2.*(4.*np.pi*self.fpi)**2) * self.clog(self.mka**2+self.a2dm,self.mu**2)
        self.ma1 = c1-c2-c3
    def l2(self):
        # line two of C4
        c1 = self.frac(self.mka**2,2.*(4.*np.pi*self.fpi)**2) * self.clog(self.mka**2,self.mu**2)
        c2A = self.frac(3.,4.*(4.*np.pi)**2)
        c2a = self.frac(self.met**2+self.a2di,self.fpi**2) * self.clog(self.met**2+self.a2di,self.mu**2)
        c2b = self.frac(self.met**2,self.fpi**2) * self.clog(self.met**2,self.mu**2)
        c2 = c2A * (c2a-c2b)
        self.ma2 = c1-c2
    def l3(self):
        # line three of C4
        c1A = self.frac(1.,2.*(4.*np.pi)**2)
        c1a = self.frac(self.mss**2+self.a2dm,self.fpi**2) * self.clog(self.mss**2+self.a2dm,self.mu**2)
        c1b = self.frac(self.mss**2,self.fpi**2) * self.clog(self.mss**2,self.mu**2)
        c1 = c1A * (c1a-c1b)
        self.ma3 = -c1
    def l4(self):
        # lines 4, 5, 6 of C4
        c1A = self.frac(self.a2di,self.fpi**2)
        c1B = self.frac(1.,12.*(4.*np.pi)**2)
        c1a = 3.
        c1b = self.frac(4.*(self.mka**2-self.mpi**2),self.mss**2-self.met**2-self.a2di)
        c1c = self.frac(3.*(self.met**2+self.a2di+self.mpi**2),self.met**2+self.a2di-self.mpi**2) * self.clog(self.mpi**2,self.mu**2)
        c1d = self.frac(2.*(3*self.mss**4-(self.met**2+self.a2di)*(3.*self.mss**2-2.*self.mka**2+2.*self.mpi**2)),(self.mss**2-self.met**2-self.a2di)**2) * self.clog(self.mss**2,self.mu**2)
        c1eA = 2.*(self.met**2+self.a2di)
        c1eB = self.clog(self.met**2+self.a2di,self.mu**2)
        c1ea = self.frac(3.*self.mss**2-3.*(self.met**2+self.a2di)+2.*self.mka**2-2.*self.mpi**2,(self.mss**2-self.met**2-self.a2di)**2)
        c1eb = self.frac(3.,self.met**2+self.a2di-self.mpi**2)
        c1 = c1A*c1B * ( c1a +c1b +c1c -c1d + c1eA*c1eB*(c1ea-c1eb) )
        self.ma4 = -c1
    def l7(self):
        c1A = (self.frac(self.a2di,self.fpi**2))**2
        c1B = self.frac(1.,12.*(4.*np.pi)**2)
        c1a = self.frac(self.fpi**2,self.met**2+self.a2di-self.mpi**2)
        c1b = self.frac(2.*self.fpi**2,self.met**2+self.a2di-self.mss**2)
        c1cA = self.clog(self.mpi**2,self.mu**2)
        c1ca = self.frac(self.fpi**2*(self.met**2+self.a2di),(self.met**2+self.a2di-self.mpi**2)**2)
        c1cb = self.frac(4.*self.fpi**2*self.mpi**2,(self.met**2+self.a2di-self.mpi**2)*(self.mss**2-self.mpi**2))
        c1d = self.clog(self.mss**2,self.mu**2) * self.frac(2.*self.fpi**2*(2.*self.mss**4-(self.met**2+self.a2di)*(self.mss**2+self.mpi**2)),(self.mss**2-self.mpi**2)*(self.met**2+self.a2di-self.mss**2)**2)
        c1eA = self.frac(self.met**2+self.a2di,self.fpi**2)
        c1eB = self.clog(self.met**2+self.a2di,self.mu**2)
        c1ea = self.frac(self.fpi**4,(self.met**2+self.a2di-self.mpi**2)**2)
        c1eb = self.frac(2.*self.fpi**4*(self.met**2+self.a2di+self.mpi**2-2.*self.mss**2),(self.met**2+self.a2di-self.mpi**2)*(self.met**2+self.a2di-self.mss**2)**2)
        c1 = c1A*c1B * (c1a - c1b + c1cA*(c1ca-c1cb) - c1d - c1eA*c1eB*(c1ea+c1eb) )
        self.ma7 = c1
