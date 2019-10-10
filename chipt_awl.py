import sys
import scipy.special as spsp
import numpy as np
import matplotlib.pyplot as plt
import gvar as gv
import lsqfit


class Fit(object):
    def __init__(self,switches,xyp_init):
        self.switches = switches
        self.x        = xyp_init['x']
        self.y        = xyp_init['y']
        self.p_init   = xyp_init['p']
        self.eft      = self.switches['ansatz']['model'].split('_')[0]
        self.order    = self.switches['ansatz']['model'].split('_')[1]
        self.fv       = 'FV' in self.switches['ansatz']['model']
        self.alphaS   = 'alphaS' in self.switches['ansatz']['model']
        self.logSq    = 'logSq' in self.switches['ansatz']['model']
        self.lec_nnlo = ['L5','L4','k_4','p_4','s_4','saS_4']
        self.lec_nnnlo= ['s_6','sk_6','sp_6','kp_6','k_6','p_6']

        if self.eft   == 'xpt':
            self.fit_function = self.xpt_nlo
        elif self.eft == 'xpt-ratio':
            self.fit_function = self.xpt_ratio_nlo
        elif self.eft == 'ma':
            self.fit_function = self.ma_nlo
        elif self.eft == 'ma-ratio':
            self.fit_function = self.ma_ratio_nlo
        elif self.eft == 'ma-longform':
            self.fit_function = self.ma_longform_nlo

        if self.fv:
            self.fv_I  = dict()
            self.fv_dI = dict()
            def k1k2k3(mL):
                cn = np.array([6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24])
                n_mag = np.sqrt(np.arange(1,len(cn)+1,1))
                k1 = np.sum(cn * spsp.kn(1,n_mag * mL) / mL / n_mag)
                k2 = np.sum(cn * spsp.kn(2,n_mag * mL) )
                k0 = np.sum(cn * spsp.kn(0,n_mag * mL) )
                return k0,k1,k2
            for e in self.switches['ensembles']:
                # mpi
                Lchi     = self.p_init[(e,'Lchi_'+self.switches['scale'])]
                esq      = self.p_init[(e,'mpi')]**2 / Lchi**2
                mL       = self.x[e]['mpiL']
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'p2')]  = 4 * k1
                self.fv_dI[(e,'p2')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                #print(e,'p2',self.fv_I[(e,'p2')],self.fv_dI[(e,'p2')])
                # mK
                esq      = self.p_init[(e,'mk')]**2 / Lchi**2
                mL       = self.x[e]['mkL']
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'k2')]  = 4 * k1
                self.fv_dI[(e,'k2')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                # meta
                esq      = 4./3 * self.p_init[(e,'mk')]**2 / Lchi**2
                esq     += -1./3 *self.p_init[(e,'mpi')]**2 / Lchi**2
                mL       = np.sqrt(4./3*self.x[e]['mkL']**2 - 1./3 * self.x[e]['mpiL']**2)
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'e2')]  = 4 * k1
                self.fv_dI[(e,'e2')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                if 'ma' in self.eft:
                    # mss
                    esq = self.p_init[(e,'mss')]**2 / Lchi**2
                    mL  = self.x[e]['mssL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'s2')]  = 4 * k1
                    self.fv_dI[(e,'s2')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                    # mju
                    esq = self.p_init[(e,'mju')]**2 / Lchi**2
                    mL  = self.x[e]['mjuL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'ju')]  = 4 * k1
                    self.fv_dI[(e,'ju')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                    # mjs
                    esq = self.p_init[(e,'mjs')]**2 / Lchi**2
                    mL  = self.x[e]['mjsL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'js')]  = 4 * k1
                    self.fv_dI[(e,'js')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                    # mru
                    esq = self.p_init[(e,'mru')]**2 / Lchi**2
                    mL  = self.x[e]['mruL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'ru')]  = 4 * k1
                    self.fv_dI[(e,'ru')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                    # mrs
                    esq = self.p_init[(e,'mrs')]**2 / Lchi**2
                    mL  = self.x[e]['mrsL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'rs')]  = 4 * k1
                    self.fv_dI[(e,'rs')] = (2*k1 - k0 - k2) / (4*np.pi)**2
                    # mX
                    e2  = 4./3 * self.p_init[(e,'mk')]**2 / Lchi**2
                    e2 += -1./3 *self.p_init[(e,'mpi')]**2 / Lchi**2
                    esq = e2 + self.p_init[(e,'a2DI')] / Lchi**2
                    mL  = np.sqrt(4./3*self.x[e]['mkL']**2 - 1./3 * self.x[e]['mpiL']**2)
                    mL  = np.sqrt(esq / e2) * mL
                    self.fv_I[(e,'x2')]  = 4 * k1
                    self.fv_dI[(e,'x2')] = (2*k1 - k0 - k2) / (4*np.pi)**2

    def prune_x(self):
        x = gv.BufferDict()
        for e in self.switches['ensembles_fit']:
            x[e] = self.x[e]
        return x

    def prune_data(self):
        y = gv.BufferDict()
        for e in self.switches['ensembles_fit']:
            y[e] = self.y[e]
        return y

    def prune_priors(self):
        p = gv.BufferDict()
        for e in self.switches['ensembles_fit']:
            Lchi = self.p_init[(e,'Lchi_'+self.switches['scale'])]
            #Lchi = self.x[e]['Lchi_'+self.switches['scale']]
            p[(e,'p2')] = self.p_init[(e,'mpi')]**2 / Lchi**2
            p[(e,'k2')] = self.p_init[(e,'mk')]**2 /  Lchi**2
            p[(e,'e2')] = 4./3 * p[(e,'k2')] - 1./3 * p[(e,'p2')]
            if 'ma' in self.eft:
                p[(e,'s2')] = self.p_init[(e,'mss')]**2 / Lchi**2
                ''' Add ability to use average mixed meson splitting? '''
                p[(e,'ju')] = self.p_init[(e,'mju')]**2 / Lchi**2
                p[(e,'js')] = self.p_init[(e,'mjs')]**2 / Lchi**2
                p[(e,'ru')] = self.p_init[(e,'mru')]**2 / Lchi**2
                p[(e,'rs')] = self.p_init[(e,'mrs')]**2 / Lchi**2
                p[(e,'x2')] = p[(e,'e2')] + self.p_init[(e,'a2DI')] / Lchi**2
                p[(e,'dju2')] = self.p_init[(e,'a2DI')] / Lchi**2
                p[(e,'drs2')] = self.p_init[(e,'a2DI')] / Lchi**2
            if self.order in ['nnlo','nnnlo']:
                p[(e,'a2')]   = self.p_init[(e,'aw0')]**2 / (4 * np.pi)
        ''' Add LECs '''
        if self.order in ['nlo','nnlo','nnnlo']:
            p['L5'] = self.p_init['L5']
            if 'ratio' in self.eft:
                p['L4'] = self.p_init['L4']
        if self.order in ['nnlo','nnnlo']:
            p['s_4'] = self.p_init['s_4']
            p['k_4'] = self.p_init['k_4']
            p['p_4'] = self.p_init['p_4']
            if self.alphaS:
                p['saS_4'] = self.p_init['saS_4']
        if self.order in ['nnnlo']:
            p['kp_6'] = self.p_init['kp_6']
            p['k_6']  = self.p_init['k_6']
            p['p_6']  = self.p_init['p_6']
            p['s_6']  = self.p_init['s_6']
            p['sk_6'] = self.p_init['sk_6']
            p['sp_6'] = self.p_init['sp_6']

        return p

    def make_x_lec(self,x,p,e):
        x_par = gv.BufferDict()
        x_par['p2'] = {'esq':p[(e,'p2')]}
        x_par['k2'] = {'esq':p[(e,'k2')]}
        x_par['e2'] = {'esq':p[(e,'e2')]}
        if self.fv:
            x_par['p2'].update({'fvI':self.fv_I[(e,'p2')], 'fvdI':self.fv_dI[(e,'p2')]})
            x_par['k2'].update({'fvI':self.fv_I[(e,'k2')], 'fvdI':self.fv_dI[(e,'k2')]})
            x_par['e2'].update({'fvI':self.fv_I[(e,'e2')], 'fvdI':self.fv_dI[(e,'e2')]})
        if 'ma' in self.eft:
            x_par['s2'] = {'esq':p[(e,'s2')]}
            x_par['ju'] = {'esq':p[(e,'ju')]}
            x_par['js'] = {'esq':p[(e,'js')]}
            x_par['ru'] = {'esq':p[(e,'ru')]}
            x_par['rs'] = {'esq':p[(e,'rs')]}
            x_par['x2'] = {'esq':p[(e,'x2')]}
            if self.fv:
                x_par['s2'].update({'fvI':self.fv_I[(e,'s2')], 'fvdI':self.fv_dI[(e,'s2')]})
                x_par['ju'].update({'fvI':self.fv_I[(e,'ju')], 'fvdI':self.fv_dI[(e,'ju')]})
                x_par['js'].update({'fvI':self.fv_I[(e,'js')], 'fvdI':self.fv_dI[(e,'js')]})
                x_par['ru'].update({'fvI':self.fv_I[(e,'ru')], 'fvdI':self.fv_dI[(e,'ru')]})
                x_par['rs'].update({'fvI':self.fv_I[(e,'rs')], 'fvdI':self.fv_dI[(e,'rs')]})
                x_par['x2'].update({'fvI':self.fv_I[(e,'x2')], 'fvdI':self.fv_dI[(e,'x2')]})

            x_par['dju2'] = p[(e,'dju2')]
            x_par['drs2'] = p[(e,'drs2')]
        if self.order in ['nnlo','nnnlo']:
            x_par['a2']   = p[(e,'a2')]
        if self.alphaS:
            x_par['alphaS'] = x[e]['alphaS']

        lec = {key:val for key,val in p.items() if isinstance(key,str)}

        return x_par,lec

    def counterterms(self,x,lec):
        ct = 0.
        # NLO terms in FK and Fpi functions
        '''
        ct = 4 * (k2 - p2) * (4*np.pi)**2 * p['L5']
        '''
        k2 = x['k2']['esq']
        p2 = x['p2']['esq']
        if self.order in ['nnlo','nnnlo']:
            a2 = x['a2']
        if self.alphaS:
            alphaS = x['alphaS']
        if self.order in ['nnlo','nnnlo']:
            ct += (k2 - p2)    * a2 * lec['s_4'] # a^2 m^2
            ct += (k2 - p2)    * k2 * lec['k_4'] # m^4
            ct += (k2 - p2)    * p2 * lec['p_4'] # m^4
            '''
            given the SU(3) flavor constraint, the following ct is redundant
                (k2 - p2)**2
            '''
            if self.alphaS:
                ct += (k2 - p2) * a2 * alphaS * lec['saS_4']
            if self.logSq:
                ct += 1./96 * (k2 - p2)* (17*k2 + 37*p2)*(np.log(np.sqrt(k2*p2)))**2
                #ct += 1./96 * (k2 - p2)* (17*k2 + 37*p2)*(np.log(p2))**2
        if self.order in ['nnnlo']:
            ct += (k2 - p2)**2 * k2 * p2 * lec['kp_6']
            ct += (k2 - p2)**2 * k2      * lec['k_6']
            ct += (k2 - p2)**2 * p2      * lec['p_6']
            ct += (k2 - p2)    * a2**2   * lec['s_6']
            ct += (k2 - p2)    * k2 * a2 * lec['sk_6']
            ct += (k2 - p2)    * p2 * a2 * lec['sp_6']
        return ct

    def I(self,x):
        r = x['esq'] * np.log(x['esq'])
        if self.fv:
            r += x['esq'] * x['fvI']
        return r
    def dI(self,x):
        r = 1 + np.log(x['esq'])
        if self.fv:
            r += x['fvdI']
        return r
    def K(self,x1,x2):
        r  = self.I(x2)
        r += -self.I(x1)
        r  = r / (x2['esq'] - x1['esq'])
        return r
    def K21(self,x1,x2):
        r  = self.K(x1,x2) / (x2['esq'] - x1['esq'])
        r += -self.dI(x1) / (x2['esq'] - x1['esq'])
        return r
    def K123(self,x1,x2,x3):
        r  = self.I(x1) / (x1['esq'] - x2['esq']) / (x1['esq'] - x3['esq'])
        r += self.I(x2) / (x2['esq'] - x1['esq']) / (x2['esq'] - x3['esq'])
        r += self.I(x3) / (x3['esq'] - x1['esq']) / (x3['esq'] - x2['esq'])
        return r

    def Fpi_xpt_nlo(self,x,lec):
        r  = -self.I(x['p2'])
        r += -0.5 * self.I(x['k2'])
        r += lec['L5'] * (4*np.pi)**2 * 4 * x['p2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*np.pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
        return r

    def FK_xpt_nlo(self,x,lec):
        r  = -3./8 * self.I(x['p2'])
        r += -3./4 * self.I(x['k2'])
        r += -3./8 * self.I(x['e2'])
        r += lec['L5'] * (4*np.pi)**2 * 4 * x['k2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*np.pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
        return r

    def Fpi_ma_nlo(self,x,lec):
        r  = -self.I(x['ju'])
        r += -0.5 * self.I(x['ru'])
        r += lec['L5'] * (4*np.pi)**2 * 4 * x['p2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*np.pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
        return r

    def FK_ma_nlo(self,x,lec):
        self.meson = 'mjuL'
        r  = -1./2 * self.I(x['ju'])
        r +=  1./8 * self.I(x['p2'])
        r += -1./4 * self.I(x['ru'])
        r += -1./2 * self.I(x['js'])
        r += -1./4 * self.I(x['rs'])
        r +=  1./4 * self.I(x['s2'])
        r += -3./8 * self.I(x['x2'])
        r += lec['L5'] * (4*np.pi)**2 * 4 * x['k2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*np.pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])

        r += x['dju2']             * -1./8  * self.dI(x['p2'])
        r += x['dju2']             *  1./4  * self.K(x['p2'],x['x2'])
        r += x['dju2']**2          * -1./24 * self.K21(x['p2'],x['x2'])
        r += x['dju2'] * x['drs2'] *  1./12 * self.K21(x['s2'],x['x2'])
        r += x['dju2'] * x['drs2'] * -1./6  * self.K123(x['p2'],x['s2'],x['x2'])
        r += x['drs2']             *  1./4  * self.K(x['s2'],x['x2'])
        r += x['drs2'] * x['k2']['esq'] * -1./6  * self.K21(x['s2'],x['x2'])
        r += x['drs2'] * x['p2']['esq'] *  1./6  * self.K21(x['s2'],x['x2'])

        return r

    def xpt_nlo(self,x,p):
        r = dict()
        for e in x:
            x_par,lec = self.make_x_lec(x,p,e)
            r[e]  = self.counterterms(x_par,lec)
            r[e] += 1.
            r[e] +=  self.FK_xpt_nlo(x_par,lec)
            r[e] += -self.Fpi_xpt_nlo(x_par,lec)
        return r
    def xpt_ratio_nlo(self,x,p):
        r = dict()
        for e in x:
            x_par,lec = self.make_x_lec(x,p,e)
            r[e]  = self.counterterms(x_par,lec)
            num   = 1. + self.FK_xpt_nlo(x_par,lec)
            den   = 1. + self.Fpi_xpt_nlo(x_par,lec)
            r[e] += num / den
        return r
    def ma_nlo(self,x,p):
        r = dict()
        for e in x:
            x_par,lec = self.make_x_lec(x,p,e)
            r[e]  = self.counterterms(x_par,lec)
            r[e] += 1.
            r[e] += self.FK_ma_nlo(x_par,lec)
            r[e] += -self.Fpi_ma_nlo(x_par,lec)
        return r
    def ma_ratio_nlo(self,x,p):
        r = dict()
        for e in x:
            x_par,lec = self.make_x_lec(x,p,e)
            r[e]  = self.counterterms(x_par,lec)
            num   = 1. + self.FK_ma_nlo(x_par,lec)
            den   = 1. + self.Fpi_ma_nlo(x_par,lec)
            r[e] += num / den
        return r
    def ma_longform_nlo(self,x,p):
        ''' for debugging - a cross check expression '''
        r = dict()
        for e in x:
            x_par,lec = self.make_x_lec(x,p,e)
            r[e]  = self.counterterms(x_par,lec)
            r[e] += 1.
            r[e] += x_par['ju'] * np.log(x_par['ju'])
            r[e] += 0.5 * x_par['ru'] * np.log(x_par['ru'])
            # (+) kaon log terms
            r[e] += -0.5 * x_par['ju'] * np.log(x_par['ju'])
            r[e] += -1./4 * x_par['ru'] * np.log(x_par['ru'])
            r[e] += -0.5*x_par['js'] * np.log(x_par['js'])
            r[e] += -1./4 * x_par['rs'] * np.log(x_par['rs'])
            r[e] += -x_par['dju2'] / 8
            r[e] +=  x_par['dju2']**2 / 24 / (x_par['x2'] - x_par['p2'])
            r[e] +=  x_par['drs2'] * (x_par['k2']-x_par['p2']) / 6 / (x_par['x2']-x_par['s2'])
            r[e] += -x_par['dju2'] * x_par['drs2'] / 12 / (x_par['x2'] - x_par['s2'])
            r[e] += np.log(x_par['p2'])/24 * (3*x_par['p2'] \
                    - 3*x_par['dju2']*(x_par['x2']+x_par['p2'])/(x_par['x2']-x_par['p2']) \
                    + x_par['dju2']**2 * x_par['x2']/(x_par['x2']-x_par['p2'])**2\
                    -4*x_par['dju2']*x_par['drs2']*x_par['p2']/(x_par['x2']-x_par['p2'])/(x_par['s2']-x_par['p2'])\
                    )
            r[e] += -x_par['x2']/24 * np.log(x_par['x2'])*(9 \
                    -6*x_par['dju2']/(x_par['x2']-x_par['p2']) \
                    + x_par['dju2']**2/(x_par['x2']-x_par['p2'])**2\
                    +x_par['drs2']*(4*(x_par['k2']-x_par['p2'])+6*(x_par['s2']-x_par['x2']))/(x_par['x2']- x_par['s2'])**2 \
                    -2*x_par['dju2']*x_par['drs2']*(2*x_par['s2']-x_par['p2']-x_par['x2'])/(x_par['x2']-x_par['s2'])**2/(x_par['x2']-x_par['p2'])\
                    )
            r[e] += np.log(x_par['s2'])/12 * (3*x_par['s2'] \
                    +x_par['drs2']*(3*x_par['s2']**2 + 2*(x_par['k2']-x_par['p2'])*x_par['x2'] -3*x_par['s2']*x_par['x2'])/(x_par['x2']-x_par['s2'])**2\
                    -x_par['dju2']*x_par['drs2']*(2*x_par['s2']**2 - x_par['x2']*(x_par['s2']+x_par['p2']))/(x_par['x2']-x_par['s2'])**2 / (x_par['s2']-x_par['p2'])\
                    )
        return r

    def fit_data(self):
        x = self.prune_x()
        y = self.prune_data()
        p = self.prune_priors()
        self.fit = lsqfit.nonlinear_fit(data=(x,y),prior=p,fcn=self.fit_function)

    def report_phys_point(self,phys_point):
        # get copy of all self attributes
        self_dict = self.__dict__.copy()
        # set physical point
        x_phys = dict()
        x_phys['phys'] = {k:np.inf for k in ['mpiL','mkL']}
        x_phys['phys']['alphaS'] = 0.

        y_phys = dict()
        y_phys['phys'] = phys_point['FK'] / phys_point['Fpi']
        y_phys['FLAG'] = phys_point['FKFPi_FLAG']

        p_phys = dict()
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p_phys[k] = 0.
        Lchi_phys = phys_point['Lchi']
        p_phys[('phys','p2')] = phys_point['mpi']**2 / Lchi_phys**2
        p_phys[('phys','k2')] = phys_point['mk']**2  / Lchi_phys**2
        p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
        p_phys[('phys','a2')] = 0.
        for k in self.fit.p:
            if isinstance(k,str):
                print(k,self.fit.p[k])
        for k in ['L5','L4','k_4','p_4','kp_6','k_6','p_6']:
            if k in self.fit.p:
                p_phys[k] = self.fit.p[k]
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio_nlo
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt_nlo
        self.fv = False

        print('chi2/dof [dof] = %.2f [%d]    Q = %.2e    logGBF = %.3f' \
            %(self.fit.chi2/self.fit.dof,self.fit.dof,self.fit.Q,self.fit.logGBF))
        print(self.fit_function(x_phys,p_phys),'\n')

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

    def check_fit(self,phys_point):
        # get copy of all self attributes
        self_dict = self.__dict__.copy()
        # set physical point
        x_phys = dict()
        x_phys['phys'] = {k:np.inf for k in ['mpiL','mkL']}
        x_phys['phys']['alphaS'] = 0.

        y_phys = dict()
        y_phys['phys'] = phys_point['FK'] / phys_point['Fpi']
        y_phys['FLAG'] = phys_point['FKFPi_FLAG']

        p_phys = dict()
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p_phys[k] = 0.
        Lchi_phys = phys_point['Lchi']
        p_phys[('phys','p2')] = phys_point['mpi']**2 / Lchi_phys**2
        p_phys[('phys','k2')] = phys_point['mk']**2  / Lchi_phys**2
        p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
        p_phys[('phys','a2')] = 0.
        #for k in self.fit.p:
            #if isinstance(k,str):
                #print(k,self.fit.p[k])
        for k in ['L5','L4','k_4','p_4','kp_6','k_6','p_6']:
            if k in self.fit.p:
                p_phys[k] = self.fit.p[k]
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio_nlo
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt_nlo
        self.fv = False

        #print('chi2/dof [dof] = %.2f [%d]    Q = %.2e    logGBF = %.3f' \
        #    %(self.fit.chi2/self.fit.dof,self.fit.dof,self.fit.Q,self.fit.logGBF))
        fkfpi = self.fit_function(x_phys,p_phys)

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return self.fit.chi2/self.fit.dof, self.fit.dof, self.fit.logGBF, self.fit.p['L5'], self.fit.p['k_4'], self.fit.p['p_4'], self.fit.p['s_4'], fkfpi['phys']


    def vs_epi(self,epi_range, ax):
        # get copy of all self attributes
        self_dict = self.__dict__.copy()
        # switch to continuum fit func
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio_nlo
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt_nlo
        self.fv = False
        # set x,y,p
        x = dict()
        y_phys = dict()
        p = dict()
        for k in self.fit.p:
            if isinstance(k,str):
                print(k,self.fit.p[k])
        for k in self.lec_nnlo + self.lec_nnnlo:
            if k in self.fit.p:
                p[k] = self.fit.p[k]

        Lchi = epi_range['Lchi']
        x_plot = []
        y_plot = {}
        y_plot['a0']  = []
        y_plot['a09'] = []
        y_plot['a12'] = []
        y_plot['a15'] = []
        for mpi in epi_range['mpi']:
            x[mpi] = {k:np.inf for k in ['mpiL','mkL']}
            x[mpi]['alphaS'] = 0.

            p[(mpi,'p2')] = mpi**2 / Lchi**2
            p[(mpi,'k2')] = epi_range['mk']**2  / Lchi**2
            p[(mpi,'e2')] = 4./3*p[(mpi,'k2')] - 1./3 * p[(mpi,'p2')]
            p[(mpi,'a2')] = 0.
            x_plot.append(p[(mpi,'p2')].mean)
            y_plot['a0'].append(self.fit_function(x,p)[mpi])
            # finite a
            p[(mpi,'a2')] = self.p_init[('a09m310','aw0')]**2 / (4 * np.pi)
            y_plot['a09'].append(self.fit_function(x,p)[mpi])
            p[(mpi,'a2')] = self.p_init[('a12m310','aw0')]**2 / (4 * np.pi)
            y_plot['a12'].append(self.fit_function(x,p)[mpi])
            p[(mpi,'a2')] = self.p_init[('a15m310','aw0')]**2 / (4 * np.pi)
            y_plot['a15'].append(self.fit_function(x,p)[mpi])

        x_plot = np.array(x_plot)
        y  = np.array([k.mean for k in y_plot['a0']])
        dy = np.array([k.sdev for k in y_plot['a0']])

        ax.fill_between(x_plot, y-dy, y+dy,color='#b36ae2',alpha=0.4)
        colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9'}
        for a in ['a15','a12','a09']:
            y = np.array([k.mean for k in y_plot[a]])
            ax.plot(x_plot, y, color=colors[a])
        epi_phys =  epi_range['mpi_phys']**2 / Lchi**2
        ax.axvline(epi_phys.mean,linestyle='--',color='#a6aaa9')
        ax.axvspan(epi_phys.mean -epi_phys.sdev, epi_phys.mean +epi_phys.sdev,
            alpha=0.4, color='#a6aaa9')
        ax.axvline(epi_phys.mean,linestyle='--',color='#a6aaa9')
        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return ax

    def shift_phys_epsK(self,phys_point):
        '''
            Compute the shift in FK/FPi from the fit, to the point where
            m_K = m_K^phys
            if MA:
                m_ju = m_pi
                m_js = m_K^phys
                m_ru = m_K^phys
                m_rs = m_ss
                dju  = 0
                drs  = 0
            and apply this shift to the numerical data, so that it can be
            plotted with a fit that only depends upon eps_pi = m_pi**2 / L_chi**2

            1: compute FK/FPi from the particular fit
                y_fit

            2: computee FK/FPi from the corresponding XPT fit function + a**2 terms
                y_shift
            3: shift all the data by
                FK/FPi -> FK/FPi + (y_shift - y_fit).mean

            return a dictionary of shifts
        '''
        # get copy of all self attributes
        self_dict = self.__dict__.copy()

        # determine fitted values of FK/FPi
        y_shift = dict()
        print('correcting FK/Fpi data for plot')
        print('for model %s' %self.switches['ansatz']['model'])
        x = self.prune_x()
        for e in x:
            Lchi = self.p_init[(e,'Lchi_'+self.switches['scale'])]
            p    = self.prune_priors()
            for k in self.lec_nnlo + self.lec_nnnlo:
                if k in self.fit.p:
                    p[k] = self.fit.p[k]
            y_shift[e] = -self.fit_function(x,p)[e]
        # switch to XPT + a**2 terms
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio_nlo
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt_nlo
        self.fv = False
        if self.switches['verbose']:
            print("%9s %12s %12s" %('ensemble','data', 'shift'))
            print('-------- | ----------- | ----------')
        for e in x:
            x_tmp = dict()
            p_tmp = dict()
            x_tmp[e] = {k:np.inf for k in ['mpiL','mkL']}
            x_tmp[e]['alphaS'] = self.x[e]['alphaS']
            for k in self.lec_nnlo + self.lec_nnnlo:
                if k in self.fit.p:
                    p_tmp[k] = self.fit.p[k]
            Lchi = self.p_init[(e,'Lchi_'+self.switches['scale'])]
            p_tmp[(e,'p2')] = self.p_init[(e,'mpi')]**2 / Lchi**2
            p_tmp[(e,'k2')] = phys_point['mk']**2 / phys_point['Lchi']**2
            p_tmp[(e,'e2')] = 4./3 *p_tmp[(e,'k2')] - 1./3 *p_tmp[(e,'p2')]
            p_tmp[(e,'a2')] = self.p_init[(e,'aw0')]**2 / (4 * np.pi)
            y_shift[e] += self.fit_function(x_tmp,p_tmp)[e]
            if self.switches['verbose']:
                print("%9s %12s %12s" %(e, self.y[e], y_shift[e].mean))

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return y_shift
