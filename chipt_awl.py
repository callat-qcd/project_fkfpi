## chipt functions
import scipy.special as spsp
import numpy as np
import gvar as gv

def fv_correction(switches,fv_params,priors):
    mpiL = fv_params['mpiL']
    mjuL = fv_params['mjuL']
    Lchi = fv_params['Lchi']
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
        dju2 = priors['a2di'] / Lchi**2
        drs2 = priors['a2di'] / Lchi**2
        fv_mns_inv  = 2 * ju * k1ju + 0.5 * p2 * k1pi
        fv_mns_inv -= dju2/8 * (2*k1pi - k0pi - k2pi)
        fv_mns_inv -= dju2 * p2 / (x2-p2) * k1pi
        fv_mns_inv += dju2**2 / 24 * (4*p2*k1pi/(x2-p2)**2 \
            +(2*k1pi - k0pi -k2pi)/(x2-p2))
        fv_mns_inv += 2./3 * dju2 * drs2 * p2 * k1pi/(x2-p2)/(s2-p2)
    return fv_mns_inv

class Fit(object):
    def __init__(self,switches,xyp_init):
        self.switches = switches
        self.x        = xyp_init['x']
        self.y        = xyp_init['y']
        self.p_init   = xyp_init['p']
        self.order    = self.switches['ansatz']['model'].split('_')[-1]
        self.eft      = self.switches['ansatz']['model'].split('_')[0]
        self.fv       = self.switches['ansatz']['FV']
        if self.fv:
            #self.fv = dict()
            cn = np.array([6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24])
            n_mag = np.sqrt(np.arange(1,len(cn)+1,1))
            #self.fv[(e,'ju',mpiL)] =

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
                p[(e,'a2')]   = self.p_init[(e,'aw0')] / (4 * np.pi)
        ''' Add LECs '''
        if self.order in ['nlo','nnlo','nnnlo']:
            p['L5'] = self.p_init['L5']
            if 'ratio' in self.eft:
                p['L4'] = self.p_init['L4']
        if self.order in ['nnlo','nnnlo']:
            p['s4'] = self.p_init['s4']
            if self.switches['ansatz']['alpha_S']:
                p['s4aS'] = self.p_init['s4aS']
            p['c4'] = self.p_init['c4']
            p['d4'] = self.p_init['d4']
            p['e4'] = self.p_init['e4']
        if self.order in ['nnnlo']:
            p['s6']  = self.p_init['s6']
            p['sc6'] = self.p_init['sc6']

        return p

    def make_x_lec(self,x,p,e):
        x_par = gv.BufferDict()
        x_par['p2']   = {'esq':p[(e,'p2')], 'mL':self.x[e]['mpiL']}
        x_par['k2']   = {'esq':p[(e,'k2')], 'mL':self.x[e]['mkL']}
        meL = np.sqrt(4./3 * self.x[e]['mkL']**2 -1./3 * self.x[e]['mpiL']**2)
        x_par['e2']   = {'esq':p[(e,'e2')], 'mL':meL}
        if 'ma' in self.eft:
            x_par['s2']   = {'esq':p[(e,'s2')], 'mL':self.x[e]['mssL']}
            x_par['ju']   = {'esq':p[(e,'ju')], 'mL':self.x[e]['mjuL']}
            x_par['js']   = {'esq':p[(e,'js')], 'mL':self.x[e]['mjsL']}
            x_par['ru']   = {'esq':p[(e,'ru')], 'mL':self.x[e]['mruL']}
            x_par['rs']   = {'esq':p[(e,'rs')], 'mL':self.x[e]['mrsL']}
            mxL = np.sqrt(p[(e,'x2')] / p[(e,'e2')]) * meL
            x_par['x2']   = {'esq':p[(e,'x2')], 'mL':mxL}
            x_par['dju2'] = p[(e,'dju2')]
            x_par['drs2'] = p[(e,'drs2')]
        if self.order in ['nnlo','nnnlo']:
            x_par['a2']   = p[(e,'a2')]
        if self.switches['ansatz']['alpha_S']:
            x_par['alpha_S'] = x[e]['alpha_S']

        lec = {key:val for key,val in p.items() if isinstance(key,str)}

        return x_par,lec

    def counterterms(self,x,lec):
        ct = 0.
        # NLO terms in FK and Fpi functions
        '''
        ct = 4 * (k2 - p2) * (4*np.pi)**2 * p['L5']
        '''
        if self.order in ['nnlo','nnnlo']:
            ct += x['a2'] * (x['k2'] -x['p2']) * lec['s4'] # a^2 m^2
            if self.switches['ansatz']['alpha_S']:
                ct += x['a2'] * x['alpha_S'] * (x['k2'] -x['p2']) * p['s4aS']
            ct += (x['k2'] -x['p2'])**2        * lec['c4'] # m^4
            ct += x['k2']*(x['k2'] -x['p2'])   * lec['d4'] # m^4
            ct += x['p2']*(x['k2'] -x['p2'])   * lec['e4'] # m^4
        if self.order in ['nnnlo']:
            ct += x['a2']**2 *(x['k2'] -x['p2'])        * lec['s6']  # a^4*m^2
            ct += x['a2'] * (x['k2'] -x['p2'])**2       * lec['sc6'] # a^2*m^4
            ct += x['a2'] * x['k2'] *(x['k2'] -x['p2']) * lec['sd6']
            ct += x['a2'] * x['p2'] *(x['k2'] -x['p2']) * lec['se6']
        return ct

    def I(self,x):
        r = x['esq'] * np.log(x['esq'])
        if self.fv:
            r += self.fv_I(x['esq'],x['mL'])
        return r
    def dI(self,x):
        r = 1 + np.log(x['esq'])
        if self.fv:
            r += self.fv_dI(x['mL'])
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

    def fit_function(self,x,p):
        r = dict()
        for e in x:

            x_par,lec = self.make_x_lec(x,p,e)
            r[e] = self.counterterms(x_par,lec)

            if self.eft == 'xpt':
                r[e] += 1.
                r[e] +=  self.FK_xpt_nlo(x_par,lec)
                r[e] += -self.Fpi_xpt_nlo(x_par,lec)

            elif self.eft == 'xpt-ratio':
                num   = 1. + self.FK_xpt_nlo(x_par,lec)
                den   = 1. + self.Fpi_xpt_nlo(x_par,lec)
                r[e] += num / den

            elif self.eft == 'ma':
                r[e] += 1.
                r[e] += self.FK_ma_nlo(x_par,lec)
                r[e] += -self.Fpi_ma_nlo(x_par,lec)

            elif self.eft == 'ma-ratio':
                num   = 1. + self.FK_ma_nlo(x_par,lec)
                den   = 1. + self.Fpi_ma_nlo(x_par,lec)
                r[e] += num / den

            elif self.eft == 'ma-longform':
                ''' for debugging - a cross check expression '''
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
