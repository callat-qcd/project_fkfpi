## chipt functions
import scipy.special as spsp
import numpy as np

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
        dju = priors['a2di'] / Lchi**2
        drs = priors['a2di'] / Lchi**2
        fv_mns_inv  = 2 * ju * k1ju + 0.5 * p2 * k1pi
        fv_mns_inv -= dju/8 * (2*k1pi - k0pi - k2pi)
        fv_mns_inv -= dju * p2 / (x2-p2) * k1pi
        fv_mns_inv += dju**2 / 24 * (4*p2*k1pi/(x2-p2)**2 \
            +(2*k1pi - k0pi -k2pi)/(x2-p2))
        fv_mns_inv += 2./3 * dju * drs * p2 * k1pi/(x2-p2)/(s2-p2)
    return fv_mns_inv

class Fit(object):
    def __init__(self,switches):
        self.switches = switches
        #self.n = switches['ansatz']['truncation']
    def counterterms(self,x,p,e):
        Lchi = x[e]['Lchi_'+self.switches['scale']]
        p2   = p[(e,'mpi')]**2 / Lchi**2
        k2   = p[(e,'mk')]**2 / Lchi**2
        e2   = 4./3 * k2 - 1./3 * p2
        a2   = p[(e,'aw0')]**2
        if self.switches['ansatz']['model'].split('_')[-1] in ['nlo','nnlo','nnnlo']:
            ct = 4 * (k2 - p2) * (4*np.pi)**2 * p['L5'] # m^2
        if self.switches['ansatz']['model'].split('_')[-1] in ['nnlo','nnnlo']:
            ct += a2 * (k2 -p2) * p['s4'] # a^2 m^2
            ct += (k2 -p2)**2   * p['c4'] # m^4
            ct += k2*(k2 -p2)   * p['d4'] # m^4
            ct += p2* (k2 -p2)  * p['e4'] # m^4
        if self.switches['ansatz']['model'].split('_')[-1] in ['nnnlo']:
            ct += a2**2 *(k2 -p2)   * p['s6']  # a^4*m^2
            ct += a2 * (k2 -p2)**2  * p['sc6'] # a^2*m^4
            ct += a2 * k2 *(k2 -p2) * p['sd6']
            ct += a2 * p2 *(k2 -p2) * p['se6']
        return ct

    def fit_function(self,x,p):
        r = dict()
        for e in x:
            Lchi = x[e]['Lchi_'+self.switches['scale']]
            p2   = p[(e,'mpi')]**2 / Lchi**2
            k2   = p[(e,'mk')]**2 / Lchi**2
            e2   = 4./3 * k2 - 1./3 * p2
            a2   = p[(e,'aw0')]**2 / 4 / np.pi
            if self.switches['debug']:
                print('DEBUG: fit_function x,y values')
                print(e,p2,k2,e2,a2)

            r[e]  = 1.
            r[e] += self.counterterms(x,p,e)

            if self.switches['ansatz']['model'].split('_')[0] == 'xpt':
                r[e] += 5./8 * p2 * np.log(p2)
                r[e] -= 1./4 * k2 * np.log(k2)
                r[e] -= 3./8 * e2 * np.log(e2)

            elif self.switches['ansatz']['model'].split('_')[0] == 'ma':
                s2 = p[(e,'mss')]**2 / Lchi**2
                if self.switches['ansatz']['a2dm'] == 'avg':
                    ju = p2 + p['a2dm'] / Lchi**2
                    sj = k2 + p['a2dm'] / Lchi**2
                    ru = k2 + p['a2dm'] / Lchi**2
                    rs = s2 + p['a2dm'] / Lchi**2
                elif self.switches['ansatz']['a2dm'] == 'individual':
                    ju = p[(e,'mju')]**2 / Lchi**2
                    sj = p[(e,'mjs')]**2 / Lchi**2
                    ru = p[(e,'mru')]**2 / Lchi**2
                    rs = p[(e,'mrs')]**2 / Lchi**2

                x2  = e2 + p[(e,'a2DI')] / Lchi**2
                dju = p[(e,'a2DI')] / Lchi**2
                drs = p[(e,'a2DI')] / Lchi**2
                # (-) pion log terms
                r[e] += ju * np.log(ju)
                r[e] += 0.5 * ru * np.log(ru)
                # (+) kaon log terms
                r[e] += -0.5 * ju * np.log(ju)
                r[e] += -1./4 * ru * np.log(ru)
                r[e] += -0.5*sj * np.log(sj)
                r[e] += -1./4 * rs * np.log(rs)
                r[e] += -dju / 8
                r[e] +=  dju**2 / 24 / (x2 - p2)
                r[e] +=  drs * (k2-p2) / 6 / (x2-s2)
                r[e] += -dju * drs / 12 / (x2 - s2)
                r[e] += np.log(p2)/24 * (3*p2 \
                        - 3*dju*(x2+p2)/(x2-p2) \
                        + dju**2 * x2/(x2-p2)**2\
                        -4*dju*drs*p2/(x2-p2)/(s2-p2)\
                        )
                r[e] += -x2/24 * np.log(x2)*(9 \
                        -6*dju/(x2-p2) \
                        + dju**2/(x2-p2)**2\
                        +drs*(4*(k2-p2)+6*(s2-x2))/(x2-s2)**2 \
                        -2*dju*drs*(2*s2-p2-x2)/(x2-s2)**2/(x2-p2)\
                        )
                r[e] += np.log(s2)/12 * (3*s2 \
                        +drs*(3*s2**2 + 2*(k2-p2)*x2 -3*s2*x2)/(x2-s2)**2\
                        -dju*drs*(2*s2**2 - x2*(s2+p2))/(x2-s2)**2 / (s2-p2)\
                        )
            #r[e] = np.array(r[e])
        return r
