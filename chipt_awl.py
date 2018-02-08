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
