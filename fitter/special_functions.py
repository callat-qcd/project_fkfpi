import gvar as gv
import numpy as np
import scipy.special as ss
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/py_chiron")
import chiron


### volume corrections
# gvar version of modified Bessel function of the second kind
def fcn_Kn(n, g):
    # input is a gvar scalar
    if isinstance(g, gv._gvarcore.GVar):
        f = ss.kn(n, gv.mean(g))
        dfdg = ss.kvp(n, gv.mean(g), 1)
        return gv.gvar_function(g, f, dfdg)

    # input is a gvar vector
    elif hasattr(g, "__len__") and isinstance(g[0], gv._gvarcore.GVar):
        f = ss.kn(n, gv.mean(g))
        dfdg = ss.kvp(n, gv.mean(g), 1)
        return [gv.gvar_function(g[j], f[j], dfdg[j]) for j in range(len(g))]

    # input is not a gvar variable
    else:
        return ss.kn(n, gv.mean(g))

# I(m) in notes
def fcn_I_m(m, L, mu, order):
    c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]

    output = (m/(4*np.pi))**2 *np.log((m/mu)**2)

    for n in range(1, np.min((order+1, 11))):
        output = output + (m**2/(4*np.pi**2)) *(c[n]/(m *L *np.sqrt(n))) *fcn_Kn(1, m *L *np.sqrt(n))
    return output

# dI(m) in notes
def fcn_dI_m(m, L, mu, order):
    c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]

    output = 1/(4 *np.pi)**2 + fcn_I_m(m, L, mu, order)/m**2

    for n in range(1, np.min((order, 11))):
        output = output + c[n]/(4 *np.pi)**2 *(
            2 *fcn_Kn(1, m *L *np.sqrt(n))/(m *L *np.sqrt(n))
            - fcn_Kn(0, m *L *np.sqrt(n))
            - fcn_Kn(2, m *L *np.sqrt(n))
        )
    return output

# K(m, M) in notes
def fcn_K_mM(mM, L, mu, order):
    m, M = mM
    output = 1.0/(M**2 - m**2) *(fcn_I_m(M, L, mu, order) - fcn_I_m(m, L, mu, order))
    return output

# K_21(m, M) in notes
def fcn_K21_mM(mM, L, mu, order):
    m, M = mM
    output = (
        1.0/(M**2 - m**2)**2 *(fcn_I_m(M, L, mu, order) - fcn_I_m(m, L, mu, order))
        - 1.0/(M**2 - m**2) *fcn_dI_m(m, L, mu, order)
    )
    return output

def fcn_K_m1m2m3(m1m2m3, L, mu, order):
    m1, m2, m3 = m1m2m3
    output = (
        1.0/((m1**2 - m2**2) *(m1**2 - m3**2)) *fcn_I_m(m1, L, mu, order)
        + 1.0/((m2**2 - m1**2) *(m2**2 - m3**2)) *fcn_I_m(m2, L, mu, order)
        + 1.0/((m3**2 - m1**2) *(m3**2 - m2**2)) *fcn_I_m(m3, L, mu, order)
    )
    return output

### functions from hep-ph/1711.11328
def fcn_Kr_j(j, mpi, mk, lam2_chi):

    c = {}
    if j == 1:
        c['kk'], c['kp'], c['pp'] = [0, 11./24, 131./192]
    elif j == 2:
        c['kk'], c['kp'], c['pp'] = [0, -41./96, -3./32]
    elif j == 3:
        c['kk'], c['kp'], c['pp'] = [0, 13./24, 59./96]
    elif j == 4:
        c['kk'], c['kp'], c['pp'] = [17./36, 7./144, 0]
    elif j == 5:
        c['kk'], c['kp'], c['pp'] = [-163./144, -67./288, 3./32]
    elif j ==6:
        c['kk'], c['kp'], c['pp'] = [241./288, -13./72, -61./192]

    output = (
        + c['kk'] *mk**2
        + c['kp'] *mk *mpi
        + c['pp'] *mpi**2
    ) / lam2_chi

    return output

def fcn_Cr_j(j, mpi, mk, lam2_chi, p):

    c = {}
    if j in [1, 2, 3]:
        if j == 1:
            c['kk'] = 0
            c['kp'] = -(+ 7./9
                        + (4 *np.pi)**2 *(
                            11./2 *p['L_5']
                            )
                        )
            c['pp'] = -(+ 113./72
                        + (4 *np.pi)**2 *(
                            4 *p['L_1'] + 10 *p['L_2'] + 13./2 *p['L_3'] - 21./2 * p['L_5']
                            )
                        )

        elif j == 2:
            c['kk'] = +(+ 53./96
                        + (4 *np.pi)**2 *(
                            4 *p['L_1'] + 10 *p['L_2'] + 5 *p['L_3'] - 5 *p['L_5']
                            )
                        )
            c['kp'] = +(+ 209./144
                        + (4 *np.pi)**2 *(
                            3 *p['L_5']
                            )
                        )
            c['pp'] = 0

        elif j == 3:
            c['kk'] = +(+ 13./8
                        + (4 *np.pi)**2 *(
                            8./3 *p['L_3'] - 2./3 *p['L_5'] - 16 *p['L_7'] - 8 *p['L_8']
                            )
                        )
            c['kp'] = -(+ 4./9
                        + (4 *np.pi)**2 *(
                            4./3 *p['L_3'] + 25./6 *p['L_5'] - 32 *p['L_7'] - 16 *p['L_8']
                            )
                        )
            c['pp'] = +(+ 19./288
                        + (4 *np.pi)**2 *(
                            1./6 *p['L_3'] + 11./6 *p['L_5'] - 16 *p['L_7'] - 8 *p['L_8']
                            )
                        )

        output = (
            + c['kk'] *mk**2
            + c['kp'] *mk *mpi
            + c['pp'] *mpi**2
        ) / lam2_chi

        return output

    elif j ==4:
        return None

    else:
        return None

def fcn_FF(x):
    stepSize = 1e-7
    # input is a vector
    if hasattr(x, "__len__"):
        if isinstance(x[0], gv.GVar):
            f = [chiron.FF(e.mean) for e in x]
            dfdx = [0.5*(chiron.FF(e.mean+stepSize) - chiron.FF(e.mean-stepSize))/stepSize for e in x]
            return [gv.gvar_function(x[j], f[j], dfdx[j]) for j in range(len(x))]
        else:
            return [chiron.FF(e) for e in x]
    # input is a scalar
    else:
        if isinstance(x, gv.GVar):
            f = chiron.FF(x.mean)
            dfdx = 0.5*(chiron.FF(x.mean+stepSize) - chiron.FF(x.mean-stepSize))/stepSize
            return gv.gvar_function(x, f, dfdx)
        else:
            return chiron.FF(x)
