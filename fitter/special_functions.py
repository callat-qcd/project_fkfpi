import gvar as gv
import numpy as np
import scipy.special as ss

# gvar version of modified Bessel function of the second kind
def fcn_Kn(n, g):
    ymean = ss.kn(n, gv.mean(g))
    
    if isinstance(g, gv._gvarcore.GVar):
        ysdev = np.abs(ss.kvp(n, gv.mean(g), 1)) * gv.sdev(g)
        return gv.gvar(ymean, ysdev)
    else:
        return ymean

# I(m) in notes
def fcn_I_m(m, mL, mu, order):
    c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]

    output = (m/(4*np.pi))**2 *np.log((m/mu)**2)

    for n in range(1, np.min((order+1, 11))):
        output = output + (m/(4*np.pi))**2 *(c[n]/(mL *np.sqrt(n))) *fcn_Kn(1, mL *np.sqrt(n))
    return output

# dI(m) in notes
def fcn_dI_m(m, mL, mu, order):
    c = [None, 6, 12, 8, 6, 24, 24, 0, 12, 30, 24]

    output = 1/(4 *np.pi)**2 + fcn_I_m(m, mL, mu, order)/m**2

    for n in range(1, np.min((order, 11))):
        output = output + c[n]/(4 *np.pi)**2 *(
            2 *fcn_Kn(1, mL *np.sqrt(n))/(mL *np.sqrt(n))
            - fcn_Kn(0, mL *np.sqrt(n))
            - fcn_Kn(2, mL *np.sqrt(n))
        )
    return output

# K(m, M) in notes
def fcn_K_mM((m, M), mL, mu, order):
    output = 1.0/(M**2 - m**2) *(fcn_I_m(M, mL, mu, order) - fcn_I_m(m, mL, mu, order))
    return output

# K_21(m, M) in notes
def fcn_K21_mM((m, M), mL, mu, order):
    output = (
        1.0/(M**2 - m**2)**2 *(fcn_I_m(M, mL, mu, order) - fcn_I_m(m, mL, mu, order))
        - 1.0/(M**2 - m**2) *fcn_dI_m(m, mL, mu, order)
    )
    return output

def fcn_K_m1m2m3((m1, m2, m3), mL, mu, order):
    output = (
        1.0/((m1**2 - m2**2) *(m1**2 - m3**2)) *fcn_I_m(m1, mL, mu, order)
        + 1.0/((m2**2 - m1**2) *(m2**2 - m3**2)) *fcn_I_m(m2, mL, mu, order)
        + 1.0/((m3**2 - m1**2) *(m3**2 - m2**2)) *fcn_I_m(m3, mL, mu, order)
    )
    return output
