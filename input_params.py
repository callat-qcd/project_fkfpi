import gvar as gv
import numpy as np

switches = dict()
# Make two sets of ensembles - so we can plot data excluded from fit
switches['ensembles'] = [
    'a15m400'  ,'a12m400' ,'a09m400',
    'a15m350'  ,'a12m350' ,'a09m350',
    'a15m310'  ,'a12m310' ,'a09m310','a06m310L',
    'a15m220'  ,'a12m220S','a12m220','a12m220L','a09m220',
    'a15m135XL','a12m130','a09m135',]
switches['ensembles_fit'] = [
    'a15m400'  ,'a12m400' ,'a09m400',
    'a15m350'  ,'a12m350' ,'a09m350',
    'a15m310'  ,'a12m310' ,'a09m310','a06m310L',
    'a15m220'  ,'a12m220' ,'a09m220','a12m220S','a12m220L',
    'a15m135XL','a12m130' ,'a09m135',
    ]

# FIT MODELS
switches['ansatz'] = dict()
switches['ansatz']['models'] = ['xpt_nnnlo_FV']
''' The full list of models can be rather long.  The sys switches help loop
    over them.  Example other base models are
        taylor_nnnlo_FV
        ma_nnnlo_FV
'''
# SYSTEMATIC SWITCHES
switches['sys'] = dict()     # these cause the fitter to loop over various options
switches['sys']['Lam_chi']   = True # FF = Fpi Fpi, Fpi FK, FK FK
switches['sys']['alphaS']    = True # include alphaS at NNLO?
switches['sys']['nnlo_ct']   = True # NNLO = full XPT or just counterterm
switches['sys']['ratio']     = True # use ratio version of NLO fit
# OLDER SYSTEMATICS - still work, but not used
switches['sys']['FV']        = False # turn on/off FV corrections
switches['sys']['logSq']     = False # only include logSq and ct (no log)
switches['sys']['a4']        = False # only include mpi^2 a^4 N3LO ct

switches['scales']           = ['PP','PK','KK'] # choices of F**2 to loop over
switches['scale']            = 'PP' # PP, PK, KK, LamChi = 4 * pi * sqrt(FA * FB)

switches['print_lattice']    = False # print data for paper - not fitting will occur

# Fitting options
switches['check_fit']        = False # print pieces of fit function - no fitting will occur
switches['bs_bias']          = True  # shift bs avg to b0?
switches['scipy']            = False # use scipy minimizer instead of gsl?
switches['print_fit']        = False # print lsqfit results?
switches['save_fits']        = True  # save fits in pickle file?
switches['prior_search']     = False # perform a crude grid search to optimize
switches['prior_verbose']    = False # NNLO and NNNLO prior widths

# Plotting options
switches['make_plots']       = True  # make plots
switches['save_figs']        = True  # save figures
switches['milc_compare']     = True  # compare with MILCs result

# DEBUGGING
switches['debug_models']     = False # print list of models being generated
switches['debug_save_fit']   = False # check pickling of fit works
switches['debug_shift']      = False # check the shifting of raw data to extrapolated points
switches['debug_bs']         = False # debug shape of bs lists

# PRIORS for fit
priors = dict()
priors['L5']   = gv.gvar(0, 0.04)
priors['L4']   = gv.gvar(0, 0.005)
# Ananthanarayan et al, 1711.11328 --> Bijnens, Ecker, 1405.6488
priors['L1']   = gv.gvar( 0.000372,0.000372)
priors['L2']   = gv.gvar( 0.000493,0.000493)
priors['L3']   = gv.gvar(-0.003070,0.003070)
priors['L6']   = gv.gvar( 0.000011,0.00011)
priors['L7']   = gv.gvar(-0.000340,0.000340)
priors['L8']   = gv.gvar( 0.000294,0.000294)

# Taylor priors - beyond NLO - use "XPT" NiLO priors
priors['c2'] = gv.gvar(0,10)
priors['t_fv'] = gv.gvar(0,100)

# Counter terms
nnlo_x = 2
nnlo_a = 2
nnlo_a = nnlo_x # from prior optimization, we found holding them the same is good
priors['k_4']   = gv.gvar(0.0, nnlo_x)
priors['p_4']   = gv.gvar(0.0, nnlo_x)
priors['s_4']   = gv.gvar(0.0, nnlo_a)
priors['saS_4'] = gv.gvar(0.0, nnlo_a)

n3lo_x = 5
n3lo_a = 5
n3lo_a = n3lo_x
priors['kp_6']  = gv.gvar(0.0, n3lo_x)
priors['k_6']   = gv.gvar(0.0, n3lo_x)
priors['p_6']   = gv.gvar(0.0, n3lo_x)
priors['s_6']   = gv.gvar(0.0, n3lo_a)
priors['sk_6']  = gv.gvar(0.0, n3lo_a)
priors['sp_6']  = gv.gvar(0.0, n3lo_a)

''' Physical point extrapolation
'''
phys_point = dict()
# Physical point values taken from FLAG
# FLAG[2019] = 1902.08191
# FLAG[2017] = 1607.00299
FPi_phys = gv.gvar(130.2/np.sqrt(2), 0.8/np.sqrt(2)) # FLAG[2019] (84)
FK_phys  = gv.gvar(155.7/np.sqrt(2), 0.7/np.sqrt(2)) # FLAG[2019] (85) - use NF=2+1 for consistency

phys_point = {
    'p':{
        'Fpi'     : FPi_phys,
        'FK'      : FK_phys,
        'Lchi_PP' : 4 * np.pi * FPi_phys,
        'Lchi_PK' : 4 * np.pi * np.sqrt(FPi_phys * FK_phys),
        'Lchi_KK' : 4 * np.pi * FK_phys,
        'mpi'     : gv.gvar(134.8, 0.3), #FLAG 2017 (16)
        #'mk'      : gv.gvar(494.2, 0.3), #FLAG 2017 (16) isospin symmetric
        'mk'      : gv.gvar(491.5, 0.7), # our estimate of MK+ without QED - see paper
        'aw0'     : gv.gvar(0,0),
        'a2DI'    : gv.gvar(0,0),
        'w0'      : gv.gvar(0.1714,0),

        'F0'      : gv.gvar(80,20), #FLAG use of F0 in SU(2) correction for FK/Fpi
        'FKFpi_FLAG' : gv.gvar(1.1932, 0.0019)
    },
    'x' : {'alphaS':0},
    'y' : {},
}

''' Values for checking fit function
'''
Fpi_check = 92.2
FK_check  = 110.
mpi_check = 135.0
mk_check  = 495.5
me_check  = np.sqrt(4./3*mk_check**2 - 1./3*mpi_check**2)
L_check   = 3.5 / mpi_check
check_fit = {
    'p':{
        'mpi'    : mpi_check,
        'mk'     : mk_check,
        'Fpi'    : Fpi_check,
        'FK'     : FK_check,
        'Lchi_PP': 4 * np.pi * Fpi_check,
        'Lchi_PK': 4 * np.pi * np.sqrt(Fpi_check * FK_check),
        'Lchi_KK': 4 * np.pi * FK_check,
        'L1'     :  0.000372,
        'L2'     :  0.000493,
        'L3'     : -0.003070,
        'L4'     :  0.000089,
        'L5'     :  0.000377,
        'L6'     :  0.000011,
        'L7'     : -0.000340,
        'L8'     :  0.000294,
        'k_4'    : -3.0,
        'p_4'    :  4.0,
        # discretization
        'aw0'    : 0.8,
        's_4'    : 1.5,
        'saS_4'  : 2.5,
        'kp_6'   : 2.1,
        'k_6'    : 2.2,
        'p_6'    : 2.3,
        's_6'    : 2.4,
        'sk_6'   : 2.5,
        'sp_6'   : 2.6,
        # mixed action
        'mss'    : 520.,
        'mju'    : 200.,
        'mjs'    : 510.,
        'mru'    : 505.,
        'mrs'    : 525.,
        'a2DI'   : 400.**2,
    },
    'x':{'alphaS':0.2, 'meL':me_check*L_check}
}
for mphi in ['mpi','mk', 'mju', 'mjs', 'mru', 'mrs', 'mss']:
    if mphi in check_fit['p']:
        check_fit['x'][mphi+'L'] = check_fit['p'][mphi] * L_check
check_fit['x']['mxL'] = L_check * np.sqrt(4./3*mk_check**2 - 1./3*mpi_check**2 + check_fit['p']['a2DI'])
