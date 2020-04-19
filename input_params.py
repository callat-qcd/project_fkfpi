import gvar as gv
import numpy as np
import prior_width_study as pw

switches = dict()
switches['ensembles'] = [
    'a15m400'  ,'a12m400' ,'a09m400',
    'a15m350'  ,'a12m350' ,'a09m350',
    'a15m310'  ,'a12m310' ,'a09m310','a06m310L',
    'a15m220'  ,'a12m220S','a12m220','a12m220L','a09m220','a09m135',
    'a15m135XL','a12m130',]
switches['ensembles_fit'] = [
    'a15m400'  ,'a12m400' ,'a09m400',
    'a15m350'  ,'a12m350' ,'a09m350',
    'a15m310'  ,'a12m310' ,'a09m310','a06m310L',
    'a15m220'  ,'a12m220S','a12m220','a12m220L','a09m220',
    'a15m135XL','a12m130','a09m135',
    ]
#switches['ensembles_fit'] = [
#    'a12m400' ,'a09m400',
#    'a12m350' ,'a09m350',
#    'a12m310' ,'a09m310','a06m310L',
#    'a12m220S','a12m220','a12m220L','a09m220',
#    'a12m130',
#    ]
#switches['ensembles_fit'] = ['a09m135']

# FIT MODELS
switches['ansatz'] = dict()
switches['ansatz']['models'] = ['xpt_nnlo','ma_nnlo']
switches['ansatz']['models'] = ['xpt_nnlo_FV_a4','ma_nnlo_FV_a4']#,'xpt_nnnlo_FV']
switches['ansatz']['models'] = ['xpt_nnlo_FV_a4', 'xpt_nnnlo_FV']#,'ma_nnlo', 'xpt_nnnlo', 'ma_nnnlo']#, 'xpt_nnnlo_FV']#,'ma_nnlo_FV_a4']#,'xpt_nnnlo_FV']
# SYSTEMATIC SWITCHES
switches['sys'] = dict()
switches['sys']['Lam_chi']   = False
switches['sys']['FV']        = False
switches['sys']['alphaS']    = False
switches['sys']['nnlo_ct']   = False
switches['sys']['logSq']     = False
switches['sys']['a4']        = False
switches['sys']['ratio']     = False

switches['ansatz']['a2dm']   = 'individual' # avg or individual
switches['scales']           = ['PP','PK','KK']#,'PK','KK']
switches['scale']            = 'PP' # PP, PK, KK, LamChi = 4 * pi * sqrt(FA * FB)
switches['do_analysis']      = True
switches['save_fits']        = False
switches['model_avg']        = True
# use optimized (True) or default (False) priors
switches['optimized_priors'] = False
# set mean to boot0 vs add boot0 to bs list
switches['bs_bias']          = True
# use scipy instead of GSL?
switches['scipy']            = False
# fit options
switches['print_fit']        = False
switches['report_fit']       = False
switches['make_plots']       = False
switches['save_figs']        = True
switches['plot_raw_data']    = False
switches['verbose']          = False
switches['milc_compare']     = False
switches['plot_asq_converg'] = True
# for tuning prior widths
switches['nnlo_priors']      = False
switches['prior_group']      = True
switches['refine_prior']     = True
# fit checks
switches['nlo_report']       = False
switches['nlo_fv_report']    = False
switches['simple_fit']       = False
#switches['check_fit']       = False
# DEBUGGING
switches['debug']            = False
switches['debug_bs']         = False
switches['debug_x']          = False
switches['debug_phys']       = False
switches['debug_nlo_check']  = False
switches['debug_nnlo_check'] = False
switches['debug_save_fit']   = False # also need 'save_fits' to work
switches['debug_models']     = False
switches['print_lattice']    = False # this will turn off all fitting - only reads data
# testing NNLO function
switches['check_fit']        = False # need a new name
switches['FF_approximate']   = False

switches['LECs'] = ['L1','L2','L3','L4','L5','L6','L7','L8','p_4','k_4','s_4','saS_4']

priors = dict()
priors['L5']   = gv.gvar(0, 0.04)
priors['L4']   = gv.gvar(0, 0.005)

priors['L1']   = gv.gvar( 0.000372,0.000372)
priors['L2']   = gv.gvar( 0.000493,0.000493)
priors['L3']   = gv.gvar(-0.003070,0.003070)
priors['L6']   = gv.gvar( 0.000011,0.00011)
priors['L7']   = gv.gvar(-0.000340,0.000340)
priors['L8']   = gv.gvar( 0.000294,0.000294)

# Taylor priors
priors['c2'] = gv.gvar(0,10)

''' default values '''
n2lo_width = 5.
#priors['s_4']   = gv.gvar(0.0, n2lo_width)
priors['s_4']   = gv.gvar(0, 10)
#priors['s_4']   = gv.gvar(0, 100)
priors['k_4']   = gv.gvar(0.0, n2lo_width)
priors['p_4']   = gv.gvar(0.0, n2lo_width)
priors['saS_4'] = gv.gvar(0.0, 10)
#priors['saS_4'] = gv.gvar(0.0, n2lo_width)

# Ananthanarayan et al, 1711.11328
priors['k_4']   = gv.gvar(2.2, 4)
priors['p_4']   = gv.gvar(7.9, 10)

'''
n3lo_width = 5.
priors['kp_6']  = gv.gvar(0.0, n3lo_width)
priors['k_6']   = gv.gvar(0.0, n3lo_width)
priors['p_6']   = gv.gvar(0.0, n3lo_width)
priors['s_6']   = gv.gvar(100, 100)
priors['sk_6']  = gv.gvar(0.0, 100)
priors['sp_6']  = gv.gvar(0.0, 100)
'''
n3lo_width = 100.
priors['kp_6']  = gv.gvar(0.0, n3lo_width)
priors['k_6']   = gv.gvar(0.0, n3lo_width)
priors['p_6']   = gv.gvar(0.0, n3lo_width)
priors['s_6']   = gv.gvar(0.0, n3lo_width)
priors['sk_6']  = gv.gvar(0.0, n3lo_width)
priors['sp_6']  = gv.gvar(0.0, n3lo_width)



''' load optimized priors  '''
switches['nnlo_priors_model'] = pw.model
priors['p_range'] = pw.p_range
priors['a_range'] = pw.a_range
nnlo_width = pw.nnlo_width

phys_point = dict()
# PDG valuese
# http://pdg.lbl.gov/2016/tables/rpp2016-tab-mesons-light.pdf
#phys_point['mpi'] = gv.gvar(139,0.3)
#phys_point['mk']  = gv.gvar(497,0.3)
#phys_point['meta'] = gv.gvar(547.862,0.017) #PDG
# http://pdg.lbl.gov/2015/reviews/rpp2015-rev-pseudoscalar-meson-decay-cons.pdf
FPi_phys = gv.gvar(130.2/np.sqrt(2), 0.8/np.sqrt(2))
FK_phys  = gv.gvar(155.5/np.sqrt(2), 0.7/np.sqrt(2))

phys_point = {
    'p':{
        'Fpi'     : FPi_phys, #PDG fpi+ eq(16)
        'FK'      : FK_phys, #PDG fK+ eq(16)
        'Lchi_PP' : 4 * np.pi * FPi_phys,
        'Lchi_PK' : 4 * np.pi * np.sqrt(FPi_phys * FK_phys),
        'Lchi_KK' : 4 * np.pi * FK_phys,
        'mpi'     : gv.gvar(134.8, 0.3), #FLAG 2017 (16)
        #'mk'      : gv.gvar(494.2, 0.3), #FLAG 2017 (16) isospin symmetric
        'mk'      : gv.gvar(491.5, 0.7), # my estimate of MK+ without QED

        'aw0'     : gv.gvar(0,0),
        'a2DI'    : gv.gvar(0,0),
        'w0'      : gv.gvar(0.1714,0),

        'F0'      : gv.gvar(80,20), #FLAG use of F0 in SU(2) correction for FK/Fpi
        'FKFpi_FLAG' : gv.gvar(1.1932, 0.0019)
    },
    'x' : {'alphaS':0},
    'y' : {},
}


check_fit = dict()
check_fit['mpi'] =  135.0
check_fit['mk']  =  495.5
check_fit['Fpi'] =  92.2
check_fit['FK']  =  110.5
check_fit['L1']  =  0.000372
check_fit['L2']  =  0.000493
check_fit['L3']  = -0.003070
check_fit['L4']  =  0.000089
check_fit['L5']  =  0.000377
check_fit['L6']  =  0.000011
check_fit['L7']  = -0.000340
check_fit['L8']  =  0.000294
check_fit['k_4'] = -3.0
check_fit['p_4'] =  4.0
