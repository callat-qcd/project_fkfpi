import gvar as gv
import numpy as np

switches = dict()
switches['ensembles'] = [
    'a15m400','a12m400','a09m400',
    'a15m350','a12m350','a09m350',
    'a15m310','a12m310','a09m310',
    'a15m220','a12m220S','a12m220','a12m220L','a09m220',
    'a12m130','a15m135XL']
switches['ensembles_fit'] = [
    'a15m400','a12m400','a09m400',
    'a15m350','a12m350','a09m350',
    'a15m310','a12m310','a09m310',
    'a15m220','a12m220S','a12m220','a12m220L','a09m220',
    'a12m130','a15m135XL']
switches['ansatz'] = dict()
switches['ansatz']['models'] = [
    'xpt_nnlo'                ,'ma_nnlo',
    'xpt_nnlo_FV'             ,'ma_nnlo_FV',
    'xpt_nnlo_FV_alphaS'      ,'ma_nnlo_FV_alphaS',
    'xpt_nnlo_FV_alphaS_logSq','ma_nnlo_FV_alphaS_logSq',
    'xpt_nnnlo',
    'xpt_nnnlo_FV'            ,'ma_nnnlo_FV',
    'xpt_nnnlo_FV_alphaS'     ,'ma_nnnlo_FV_alphaS',
    'xpt-ratio_nnlo'                ,'ma-ratio_nnlo',
    'xpt-ratio_nnlo_FV'             ,'ma-ratio_nnlo_FV',
    'xpt-ratio_nnlo_FV_alphaS'      ,'ma-ratio_nnlo_FV_alphaS',
    'xpt-ratio_nnlo_FV_alphaS_logSq','ma-ratio_nnlo_FV_alphaS_logSq'
]
#switches['ansatz']['models'] = ['xpt_nnlo_FV_alphaS_logSq','xpt_nnlo_FV_logSq']
switches['ansatz']['a2dm']  = 'individual' # avg or individual
switches['scale']           = 'PK' # PP, PK, KK, LamChi = 4 * pi * sqrt(FA * FB)
switches['do_analysis']     = False
switches['print_fit']       = False
switches['debug']           = False
switches['nlo_report']      = False # Do pure NLO fits for L5 test?
switches['nlo_fv_report']   = False # Do pure NLO fits for L5 with FV test?
switches['make_plots']      = False
switches['nnlo_priors']     = False
switches['prior_group']     = True
switches['verbose']         = True
switches['check_fit']       = True

priors = dict()
priors['L5']   = gv.gvar(0, 0.005)
#priors['L5']   = gv.gvar(0,0.005)
priors['L4']   = gv.gvar(0,0.005)

''' default values '''
switches['default_priors'] = True
priors['s_4']   = gv.gvar(0.0, 5)
priors['k_4']   = gv.gvar(0.0, 5)
priors['p_4']   = gv.gvar(0.0, 5)
priors['saS_4'] = gv.gvar(0.0, 5)

''' optimized values '''
nnlo_width = dict()
nnlo_width['xpt_nnlo_FV_alphaS'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS']['PK'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['p_4']   = 3.5
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['k_4']   = 3.5
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['s_4']   = 2.7
nnlo_width['xpt_nnlo_FV_alphaS']['PK']['saS_4'] = 2.7

nnlo_width['xpt_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 1.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 1.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 2.9
nnlo_width['xpt_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 2.9

nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 9.8
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 9.8
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 3.3
nnlo_width['xpt-ratio_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 3.3

nnlo_width['ma_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 0.9
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 0.9
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 0.8
nnlo_width['ma_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 0.8

nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq'] = dict()
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK'] = dict()
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['p_4']   = 0.9
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['k_4']   = 0.9
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['s_4']   = 0.1
nnlo_width['ma-ratio_nnlo_FV_alphaS_logSq']['PK']['saS_4'] = 0.1


n3lo_width = 5.

priors['kp_6']  = gv.gvar(0.0, n3lo_width)
priors['k_6']   = gv.gvar(0.0, n3lo_width)
priors['p_6']   = gv.gvar(0.0, n3lo_width)
priors['s_6']   = gv.gvar(0.0, n3lo_width)
priors['sk_6']  = gv.gvar(0.0, n3lo_width)
priors['sp_6']  = gv.gvar(0.0, n3lo_width)

phys_point = dict()
# http://pdg.lbl.gov/2016/tables/rpp2016-tab-mesons-light.pdf
phys_point['mpi'] = gv.gvar(134.8, 0.3) #FLAG 2017 (16)
phys_point['mk']  = gv.gvar(494.2 , 0.3) #FLAG 2017 (16)

#phys_point['mpi'] = gv.gvar(139,0.3)
#phys_point['mk']  = gv.gvar(497,0.3)

phys_point['F0']  = gv.gvar(80,20) #FLAG use of F0 in SU(2) correction for FK/Fpi
phys_point['meta'] = gv.gvar(547.862,0.017) #PDG
# http://pdg.lbl.gov/2015/reviews/rpp2015-rev-pseudoscalar-meson-decay-cons.pdf
phys_point['Fpi'] = gv.gvar(130.2/np.sqrt(2), 1.7/np.sqrt(2)) #PDG fpi+ eq(16)
phys_point['FK']  = gv.gvar(155.6/np.sqrt(2), 0.4/np.sqrt(2)) #PDG fK++ eq(16)
if switches['scale'] == 'PK':
    phys_point['Lchi'] = 4*np.pi*np.sqrt(phys_point['Fpi']*phys_point['FK'])
elif switches['scale'] == 'PP':
    phys_point['Lchi'] = 4*np.pi*phys_point['Fpi']
elif switches['scale'] == 'KK':
    phys_point['Lchi'] = 4*np.pi*phys_point['FK']
phys_point['aw0'] = 0
phys_point['FKFPi_FLAG'] = gv.gvar(1.1933, 0.0029)
