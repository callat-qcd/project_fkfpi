import gvar as gv
import numpy as np

switches = dict()
switches['ensembles'] = [
    'a15m400','a12m400','a09m400',
    'a15m350','a12m350','a09m350',
    'a15m310','a12m310','a09m310',
    'a15m220','a12m220S','a12m220','a12m220L','a09m220',
    'a12m130','a15m130','a15m135XL']
switches['ensembles_fit'] = [
    'a15m400','a12m400','a09m400',
    'a15m350','a12m350','a09m350',
    'a15m310','a12m310','a09m310',
    'a15m220','a12m220S','a12m220','a12m220L','a09m220',
    'a12m130','a15m135XL']
switches['ansatz'] = dict()
switches['ansatz']['model'] = 'xpt_nlo_FV' # Type of fit: 'xpt' or 'MA'
switches['ansatz']['a2dm']  = 'individual' # avg or individual
switches['scale']           = 'PP' # PP, PK, KK, LamChi = 4 * pi * sqrt(FA * FB)
switches['do_analysis']     = False
switches['debug']           = False
switches['nlo_report']      = False # Do pure NLO fits for L5 test?
switches['nlo_fv_report']   = True # Do pure NLO fits for L5 with FV test?

priors = dict()
priors['L5']   = gv.gvar(1e-3, 5.e-3)
#priors['L5']   = gv.gvar(0,0.005)
priors['L4']   = gv.gvar(0,0.005)
nnlo_wdith = 2.
n3lo_width = 1.
priors['s_4']   = gv.gvar(0.0, nnlo_wdith)
priors['k_4']   = gv.gvar(0.0, nnlo_wdith)
priors['p_4']   = gv.gvar(0.0, nnlo_wdith)
priors['saS_4'] = gv.gvar(0.0, nnlo_wdith)
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
