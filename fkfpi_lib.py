from __future__ import print_function
import os, sys, shutil
import numpy as np
import scipy.stats as stats
import tables as h5
import pandas as pd
import pickle
import lsqfit
import gvar as gv
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_style("ticks")
import yaml
import chipt as xpt

def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

def main():
    import input_params as ip
    print("python     version:", sys.version)
    print("numpy      version:", np.__version__)
    print("pandas     version:", pd.__version__)
    print("matplotlib version:", matplotlib.__version__)
    print("gvar       version:", gv.__version__)
    print("lsqfit     version:", lsqfit.__version__)
    print("yaml       version:", yaml.__version__)
    print('')

    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.ion()

    # Load input params
    switches   = ip.switches
    priors     = ip.priors
    phys_point = ip.phys_point
    check_fit  = ip.check_fit

    # check fit
    if switches['check_fit']:
        check_fit_function(switches,check_fit)
        sys.exit()

    # Load data
    data    = h5.open_file('data/FK_Fpi_data.h5','r')
    gv_data = format_h5_data(switches,data)
    data.close()

    if switches['simple_fit']:
        import simple_fit as sf
        d_e = dict()
        d_e['x'] = {k:gv_data['x'][k] for k in switches['ensembles']}
        d_e['y'] = {k:gv_data['y'][k] for k in switches['ensembles']}
        d_e['p'] = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
        for f_model in ['simple','simple_su3']:
            print('\n=================================================================')
            print(f_model)
            print('=================================================================')
            switches['ansatz']['model'] = f_model
            fit = sf.SimpleFit(switches,xyp=d_e)
            fit.fit_data()
            print(fit.fit.format(maxline=True))

            ea_range = dict()
            ea_range['Lchi'] = phys_point['Lchi_'+switches['scale']]
            ea_range['mpi']  = phys_point['mpi']
            ea_range['mk']   = phys_point['mk']
            ea_range['a']    = np.sqrt(np.arange(0,.16**2,.16**2/50))
            ea_range['w0']   = phys_point['w0']
            fig_vs_ea = plt.figure('FKFpi_vs_ea_'+f_model)
            ax = plt.axes([.12, .12, .85, .85])
            ax = fit.plot_fit(ea_range,ax)
            ax = fit.plot_data(ax)

            ax.legend()
            ax.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=16)
            ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
            ax.set_xlim(0,.065)
            ax.set_ylim(1.09, 1.13)

            plt.savefig('figures/FKFpi_vs_ea_'+f_model+'.pdf',transparent=True)


    if switches['do_analysis']:
        fit_results = perform_analysis(switches,gv_data,priors,phys_point)
        if switches['model_avg']:
            model_average(fit_results,switches,phys_point)
        if switches['make_plots']:
            for model in switches['ansatz']['models']:
                if model in ['xpt_nnlo_FV_a4_PK','xpt_nnlo_FV_a4_PP','xpt_nnlo_FV_PP']:
                    fit_results[model].vs_ea()

    if switches['check_fit']:
        fit_checker(switches,priors)

    if switches['nnlo_priors']:
        nnlo_prior_scan(switches,priors)

    if switches['nlo_fv_report']:
        nlo_report(switches,priors,FV=True)
    if switches['nlo_report']:
        nlo_report(switches,priors,FV=False)

    plt.ioff()
    if run_from_ipython():
        plt.show(block=False)
    else:
        plt.show()

ens_long = {
    'a15m400'  :'l1648f211b580m0217m065m838',
    'a15m350'  :'l1648f211b580m0166m065m838',
    'a15m310'  :'l1648f211b580m013m065m838',
    'a15m220'  :'l2448f211b580m0064m0640m828',
    'a15m130'  :'l3248f211b580m00235m0647m831',
    'a15m135XL':'l4864f211b580m002426m06730m8447',
    'a12m400'  :'l2464f211b600m0170m0509m635',
    'a12m350'  :'l2464f211b600m0130m0509m635',
    'a12m310'  :'l2464f211b600m0102m0509m635',
    'a12m220'  :'l3264f211b600m00507m0507m628',
    'a12m220S' :'l2464f211b600m00507m0507m628',
    'a12m220L' :'l4064f211b600m00507m0507m628',
    'a12m130'  :'l4864f211b600m00184m0507m628',
    'a09m400'  :'l3264f211b630m0124m037m440',
    'a09m350'  :'l3264f211b630m00945m037m440',
    'a09m310'  :'l3296f211b630m0074m037m440',
    'a09m220'  :'l4896f211b630m00363m0363m430',
}
r_a = {
    'a15m400'  :gv.gvar(2.059,0.023),
    'a15m350'  :gv.gvar(2.059,0.023),
    'a15m310'  :gv.gvar(2.059,0.023),
    'a15m220'  :gv.gvar(2.073,0.013),
    'a15m130'  :gv.gvar(2.089,0.008),
    'a15m135XL':gv.gvar(2.089,0.008),
    'a12m400'  :gv.gvar(2.575,0.017),
    'a12m350'  :gv.gvar(2.575,0.017),
    'a12m310'  :gv.gvar(2.575,0.017),
    'a12m220S' :gv.gvar(2.585,0.019),
    'a12m220'  :gv.gvar(2.626,0.013),
    'a12m220L' :gv.gvar(2.614,0.009),
    'a12m130'  :gv.gvar(2.608,0.008),
    'a09m400'  :gv.gvar(3.499,0.024),
    'a09m350'  :gv.gvar(3.499,0.024),
    'a09m310'  :gv.gvar(3.499,0.024),
    'a09m220'  :gv.gvar(3.566,0.014)
    }

a_ens = {
    'a15':gv.gvar(0.1509,0.0013),
    'a12':gv.gvar(0.1206,0.0010),
    'a09':gv.gvar(0.0875,0.0008),
    'a06':gv.gvar(0.0566,0.0005),
}

def format_h5_data(switches,data):
    x = dict()
    y = dict()
    p = dict()
    if switches['bs_bias']:
        print('Shifting BS data to boot0')
    if switches['print_lattice']:
        lattice_fits = []
        mixed_fits   = []

    print('%9s FK/Fpi' %'ensemble')
    print('-----------------------------------------------------------------')
    for ens in switches['ensembles']:
        x[ens] = dict()
        data_dict = dict()
        for m in ['mpi','mk','mss','mju','mjs','mrs','mru']:
            data_dict[m] = data.get_node('/'+ens+'/'+m).read()
        for f in ['FK','Fpi']:
            data_dict[f] = data.get_node('/'+ens+'/'+f).read()
        for m in ['mres_l','mres_s']:
            data_dict[m] = data.get_node('/'+ens+'/'+m).read()
        if switches['debug']:
            fkfpi = data_dict['FK']/data_dict['Fpi']
            print('  FK/Fpi = %.5f +- %.5f' %(fkfpi.mean(),fkfpi.std()))
        # HACK
        if True:
            if ens == 'a06m310L':
                for f in ['FK','Fpi']:
                    ff = data_dict[f]
                    df = ff - ff.mean()
                    df = df / 2
                    ff = ff.mean() + df
                    data_dict[f] = ff

        if switches['bs_bias']:
            data_bs = dict()
            for d in data_dict:
                data_bs[d] = data_dict[d][1:]
                if switches['debug']:
                    print(d,data_bs[d].shape)
            gvdata = gv.dataset.avg_data(data_bs,bstrap=True)
            if switches['debug']:
                gvdata_copy = dict(gvdata)
            for d in data_bs:
                gvdata[d] = gvdata[d] + (data_dict[d][0] - gvdata[d].mean)
                if switches['debug']:
                    print(d,gvdata[d].mean,gvdata_copy[d].mean)
        else:
            gvdata = gv.dataset.avg_data(data_dict,bstrap=True)
        if switches['debug']:
            print(gvdata)
        if switches['debug']:
            print('data_dict')
            print(data_dict)

        y[ens] = gvdata['FK']/gvdata['Fpi']
        print("%9s %s" %(ens,y[ens]))

        x[ens]['mpiL'] = gvdata['mpi'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mkL']  = gvdata['mk'].mean  * data.get_node('/'+ens+'/L').read()
        x[ens]['mssL'] = gvdata['mss'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mjuL'] = gvdata['mju'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mjsL'] = gvdata['mjs'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mruL'] = gvdata['mru'].mean * data.get_node('/'+ens+'/L').read()
        x[ens]['mrsL'] = gvdata['mrs'].mean * data.get_node('/'+ens+'/L').read()

        # MASSES
        p[(ens,'mpi')] = gvdata['mpi']
        p[(ens,'mk')]  = gvdata['mk']
        p[(ens,'mss')] = gvdata['mss']
        p[(ens,'mju')] = gvdata['mju']
        p[(ens,'mjs')] = gvdata['mjs']
        p[(ens,'mru')] = gvdata['mru']
        p[(ens,'mrs')] = gvdata['mrs']
        p[(ens,'Lchi_PP')] = 4 * np.pi * gvdata['Fpi']
        p[(ens,'Lchi_PK')] = 4 * np.pi * np.sqrt(gvdata['FK'] * gvdata['Fpi'])
        p[(ens,'Lchi_KK')] = 4 * np.pi * gvdata['FK']
        p[(ens,'Lchi_F0')] = 4 * np.pi * gv.gvar(80,1) / 197.3 * a_ens[ens[0:3]]
        # HISQ params
        aw0 = data.get_node('/'+ens+'/aw0').read()
        p[(ens,'aw0')] = gv.gvar(aw0[0],aw0[1])
        a2di = data.get_node('/'+ens+'/a2DI').read()
        p[(ens,'a2DI')] = gv.gvar(a2di[0],a2di[1])
        x[ens]['alphaS'] = data.get_node('/'+ens+'/alpha_s').read()

        if switches['debug_x']:
            eps_pi = gvdata['mpi'] / 4 / np.pi / gvdata['Fpi']
            eps_k  = gvdata['mk']  / 4 / np.pi / gvdata['Fpi']
            eps_a  = p[(ens,'aw0')] / np.sqrt(4*np.pi)
            fkfpi  = gvdata['FK']/gvdata['Fpi']
            print('%9s: eps_pi = %s, eps_k = %s, eps_a = %s, FK/Fpi = %s' %(ens, eps_pi, eps_k, eps_a, fkfpi))
        if switches['print_lattice']:
            lattice_fits.append('%9s& %s& %s& %s& %s& %s& %s& %s\\\\' \
                %(ens,gvdata['mres_l'],gvdata['mres_s'],gvdata['mpi'],gvdata['mk'],gvdata['Fpi'],gvdata['FK'],gvdata['FK']/gvdata['Fpi']))
            dju = p[(ens,'aw0')]**(-2) * (gvdata['mju']**2 - gvdata['mpi']**2)
            djs = p[(ens,'aw0')]**(-2) * (gvdata['mjs']**2 - gvdata['mk']**2)
            dru = p[(ens,'aw0')]**(-2) * (gvdata['mru']**2 - gvdata['mk']**2)
            drs = p[(ens,'aw0')]**(-2) * (gvdata['mrs']**2 - gvdata['mss']**2)
            mixed_fits.append('%9s& %s& %s& %s& %s& %s& %s& %s& %s\\\\' \
                %(ens,gvdata['mju'],gvdata['mjs'],gvdata['mru'],gvdata['mrs'],dju,djs,dru,drs))
    if switches['print_lattice']:
        print(r'ensemble& $am^{\rm res}_l$& $am^{\rm res}_s$& $am_\pi$& $am_K$& $aF_\pi$& $aF_K$& $F_K / F_\pi$\\')
        for l in lattice_fits:
            print(l)
        print('')
        print(r'ensemble& $am_{ju}$& $am_{js}$& $am_{ru}$& $am_{rs}$& $w_0^2(m_{ju}^2 -m_\pi^2)$& $w_0^2(m_{js}^2-m_K^2)$& $w_0^2(m_{ru}^2-m_K^2)$& $w_0^2(m_{rs}^2 -m_{ss}^2)$\\')
        print(r'\hline')
        for l in mixed_fits:
            print(l)
        data.close()
        sys.exit()
    return {'x':x, 'y':y, 'p':p}

def set_priors(p_init,priors,optimized_priors=None):
    for k in priors:
        p_init[k] = priors[k]
    if optimized_priors != None:
        df = pd.DataFrame(optimized_priors)
        sp_max = df.stack().idxmax()
        print('      setting prior widths:')
        print('      (s_4, saS_4) = %s; (p_4, k_4) = %s' %(sp_max[0],sp_max[1]))
        for key in ['s_4','saS_4']:
            if key in p_init:
                p_init[key] = gv.gvar(0,float(sp_max[0]))
        for key in ['p_4','k_4']:
            if key in p_init:
                p_init[key] = gv.gvar(0,float(sp_max[1]))
    else:
        print('      using default prior widths')
    return p_init

def pickle_fit(model,fit_result,phys_point):
    fit = gv.BufferDict()
    fit['phys_extrap'] = fit_result.report_phys_point()['phys']
    for k in fit_result.fit.prior:
        if type(k) is tuple:
            fit['prior_'+k[0]+'_'+k[1]] = fit_result.fit.prior[k]
        else:
            fit['prior_'+k] = fit_result.fit.prior[k]
    for k in fit_result.fit.p:
        if type(k) is tuple:
            fit['p_'+k[0]+'_'+k[1]] = fit_result.fit.p[k]
        else:
            fit['p_'+k] = fit_result.fit.p[k]
    # we can't save the fit.cov easily, as it is an np.array
    # but, they can be reconstructed after reading in all the fit.p values
    fit['chi2']   = gv.gvar(fit_result.fit.chi2)
    fit['dof']    = gv.gvar(fit_result.fit.dof)
    fit['logGBF'] = gv.gvar(fit_result.fit.logGBF)
    fit['Q']      = gv.gvar(fit_result.fit.Q)
    for k in fit_result.fit.y:
        fit['data_'+k] = fit_result.fit.y[k]

    if not os.path.exists('pickled_fits'):
        os.makedirs('pickled_fits')
    gv.dump(fit,'pickled_fits/'+model+'.p', add_dependencies=True)

def read_fit(model):
    fit = gv.load('pickled_fits/'+model+'.p')
    fit_processed = dict()
    fit_processed['phys_extrap'] = fit['phys_extrap']
    fit_processed['prior'] = dict()
    fit_processed['p']     = dict()
    fit_processed['stat']  = dict()
    fit_processed['data']  = dict()
    for k in fit:
        if 'data' in k:
            fit_processed['data'][k.split('_')[1]] = fit[k]
        if k in ['chi2','dof','logGBF','Q']:
            fit_processed['stat'][k] = fit[k]
        if 'p_' in k and 'prior' not in k:
            if any(part in k for part in ['a06','a09','a12','a15']):
                tmp,ens,kk = k.split('_')
                fit_processed['p'][(ens,kk)] = fit[k]
            elif len(k.split('p_')) == 3:
                fit_processed['p']['p_'+k.split('p_')[2]] = fit[k]
            elif len(k.split('p_')) == 2:
                fit_processed['p'][k.split('p_')[1]] = fit[k]
            else:
                print('what to do?',k)
        if 'prior' in k:
            if any(part in k for part in ['a06','a09','a12','a15']):
                tmp,ens,kk = k.split('_')
                fit_processed['prior'][(ens,kk)] = fit[k]
            else:
                fit_processed['prior'][k.split('prior_')[1]] = fit[k]

    return fit_processed

class PickledFit():
    def __init__(self,pickled_fit):
        self.p      = pickled_fit['p']
        self.prior  = pickled_fit['prior']
        self.y      = pickled_fit['data']
        self.dof    = pickled_fit['stat']['dof'].mean
        self.chi2   = pickled_fit['stat']['chi2'].mean
        self.logGBF = pickled_fit['stat']['logGBF'].mean
        self.Q      = pickled_fit['stat']['Q'].mean
        self.cov    = gv.evalcov(self.p)
        self.corr   = gv.evalcorr(self.p)

def fkfpi_phys(x_phys,fit):
    # use metaSq = 4/3 mK**2 - 1/3 mpi**2
    print('prediction from LQCD')
    #meta = np.sqrt(4./3 * x_phys['mka']**2 -1./3 * x_phys['mpi']**2)
    #fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
    #    meta=meta,Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    Fitc = xpt.Fit(switches={'ansatz':{'type':'xpt'}},phys=phys_point)
    for k in fit.p:
        if k in x_phys:
            pass
        else:
            x_phys[k] = fit.p[k]
    fkp = Fitc.fit_function({'Lchi':x_phys['Lchi']},x_phys)
    su2 = dsu2(FKpi=fkp,mpi=x_phys['mpi'],mk=x_phys['mka'],F0=x_phys['F0'])
    print('FK / Fpi = ',fkp)
    print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))
    # use meta = meta_pdg
    #print('prediction from LQCD + meta_PDG')
    #fkp = fkfpi_nlo(mpi=x_phys['mpi'],mka=x_phys['mka'],\
    #    meta=x_phys['meta'],Lchi=x_phys['Lchi'],L5=fit.p['L5'])
    #print('FK / Fpi = ',fkp)
    #print('FK+/Fpi+ = ',fkp*np.sqrt(1+su2))
    return fkp*np.sqrt(1+su2)

def dsu2(FKpi,mpi,mk,F0):
    R = gv.gvar(35.7,np.sqrt(1.9**2 + 1.8**2))
    eps_su2 = np.sqrt(3) / 4 / R
    d = np.sqrt(3) * eps_su2 * (\
        -4./3 * (FKpi-1) \
        + 2. * (mk**2 -mpi**2 -mpi**2*np.log(mk**2/mpi**2)) / (6*(4*np.pi)**2*F0**2)
        )
    return d

def checkIfDuplicates_1(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def sys_models(switches):
    def check_model(sys_val,models,nnlo=False,nnnlo=False):
        for model in models:
            new_model = model+sys_val
            if nnlo:
                if (sys_val not in model) and (new_model not in models):
                    if nnnlo:
                        if 'nnlo' in model and 'nnnlo' not in model:
                            models.append(new_model)
                    else:
                        if 'nnlo' in model:
                            models.append(new_model)
            else:
                if (sys_val not in model) and (new_model not in models):
                    models.append(new_model)
        #return models
    models = switches['ansatz']['models'].copy()
    #print(len(models),'models')
    if switches['sys']['FV']:
        check_model('_FV',models)
    if switches['sys']['alphaS']:
        check_model('_alphaS',models,nnlo=True)
    if switches['sys']['nnlo_ct']:
        check_model('_ct',models,nnlo=True)
    if switches['sys']['logSq']:
        check_model('_logSq',models,nnlo=True)
    if switches['sys']['a4']:
        check_model('_a4',models,nnlo=True,nnnlo=True)
    if switches['sys']['ratio']:
        for model in models:
            model_ratio = model.replace('xpt','xpt-ratio').replace('ma','ma-ratio')
            if '-ratio' not in model and model_ratio not in models:
                models.append(model_ratio)
    models_FPK = []
    for model in models:
        if switches['sys']['Lam_chi']:
            for FPK in switches['scales']:
                models_FPK.append(model+'_'+FPK)
        else:
            models_FPK.append(model+'_'+switches['scale'])
    if switches['debug_models']:
        for model in models_FPK:
            print(model)
    print(len(models_FPK),'models')
    print('Duplicate models?',checkIfDuplicates_1(models_FPK))
    return models_FPK

def perform_analysis(switches,gv_data,priors,phys_point):
    fit_results = dict()
    print('\nSetting up all models')
    models = sys_models(switches)
    switches['ansatz']['models'] = models
    for i_m,model in enumerate(models):
        FPK = model.split('_')[-1]
        switches['scale'] = FPK
        switches['ansatz']['model'] = model

        if switches['optimized_priors']:
            if not os.path.exists('data/saved_prior_search.yaml'):
                print('ERROR: you asked for optimized priors but')
                print('data/saved_prior_search.yaml')
                print('does not exist - exiting')
                sys.exit()
            with open('data/saved_prior_search.yaml','r') as fin:
                prior_grid = yaml.safe_load(fin.read())

        print('%3s of %3s, EFT: %s' %(i_m+1,len(models),model))
        x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
        y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
        p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in
        switches['ensembles']}
        if switches['optimized_priors']:
            p_e = set_priors(p_e,priors,prior_grid[model])
        else:
            p_e = set_priors(p_e,priors)

        do_fit = False
        if switches['optimized_priors']:
            p_fit = model+'_optimized_priors'
        else:
            p_fit = model
        if switches['save_fits']:
            if os.path.exists('pickled_fits/'+p_fit+'.p'):
                print('reading pickled_fits/'+p_fit+'.p')
                pickled_fit = read_fit(p_fit)
                fit_p = xpt.Fit(switches,xyp_init={'x':x_e,'y':y_e,'p':p_e},phys=phys_point)
                fit_p.fit = PickledFit(pickled_fit)
                fit_p.report_phys_point()
                fit_results[model] = fit_p

            else:
                do_fit = True
        else:
            do_fit = True
        if do_fit or switches['debug_save_fit']:
            '''
            for key in priors:
                p_e[key] = priors[key]
            if switches['optimized_priors']:
                df = pd.DataFrame(prior_grid[model])
                sp_max = df.stack().idxmax()
                print('      setting prior widths:')
                print('      (s_4, saS_4) = %s; (p_4, k_4) = %s' %(sp_max[0],sp_max[1]))
                for key in ['s_4','saS_4']:
                    if key in p_e:
                        p_e[key] = gv.gvar(0,float(sp_max[0]))
                for key in ['p_4','k_4']:
                    if key in p_e:
                        p_e[key] = gv.gvar(0,float(sp_max[1]))
            else:
                print('      using default prior widths')
            '''
            d_e = dict()
            d_e['x'] = x_e
            d_e['y'] = y_e
            d_e['p'] = p_e
            fit_e = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
            fit_e.fit_data()
            fit_results[model] = fit_e
            if switches['save_fits']:
                pickle_fit(p_fit,fit_e,phys_point)
            else:
                fit_p = fit_e

        if switches['debug_save_fit'] and not do_fit:
            print('=================================================================')
            print('Comparing Live fit to Pickled Fit')
            print('=================================================================')
            print('live fit   : FK/Fpi             = ',fit_e.report_phys_point()['phys'])
            print('error breakdown: y, L_i')
            print(fit_e.report_phys_point()['phys'].partialsdev(fit_e.fit.y))
            fit_priors = dict()
            for k in fit_e.lec_nnlo + fit_e.lec_nnnlo:
                if k in fit_e.fit.prior:
                    fit_priors[k] = fit_e.fit.prior[k]
            print(fit_e.report_phys_point()['phys'].partialsdev(fit_priors))

            print('pickled fit: FK/Fpi.report_phys = ',fit_p.report_phys_point()['phys'])
            print('pickled fit: FK/Fpi.saved_fit   = ',fit_p.phys_extrap)
            print('error breakdown: y, L_i')
            print(fit_p.report_phys_point()['phys'].partialsdev(fit_p.fit.y))
            fit_priors = dict()
            for k in fit_p.lec_nnlo + fit_p.lec_nnnlo:
                if k in fit_p.fit.prior:
                    fit_priors[k] = fit_p.fit.prior[k]
            print(fit_p.report_phys_point()['phys'].partialsdev(fit_priors))
            print('=================================================================')

        if not do_fit:
            fit_e = fit_p
        print('FK/Fpi = ',fit_e.report_phys_point()['phys'])
        print('do fit',do_fit)
        if switches['print_fit'] and do_fit:
            print(fit_e.fit.format(maxline=True))
        if switches['debug_phys']:
            fit_e.report_phys_point()['phys']

    return fit_results

def plot_continuum(fit_results,switches,phys_point):
    for model in switches['ansatz']['models']:
        FPK = model.split('_')[-1]
        base_model = model.strip('_'+FPK)
        switches['scale'] = FPK
        switches['ansatz']['model'] = model
        if model in ['xpt_nnnlo_FV_a4_PK','xpt_nnlo_FV_a4_PK']:
            fit = fit_results[model]
            ea_range = dict()
            ea_range['Lchi'] = phys_point['Lchi_'+FPK]
            ea_range['mpi']  = phys_point['mpi']
            ea_range['mk']   = phys_point['mk']
            ea_range['a']    = np.sqrt(np.arange(0,.16**2,.16**2/50))
            ea_range['w0']   = phys_point['w0']

            fig_vs_ea = plt.figure('FKFpi_vs_ea_'+model)
            ax = plt.axes([.12, .12, .85, .85])
            #print(fit_e.fit.format(maxline=True))
            ax = fit.vs_ea(ea_range, ax)
            # add data
            y_shift = fit.shift_phys_mass(phys_point)
            for e in switches['ensembles_fit']:
                x = fit.fit.p[(e,'a2')]
                y = fit.y[e] + y_shift[e]
                #print(e,y_e[e],y_shift[e])
                ax = plot_data(ax, e, x, y, offset='cont')
                if switches['plot_raw_data']:
                    ax = plot_data(ax, e, x, fit.y[e], raw=True)
            handles, labels = ax.get_legend_handles_labels()
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
            ax.legend(handles, labels, ncol=4)
            ax.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=16)
            ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
            ax.set_xlim(0,.065)
            if switches['plot_raw_data']:
                ax.set_ylim(1.09, 1.225)
            else:
                ax.set_ylim(1.135, 1.225)
            plt.savefig('figures/vs_epasq_'+model+'.pdf',transparent=True)


def model_average(fit_results,switches,phys_point):
    model_avg,model_var = bayes_model_avg(switches,fit_results,phys_point)
    print('SU(2) correction',dsu2(model_avg,phys_point['mpi'],phys_point['mk'],phys_point['F0']))

''' this is make plots and bma - needs own functions
            if switches['make_plots']:
                if base_model in ['xpt_nnlo_FV','xpt_nnlo_FV_a4']:#, 'xpt_nnlo_FV_logSq', 'xpt-ratio_nnlo_FV', 'xpt-ratio_nnlo_FV_logSq']:
                    ea_range = dict()
                    ea_range['Lchi'] = phys_point['Lchi_'+FPK]
                    ea_range['mpi']  = phys_point['mpi']
                    ea_range['mk']   = phys_point['mk']
                    ea_range['a']    = np.sqrt(np.arange(0,.16**2,.16**2/100))
                    ea_range['w0']   = phys_point['w0']

                    fig_vs_ea = plt.figure('FKFpi_vs_ea_'+model)
                    ax = plt.axes([.12, .12, .85, .85])
                    #print(fit_e.fit.format(maxline=True))
                    ax = fit_e.vs_ea(ea_range, ax)
                    # add data
                    y_shift = fit_e.shift_phys_mass(phys_point)
                    for e in switches['ensembles_fit']:
                        x = fit_e.fit.p[(e,'a2')]
                        y = y_e[e] + y_shift[e]
                        #print(e,y_e[e],y_shift[e])
                        ax = plot_data(ax, e, x, y)
                    ax.legend(ncol=3)
                    ax.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=16)
                    ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
                    ax.set_xlim(0,.065)
                    ax.set_ylim(1.135, 1.205)
                    plt.savefig('figures/vs_epasq_'+model+'.pdf',transparent=True)
                if False:
                    epi_range = dict()
                    epi_range['Lchi'] = phys_point['Lchi']
                    epi_range['mk']   = phys_point['mk']
                    epi_range['mpi']  = np.sqrt(np.arange(1,400**2,400**2 // 100))
                    epi_range['mpi_phys'] = phys_point['mpi']

                    fig_vs_epi = plt.figure('FKFpi_vs_epi_'+model)
                    ax = plt.axes([.12, .12, .85, .85])
                    ax = fit_e.vs_epi(epi_range, ax)
                    y_shift = fit_e.shift_phys_epsK(phys_point)
                    for e in switches['ensembles_fit']:
                        x = fit_e.fit.p[(e,'p2')]
                        y = y_e[e] + y_shift[e]
                        ax = plot_data(ax, e, x, y)
                    ax.legend()
                    ax.set_xlabel(r'$\epsilon_\pi^2$',fontsize=16)
                    ax.set_ylabel(r'$F_K / F_\pi$',fontsize=16)
                    ax.set_xlim(0,.1)
                    ax.set_ylim(1.05, 1.25)
                    plt.savefig('figures/vs_episq_'+model+'.pdf',transparent=True)
'''

def bayes_model_avg(switches,results,phys_point):
    logGBF_list = []
    r_list      = []
    w_list      = []
    models      = []
    FKFpi       = 0.
    pdf         = 0
    cdf         = 0
    pdf_pp      = 0
    pdf_pk      = 0
    pdf_kk      = 0
    x = np.arange(1.15,1.2301,.0001)
    for model in results:
        logGBF_list.append(results[model].fit.logGBF)
        FPK=model.split('_')[-1]
        results[model].switches['scale'] = FPK
        r_list.append(results[model].report_phys_point()['phys'])
        if switches['debug']:
            print('\nDEBUG:',model,results[model].switches['scale'])
    print('\nAll Models')
    print("%33s %5s %8s %8s %s" %('model','Q','logGBF','w','FK/Fpi'))
    print('-------------------------------------------------------------------')
    for i_m,model in enumerate(results):
        models.append(model)
        r = r_list[i_m]
        w = np.exp(logGBF_list[i_m])
        w = w / np.sum(np.exp(np.array(logGBF_list)))
        w_list.append(w)
        print("%33s %.3f %.3f %.2e %s" %(model,results[model].fit.Q,results[model].fit.logGBF,w,r))
        FKFpi += gv.gvar(w*r.mean, np.sqrt(w)*r.sdev)
        p = stats.norm.pdf(x,r.mean,r.sdev)
        pdf += w * p
        if '_PP' in model:
            pdf_pp += w * p
        elif '_PK' in model:
            pdf_pk += w * p
        elif '_KK' in model:
            pdf_kk += w * p
        c = stats.norm.cdf(x,r.mean,r.sdev)
        cdf += w * c
    w_list = np.array(w_list)
    r_list = np.array(r_list)
    model_var = np.sum(w_list * np.array([ r.mean**2 for r in r_list])) - FKFpi.mean**2

    print('\nModels with Q > 0.05')
    print("%33s %5s %6s %5s %s" %('model','Q','logGBF','w','FK/Fpi'))
    print('-------------------------------------------------------------------')
    for i_m,model in enumerate(models):
        if results[model].fit.Q > 0.05:
            print("%33s %.3f %.3f %.3f %s %f %f" %(model,results[model].fit.Q,results[model].fit.logGBF,w_list[i_m],r_list[i_m],r_list[i_m].mean,r_list[i_m].sdev))
        else:
            if w_list[i_m] > 0.05:
                print("%33s %.2f %.3f %s" %(model,results[model].fit.Q,w_list[i_m],'NEGLECTED'))

    print('\nFull average with %d Models' %len(results))
    print("%s +- %.4f" %(FKFpi,np.sqrt(model_var)))

    fig = plt.figure('hist')
    ax = plt.axes([0.11,0.12,0.87,0.87])
    x0 = FKFpi.mean
    dx = np.sqrt(model_var + FKFpi.sdev**2)
    #ax.axvspan(x0-dx,x0+dx,color='k',alpha=0.05)
    #ax.axvline(x0-dx,color='k',linestyle='-.')
    #ax.axvline(x0+dx,color='k',linestyle='-.')
    #dx = FKFpi.sdev
    #ax.axvspan(x0-dx,x0+dx,color='k',alpha=0.2)
    #ax.axvline(x0,color='k')
    #ax.axvline(x0-dx,color='k',linestyle='--')
    #ax.axvline(x0+dx,color='k',linestyle='--')
    ax.plot(x,pdf,color='k')
    ax.fill_between(x=x,y1=pdf,color='k',alpha=0.1)
    # 95%
    lidx95 = abs(cdf-0.025).argmin()
    uidx95 = abs(cdf-0.975).argmin()
    ax.fill_between(x=x[lidx95:uidx95],y1=pdf[lidx95:uidx95],facecolor='k',edgecolor='k',alpha=0.1)
    #68%
    lidx68 = abs(cdf-0.158655254).argmin()
    uidx68 = abs(cdf-0.841344746).argmin()
    ax.fill_between(x=x[lidx68:uidx68],y1=pdf[lidx68:uidx68],facecolor='k',edgecolor='k',alpha=0.1)
    # black lines
    ax.errorbar(x=[x[lidx95],x[lidx95]],y=[0,pdf[lidx95]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[uidx95],x[uidx95]],y=[0,pdf[uidx95]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[lidx68],x[lidx68]],y=[0,pdf[lidx68]],color='k')#,lw=0.5)
    ax.errorbar(x=[x[uidx68],x[uidx68]],y=[0,pdf[uidx68]],color='k')#,lw=0.5)

    # split Fpi, FK and Fpi FK results
    ax.fill_between(x=x,y1=pdf_kk,color='b',alpha=0.6,label=r'$F^2 \rightarrow F_K^2$')
    ax.fill_between(x=x,y1=pdf_pk,color='g',alpha=0.6,label=r'$F^2 \rightarrow F_\pi F_K$')
    ax.fill_between(x=x,y1=pdf_pp,color='r',alpha=0.6,label=r'$F^2 \rightarrow F_\pi^2$')
    ax.set_xlim([1.1675,1.2275])
    ax.set_ylim(ymin=0)
    ax.set_xlabel(r'$F_K / F_\pi$',fontsize=16)
    ax.legend()
    ax.set_ylabel(r'Bayes Model Avg PDF',fontsize=16)
    if switches['optimized_priors']:
        ax.text(0.05,0.95, 'optimized prior widths',verticalalignment='center',horizontalalignment='left',transform=ax.transAxes,fontsize=16)
        plt.savefig('figures/model_avg_hist_logGBF_optimal_priors.pdf',transpareent=True)
    else:
        ax.text(0.05,0.95, 'default prior widths',verticalalignment='center',horizontalalignment='left',transform=ax.transAxes,fontsize=16)
        plt.savefig('figures/model_avg_hist_default_priors.pdf',transpareent=True)

    return FKFpi, model_var

def bma(switches,result,isospin):
    # read Bayes Factors
    logGBF_list = []
    for a in switches['ansatz']['type']:
        logGBF_list.append(result[a]['fit'].logGBF)
    # initiate a bunch of parameters
    # gA
    gA = 0
    gA_lst = []
    gA_dict = dict()
    # weights
    w_lst = []
    wd = dict()
    # p. dist. fcn
    pdf = 0
    pdfdict = dict()
    # c. dist. fcn.
    cdf = 0
    cdfdict = dict()
    # for plotting
    x = np.linspace(1.222,1.352,13000)
    # error breakdown for each model
    model_error = error_budget(switches,result)
    model_budget = {k:0 for k in model_error[list(model_error.keys())[0]]['std'].keys()}
    for a in switches['ansatz']['type']:
        r = result[a]['phys']['result']
        gA_dict[a] = r
        w = 1/sum(np.exp(np.array(logGBF_list)-result[a]['fit'].logGBF))
        sqrtw = np.sqrt(w) # sqrt scales the std dev correctly
        wd[a] = w
        w_lst.append(w)
        gA += gv.gvar(w*r.mean,sqrtw*r.sdev)
        gA_lst.append(r.mean)
        p = stats.norm.pdf(x,r.mean,r.sdev)
        pdf += w*p
        pdfdict[a] = w*p
        c = stats.norm.cdf(x,r.mean,r.sdev)
        cdf += w*c
        cdfdict[a] = w*c
        # error breakdown
        model_std = model_error[a]['std']
        model_budget = {k:model_budget[k]+w*model_std[k]**2 for k in model_std} # variance breakdown of model average
    gA_lst = np.array(gA_lst)
    w_lst = np.array(w_lst)
    model_var = np.sum(w_lst*gA_lst**2) - gA.mean**2
    final_error = np.sqrt(gA.sdev**2 + isospin**2)
    model_budget['isospin'] = isospin**2
    model_budget['model'] = model_var
    model_budget['total'] = model_budget['total']+model_budget['isospin']+model_budget['model'] # add in quadrature isospin and model variance
    pct_budget = {k:[np.sqrt(model_budget[k])/gA.mean*100] for k in model_budget} # percent breakdown of model average
    error = {'E(gA)': gA.mean, 's(gA)': final_error, 's(Mk)': np.sqrt(model_var), 'weights': wd, 'error_budget': model_budget, 'pct_budget': pct_budget, 'gA_dict':gA_dict}
    plot_params = {'x':x, 'pdf':pdf, 'pdfdict':pdfdict, 'cdf':cdf, 'cdfdict':cdfdict}
    return error, plot_params

def nnlo_prior_scan(switches,priors):
    for model in switches['ansatz']['models']:
        FPK = model.split('_')[-1]
        if os.path.exists('data/saved_prior_search.yaml'):
            with open('data/saved_prior_search.yaml','r') as fin:
                prior_grid = yaml.safe_load(fin.read())
            shutil.copyfile('data/saved_prior_search.yaml','data/saved_prior_search.yaml.bak')
        else:
            prior_grid = dict()
        switches['scale'] = FPK
        switches['ansatz']['model'] = model
        print('Prior width study: ',model)

        logGBF_array = []

        if model not in prior_grid:
            prior_grid[model] = dict()
        x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
        y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
        p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
        for key in priors:
            p_e[key] = priors[key]
        d_e = dict()
        d_e['x'] = x_e
        d_e['y'] = y_e

        if switches['prior_group']:
            p_range = priors['p_range']
            a_range = priors['a_range']
            if switches['refine_prior'] and len(prior_grid[model]) > 0:
                df = pd.DataFrame(prior_grid[model])
                print('    maximum logGBF = %.4f' %(df.stack().max()))
                sp_max = df.stack().idxmax()
                print('    (s_4, p_4) = ', sp_max[0],sp_max[1])
                new_s = []
                sf = float("%.1f" %float(sp_max[0]))
                for i in np.arange(sf-0.5,sf+0.5,.1):
                    ii = float("%.1f" %i)
                    if ii not in a_range and ii > 0.:
                        new_s.append(ii)
                new_s = np.array(new_s)
                a_range = np.concatenate((a_range,new_s))
                new_p = []
                pf = float(sp_max[1])
                for i in np.arange(pf-0.5,pf+0.5,.1):
                    ii = float("%.1f" %i)
                    if ii not in p_range and ii > 0.:
                        new_p.append(ii)
                new_p = np.array(new_p)
                p_range = np.concatenate((p_range,new_p))
                print(a_range)
                print(p_range)
                #sys.exit()
            z = np.zeros([len(a_range),len(p_range)])
            tot = len(a_range) * len(p_range)
            i_t = 0
            for i_s,s4 in enumerate(a_range):
                s4_s = str(s4)
                if s4_s not in prior_grid[model]:
                    prior_grid[model][s4_s] = dict()
                p_e['s_4']   = gv.gvar(0,s4)
                if 'alphaS' in model:
                    p_e['saS_4'] = gv.gvar(0,s4)
                for i_p,p4 in enumerate(p_range):
                    p4_s = str(p4)
                    sys.stdout.write('%4d out of %d, s_4 = %s p_4 = %s\r' %(i_t,tot,s4_s,p4_s))
                    sys.stdout.flush()
                    if p4_s not in prior_grid[model][s4_s]:
                        p_e['p_4'] = gv.gvar(0,p4)
                        p_e['k_4'] = gv.gvar(0,p4)
                        d_e['p'] = p_e
                        fit_e = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
                        fit_e.fit_data()
                        tmp = [ p_e[k].sdev for k in ['p_4','s_4'] ]
                        tmp.append(fit_e.fit.logGBF)
                        logGBF_array.append(tmp)
                        prior_grid[model][s4_s][p4_s] = float(fit_e.fit.logGBF)
                        i_t += 1
                        z[i_s,i_p] = fit_e.fit.logGBF
                    else:
                        tmp = [ p4, s4 ]
                        tmp.append(prior_grid[model][s4_s][p4_s])
                        logGBF_array.append(tmp)
                        i_t += 1
                        z[i_s,i_p] = prior_grid[model][s4_s][p4_s]
        else:
            print('individual prior width study not supported [yet]')

        df = pd.DataFrame(prior_grid[model])
        print('    maximum logGBF', df.stack().max())
        sp_max = df.stack().idxmax()
        print('    (s_4, p_4) = ', sp_max[0],sp_max[1])
        logGBF_array = np.array(logGBF_array)
        logGBF_max = np.argmax(logGBF_array[:,-1])
        print('optimal prior widths for %s' %model)
        logGBF_optimal = logGBF_array[logGBF_max]
        tmp = ''
        for i_k,k in enumerate(['p_4','s_4']):
            tmp += '%s = %.2f ' %(k,logGBF_optimal[i_k])
        print(tmp,'logGBF = %.4f\n' %(logGBF_optimal[-1]))

        lgbf = logGBF_array[:,-1]
        w = np.exp(lgbf - logGBF_optimal[-1])
        #w = w / w.sum()
        logGBF_w = np.copy(logGBF_array)
        logGBF_w[:,-1] = w

        prior_file = open('data/saved_prior_search.yaml', 'w')
        yaml.dump(prior_grid, prior_file)
        prior_file.close()


    sys.exit()

    if switches['prior_group']:
        z = np.exp(z - logGBF_optimal[-1])
        cmap_old = sns.cubehelix_palette(8, as_cmap=True)
        levels = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
        from matplotlib.colors import ListedColormap
        cmap = ListedColormap(sns.color_palette("BrBG_r", len(levels)).as_hex())
        fig, ax = plt.subplots(constrained_layout=True)
        CS = ax.contour(p_range,a_range,z,levels=levels, cmap=cmap)
        ax.clabel(CS, CS.levels[0::2], fmt='%2.1f', colors='k', fontsize=12)
        ax.plot(logGBF_optimal[0], logGBF_optimal[1],marker='X',color='r',markersize=20)
        ax.set_xlabel(r'$\sigma_{nnlo,\chi{\rm PT}}$',fontsize=16)
        ax.set_ylabel(r'$\sigma_{nnlo,a^2}$',fontsize=16)
        fig_name = 'prior_width_'+switches['ansatz']['model']+'_'+switches['scale']
        if not os.path.exists('figures'):
            os.makedirs('figures')
        plt.savefig('figures/'+fig_name+'.pdf',transparent=True)

def nlo_report(switches,priors,FV=True):
    latex = {
        'ma_nlo':'ma nlo',   'ma_nlo_FV':'ma nlo w/FV',
        'xpt_nlo':'xpt nlo', 'xpt_nlo_FV':'xpt nlo w/FV'}
    marker = {
        'ma_nlo':'s',  'ma_nlo_FV':'o',
        'xpt_nlo':'*', 'xpt_nlo_FV':'d'}
    color = {
        'ma_nlo':'r',  'ma_nlo_FV':'b',
        'xpt_nlo':'m', 'xpt_nlo_FV':'g'}
    if FV:
        models = ['ma_nlo','ma_nlo_FV','xpt_nlo','xpt_nlo_FV']
    else:
        models = ['ma_nlo','xpt_nlo']

    fit_results = dict()
    for model in models:
        switches['ansatz']['model'] = model
        print('EFT: ',model)
        model_result = dict()
        for e in switches['ensembles']:
            switches['ensembles_fit'] = [e]
            x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
            y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
            p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
            p_e['L5'] = priors['L5']
            d_e = dict()
            d_e['x'] = x_e
            d_e['y'] = y_e
            d_e['p'] = p_e
            fit_e = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
            fit_e.fit_data()
            model_result[e] = fit_e
        fit_results[model] = model_result
    fig_name = 'nlo_report'
    if FV:
        fig_name += '_FV'

    fig = plt.figure(fig_name)
    ax  = plt.axes([.13, .13, .84, .84])

    if len(models) == 2:
        print("%8s & %13s & %13s & \\\\" \
            %('ensemble',latex[models[0]], latex[models[1]]))
    else:
        print("%8s & %13s & %13s & %13s & %13s\\\\" \
            %('ensemble',latex[models[0]], latex[models[1]],latex[models[2]], latex[models[3]]))
    print("\\hline")
    for i_e,e in enumerate(switches['ensembles']):
        s = "%8s" %e
        for model in models:
            s += " & %13s" %str(fit_results[model][e].fit.p['L5'])
        s += "\\\\"
        print(s)

        for i_m,model in enumerate(models):
            if i_e == 0:
                label=latex[models[i_m]]
            else:
                label=''
            y = len(switches['ensembles']) - i_e + 0.1*i_m
            ax.errorbar(x=fit_results[model][e].fit.p['L5'].mean,y=y,
                xerr=fit_results[model][e].fit.p['L5'].sdev,linestyle='None',
                marker=marker[model],mfc='None',color=color[model],label=label)
    plt.yticks(np.arange(len(switches['ensembles']),0,-1),tuple(switches['ensembles']))
    ax.set_xlabel(r'$L_5$',fontsize=16)
    ax.legend(loc=1,fontsize=16)
    ax.set_xlim(-0.0002,0.001)
    if not os.path.exists('figures'):
        os.makedirs('figures')
    plt.savefig('figures/'+fig_name+'.pdf',transparent=True)

def fit_checker(switches,priors):
    fit_results = dict()
    for model in switches['ansatz']['models']:
        switches['ansatz']['model'] = model
        print('EFT: ',model)
        x_e = {k:gv_data['x'][k] for k in switches['ensembles']}
        y_e = {k:gv_data['y'][k] for k in switches['ensembles']}
        p_e = {k: gv_data['p'][k] for k in gv_data['p'] if k[0] in switches['ensembles']}
        for key in priors:
            p_e[key] = priors[key]
        if switches['optimized_priors']:
            try:
                for key in ['p_4','k_4','s_4','saS_4']:
                    p_e[key] = gv.gvar(0,ip.nnlo_width[model][switches['scale']][key])
            except Exception as e:
                print('WARNING: non optimized NNLO prior widths')
                print(str(e))
        else:
            print('using default prior widths')
        d_e = dict()
        d_e['x'] = x_e
        d_e['y'] = y_e
        d_e['p'] = p_e
        fit_e = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
        fit_e.fit_data()
        fit_results[model] = fit_e
        if switches['print_fit']:
            print(fit_e.fit.format(maxline=True))

    print('model & PK& chisq/dof& logGBF& $L_5$& $k_4$& $p_4$& $s_4$& $F_K/F_\pi$\\\\')
    for model in switches['ansatz']['models']:
        chi2dof,dof,logGBF,L5,k4,p4,s4,FKFpi = fit_results[model].check_fit(phys_point)
        print('%34s& %s& %.3f [%2d]& %.2f& %14s& %8s& %8s& %8s& %10s\\\\'\
            %(model.replace('_','\_'),switches['scale'],chi2dof,dof,logGBF,L5,k4,p4,s4,FKFpi))

def check_fit_function(switches,check_point):
    print('Performing check of fit function')
    for k in check_point:
        print('%s = %f' %(k,check_point[k]))
    switches['ensembles_fit'] = ['check']

    x_check = dict()
    x_check['check'] = {k:np.inf for k in ['mpiL','mkL']}
    x_check['check']['alphaS'] = 0.

    y_check = dict()
    y_check['check'] = 1.19897

    p_check = dict()
    for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
        p_check[k] = 0.
    for key in ['L1','L2','L3','L4','L5','L6','L7','L8','k_4','p_4']:
        p_check[key] = check_point[key]

    d_e = dict()
    d_e['x'] = x_check
    d_e['y'] = y_check

    for FPK in switches['scales']:
        model = 'xpt_nnlo_'+FPK
        switches['scale'] = FPK
        switches['ansatz']['model'] = model
        if FPK == 'PP':
            Lchi_check = 4*np.pi * check_point['Fpi']
        elif FPK == 'PK':
            Lchi_check = 4*np.pi * np.sqrt(check_point['Fpi'] * check_point['FK'])
        elif FPK == 'KK':
            Lchi_check = 4*np.pi * check_point['FK']
        p_check[('check','Lchi_'+FPK)] = Lchi_check
        d_e['p'] = p_check
        #print('debug p',d_e['p'])

        fit_check = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
        fit_check.check_fit_function(check_point)

        # ratio
        model = 'xpt-ratio_nnlo_'+FPK
        switches['ansatz']['model'] = model
        fit_check = xpt.Fit(switches,xyp_init=d_e,phys=phys_point)
        fit_check.check_fit_function(check_point)


if __name__ == "__main__":
    main()
