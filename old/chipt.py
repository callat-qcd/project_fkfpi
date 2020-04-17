import sys
import copy
import scipy.special as spsp
import numpy as np
import matplotlib.pyplot as plt
import gvar as gv
import lsqfit

# Figure formatting
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
fig_size2 = (fig_width, fig_width * 1.6)
plt_axes  = [0.14,0.14,0.855,0.855]
fs_text   = 20 # font size of text
fs_leg    = 16 # legend font size
mrk_size  = '5' # marker size
tick_size = 16 # tick size
lw        = 1 # line width

pi = np.pi
# FF function defined in 1711.11328
sys.path.append('py_chiron')
import chiron
def FF(x):
    if isinstance(x, gv.GVar):
        f = chiron.FF(x.mean)
        stepSize = 1e-7# * chiron.FF(x.mean)
        dfdx = 0.5*(chiron.FF(x.mean+stepSize) - chiron.FF(x.mean-stepSize))/stepSize
        return gv.gvar_function(x, f, dfdx)
    else:
        return chiron.FF(x)
# approximate version - ai in file from Bijnens, some of them do not match what
#                       is provided in 1711.11328: a1, a2, a5, a8, a11
#                       but these values are necessary to reproduce the approximate
#                       value of FK/Fpi provided in the paper
ai = dict()
ai['a1']  =  0.8740
ai['a2']  = -2.172
ai['a3']  =  0.8294
ai['a4']  = -0.4583
ai['a5']  =  3.716
ai['a6']  = -0.1113
ai['a7']  =  0.8776
ai['a8']  = -1.635
ai['a9']  =  1.4697
ai['a10'] = -0.1406
ai['a11'] = -1.343
ai['a12'] =  0.2731
ai['a13'] = -0.2109
def FF_approximate(x):
    ff  = ai['a1']
    ff += (ai['a2'] + ai['a3'] * np.log(x) + ai['a4'] * np.log(x)**2) * x
    ff += (ai['a5'] + ai['a6'] * np.log(x) + ai['a7'] * np.log(x)**2) * x**2
    ff += (ai['a8'] + ai['a9'] * np.log(x) + ai['a10']* np.log(x)**2) * x**3
    ff += (ai['a11']+ ai['a12']* np.log(x) + ai['a13']* np.log(x)**2) * x**4
    return ff

class Fit(object):
    def __init__(self,switches,xyp_init,phys, name=None):
        ''' In the init, we will use parameters specified for the fit to declare
            the precise fit function used, such that the minimizer has no
            if/else decisions to make
        '''
        # collect all parameters for a given "model"
        self.phys      = phys
        self.switches  = switches
        self.x         = xyp_init['x']
        self.y         = xyp_init['y']
        self.p_init    = xyp_init['p']
        self.model     = self.switches['ansatz']['model']
        self.eft       = self.switches['ansatz']['model'].split('_')[0]
        self.order     = self.switches['ansatz']['model'].split('_')[1]
        self.FF        = self.switches['ansatz']['model'].split('_')[-1]
        self.fv        = 'FV'     in self.switches['ansatz']['model']
        self.alphaS    = 'alphaS' in self.switches['ansatz']['model']
        self.logSq     = 'logSq'  in self.switches['ansatz']['model']
        self.ct        = 'ct'     in self.switches['ansatz']['model']
        if self.logSq and self.ct:
            sys.exit('Setting logSq and ct = True == nnlo_full: redundant, aborting')
        self.nnlo_full = not any(md in self.switches['ansatz']['model'] for md in ['logSq','ct'])
        self.a4        = 'a4'     in self.switches['ansatz']['model']
        self.lec_nnlo  = ['L5','L4','k_4','p_4','s_4','saS_4']
        if self.nnlo_full:
            self.lec_nnlo += ['L1','L2','L3','L6','L7','L8']
        self.lec_nnnlo = ['s_6','sk_6','sp_6','kp_6','k_6','p_6']

        # create dictionary of function pieces
        self.fit_func_elementals = {# the NLO terms are LO + NLO, while the rest ar NnLO
            'xpt_nlo'         : self.xpt_nlo,
            'xpt-ratio_nlo'   : self.xpt_ratio_nlo,
            'ma_nlo'          : self.ma_nlo,
            'ma-longform_nlo' : self.ma_longform_nlo,
            'ma-ratio_nlo'    : self.ma_ratio_nlo,
            'taylor_nlo'      : self.taylor_nlo,
            # NNLO terms
            'nnlo_ct'         : self.nnlo_ct,
            'xpt_nnlo_logSq'  : self.xpt_nnlo_logSq,
            'xpt_nnlo_log'    : self.xpt_nnlo_log,
            'xpt_ratio_nnlo'  : self.xpt_ratio_nnlo,
            'nnlo_alphaS'     : self.nnlo_alphaS,
            # NNLO scale and F=Fpi, FK fix
            'xpt_nnlo_FF_PP'  : self.xpt_nnlo_FF_PP,
            #'xpt_nnlo_FF_PK'  : self.nnlo_FF_PK,
            #'xpt_nnlo_FF_KK'  : self.nnlo_FF_KK,
            # NNNLO terms
            #'nnnlo_ct'        : self.nnnlo_ct, # all NNNLO terms
            'nnnlo_a4'        : self.nnnlo_a4,
        }

        # fit function
        self.func_elements = self.make_function_elementals()
        #self.fit_function = self.make_fit_function

        # Init the finite volume stuff
        if not self.fv:
            self.I  = self.I_IV
            self.dI = self.dI_IV
        else:
            # declare tadpole functions
            self.I  = self.I_FV
            self.dI = self.dI_FV
            # make finite volume corrections
            self.fv_I  = dict()
            self.fv_dI = dict()
            def k1k2k3(mL):
                cn = np.array([6,12,8,6,24,24,0,12,30,24,24,8,24,48,0,6,48,36,24,24])
                n_mag = np.sqrt(np.arange(1,len(cn)+1,1))
                k1 = np.sum(cn * spsp.kn(1,n_mag * mL) / mL / n_mag)
                k2 = np.sum(cn * spsp.kn(2,n_mag * mL) )
                k0 = np.sum(cn * spsp.kn(0,n_mag * mL) )
                return k0,k1,k2
            for e in self.switches['ensembles_fit']:
                # mpi
                Lchi     = self.p_init[(e,'Lchi_'+self.FF)]
                esq      = self.p_init[(e,'mpi')]**2 / Lchi**2
                mL       = self.x[e]['mpiL']
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'p2')]  = 4 * k1
                self.fv_dI[(e,'p2')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                #print(e,'p2',self.fv_I[(e,'p2')],self.fv_dI[(e,'p2')])
                # mK
                esq      = self.p_init[(e,'mk')]**2 / Lchi**2
                mL       = self.x[e]['mkL']
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'k2')]  = 4 * k1
                self.fv_dI[(e,'k2')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                # meta
                esq      = 4./3 * self.p_init[(e,'mk')]**2 / Lchi**2
                esq     += -1./3 *self.p_init[(e,'mpi')]**2 / Lchi**2
                mL       = np.sqrt(4./3*self.x[e]['mkL']**2 - 1./3 * self.x[e]['mpiL']**2)
                k0,k1,k2 = k1k2k3(mL)
                self.fv_I[(e,'e2')]  = 4 * k1
                self.fv_dI[(e,'e2')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                if 'ma' in self.eft:
                    # mss
                    esq = self.p_init[(e,'mss')]**2 / Lchi**2
                    mL  = self.x[e]['mssL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'s2')]  = 4 * k1
                    self.fv_dI[(e,'s2')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                    # mju
                    esq = self.p_init[(e,'mju')]**2 / Lchi**2
                    mL  = self.x[e]['mjuL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'ju')]  = 4 * k1
                    self.fv_dI[(e,'ju')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                    # mjs
                    esq = self.p_init[(e,'mjs')]**2 / Lchi**2
                    mL  = self.x[e]['mjsL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'js')]  = 4 * k1
                    self.fv_dI[(e,'js')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                    # mru
                    esq = self.p_init[(e,'mru')]**2 / Lchi**2
                    mL  = self.x[e]['mruL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'ru')]  = 4 * k1
                    self.fv_dI[(e,'ru')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                    # mrs
                    esq = self.p_init[(e,'mrs')]**2 / Lchi**2
                    mL  = self.x[e]['mrsL']
                    k0,k1,k2 = k1k2k3(mL)
                    self.fv_I[(e,'rs')]  = 4 * k1
                    self.fv_dI[(e,'rs')] = (2*k1 - k0 - k2) #/ (4*pi)**2
                    # mX
                    e2  = 4./3 * self.p_init[(e,'mk')]**2 / Lchi**2
                    e2 += -1./3 *self.p_init[(e,'mpi')]**2 / Lchi**2
                    esq = e2 + self.p_init[(e,'a2DI')] / Lchi**2
                    mL  = np.sqrt(4./3*self.x[e]['mkL']**2 - 1./3 * self.x[e]['mpiL']**2)
                    mL  = np.sqrt(esq / e2) * mL
                    self.fv_I[(e,'x2')]  = 4 * k1
                    self.fv_dI[(e,'x2')] = (2*k1 - k0 - k2) #/ (4*pi)**2

    def make_function_elementals(self):
        func_elements = [self.eft+'_nlo']

        if self.eft == 'taylor':
            if self.order in ['nnlo','nnnlo']:
                func_elements += ['nnlo_ct']
            if self.order in ['nnnlo']:
                func_elements += ['nnnlo_ct']
        else:
            if self.order in ['nnlo','nnnlo']:
                if self.alphaS:
                    func_elements += ['nnlo_alphaS']
                if self.nnlo_full:
                    func_elements += ['nnlo_ct', 'xpt_nnlo_logSq', 'xpt_nnlo_log']
                if 'ratio' in self.eft:
                    func_elements += ['xpt_ratio_nnlo']
                else:
                    if self.ct:
                        func_elements += ['nnlo_ct']
                    if self.logSq:
                        func_elements += ['xpt_nnlo_logSq']
                # FF choice
                if self.FF not in ['PP','PK','KK']:
                    sys.exit('unrecognized FF choice [PP, PK, KK]: '+self.FF)
                else:
                    func_elements += ['xpt_nnlo_FF_'+self.FF]

                # a^4?
                if self.order == 'nnlo' and self.a4:
                    func_elements += ['nnnlo_a4']
                if self.order == 'nnnlo': # this gets all fits
                    func_elements += ['nnnlo_ct']
        return func_elements

    #def __reduce__(self):
    #    return (self.__class__, (self.make_fit_function, ))

    """
    def __getstate__(self):
    #    return dict(_registry=self._registry, _obj=self._obj)
        print('DEBUG self')
        for k in self.instances:
            print(k.name)
        state = self.__dict__.copy()
        # gvar dump can not handle this function
        #print('DEBUG: state keys')
        #for k in state:
        #    print(k)
        print('===================================================')
        print('===================================================')
        print('===================================================')
        print('===================================================')
        #del state['make_fit_function']
        return state
    """
    def make_fit_function(self,term_list):
        def fit_function(*args):
            results = dict()
            for elem in term_list:
                tmp_results = self.fit_func_elementals[elem](*args)
                for k,v in tmp_results.items():
                    results[k] = results.get(k, 0.) + v
            return results
        return fit_function

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
        if self.switches['debug']:
            print('DEBUG:',"%9s" %'ens','  ',"%11s" %'Lchi','eps_pi**2')
        for e in self.switches['ensembles_fit']:
            Lchi = self.p_init[(e,'Lchi_'+self.FF)]
            #Lchi = self.x[e]['Lchi_'+self.FF]
            p[(e,'p2')] = self.p_init[(e,'mpi')]**2 / Lchi**2
            p[(e,'k2')] = self.p_init[(e,'mk')]**2 /  Lchi**2
            p[(e,'e2')] = 4./3 * p[(e,'k2')] - 1./3 * p[(e,'p2')]
            if self.switches['debug']:
                print('DEBUG:',"%9s" %e,self.FF,"%11s" %Lchi,p[(e,'p2')])
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
            if self.order in ['nnlo','nnnlo'] or (self.order == 'nlo' and self.a4):
                p[(e,'a2')]   = self.p_init[(e,'aw0')]**2 / (4 * pi)
        ''' Add LECs '''
        if self.order in ['nlo','nnlo','nnnlo']:
            p['L5'] = self.p_init['L5']
        if self.order in ['nnlo','nnnlo']:
            p['L4']  = self.p_init['L4'] #L4 appears in the mu correction
            p['s_4'] = self.p_init['s_4']
            p['k_4'] = self.p_init['k_4']
            p['p_4'] = self.p_init['p_4']
            if self.nnlo_full:
                p['L1'] = self.p_init['L1']
                p['L2'] = self.p_init['L2']
                p['L3'] = self.p_init['L3']
                p['L6'] = self.p_init['L6']
                p['L7'] = self.p_init['L7']
                p['L8'] = self.p_init['L8']
            if self.alphaS:
                p['saS_4'] = self.p_init['saS_4']
            # enable test of a**4 term without any other N3LO
            if self.a4 and self.order not in ['nnnlo']:
                p['s_6'] = self.p_init['s_6']
                #p['sk_6'] = self.p_init['sk_6']
                #p['sp_6'] = self.p_init['sp_6']
        if self.order in ['nnnlo']:
            p['kp_6'] = self.p_init['kp_6']
            p['k_6']  = self.p_init['k_6']
            p['p_6']  = self.p_init['p_6']
            p['s_6']  = self.p_init['s_6']
            p['sk_6'] = self.p_init['sk_6']
            p['sp_6'] = self.p_init['sp_6']
        if self.order == 'nlo' and self.a4:
            p['s_4'] = self.p_init['s_4']
            p['s_6'] = self.p_init['s_6']
            if self.alphaS:
                p['saS_4'] = self.p_init['saS_4']

        return p

    def make_x_lec(self,x,p):
        x_e = dict()
        lec = {key:val for key,val in p.items() if isinstance(key,str)}
        for e in x:
            x_par = gv.BufferDict()
            x_par['p2'] = {'esq':p[(e,'p2')]}
            x_par['k2'] = {'esq':p[(e,'k2')]}
            x_par['e2'] = {'esq':p[(e,'e2')]}
            if self.fv:
                x_par['p2'].update({'fvI':self.fv_I[(e,'p2')], 'fvdI':self.fv_dI[(e,'p2')]})
                x_par['k2'].update({'fvI':self.fv_I[(e,'k2')], 'fvdI':self.fv_dI[(e,'k2')]})
                x_par['e2'].update({'fvI':self.fv_I[(e,'e2')], 'fvdI':self.fv_dI[(e,'e2')]})
            if 'ma' in self.eft:
                x_par['s2'] = {'esq':p[(e,'s2')]}
                x_par['ju'] = {'esq':p[(e,'ju')]}
                x_par['js'] = {'esq':p[(e,'js')]}
                x_par['ru'] = {'esq':p[(e,'ru')]}
                x_par['rs'] = {'esq':p[(e,'rs')]}
                x_par['x2'] = {'esq':p[(e,'x2')]}
                if self.fv:
                    x_par['s2'].update({'fvI':self.fv_I[(e,'s2')], 'fvdI':self.fv_dI[(e,'s2')]})
                    x_par['ju'].update({'fvI':self.fv_I[(e,'ju')], 'fvdI':self.fv_dI[(e,'ju')]})
                    x_par['js'].update({'fvI':self.fv_I[(e,'js')], 'fvdI':self.fv_dI[(e,'js')]})
                    x_par['ru'].update({'fvI':self.fv_I[(e,'ru')], 'fvdI':self.fv_dI[(e,'ru')]})
                    x_par['rs'].update({'fvI':self.fv_I[(e,'rs')], 'fvdI':self.fv_dI[(e,'rs')]})
                    x_par['x2'].update({'fvI':self.fv_I[(e,'x2')], 'fvdI':self.fv_dI[(e,'x2')]})

                x_par['dju2'] = p[(e,'dju2')]
                x_par['drs2'] = p[(e,'drs2')]
            if self.order in ['nnlo','nnnlo'] or (self.order == 'nlo' and self.a4):
                x_par['a2']   = p[(e,'a2')]
            #if self.alphaS:
            x_par['alphaS'] = x[e]['alphaS']

            x_e[e] = x_par

        return x_e,lec

    # Fit Functions
    # NLO terms
    def xpt_nlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            r[e]  = 1.
            r[e] +=  self.FK_xpt_nlo(x_e[e],lec)
            r[e] += -self.Fpi_xpt_nlo(x_e[e],lec)
        return r

    def xpt_ratio_nlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            num  = 1. + self.FK_xpt_nlo(x_e[e],lec)
            den  = 1. + self.Fpi_xpt_nlo(x_e[e],lec)
            r[e] = num / den
        return r

    def Fpi_xpt_nlo(self,x,lec):
        r  = -self.I(x['p2'])
        r += -0.5 * self.I(x['k2'])
        r += lec['L5'] * (4*pi)**2 * 4 * x['p2']['esq']
        if self.switches['debug_nlo_check']:
            print('    -eps2_p log(eps2_p) = %s' %(-self.I(x['p2'])))
            print('-1/2 eps2_k log(eps2_k) = %s' %(-0.5 * self.I(x['k2'])))
            print('   4 eps2_p L5bar      = %s' %(lec['L5'] * (4*pi)**2 * 4 * x['p2']['esq']))
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
        return r

    def FK_xpt_nlo(self,x,lec):
        r  = -3./8 * self.I(x['p2'])
        r += -3./4 * self.I(x['k2'])
        r += -3./8 * self.I(x['e2'])
        r += lec['L5'] * (4*pi)**2 * 4 * x['k2']['esq']
        if self.switches['debug_nlo_check']:
            print('-3/8 eps2_p log(eps2_p) = %s' %(-3./8 * self.I(x['p2'])))
            print('-3/4 eps2_k log(eps2_k) = %s' %(-3./4 * self.I(x['k2'])))
            print('-3/8 eps2_e log(eps2_e) = %s' %(-3./8 * self.I(x['e2'])))
            print('-3/8 eps2_k L5bar      = %s' %(lec['L5'] * (4*pi)**2 * 4 * x['k2']['esq']))
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
        return r

    # MA NLO - the long form is below with the Tadpole Integrals
    def ma_nlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            r[e]  = 1.
            r[e] += self.FK_ma_nlo(x_e[e],lec)
            r[e] += -self.Fpi_ma_nlo(x_e[e],lec)
        return r

    def ma_ratio_nlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            num  = 1. + self.FK_ma_nlo(x_e[e],lec)
            den  = 1. + self.Fpi_ma_nlo(x_e[e],lec)
            r[e] = num / den
        return r

    # Taylor Expansion
    def taylor_nlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            r[e]  = 1.
            r[e] +=  lec['L5'] * (4*pi)**2 * 4 * x['k2']['esq']
            r[e] += -lec['L5'] * (4*pi)**2 * 4 * x['p2']['esq']
        return r

    # NNLO Terms
    def nnlo_ct(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            r[e]  = (k2 - p2) * k2 * lec['k_4']
            r[e] += (k2 - p2) * p2 * lec['p_4']
            a2 = x_e[e]['a2']
            r[e] += (k2 - p2) * a2 * lec['s_4']
            if self.switches['debug_nnlo_check']:
                print('\nct        = %f' %r[e])
        return r

    def nnlo_alphaS(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            a2 = x_e[e]['a2']
            r[e] = (k2 - p2) * a2 * x_e[e]['alphaS'] * lec['saS_4']
        return r

    def xpt_nnlo_logSq(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            lk = np.log(k2)
            lp = np.log(p2)
            e2 = 4./3 * k2 - 1./3 * p2
            le = np.log(e2)
            r[e]  = lp * lp   * ( 11./24   * p2 *k2  -131./192 * p2 *p2)
            r[e] += lp * lk   * (-41./96   * p2 *k2  -3./32    * p2 *p2)
            r[e] += lp * le   * ( 13./24   * p2 *k2  +59./96   * p2 *p2)
            r[e] += lk * lk   * ( 17./36   * k2 *k2  +7./144   * p2 *k2)
            r[e] += lk * le   * (-163./144 * k2 *k2  -67./288  * p2 *k2 + 3./32   * p2 *p2)
            r[e] += le * le   * ( 241./288 * k2 *k2  -13./72   * p2 *k2 - 61./192 * p2 *p2)
            tmp = gv.gvar(r[e])
            r[e] += k2**2 * FF(p2/k2)
            if self.switches['debug_nnlo_check']:
                print('logSq     = %f' %(tmp))
                print('FF        = %f' %(k2**2 * FF(p2/k2)))
                FF_approx = k2**2 * FF_approximate(p2/k2)
                print('FF_approx = %f' %(k2**2 * FF_approximate(p2/k2)))
        return r

    def xpt_nnlo_log(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        L1 = lec['L1']
        L2 = lec['L2']
        L3 = lec['L3']
        L4 = lec['L4']
        L5 = lec['L5']
        L6 = lec['L6']
        L7 = lec['L7']
        L8 = lec['L8']
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            lk = np.log(k2)
            lp = np.log(p2)
            e2 = 4./3 * k2 - 1./3 * p2
            le = np.log(e2)
            # extra counter terms
            tk  = 8*(4*pi)**2 *L5 *(8*L4 +3*L5 -16*L6 -8*L8) -2*L1 -L2 -1./18*L3 +4./3*L5 -16*L7 -8*L8
            tp  = 8*(4*pi)**2 *L5 *(4*L4 +5*L5  -8*L6 -8*L8) -2*L1 -L2 -5./18*L3 -4./3*L5 +16*L7 +8*L8
            ct  = (4*pi)**2 * (k2 - p2) * tk * k2
            ct += (4*pi)**2 * (k2 - p2) *tp * p2
            r[e] = ct

            # single log terms
            C1  = -(7./9     + (4*pi)**2 * 11./2 *L5) *p2 *k2
            C1 += -(113./72  + (4*pi)**2 *(4*L1 +10*L2 +13./2*L3 -21./2*L5))*p2**2

            C2  =  (209./144 + (4*pi)**2 *3*L5)*p2*k2
            C2 +=  (53./96   + (4*pi)**2 *(4*L1 +10*L2 +5*L3 -5*L5))*k2**2

            C3  =  (13./18   + (4*pi)**2 *(8./3*L3 -2./3 *L5 -16*L7  -8*L8))* k2*k2
            C3 += -(4./9     + (4*pi)**2 *(4./3*L3 +25./6*L5 -32*L7 -16*L8))* p2*k2
            C3 +=  (19./288  + (4*pi)**2 *(1./6*L3 +11./6*L5 -16*L7  -8*L8))* p2*p2

            r[e] += C1 * lp  + C2 * lk  + C3 * le
            if self.switches['debug_nnlo_check']:
                print('extra ct  = %f' %ct)
                print('C1 = %s' %C1)
                print('C2 = %s' %C2)
                print('C3 = %s' %C3)
                print('log terms = %f' %(r[e] - ct))
        return r

    def xpt_ratio_nnlo(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        L4 = lec['L4']
        L5 = lec['L5']
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            lk = np.log(k2)
            lp = np.log(p2)
            e2 = 4./3 * k2 - 1./3 * p2
            le = np.log(e2)

            r[e]  = (p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4)) \
                * (3./8*p2*lp +3./4*k2*lk +3./8*e2*le +4*(4*pi)**2 *(p2*L4 +k2*(2*L4+L5)))
            #r[e] += -1.0*( p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4) )**2
            r[e] += 0.5*( p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4) )**2
            if self.switches['debug_nnlo_check']:
                print('ratio fix = %f\n' %(r[e]))
        return r

    def xpt_nnlo_FF_PP(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        L4 = lec['L4']
        L5 = lec['L5']
        for e in x_e:
            k2 = x_e[e]['k2']['esq']
            p2 = x_e[e]['p2']['esq']
            lk = np.log(k2)
            lp = np.log(p2)
            r[e] = 3./2 *(k2-p2) *(p2*lp + 0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4))
            if self.switches['debug_nnlo_check']:
                print('mu fix    = %f\n' %(r[e]))
        return r

    # NNNLO terms
    def nnnlo_a4(self,x,p):
        r = dict()
        x_e, lec = self.make_x_lec(x,p)
        for e in x_e:
            r[e] = (x_e[e]['k2']['esq'] - x_e[e]['p2']['esq']) * x_e[e]['a2']**2 * lec['s_6']
        return r

    # MA NLO Long Form for debugging
    def ma_longform_nlo(self,x,p):
        ''' for debugging - a cross check expression '''
        r = dict()
        for e in x:
            x_e[e] = self.x_p_ens[e]
            lec   = self.lec_ens[e]
            r[e]  = self.counterterms(x_e[e],lec)
            r[e] += 1.
            r[e] += x_e[e]['ju'] * np.log(x_e[e]['ju'])
            r[e] += 0.5 * x_e[e]['ru'] * np.log(x_e[e]['ru'])
            # (+) kaon log terms
            r[e] += -0.5 * x_e[e]['ju'] * np.log(x_e[e]['ju'])
            r[e] += -1./4 * x_e[e]['ru'] * np.log(x_e[e]['ru'])
            r[e] += -0.5*x_e[e]['js'] * np.log(x_e[e]['js'])
            r[e] += -1./4 * x_e[e]['rs'] * np.log(x_e[e]['rs'])
            r[e] += -x_e[e]['dju2'] / 8
            r[e] +=  x_e[e]['dju2']**2 / 24 / (x_e[e]['x2'] - x_e[e]['p2'])
            r[e] +=  x_e[e]['drs2'] * (x_e[e]['k2']-x_e[e]['p2']) / 6 / (x_e[e]['x2']-x_e[e]['s2'])
            r[e] += -x_e[e]['dju2'] * x_e[e]['drs2'] / 12 / (x_e[e]['x2'] - x_e[e]['s2'])
            r[e] += np.log(x_e[e]['p2'])/24 * (3*x_e[e]['p2'] \
                    - 3*x_e[e]['dju2']*(x_e[e]['x2']+x_e[e]['p2'])/(x_e[e]['x2']-x_e[e]['p2']) \
                    + x_e[e]['dju2']**2 * x_e[e]['x2']/(x_e[e]['x2']-x_e[e]['p2'])**2\
                    -4*x_e[e]['dju2']*x_e[e]['drs2']*x_e[e]['p2']/(x_e[e]['x2']-x_e[e]['p2'])/(x_e[e]['s2']-x_e[e]['p2'])\
                    )
            r[e] += -x_e[e]['x2']/24 * np.log(x_e[e]['x2'])*(9 \
                    -6*x_e[e]['dju2']/(x_e[e]['x2']-x_e[e]['p2']) \
                    + x_e[e]['dju2']**2/(x_e[e]['x2']-x_e[e]['p2'])**2\
                    +x_e[e]['drs2']*(4*(x_e[e]['k2']-x_e[e]['p2'])+6*(x_e[e]['s2']-x_e[e]['x2']))/(x_e[e]['x2']- x_e[e]['s2'])**2 \
                    -2*x_e[e]['dju2']*x_e[e]['drs2']*(2*x_e[e]['s2']-x_e[e]['p2']-x_e[e]['x2'])/(x_e[e]['x2']-x_e[e]['s2'])**2/(x_e[e]['x2']-x_e[e]['p2'])\
                    )
            r[e] += np.log(x_e[e]['s2'])/12 * (3*x_e[e]['s2'] \
                    +x_e[e]['drs2']*(3*x_e[e]['s2']**2 + 2*(x_e[e]['k2']-x_e[e]['p2'])*x_e[e]['x2'] -3*x_e[e]['s2']*x_e[e]['x2'])/(x_e[e]['x2']-x_e[e]['s2'])**2\
                    -x_e[e]['dju2']*x_e[e]['drs2']*(2*x_e[e]['s2']**2 - x_e[e]['x2']*(x_e[e]['s2']+x_e[e]['p2']))/(x_e[e]['x2']-x_e[e]['s2'])**2 / (x_e[e]['s2']-x_e[e]['p2'])\
                    )
        return r

    # Tadpole Integrals
    def I_IV(self,x):
        return x['esq'] * np.log(x['esq'])
    def I_FV(self,x):
        r  = x['esq'] * np.log(x['esq'])
        r += x['esq'] * x['fvI']
        return r
    def dI_IV(self,x):
        return 1 + np.log(x['esq'])
    def dI_FV(self,x):
        r = 1 + np.log(x['esq'])
        r += x['fvdI']
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


    def Fpi_ma_nlo(self,x,lec):
        r  = -self.I(x['ju'])
        r += -0.5 * self.I(x['ru'])
        r += lec['L5'] * (4*pi)**2 * 4 * x['p2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])
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
        r += lec['L5'] * (4*pi)**2 * 4 * x['k2']['esq']
        if 'ratio' in self.eft:
            r += lec['L4'] * (4*pi)**2 * (4 * x['p2']['esq'] + 8 * x['k2']['esq'])

        r += x['dju2']             * -1./8  * self.dI(x['p2'])
        r += x['dju2']             *  1./4  * self.K(x['p2'],x['x2'])
        r += x['dju2']**2          * -1./24 * self.K21(x['p2'],x['x2'])
        r += x['dju2'] * x['drs2'] *  1./12 * self.K21(x['s2'],x['x2'])
        r += x['dju2'] * x['drs2'] * -1./6  * self.K123(x['p2'],x['s2'],x['x2'])
        r += x['drs2']             *  1./4  * self.K(x['s2'],x['x2'])
        r += x['drs2'] * x['k2']['esq'] * -1./6  * self.K21(x['s2'],x['x2'])
        r += x['drs2'] * x['p2']['esq'] *  1./6  * self.K21(x['s2'],x['x2'])

        return r

    def fkfpi_nnlo(self,x,lec):
        # parameters
        k2 = x['k2']['esq']
        p2 = x['p2']['esq']
        lk = np.log(k2)
        lp = np.log(p2)
        e2 = 4./3 * k2 - 1./3 * p2
        le = np.log(e2)

        # LECs
        L4  = lec['L4']
        L5  = lec['L5']
        k_4 = lec['k_4']
        p_4 = lec['p_4']
        if self.nnlo_full:
            L1 = lec['L1']
            L2 = lec['L2']
            L3 = lec['L3']
            L6 = lec['L6']
            L7 = lec['L7']
            L8 = lec['L8']

        # discretization terms
        a2 = x['a2']
        s_4 = lec['s_4']
        if self.alphaS:
            alphaS = x['alphaS']
            saS_4  = lec['saS_4']

        # nnlo result (r)
        r = 0.
        # counter terms
        r += (k2 - p2) * k2 * k_4
        r += (k2 - p2) * p2 * p_4
        '''
        SU(3) flavor makes (k2-p2)**2 redundant
        '''
        r += (k2 -p2) * a2 * s_4
        if self.alphaS:
            r += (k2 -p2) * a2 * alphaS * saS_4
        if self.switches['debug_nnlo_check']:
            tmp = (k2 - p2) * k2 * k_4 + (k2 - p2) * p2 * p_4
            print('\nct        = %f' %tmp)
        if not self.ct:
            # NNLO terms from 1711.11328, B. Ananthanarayan, Johan Bijnens, Samuel Friot, Shayan Ghosh
            # logSq terms
            if self.switches['debug_nnlo_check']: tmp = float(r)
            r += lp * lp   * ( 11./24   * p2 *k2  -131./192 * p2 *p2)
            r += lp * lk   * (-41./96   * p2 *k2  -3./32    * p2 *p2)
            r += lp * le   * ( 13./24   * p2 *k2  +59./96   * p2 *p2)
            r += lk * lk   * ( 17./36   * k2 *k2  +7./144   * p2 *k2)
            r += lk * le   * (-163./144 * k2 *k2  -67./288  * p2 *k2 + 3./32   * p2 *p2)
            r += le * le   * ( 241./288 * k2 *k2  -13./72   * p2 *k2 - 61./192 * p2 *p2)
            if self.switches['debug_nnlo_check']:
                print('logSq     = %f' %(r -tmp))
            # FF(mpi**2 / mK**2) term
            if self.switches['FF_approximate']:
                r += k2**2 * FF_approximate(p2/k2)
                if self.switches['debug_nnlo_check']:
                    print('FF_approx = %f' %(k2**2 * FF_approximate(p2/k2)))
            else:
                r += k2**2 * FF(p2/k2)
                if self.switches['debug_nnlo_check']:
                    print('FF        = %f' %(k2**2 * FF(p2/k2)))

        if self.nnlo_full:
            # counter term pieces from Li
            tk = 8*(4*pi)**2 *L5 *(8*L4 +3*L5 -16*L6 -8*L8) -2*L1 -L2 -1./18*L3 +4./3*L5 -16*L7 -8*L8
            tp = 8*(4*pi)**2 *L5 *(4*L4 +5*L5  -8*L6 -8*L8) -2*L1 -L2 -5./18*L3 -4./3*L5 +16*L7 +8*L8
            r += (4*pi)**2 *(k2-p2) *tk * k2
            r += (4*pi)**2 *(k2-p2) *tp * p2
            if self.switches['debug_nnlo_check']:
                tmp = (4*pi)**2 *(k2-p2) *tk * k2 + (4*pi)**2 *(k2-p2) *tp * p2
                print('extra ct  = %f' %tmp)

            # log terms
            C1  = -(7./9     + (4*pi)**2 * 11./2 *L5) *p2 *k2
            C1 += -(113./72  + (4*pi)**2 *(4*L1 +10*L2 +13./2*L3 -21./2*L5))*p2**2

            C2  =  (209./144 + (4*pi)**2 *3*L5)*p2*k2
            C2 +=  (53./96   + (4*pi)**2 *(4*L1 +10*L2 +5*L3 -5*L5))*k2**2

            C3  =  (13./18   + (4*pi)**2 *(8./3*L3 -2./3 *L5 -16*L7  -8*L8))* k2*k2
            C3 += -(4./9     + (4*pi)**2 *(4./3*L3 +25./6*L5 -32*L7 -16*L8))* p2*k2
            C3 +=  (19./288  + (4*pi)**2 *(1./6*L3 +11./6*L5 -16*L7  -8*L8))* p2*p2

            r  += C1 * lp
            r  += C2 * lk
            r  += C3 * le

            if self.switches['debug_nnlo_check']:
                tmp = C1 * lp + C2 * lk + C3 * le
                print('C1 = %s' %C1)
                print('C2 = %s' %C2)
                print('C3 = %s' %C3)
                print('log terms = %f' %tmp)


        # correct for renormalization scale defined by F_latt
        #         and for changing Fpi -> sqrt(Fpi FK) or FK
        if self.FF == 'PP':
            if self.switches['debug_nnlo_check']: tmp = float(r)
            r += 3./2 *(k2-p2) *(p2*lp + 0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4))
            if self.switches['debug_nnlo_check']:
                if 'ratio' in self.eft:
                    print('mu fix    = %f' %(r-tmp))
                else:
                    print('mu fix    = %f\n' %(r-tmp))
        else:
            dFKFpi_nlo  = 5./8 *p2*lp -1./4 *k2*lk -3./8 *e2*le +4*(4*pi)**2 * L5 * (k2-p2)
            if self.FF == 'PK':
                if self.switches['debug_nnlo_check']: tmp = float(r)
                r += 3./4 *(k2-p2) * (11./8 *p2*lp +5./4 *k2*lk +3./8 *e2*le )
                r += 3./4 *(k2-p2) * (-4)*(4*pi)**2* (p2*(2*L4 +L5) + k2*(4*L4 +L5))
                if self.switches['debug_nnlo_check']:
                    print('mu fix    = %f' %(r-tmp))
                    tmp = float(r)
                r += dFKFpi_nlo**2
                if self.switches['debug_nnlo_check']:
                    print('F fix     = %f\n' %(r-tmp))
            elif self.FF == 'KK':
                if self.switches['debug_nnlo_check']: tmp = float(r)
                r += 3./2 *(k2-p2) * (3./8 *p2*lp  +3./4 *k2*lk +3./8 *e2*le )
                r += 3./2 *(k2-p2) * (-4)*(4*pi)**2* (p2*L4 +k2*(2*L4 +L5))
                if self.switches['debug_nnlo_check']:
                    print('mu fix    = %f' %(r-tmp))
                    tmp = float(r)
                r += 2 * dFKFpi_nlo**2
                if self.switches['debug_nnlo_check']:
                    if 'ratio' in self.eft:
                        print('F fix     = %f' %(r-tmp))
                    else:
                        print('F fix     = %f\n' %(r-tmp))
            else:
                sys.exit('unrecognized F=F_{pi,K} option: '+self.FF)

        # correct if ratio functions are used
        if 'ratio' in self.eft:
            if self.switches['debug_nnlo_check']:
                tmp = float(r)
            r += (p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4)) \
                * (3./8*p2*lp +3./4*k2*lk +3./8*e2*le +4*(4*pi)**2 *(p2*L4 +k2*(2*L4+L5)))
            #r += -1.0*( p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4) )**2
            r += 0.5*( p2*lp +0.5*k2*lk -4*(4*pi)**2 *(p2*(L4+L5) +2*k2*L4) )**2
            if self.switches['debug_nnlo_check']:
                print('ratio fix = %f\n' %(r-tmp))

        return r

    def nnnlo(self,x,lec):
        k2 = x['k2']['esq']
        p2 = x['p2']['esq']
        a2 = x['a2']
        r  = (k2 - p2) * a2**2 * lec['s_6']
        #r += (k2 - p2)    * k2 * a2 * lec['sk_6']
        #r += (k2 - p2)    * p2 * a2 * lec['sp_6']
        if not self.a4:
            r += (k2 - p2)    * k2 * p2 * lec['kp_6']
            r += (k2 - p2)**2 * k2      * lec['k_6']
            r += (k2 - p2)**2 * p2      * lec['p_6']
            r += (k2 - p2)    * k2 * a2 * lec['sk_6']
            r += (k2 - p2)    * p2 * a2 * lec['sp_6']

        return r

    def nlo_a4(self,x,lec):
        k2 = x['k2']['esq']
        p2 = x['p2']['esq']
        a2 = x['a2']
        r  = (k2 - p2) * a2 * lec['s_4']
        r += (k2 - p2) * a2**2 * lec['s_6']
        if self.alphaS:
            r += (k2 - p2) * a2 * x['alphaS'] * lec['saS_4']
        return r

    def xpt(self,x,p,debug=False):
        r = dict()
        for e in x:
            x_par = {key:val for key,val in p[e].items if not isinstance(key,str)}
            lec   = {key:val for key,val in p[e].items if isinstance(key,str)}
            if debug:
                print(e,x_par)
                print(lec)
            r[e]  = 1.
            r[e] +=  self.FK_xpt_nlo(x_par,lec)
            r[e] += -self.Fpi_xpt_nlo(x_par,lec)
            if self.order in ['nnlo','nnnlo']:
                r[e] += self.fkfpi_nnlo(x_par,lec)
            if self.order in ['nnnlo'] or (self.a4 and self.order != 'nlo'):
                r[e] += self.nnnlo(x_par,lec)
            if self.order == 'nlo' and self.a4:
                r[e] += self.nlo_a4(x_par,lec)
        return r
    def xpt_ratio(self,x,p):
        r = dict()
        for e in x:
            x_par = self.x_p_ens[e]
            lec   = self.lec_ens[e]
            num  = 1. + self.FK_xpt_nlo(x_par,lec)
            den  = 1. + self.Fpi_xpt_nlo(x_par,lec)
            r[e] = num / den
            if self.order in ['nnlo','nnnlo']:
                r[e] += self.fkfpi_nnlo(x_par,lec)
            if self.order in ['nnnlo'] or self.a4:
                r[e] += self.nnnlo(x_par,lec)
        return r

    def fit_data(self):
        x = self.prune_x()
        y = self.prune_data()
        p = self.prune_priors()
        if self.switches['scipy']:
            fitter='scipy_least_squares'
        else:
            fitter='gsl_multifit'
        fit_func = self.make_fit_function(self.func_elements)
        self.fit = lsqfit.nonlinear_fit(data=(x,y),prior=p,fcn=fit_func,fitter=fitter, debug=True)

    def report_phys_point(self,debug=False):
        # copy the class objects we will change for constructing the physical point
        fit_ensembles = list(self.switches['ensembles_fit'])
        self.switches['ensembles_fit'] = ['phys']
        func_elements = list(self.func_elements)
        if 'ma' in self.func_elements[0]:
            func_elements[0] = func_elements[0].replace('ma','xpt')
            eft = self.eft
            self.eft = eft.replace('ma','xpt')
        fit_fv = self.fv
        if fit_fv:
            self.fv = False
            self.I  = self.I_IV
            self.dI = self.dI_IV
        # set physical point
        x_phys = dict()
        x_phys['phys'] = {k:np.inf for k in ['mpiL','mkL']}
        x_phys['phys']['alphaS'] = 0.

        if debug:
            print('DEBUG: self.phys')
            print(self.phys)
        y_phys = dict()
        y_phys['phys'] = self.phys['FK'] / self.phys['Fpi']
        y_phys['FLAG'] = self.phys['FKFPi_FLAG']

        p_phys = dict()
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p_phys[k] = 0.
        Lchi_phys = self.phys['Lchi_'+self.FF]
        p_phys[('phys','p2')] = self.phys['mpi']**2 / Lchi_phys**2
        p_phys[('phys','k2')] = self.phys['mk']**2  / Lchi_phys**2
        p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
        p_phys[('phys','a2')] = 0.

        for k in ['L5','L4','k_4','p_4','kp_6','k_6','p_6','L1','L2','L3','L6','L7','L8']:
            if k in self.fit.p:
                p_phys[k] = self.fit.p[k]
        # create new function elementals if need be
        fit_function_phys = self.make_fit_function(func_elements)
        self.phys_extrap = fit_function_phys(x_phys,p_phys)
        if self.switches['debug_phys']:
            print('DEBUG: SCALE', self.FF)
            print('DEBUG: Lchi =',Lchi_phys)
            print('DEBUG: mpi =',self.phys['mpi'])
            print('DEBUG: mK  =',self.phys['mk'])
            for k in p_phys:
                print('DEBUG:',k,p_phys[k])
        if self.switches['report_fit']:
            print(self.model)
            print('chi2/dof [dof] = %.2f [%d]    Q = %.2e    logGBF = %.3f' \
                %(self.fit.chi2/self.fit.dof,self.fit.dof,self.fit.Q,self.fit.logGBF))
            for k in self.fit.p:
                if isinstance(k,str):
                    print(k,self.fit.p[k])
            print(self.phys_extrap['phys'],'\n')
        # restore objects to original
        if 'ma' in self.func_elements[0]:
            self.eft = eft
        self.switches['ensembles_fit'] = fit_ensembles
        if fit_fv:
            self.fv = True
            self.I  = self.I_FV
            self.dI = self.dI_FV
        #self.func_elements = self.make_function_elementals()
        #self.fit_function = self.make_fit_function
        #for key,val in self_dict.items():
        #    setattr(self, key, val)
        return self.phys_extrap

    def check_fit(self,phys_point):
        # get copy of all self attributes
        self_dict = copy.deepcopy(self.__dict__)
        # turn off FV
        self.fv = False
        self.I  = self.I_IV
        self.dI = self.dI_IV
        # set physical point
        x_phys = dict()
        x_phys['phys'] = {k:np.inf for k in ['mpiL','mkL']}
        x_phys['phys']['alphaS'] = 0.

        y_phys = dict()
        y_phys['phys'] = phys_point['FK'] / phys_point['Fpi']
        y_phys['FLAG'] = phys_point['FKFPi_FLAG']

        p_phys = dict()
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p_phys[k] = 0.
        Lchi_phys = phys_point['Lchi']
        p_phys[('phys','p2')] = phys_point['mpi']**2 / Lchi_phys**2
        p_phys[('phys','k2')] = phys_point['mk']**2  / Lchi_phys**2
        p_phys[('phys','e2')] = 4./3*p_phys[('phys','k2')] - 1./3 * p_phys[('phys','p2')]
        p_phys[('phys','a2')] = 0.
        #for k in self.fit.p:
            #if isinstance(k,str):
                #print(k,self.fit.p[k])
        for k in ['L5','L4','k_4','p_4','kp_6','k_6','p_6']:
            if k in self.fit.p:
                p_phys[k] = self.fit.p[k]
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt

        #print('chi2/dof [dof] = %.2f [%d]    Q = %.2e    logGBF = %.3f' \
        #    %(self.fit.chi2/self.fit.dof,self.fit.dof,self.fit.Q,self.fit.logGBF))
        fkfpi = self.fit_function(x_phys,p_phys)

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return self.fit.chi2/self.fit.dof, self.fit.dof, self.fit.logGBF, self.fit.p['L5'], self.fit.p['k_4'], self.fit.p['p_4'], self.fit.p['s_4'], fkfpi['phys']

    def check_fit_function(self,check_point):
        x = dict()
        x['check'] = {k:np.inf for k in ['mpiL','mkL']}
        y = dict()
        y['check'] = self.y['check']
        p = dict()
        for k in self.p_init:
            p[k] = self.p_init[k]
        for k in ['s_4','saS_4','s_6','sk_6','sp_6']:
            p[k] = 0.
        Lchi = self.p_init[('check','Lchi_'+self.FF)]
        p[('check','p2')] = check_point['mpi']**2 / Lchi**2
        p[('check','k2')] = check_point['mk']**2  / Lchi**2
        p[('check','e2')] = 4.*p[('check','k2')]/3 - p[('check','p2')]/3
        p[('check','a2')] = 0.
        x_in,lec = self.make_x_lec(x,p,'check')
        print('\n----------------------------------\n'+self.model)
        if self.switches['FF_approximate']:
            print('  APPROXIMATE FF(mpi**2/mK**2)')
        #for k in p:
        #    print(k,p[k])
        print('eps2_p = %s' %x_in['p2']['esq'])
        print('eps2_k = %s' %x_in['k2']['esq'])
        print('eps2_e = %s' %x_in['e2']['esq'])
        if self.switches['debug_nlo_check']:
            print('FK_NLO terms')
        nlo_K = self.FK_xpt_nlo(x_in,lec)
        if self.switches['debug_nlo_check']:
            print('               -------------------------------')
            print('               FK_NLO = %s' %nlo_K)
            print('FP_NLO terms')
        nlo_P = self.Fpi_xpt_nlo(x_in,lec)
        if self.switches['debug_nlo_check']:
            print('               -------------------------------')
            print('               FP_NLO = %s' %nlo_P)
        if 'ratio' in self.eft:
            nlo = (1 + nlo_K) / (1 + nlo_P)
        else:
            nlo = 1 + nlo_K - nlo_P
        print('NNLO\n----------------------------------------------')
        nnlo = self.fkfpi_nnlo(x_in,lec)
        print('FK - Fpi | 1 + NLO = %f' %nlo)
        print('FK - Fpi | NNLO    = %f' %nnlo)
        print('FK/Fpi = %f' %(nlo + nnlo))

    def vs_epi(self,raw_data=False):
        print('making chiral plot')
        # get copy of all self attributes
        self_dict = copy.deepcopy(self.__dict__)
        # turn off FV
        self.fv = False
        self.I  = self.I_IV
        self.dI = self.dI_IV
        # switch to continuum fit func
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt
        # set x,y,p

        x = dict()
        p = dict()
        for k in self.lec_nnlo + self.lec_nnnlo:
            if k in self.fit.p:
                p[k] = self.fit.p[k]

        x_plot = []
        y_plot = {}
        y_plot['a00'] = []
        y_plot['a06'] = []
        y_plot['a09'] = []
        y_plot['a12'] = []
        y_plot['a15'] = []
        Lchi = self.phys['Lchi_'+self.FF]
        k2 = self.phys['mk']**2 / Lchi**2
        mpi_range = np.sqrt(np.arange(1,401**2,400**2/50))
        for mpi in mpi_range:
            print(mpi)
            x[mpi] = {k:np.inf for k in ['mpiL','mkL']}
            x[mpi]['alphaS'] = 0.
            p2 = mpi**2 / Lchi**2
            p[(mpi,'p2')] = p2
            p[(mpi,'k2')] = k2
            p[(mpi,'e2')] = 4./3*k2 - 1./3*p2
            p[(mpi,'a2')] = 0.
            x_plot.append(p[(mpi,'p2')].mean)
            y_plot['a00'].append(self.fit_function(x,p)[mpi])
            # finite a
            p[(mpi,'a2')] = self.p_init[('a06m310L','aw0')]**2 / (4 * pi)
            y_plot['a06'].append(self.fit_function(x,p)[mpi])
            p[(mpi,'a2')] = self.p_init[('a09m310','aw0')]**2 / (4 * pi)
            y_plot['a09'].append(self.fit_function(x,p)[mpi])
            p[(mpi,'a2')] = self.p_init[('a12m310','aw0')]**2 / (4 * pi)
            y_plot['a12'].append(self.fit_function(x,p)[mpi])
            p[(mpi,'a2')] = self.p_init[('a15m310','aw0')]**2 / (4 * pi)
            y_plot['a15'].append(self.fit_function(x,p)[mpi])

        x_plot = np.array(x_plot)
        y  = np.array([k.mean for k in y_plot['a00']])
        dy = np.array([k.sdev for k in y_plot['a00']])

        fig = plt.figure('FKFpi_vs_epi_'+self.model,figsize=figsize)
        self.ax_chiral = plt.axes(plt_axes)
        self.ax_chiral.fill_between(x_plot, y-dy, y+dy,color='#b36ae2',alpha=0.4)
        colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9', 'a06':'#00FFFF'}
        for a in ['a15','a12','a09','a06']:
            y = np.array([k.mean for k in y_plot[a]])
            ax.plot(x_plot, y, color=colors[a])
        # plot data
        self.plot_data('ea',offset='cont')
        if self.switches['plot_raw_data']:
            self.plot_data('ea', offset='cont', raw=True)
        # plot phys point
        epi_phys =  epi_range['mpi_phys']**2 / Lchi**2
        self.ax_chiral.axvline(epi_phys.mean,linestyle='--',color='#a6aaa9')
        self.ax_chiral.axvspan(epi_phys.mean -epi_phys.sdev, epi_phys.mean +epi_phys.sdev,
            alpha=0.4, color='#a6aaa9')
        self.ax_chiral.axvline(epi_phys.mean,linestyle='--',color='#a6aaa9')
        # axis labels, legends, etc
        self.ax_chiral.set_xlabel(r'$\epsilon_\pi^2 = (m_\pi / 4\pi F_\pi)^2$',fontsize=fs_text)
        self.ax_chiral.set_ylabel(r'$F_K / F_\pi$',fontsize=fs_text)
        self.ax_chiral.set_xlim(0,.1)
        if self.switches['plot_raw_data']:
            self.ax_chiral.set_ylim(1.06, 1.228)
        else:
            self.ax_chiral.set_ylim(1.14, 1.228)

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return ax

    def vs_ea(self,raw_data=False):
        print('making continuum plot')
        # get copy of all self attributes
        self_dict = copy.deepcopy(self.__dict__)
        # turn off FV
        self.fv = False
        self.I  = self.I_IV
        self.dI = self.dI_IV
        # switch to continuum fit func
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt

        # set x,y,p
        x = dict()
        p = dict()
        for k in self.lec_nnlo + self.lec_nnnlo:
            if k in self.fit.p:
                p[k] = self.fit.p[k]

        x_plot = []
        y_plot = dict() # we will make 1 or 3 different y_plots
        y_plot['m135'] = [] # physical pion mass
        if 's_6' in p: #
            y_plot['a2']   = [] # the a**2 terms
            y_plot['a4']   = [] # the a**4 terms
            p_a2 = dict(p)
            p_a4 = dict(p)
            p_a2['s_6'] = gv.gvar(0,0)
            p_a4['s_4'] = gv.gvar(0,0)
        p2 = self.phys['mpi']**2 / self.phys['Lchi_'+self.FF]**2
        k2 = self.phys['mk']**2  / self.phys['Lchi_'+self.FF]**2
        # make equal spaced a in a**2 range
        a_range = np.sqrt(np.arange(0, .16**2, .16**2 / 50))
        for a in a_range:
            x[a] = {k:np.inf for k in ['mpiL','mkL']}
            p[(a,'p2')] = p2
            p[(a,'k2')] = k2
            p[(a,'e2')] = 4./3*k2 - 1./3*p2
            p[(a,'a2')] = (a / self.phys['w0'])**2 / 4 / pi
            x_plot.append(p[(a,'a2')])
            y_plot['m135'].append(self.fit_function(x,p)[a])
            if 's_6' in p:
                p_a2[(a,'p2')] = p2
                p_a2[(a,'k2')] = k2
                p_a2[(a,'e2')] = 4./3*k2 - 1./3*p2
                p_a2[(a,'a2')] = (a / self.phys['w0'])**2 / 4 / pi
                p_a4[(a,'p2')] = p2
                p_a4[(a,'k2')] = k2
                p_a4[(a,'e2')] = 4./3*k2 - 1./3*p2
                p_a4[(a,'a2')] = (a / self.phys['w0'])**2 / 4 / pi
                # subtract a**4 to get c + a**2
                y_plot['a2'].append(y_plot['m135'][-1] - self.fit_function(x,p_a4)[a])
                # subtract a**2 to get c + a**4
                y_plot['a4'].append(y_plot['m135'][-1] - self.fit_function(x,p_a2)[a])
        x_plot = np.array(x_plot)
        y  = np.array([k.mean for k in y_plot['m135']])
        dy = np.array([k.sdev for k in y_plot['m135']])

        # continuum extrapolation figure
        if self.switches['milc_compare']:
            figsize = fig_size2
        else:
            figsize = fig_size
        fig = plt.figure('FKFpi_vs_ea_'+self.model,figsize=figsize)
        self.ax_cont = plt.axes(plt_axes)
        self.ax_cont.fill_between(x_plot, y-dy, y+dy, color='#b36ae2', alpha=0.4)
        # plot data
        self.plot_data('ea',offset='cont')
        if self.switches['plot_raw_data']:
            self.plot_data('ea', offset='cont', raw=True)
        handles, labels = self.ax_cont.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        self.ax_cont.legend(handles, labels, ncol=4, fontsize=fs_leg)
        self.ax_cont.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=fs_text)
        self.ax_cont.set_ylabel(r'$F_K / F_\pi$',fontsize=fs_text)
        self.ax_cont.set_xlim(0,.065)
        if self.switches['plot_raw_data']:
            self.ax_cont.set_ylim(1.06, 1.228)
        else:
            self.ax_cont.set_ylim(1.14, 1.228)
        if self.switches['save_figs']:
            plt.savefig('figures/vs_epasq_'+self.switches['ansatz']['model']+'.pdf',transparent=True)

        if 's_6' in p and self.switches['plot_asq_converg']:
            y_plot['a2'] = np.array(y_plot['a2'])
            y_plot['a4'] = np.array(y_plot['a4'])
            a2_vs_a0 = y_plot['a2'] / y_plot['m135'][0]
            a4_vs_a2 = y_plot['a4'] / y_plot['m135'][0]
            fig = plt.figure('FKFpi_asq_convergence_'+self.model,figsize=fig_size)
            self.ax_convergence = plt.axes(plt_axes)
            # N2LO disco
            y_tmp  = np.array([k.mean for k in a2_vs_a0])
            dy_tmp = np.array([k.sdev for k in a2_vs_a0])
            self.ax_convergence.fill_between(x_plot, y_tmp-dy_tmp, y_tmp+dy_tmp, color='k', alpha=0.2, label=r'N$^2$LO $a^2$ corrections')
            # N3LO disco
            #y_tmp  = np.array([k.mean for k in a4_vs_a2])
            #dy_tmp = np.array([k.sdev for k in a4_vs_a2])
            #self.ax_convergence.fill_between(x_plot, y_tmp-dy_tmp, y_tmp+dy_tmp, color='r', alpha=0.5, label=r'N$^3$LO $a^2 + a^4$ corrections')
            # N2LO + N3LO
            y_tmp  = np.array([k.mean for k in a2_vs_a0 + a4_vs_a2])
            dy_tmp = np.array([k.sdev for k in a2_vs_a0 + a4_vs_a2])
            self.ax_convergence.fill_between(x_plot, y_tmp-dy_tmp, y_tmp+dy_tmp, color='#b36ae2', alpha=0.4, label=r'N$^2$LO + N$^3$LO')
            self.ax_convergence.legend()
            self.ax_convergence.set_xlabel(r'$\epsilon_a^2 = a^2 / (4\pi w_0^2)$',fontsize=fs_text)
            self.ax_convergence.set_ylabel(r'ratio N$^n$LO / continuum',fontsize=fs_text)
            self.ax_convergence.set_xlim(0,.065)
            if self.switches['save_figs']:
                plt.savefig('figures/asq_convergence_'+self.switches['ansatz']['model']+'.pdf',transparent=True)
        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

    def plot_data(self, p_type, offset=False, raw=False):
        # p_type: ea  = continuum extrapolation
        #         epi = pion mass extrapolation
        colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9', 'a06':'#00FFFF'}
        shapes = {'m400':'h', 'm350':'p', 'm310':'s', 'm220':'^', 'm130':'o', 'm135':'*'}
        labels = {
            'a15m400':'', 'a15m350':'', 'a15m310':'a15', 'a15m220':'','a15m135XL':'',
            'a12m400':'', 'a12m350':'', 'a12m310':'a12', 'a12m220':'', 'a12m130':'',
            'a12m220L':'', 'a12m220S':'',
            'a09m400':'', 'a09m350':'', 'a09m310':'a09', 'a09m220':'', 'a09m135':'',
            'a06m310L':'a06',
            }
        dx_cont = {
            'a15m400'  :0.0016, 'a12m400':0.0016, 'a09m400':0.0016,
            'a15m350'  :0.0008, 'a12m350':0.0008, 'a09m350':0.0008,
            'a15m310'  :0.,     'a12m310':0.,     'a09m310':0.,     'a06m310L':0.,
            'a15m220'  :-0.0008,'a12m220':-0.0008,'a09m220':-0.0008,
            'a15m135XL':-.0016, 'a12m130':-0.0016,'a09m135':-0.0016,
            'a12m220L' :-0.0012,'a12m220S':-0.0004,
        }

        if p_type == 'ea':
            y_shift = self.shift_phys_mass()
        elif p_type == 'epi':
            y_shift = {e:0 for e in self.switches['ensembles_fit']}

        for e in self.switches['ensembles_fit']:
            c = colors[e.split('m')[0]]
            s = shapes['m'+e.split('m')[1][0:3]]
            x = self.fit.p[(e,'a2')]
            if raw:
                y = self.fit.y[e]
                mfc = 'None'
                label = None
            else:
                y = self.fit.y[e] + y_shift[e]
                mfc = c
                label = labels[e]
            if offset=='cont':
                dx = dx_cont[e]
            else:
                dx = 0
            self.ax_cont.errorbar(x=x.mean+dx, y=y.mean, xerr=x.sdev, yerr=y.sdev,
                marker=s,color=c,mfc=mfc,linestyle='None',label=label)

    def shift_phys_epsK(self,phys_point):
        '''
            Compute the shift in FK/FPi from the fit, to the point where
            m_K = m_K^phys
            if MA:
                m_ju = m_pi
                m_js = m_K^phys
                m_ru = m_K^phys
                m_rs = m_ss
                dju  = 0
                drs  = 0
            and apply this shift to the numerical data, so that it can be
            plotted with a fit that only depends upon eps_pi = m_pi**2 / L_chi**2

            1: compute FK/FPi from the particular fit
                y_fit

            2: computee FK/FPi from the corresponding XPT fit function + a**2 terms
                y_shift
            3: shift all the data by
                FK/FPi -> FK/FPi + (y_shift - y_fit).mean

            return a dictionary of shifts
        '''
        # get copy of all self attributes
        self_dict = copy.deepcopy(self.__dict__)
        # turn off FV
        self.fv = False
        self.I  = self.I_IV
        self.dI = self.dI_IV

        # determine fitted values of FK/FPi
        y_shift = dict()
        print('correcting FK/Fpi data for plot')
        print('for model %s' %self.switches['ansatz']['model'])
        x = self.prune_x()
        for e in x:
            Lchi = self.p_init[(e,'Lchi_'+self.FF)]
            p    = self.prune_priors()
            for k in self.lec_nnlo + self.lec_nnnlo:
                if k in self.fit.p:
                    p[k] = self.fit.p[k]
            y_shift[e] = -self.fit_function(x,p)[e]
        # switch to XPT + a**2 terms
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt
        if self.switches['verbose']:
            print("%9s %12s %12s" %('ensemble','data', 'shift'))
            print('-------- | ----------- | ----------')
        for e in x:
            x_tmp = dict()
            p_tmp = dict()
            x_tmp[e] = {k:np.inf for k in ['mpiL','mkL']}
            x_tmp[e]['alphaS'] = self.x[e]['alphaS']
            for k in self.lec_nnlo + self.lec_nnnlo:
                if k in self.fit.p:
                    p_tmp[k] = self.fit.p[k]
            Lchi = self.p_init[(e,'Lchi_'+self.FF)]
            p_tmp[(e,'p2')] = self.p_init[(e,'mpi')]**2 / Lchi**2
            p_tmp[(e,'k2')] = phys_point['mk']**2 / phys_point['Lchi']**2
            p_tmp[(e,'e2')] = 4./3 *p_tmp[(e,'k2')] - 1./3 *p_tmp[(e,'p2')]
            p_tmp[(e,'a2')] = self.p_init[(e,'aw0')]**2 / (4 * pi)
            y_shift[e] += self.fit_function(x_tmp,p_tmp)[e]
            if self.switches['verbose']:
                print("%9s %12s %12s" %(e, self.y[e], y_shift[e].mean))

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return y_shift

    def shift_phys_mass(self):
        '''
            Compute the shift in FK/FPi from the fit, to the point where
            m_K = m_K^phys
            if MA:
                m_ju = m_pi
                m_js = m_K^phys
                m_ru = m_K^phys
                m_rs = m_ss
                dju  = 0
                drs  = 0
            and apply this shift to the numerical data, so that it can be
            plotted with a fit that only depends upon eps_pi = m_pi**2 / L_chi**2

            1: compute FK/FPi from the particular fit
                y_fit

            2: computee FK/FPi from the corresponding XPT fit function + a**2 terms
                y_shift
            3: shift all the data by
                FK/FPi -> FK/FPi + (y_shift - y_fit).mean

            return a dictionary of shifts
        '''
        # get copy of all self attributes
        self_dict = copy.deepcopy(self.__dict__)

        # determine fitted values of FK/FPi
        y_shift = dict()
        print('correcting FK/Fpi data for plot')
        print('for model %s' %self.switches['ansatz']['model'])
        x = self.prune_x()
        # turn FV corrections back on temporarily
        self.fv = 'FV' in self.switches['ansatz']['model']
        if self.fv:
            self.I  = self.I_FV
            self.dI = self.dI_FV
        else:
            self.I  = self.I_IV
            self.dI = self.dI_IV

        PK = self.FF
        p    = self.prune_priors()
        for k in self.lec_nnlo + self.lec_nnnlo:
            if k in self.fit.p:
                p[k] = self.fit.p[k]
        for e in x:
            y_shift[e] = -self.fit_function(x,p)[e]
        # switch to XPT + a**2 terms
        if 'ratio' in self.eft:
            self.eft          = 'xpt-ratio'
            self.fit_function = self.xpt_ratio
        else:
            self.eft          = 'xpt'
            self.fit_function = self.xpt
        self.fv = False
        if self.switches['verbose']:
            print("%9s %12s %12s" %('ensemble','data', 'shift'))
            print('-------- | ----------- | ----------')
        for e in x:
            x_tmp = dict()
            p_tmp = dict()
            x_tmp[e] = {k:np.inf for k in ['mpiL','mkL']}
            x_tmp[e]['alphaS'] = self.x[e]['alphaS']
            for k in self.lec_nnlo + self.lec_nnnlo:
                if k in self.fit.p:
                    p_tmp[k] = self.fit.p[k]
            Lchi = self.p_init[(e,'Lchi_'+PK)]
            p_tmp[(e,'p2')] = self.phys['mpi']**2 / self.phys['Lchi_'+PK]**2
            p_tmp[(e,'k2')] = self.phys['mk']**2 / self.phys['Lchi_'+PK]**2
            p_tmp[(e,'e2')] = 4./3 *p_tmp[(e,'k2')] - 1./3 *p_tmp[(e,'p2')]
            p_tmp[(e,'a2')] = self.p_init[(e,'aw0')]**2 / (4 * pi)
            y_shift[e] += self.fit_function(x_tmp,p_tmp)[e]

            if self.switches['verbose']:
                print("%9s %12s %12s" %(e, self.y[e], y_shift[e].mean))

        # restore original self attributes
        for key,val in self_dict.items():
            setattr(self, key, val)

        return y_shift
