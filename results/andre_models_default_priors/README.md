## Model Description
Using default priors w/ Bijnen's central values:
```python
prior = {}
prior.update(prior_bijnen_central) # as determined by Bijnen et al, using central value (ie, not 0)

# nlo Gasser-Leutwyler LECs
prior['L_4'] = gv.gvar('0.0(0.005)')
prior['L_5'] = gv.gvar('0.0(0.005)')

# nlo polynomial
prior['A_x'] = gv.gvar('0.0(1.0)')

# Lattice spacing terms
prior['A_a']    = gv.gvar('0.0(100.0)')
prior['A_loga'] = gv.gvar('0.0(100.0)')
prior['A_aa']   = gv.gvar('0.0(1000.0)')
prior['A_aaa']  = gv.gvar('0.0(10000.0)')

# n2lo terms
prior['A_k'] = gv.gvar('0.0(10.0)')#gv.gvar('0.0(5.0)')
prior['A_p'] = gv.gvar('0.0(10.0)')#gv.gvar('0.0(5.0)')

# n3lo terms
prior['A_ak'] = gv.gvar('0.0(100.0)')
prior['A_ap'] = gv.gvar('0.0(100.0)')
prior['A_kk'] = gv.gvar('0.0(100.0)')
prior['A_kp'] = gv.gvar('0.0(100.0)')
prior['A_pp'] = gv.gvar('0.0(100.0)')
```

Models generated by iterating over the following choices:
```python
choices = []

# n2lo fits
choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n2lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False, True],
    'semi-n2lo_corrections' : ['n2lo-ct', 'n2lo-full_bijnens'],

    # n3lo corrections
    'include_latt_n3lo' : [False, True],
    'include_latt_n4lo' : [False],
})

#n3lo fits w/o n4lo latt spacing term
choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n3lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False, True],
    'semi-n2lo_corrections' : ['n2lo-full_bijnens'],

    # n3lo corrections
    'include_latt_n3lo' : [True],
    'include_latt_n4lo' : [False],
})

# n3lo fits w/ n4lo latt spacing term
choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n3lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False],
    'semi-n2lo_corrections' : ['n2lo-full_bijnens'],

    # n3lo corrections
    'include_latt_n3lo' : [False],
    'include_latt_n4lo' : [True],
})

```

## Model Average
```yaml
FK/Fpi_pm: 1.1982(78)
[FLAG:     1.1932(19)]

---
FK/Fpi:        1.2009(80) 
delta_su(2):  -0.00439(73) 

---
Uncertainty: 
   Unexplained:  0.00474 
   Explained:    0.00639 

---
Error Budget: 
   Chiral:       0.00059 
   Phys Point:   0.00141 
   Statistical:  0.00415 

---
Highest Weight: 
   0.206:  xpt_FpiFpi_n2lo_log_logSq_sunset_a4_FV_bijnens
   0.193:  xpt_FKFpi_n2lo_log_logSq_sunset_a4_FV_bijnens
   0.160:  xpt_FKFK_n2lo_a4_FV
   0.073:  xpt_FKFpi_n2lo_log_logSq_sunset_alphaS_FV_bijnens
   0.073:  xpt_FpiFpi_n2lo_log_logSq_sunset_alphaS_FV_bijnens

---

Model: xpt_FpiFpi_n2lo_log_logSq_sunset_a4_FV_bijnens

Fitted/[FLAG] values at physical point:
	F_K / F_pi = 1.2028(33) [1.1959(20)]	(delta_su2 = -0.00444(71))

Parameters:
            L_4    -0.0002 (29)      [   0.0000 (50) ]  
            L_5   -0.00042 (53)      [   0.0000 (50) ]  
            L_2    0.00032 (16)      [  0.00033 (16) ]  
            L_3    -0.0036 (15)      [  -0.0031 (15) ]  
            L_7   -0.00036 (17)      [ -0.00034 (17) ]  
            A_k      -0.9 (7.3)      [        0 (10) ]  
            A_p       2.3 (3.8)      [        0 (10) ]  
            A_a      -8.8 (1.1)      [       0 (100) ]  
           A_aa         69 (14)      [      0 (1000) ]  

Least Square Fit:
  chi2/dof [dof] = 0.62 [18]    Q = 0.89    logGBF = 68.687

Settings:
  svdcut/n = 1e-12/0    tol = (1e-08*,1e-10,1e-10)    (itns/time = 48/24.7)

Error Budget (relative error):
  disc       3.8e-05
  chiral     3.4e-04
  pp_input   1.1e-03
  stat       2.5e-03


---

Model: xpt_FKFpi_n2lo_log_logSq_sunset_a4_FV_bijnens

Fitted/[FLAG] values at physical point:
	F_K / F_pi = 1.1984(31) [1.1958(19)]	(delta_su2 = -0.00432(70))

Parameters:
            L_4    -0.0005 (21)      [   0.0000 (50) ]  
            L_5   -0.00009 (54)      [   0.0000 (50) ]  
            L_1    0.00023 (12)      [  0.00024 (12) ]  
            L_3    -0.0035 (15)      [  -0.0031 (15) ]  
            L_7   -0.00035 (17)      [ -0.00034 (17) ]  
            A_k      -0.5 (8.2)      [        0 (10) ]  
            A_p       1.9 (4.1)      [        0 (10) ]  
            A_a      -9.7 (1.2)      [       0 (100) ]  
           A_aa         74 (16)      [      0 (1000) ]  

Least Square Fit:
  chi2/dof [dof] = 0.64 [18]    Q = 0.87    logGBF = 68.621

Settings:
  svdcut/n = 1e-12/0    tol = (1e-08*,1e-10,1e-10)    (itns/time = 30/16.6)

Error Budget (relative error):
  disc       4.2e-05
  chiral     4.7e-04
  pp_input   9.1e-04
  stat       2.4e-03


---

Model: xpt_FKFK_n2lo_a4_FV

Fitted/[FLAG] values at physical point:
	F_K / F_pi = 1.1936(38) [1.1957(19)]	(delta_su2 = -0.00418(70))

Parameters:
            L_4    -0.0013 (35)      [   0.0000 (50) ]  
            L_5    0.00027 (86)      [   0.0000 (50) ]  
            A_k       2.5 (6.9)      [        0 (10) ]  
            A_p       0.6 (3.7)      [        0 (10) ]  
            A_a     -11.1 (1.5)      [       0 (100) ]  
           A_aa         81 (23)      [      0 (1000) ]  

Least Square Fit:
  chi2/dof [dof] = 0.84 [18]    Q = 0.66    logGBF = 68.433

Settings:
  svdcut/n = 1e-12/0    tol = (1e-08*,1e-10,1e-10)    (itns/time = 21/1.3)

Error Budget (relative error):
  disc       7.7e-05
  chiral     6.0e-04
  pp_input   1.2e-03
  stat       2.9e-03


---

Model: xpt_FKFpi_n2lo_log_logSq_sunset_alphaS_FV_bijnens

Fitted/[FLAG] values at physical point:
	F_K / F_pi = 1.2091(57) [1.1960(20)]	(delta_su2 = -0.00462(72))

Parameters:
         A_loga      36.3 (8.5)      [       0 (100) ]  
            L_4    -0.0006 (19)      [   0.0000 (50) ]  
            L_5   -0.00037 (52)      [   0.0000 (50) ]  
            L_1    0.00023 (12)      [  0.00024 (12) ]  
            L_2    0.00021 (11)      [  0.00022 (11) ]  
            L_3    -0.0038 (15)      [  -0.0031 (15) ]  
            L_7   -0.00037 (17)      [ -0.00034 (17) ]  
            A_k    -0.08 (8.29)      [        0 (10) ]  
            A_p       1.3 (4.2)      [        0 (10) ]  
            A_a     -27.9 (5.6)      [       0 (100) ]  

Least Square Fit:
  chi2/dof [dof] = 0.95 [18]    Q = 0.51    logGBF = 67.657

Settings:
  svdcut/n = 1e-12/0    tol = (1e-08*,1e-10,1e-10)    (itns/time = 27/14.9)

Error Budget (relative error):
  disc       4.4e-04
  chiral     6.7e-04
  pp_input   1.2e-03
  stat       4.5e-03


---

Model: xpt_FpiFpi_n2lo_log_logSq_sunset_alphaS_FV_bijnens

Fitted/[FLAG] values at physical point:
	F_K / F_pi = 1.2149(59) [1.1961(20)]	(delta_su2 = -0.00478(73))

Parameters:
         A_loga      33.6 (7.3)      [       0 (100) ]  
            L_4    -0.0003 (26)      [   0.0000 (50) ]  
            L_5   -0.00069 (52)      [   0.0000 (50) ]  
            L_1    0.00028 (14)      [  0.00029 (14) ]  
            L_2    0.00031 (16)      [  0.00033 (16) ]  
            L_3    -0.0040 (15)      [  -0.0031 (15) ]  
            L_7   -0.00038 (17)      [ -0.00034 (17) ]  
            L_8    0.00019 (10)      [  0.00020 (10) ]  
            A_k      -0.5 (7.7)      [        0 (10) ]  
            A_p       1.8 (4.0)      [        0 (10) ]  
            A_a     -25.5 (4.8)      [       0 (100) ]  

Least Square Fit:
  chi2/dof [dof] = 0.98 [18]    Q = 0.49    logGBF = 67.648

Settings:
  svdcut/n = 1e-12/0    tol = (1e-08*,1e-10,1e-10)    (itns/time = 42/22.4)

Error Budget (relative error):
  disc       3.9e-04
  chiral     5.6e-04
  pp_input   1.6e-03
  stat       4.5e-03

```