## Model Description
These fits are missing `a09m135`.

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
FK/Fpi_pm: 1.1955(80)
[FLAG:     1.1932(19)]

---
FK/Fpi:        1.1981(81)
delta_su(2):  -0.00431(73)

---
Uncertainty:
   Unexplained: 0.00534
   Explained:   0.00605

---
Error Budget:
   Chiral:      0.000754
   Phys Point:  0.00156
   Statistical: 0.00473

---
Highest Weight:
  0.241: xpt_FKFK_n2lo_a4_FV
  0.225: xpt_FKFK_n2lo_alphaS_FV
  0.092: xpt_FKFpi_n2lo_log_logSq_sunset_a4_FV_bijnens
  0.075: xpt_FKFK_n2lo_log_logSq_sunset_alphaS_FV_bijnens
  0.074: xpt_FpiFpi_n2lo_log_logSq_sunset_a4_FV_bijnens
```