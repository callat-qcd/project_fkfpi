# project_fkfpi

This repo performs the chiral, continuum, infinite volume extrapolation of the FK/Fpi results generated by CalLat using [MDWF valence fermions solved on gradient-flowed gauge ensembles with Nf=2+1+1 dynamical HISQ fermions](https://arxiv.org/abs/1701.07559).

While most of the code in `master` branch was written by André Walker-Loud ([walkloud](https://github.com/walkloud)), the design of the code is from Ben Hoerz ([ebatz](https://github.com/ebatz)).  Nolan Miller ([millernb](https://github.com/millernb)) wrote a completely independent fitting code (`nolan` branch) and the two were cross checked.  The code was derived from an original implementation by Jason Chang ([cchang5](https://github.com/cchang5)).

After [Installation](#installing), in principle, the only file the user should have to modify is `input_params.py` where many options are controlled, such as what fit model to use, which data to include, whether to make plots, whether to save and use fits, etc.  To run the fit, simply run

```
python fit_fkfpi.py
```
The loop over various choices of the fit "model" is controlled with the `switches['sys']` dictionary.  For example, to loop over the choice of using F=Fpi, F=FK or F^2 = Fpi FK, set
```
switches['sys']['Lam_chi'] = True
```
Other options for printing the lattice results in LaTeX cp/paste format, printing figures etc. are also controlled with the `switches` dictionary.

Further down the `input_params.py` is the setting of the priors used in the extrapolation analysis, followed by a definition of the physical point, `phys_point` for setting the final extrapolation value.  Also, values of the masses and decay constants and LECs are set to perform a fit-function check with `check_fit`.






`
NOTE: if your `lsqfit` is not compiled against GSL, you have to set
```
switches['scipy']            = True
```
We installed with GSL support by specifying
```
pip install --no-cache-dir --global-option=build_ext --global-option="-lgsl" --global-option="-I/usr/local/include/gsl" --global-option="-L/usr/local/lib" gvar==11.2
pip install --no-cache-dir --global-option=build_ext --global-option="-lgsl" --global-option="-I/usr/local/include/gsl" --global-option="-L/usr/local/lib" lsqfit==11.5.1
```



## Installation

### Required Packages
- standard scientific Python libraries (numpy, scipy, PyTables, matplotlib, yaml, functools)
- gcc/clang, we used
```
Apple clang version 11.0.0 (clang-1100.0.33.12)
Target: x86_64-apple-darwin19.0.0
```
- pybind11
- [gvar](https://github.com/gplepage/gvar) v >= 11.2: pip installable.
- [lsqfit](https://github.com/gplepage/lsqfit) v >= 11.5.1: pip installable
- CHIRON [http://home.thep.lu.se/~bijnens/chiron/](http://home.thep.lu.se/~bijnens/chiron/)

### Installing
From the root directory of project_fkfpi:
- [Install pybind11](https://anaconda.org/conda-forge/pybind11).  We used

```conda install -c conda-forge pybind11```

- Download and compile CHIRON.  We used `chiron.v0.54.tar.gz`

  - Follow the instructions in py_chiron/README for editing the Makefile of `chiron`,

  - then
  ```
  cd chiron.v0.54
  make [-j N]
  cd ../py_chiron
  ```
- Build the Python Binding to `CHIRON`
  - edit the `build_<VERION>.sh` script to point to necessary directories and then
  ```
  ./build_<VERSION>.sh
  ```
  If successful, this will have created a library such as `chiron.cpython-37m-darwin.so`.

- Test the binding
```
./test_binding.py
FF(0.7**2) =  0.01841657802571616
FF(0.2**2) =  0.5084416863714661
FF(gvar(0.2,0.1)**2) =  0.51(20)
FF(gvar(0.8,0.1)**2) =  0.0082(69)
```

# Copyright Notice
project_fkfpi Copyright (c) 2020, The 
Regents of the University of California, through Lawrence Berkeley 
National Laboratory (subject to receipt of any required approvals from 
the U.S. Dept. of Energy). All rights reserved.
If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department 
of Energy and the U.S. Government consequently retains certain rights.  As 
such, the U.S. Government has been granted for itself and others acting on 
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3979560.svg)](https://doi.org/10.5281/zenodo.3979560)
