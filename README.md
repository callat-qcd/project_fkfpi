# project_fkfpi

This repo performs the chiral, continuum, infinite volume extrapolation of the FK/Fpi results generated by CalLat.

While most of the code in `master` branch was written by André Walker-Loud, the design of the code is from Ben Hoerz.  Nolan Miller wrote a completely independent fitting code (`nolan` branch) and the two were cross checked.

After [Installation](##Installing), in principle, the only file the user should have to modify is `input_params.py` where many options are controlled, such as what fit model to use, which data to include, whether to make plots, whether to save and use fits, etc.  To run the fit, simply run

```
python fit_fkfpi.py
```
NOTE: if your `lsqfit` is not compiled against GSL, you have to set
```
switches['scipy']            = True
```
We installed with GSL support by specifying
```
pip install --no-cache-dir --global-option=build_ext --global-option="-lgsl" --global-option="-I/usr/local/include/gsl" --global-option="-L/usr/local/lib" gvar=11.2
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
