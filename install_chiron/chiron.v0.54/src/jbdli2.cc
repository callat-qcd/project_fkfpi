// jbdli2.cc is part of the numerical library jbnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Calculates the Li_2 function in the complex plane
// using the Bernouilly series see 't Hooft Veltman, Nucl. Phys. B153 (1979) 365
// uses up to $B_{28}$
#include <cmath>
#include <complex>

#include "jbnumlib.h"

typedef std::complex<double> dcomplex;

dcomplex jbdli2(const dcomplex xxx){
  dcomplex z,z2,terms[15],dli2te,dli2t;
  const dcomplex one = dcomplex(1.); 
  static double bern[15]; // takes bernouilly number b(2)->b(28):B_n/(n+1)!
  static int ibegin=0;
  static double bern0,bern1,pi2,fact;
  if (ibegin == 0){
    ibegin=1;
    bern0 = 1.;
    bern1 = -1./2./2.; // (last 2 is the 2!)
    bern[1] = 1./6.;
    bern[2] = -1./30.;
    bern[3] = 1./42.;
    bern[4] = -1./30.;
    bern[5] = 5./66.;
    bern[6] = -691./2730.;
    bern[7] = 7./6.;
    bern[8] = -3617./510.;
    bern[9] = 43867./798.;
    bern[10] = -174611./330.;
    bern[11] = 854513./138.;
    bern[12] = -236364091./2730.;
    bern[13] = 8553103./6.;
    bern[14] = -23749461029./870.;
    pi2 = M_PI*M_PI;
    fact = 1.;
    for(int i=1; i<= 14; ++i){
      fact = fact*double(2*i)*double(2*i+1);
      bern[i] = bern[i]/fact;}
  }  int ifac;
  if (abs(xxx-one)<=0.3){ //        ! x->1-x
    ifac=-1;
    dli2te = pi2/6.-log(xxx)*log(one-xxx);
    z = -log(xxx);}
  else{
    if (abs(xxx) > 2.){  //        x-> 1/x
      ifac=-1;
      dli2te = -pi2/6.-0.5*pow(log(-xxx),2);
      z = -log(one-one/xxx);}
    else{
      ifac=1;
      dli2te = dcomplex(0.);
      z = -log(one-xxx);}
  }
  z2=z*z;
  terms[1] = z2*z;
  for(int i=2;i<=14;++i){
    terms[i]=terms[i-1]*z2;}
  dli2t=bern[14]*terms[14];
  for (int i=13;i>=1;--i){
    dli2t = dli2t+bern[i]*terms[i];}
  dli2t = dli2t + bern1*z2 + bern0*z;
  return double(ifac)*dli2t+dli2te;
}
