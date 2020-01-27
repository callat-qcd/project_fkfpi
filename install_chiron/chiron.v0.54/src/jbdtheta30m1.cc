// jbdtheta30m1.cc is part of the numerical library jbnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


/* uses idea behind cernlib dtheta.for, first done
    J. Bijnens 7 September 2012
   the one is subtracted here to get more accuracy for small input values */
#include <cmath>
#include <iostream>
#include "jbnumlib.h"

double jbdtheta30m1(const double q)
  /*
     computes 2\sum_{n=1,\infty} q^(n^2)

     for q<= 0.5 simply sums the above series, up to n=10 terms

     for q >= 0.5 uses
        sqrt(lambda/pi) (1+2\sum_{n=1,\infty} exp(-lambda n^2))
        up to n=2

     precision works to machine precision (tested by comparing it against
     running with much higher nmaxu and nmaxl and comparing the two
     algorithms against each other
  */
{
  const int nmaxl = 10, nmaxu = 2; 
  double a[nmaxl+1];// should be the larger of nmaxl/nmaxu
  if( (q < 0.) || (q >=1.)){
    std::cout << "not allowed input value in jbdtheta30" << std::endl;
    return 0.;}
  if( q < 0.5 ){
    double qsq = q*q;
    double atemp = 1.;
    double p = q;
    for( int i = 0; i<= nmaxl; i++){
      a[i] = atemp;
      atemp *= p;
      p *= qsq;
    }
    double result = 0.;
    for (int i=nmaxl; i>0;i--){
      result += a[i];
    }
    return 2.*result;
  }
  // case >= 0.5
  double lam1 = M_PI/fabs(log(q));
  double lam2 = exp(-M_PI*lam1);
  double lamsq = lam2*lam2;
    double atemp = 1.;
    double p = lam2;
    for( int i = 0; i<= nmaxu; i++){
      a[i] = atemp;
      atemp *= p;
      p *= lamsq;
    }
    double result = 0.;
    for (int i=nmaxu; i>0;i--){
      result += a[i];
    }
    return sqrt(lam1)*(1.+2.*result)-1.;
}
