// jbderiv2utheta3.cc is part of the numerical library jbnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

double jbderiv2utheta3(const double u,const double q)
  /*
     computes -8pi^2\sum_{n=1,\infty} n^2 q^(n^2) cos(2n\pi u)

     for q<= 0.5 simply sums the above series, up to n=10 terms

     for q >= 0.5 uses the deivative w.r.t. u of
                theta3(u,q) = sqrt(pi/|lnq|)*exp(-pi^2 u^2/|lnq|)
	*theta3(-iu pi/|lnq|,exp(-pi^2/|lnq|))

     precision works to machine precision (tested by comparing it against
     running with much higher nmaxu and nmaxl and comparing the two
     algorithms against each other
  */
{
  const int nmaxl = 10, nmaxu = 2; 
  double a[nmaxl+1];// should be the larger of nmaxl/nmaxu
  if( (q < 0.) || (q >=1.)){
    std::cout << "not allowed input value in jbderivtheta3" << std::endl;
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
    // Clenshaw resummation not used u small happens too often
    for (int i=nmaxl; i>0;i--){
      result += a[i]*double(i)*double(i)*cos(2.*M_PI*u*double(i));
    }
    return -8.*M_PI*M_PI*result;
  }
  // case >= 0.5
  double lam1 = M_PI/fabs(log(q));
  // everything put in the exponentials otherwise overflow
  double result = 0.;
  for (int i=nmaxu; i>0;i--){
    result += 
      ( 2.*(-1.+2.*pow(-u+double(i),2)*M_PI*lam1)*
	 exp(M_PI*lam1*( 2.*double(i)*u-u*u-double(i)*double(i)))
	+2.*(-1.+2.*pow(-u-double(i),2)*M_PI*lam1)*
	 exp(M_PI*lam1*(-2.*double(i)*u-u*u-double(i)*double(i))));
  }
  return sqrt(lam1)*lam1*M_PI*((-2.+4.*u*u*M_PI*lam1)*exp(-M_PI*u*u*lam1)
			       +result);// first term is n=0
}
