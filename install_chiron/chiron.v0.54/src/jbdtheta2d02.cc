// jbdtheta2d02.cc is part of the numerical library jbnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


/* uses modular invariance to speed up matters, J. Bijnens 7 September 2012 */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

const double pi2 = M_PI*M_PI;

double jbdtheta2d02(const double al, const double be, const double ga){
// Theta = sum_(n1,n2) n1^2 exp(- al n1^2-be n2^2 -ga (n1-n2)^2)
  if((al <= 0.) || (be <= 0.) || (ga <= 0.)){
    std::cout << " RiemannTheta2Real called with bad arguments\n";
  }
  double al1 = al;
  double be1=be;
  double ga1 = ga;
  double det = al1*be1+be1*ga1+ga1*al1;
  const double acc = 1e-18;// exponential should be larger than this
  const double acclog = -log(acc);// exponent should be smaller than this
  if(det >= pi2){
    // eigenvalues, we only need the smallest
    double eig2 = 0.5*(al1+be1+2.*ga1-sqrt(al1*al1+be1*be1+4.*ga1*ga1-2.*al1*be1));

    int nmax = sqrt(acclog/eig2);
    if (nmax==0) nmax = 1;// make sure at least a bit calculated
    //std::cout<<"nmax: "<<nmax<<std::endl;
    // now we simply run over the rectangle, even though
    // this wastes a lot (case with both n_i signs flipped is the factor of 2)
    al1 = -al1;
    be1 = -be1;
    ga1 = -ga1;
    double result = 0.;
    // neither zero
    for(int n1= -nmax; n1 < 0; n1++){
      for(int n2= -nmax; n2 < 0; n2++){
	double expo1 = n1*n1*al1+n2*n2*be1+(n1-n2)*(n1-n2)*ga1;
	double expo2 = n1*n1*al1+n2*n2*be1+(n1+n2)*(n1+n2)*ga1;
	result += 2.*n1*n1*(exp(expo1)+exp(expo2));
      }}
    // one zero
    for(int n1= -nmax; n1 < 0; n1++){
      double expo1 = n1*n1*(al1+ga1);
      result += 2.*n1*n1*exp(expo1);
    }
    // both zero not needed
    return result;
  }
 double ovfac = -1./(M_PI*sqrt(det));
 al1 = pi2/det*be;
 be1 = pi2/det*al;
 ga1 = pi2/det*ga;
 //std::cout << "went via inverse\n";
 
 // eigenvalues, we only need the smallest
 double eig2 = 0.5*(al1+be1+2.*ga1-sqrt(al1*al1+be1*be1+4.*ga1*ga1-2.*al1*be1));
 int nmax = sqrt(acclog/eig2);
 if (nmax == 0) nmax = 1;
 // std::cout<<"nmax: "<<nmax<<std::endl;
 // now we simply run over the rectangle, even though
 // this wastes a lot (case with both n_i signs flipped is the factor of 2)
 al1 = -al1;
 be1 = -be1;
 ga1 = -ga1;
 double bega1 = al1+ga1;
 double result = 0.;
 // neither zero
 for(int n1= -nmax; n1 < 0; n1++){
 for(int n2= -nmax; n2 < 0; n2++){
   double expo1 = n1*n1*al1+n2*n2*be1+(n1-n2)*(n1-n2)*ga1;
   double expo2 = n1*n1*al1+n2*n2*be1+(n1+n2)*(n1+n2)*ga1;
   result += 2.*( (0.5*bega1+bega1*expo1-pi2*pi2/det*n2*n2)*exp(expo1)
		 +(0.5*bega1+bega1*expo2-pi2*pi2/det*n2*n2)*exp(expo2));

   }}
  // one zero
  for(int n1= -nmax; n1 < 0; n1++){
    double expo1 = n1*n1*(al1+ga1);
    double expo2 = n1*n1*(be1+ga1);
    result += 2.*( (0.5*bega1+bega1*expo1)*exp(expo1)
		  +(0.5*bega1+bega1*expo2-pi2*pi2/det*n1*n1)*exp(expo2));
  }
  // both zero
  result += (0.5*bega1);

  return result*ovfac;
}
