// jbdcauch.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* adaption of CERNLIB dcauch fortran program to C++ by J.~Bijnens, 
   first done 21/2/2006 */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

/* integration around a singularity, simple adaptation of the corresponding
   dcauch.for routine of cernlib, a,b boundaries, s place of the singularity,
   checks the naive cancellation around the singularity */

double jbdcauch(double (*f)(double),const double a,const double b,
		const double s,const double eps){
 
  const double x[13]={0.0,
              9.6028985649753623e-1,
              7.9666647741362674e-1,
              5.2553240991632899e-1,
              1.8343464249564980e-1,
              9.8940093499164993e-1,
              9.4457502307323258e-1,
              8.6563120238783174e-1,
              7.5540440835500303e-1,
              6.1787624440264375e-1,
              4.5801677765722739e-1,
              2.8160355077925891e-1,
		9.5012509837637440e-2};

  const double w[13]={0.0,
              1.0122853629037626e-1,
              2.2238103445337447e-1,
              3.1370664587788729e-1,
              3.6268378337836198e-1,
              2.7152459411754095e-2,
              6.2253523938647893e-2,
              9.5158511682492785e-2,
              1.2462897125553387e-1,
              1.4959598881657673e-1,
              1.6915651939500254e-1,
              1.8260341504492359e-1,
		1.8945061045506850e-1};
 
  const double z1 = 1., hf = z1/2., cst = 5.*z1/1000.;
  double b0,c,aa,bb,c1,c2,c3,c4,s8,s16,u;
  double h=0.;
  if( (s==a) || (s==b)){
      h=0.;
      std::cout << "singularity at border in jbdcauch" << std::endl;
      return h;}
  else {
    if( (s< fmin(a,b)) || (s > fmax(a,b))){
      h= jbdgauss(f,a,b,eps);
    return h;}
    else{
      if(2*s<(a+b)) {
        h= jbdgauss(f,2.*s-a,b,eps);
        b0=s-a;}
      else{
        h= jbdgauss(f,a,2.*s-b,eps);
        b0=b-s;}
      c=cst/b0;
      bb=0;
      
    label1:
      aa=bb;
      bb=b0;
    label2:
      c1=hf*(bb+aa);
      c2=hf*(bb-aa);
      c3=s+c1;
      c4=s-c1;
      s8=0;
      for (int i=0; i<=4;i++){
	u=c2*x[i];
	s8=s8+w[i]*((f(c3+u)+f(c4-u))+(f(c3-u)+f(c4+u)));
      }
      s8=c2*s8;
      s16=0.;
      for (int i=5; i<=12;i++){
	u=c2*x[i];
	s16=s16+w[i]*((f(c3+u)+f(c4-u))+(f(c3-u)+f(c4+u)));
      }
      s16=c2*s16;
      if(fabs(s16-s8)<(eps*(1.+fabs(s16)))) goto label5;
      bb=c1;
      if(1.+fabs(c*c2) != 1) goto label2;
      h=0.;
      std::cout << "jbdcauch too high accuracy required" << std::endl;
      goto label9;
    label5:
      h=h+s16;
      if(bb !=  b0) goto label1;
    }}
 label9:
  return h;
}
  
