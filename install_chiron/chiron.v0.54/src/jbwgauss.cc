// jbwgauss.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Adaption of the CERNLIB routine  WGAUSS to C++, first done 25/4/2002

#include <cmath>
#include <complex>
#include <iostream>

typedef std::complex<double> dcomplex;

#include "jbnumlib.h"

dcomplex jbwgauss(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b, const double eps){
  const double x[13]={0.0
		  ,9.6028985649753623e-1
		  ,7.9666647741362674e-1
		  ,5.2553240991632899e-1
		  ,1.8343464249564980e-1
		  ,9.8940093499164993e-1
		  ,9.4457502307323258e-1
		  ,8.6563120238783174e-1
		  ,7.5540440835500303e-1
		  ,6.1787624440264375e-1
		  ,4.5801677765722739e-1
		  ,2.8160355077925891e-1
		  ,9.5012509837637440e-2};
  const double w[13]={0.0 
               ,1.0122853629037626e-1
               ,2.2238103445337447e-1
               ,3.1370664587788729e-1
               ,3.6268378337836198e-1
               ,2.7152459411754095e-2
               ,6.2253523938647893e-2
               ,9.5158511682492785e-2
               ,1.2462897125553387e-1
               ,1.4959598881657673e-1
               ,1.6915651939500254e-1
               ,1.8260341504492359e-1
	       ,1.8945061045506850e-1};
  dcomplex jbwgausst;
  dcomplex aa,bb,u,c1,c2,s8,s16;
  int i;
  /* start */
  jbwgausst = dcomplex(0.,0.);
  if (b == a) return jbwgausst;
  double cconst = 0.005/abs(b-a);
  bb=a;
  /* computational loop */
label1:
  aa=bb;
  bb=b;
label2:
  c1=0.5*(bb+aa);
  c2=0.5*(bb-aa);
  s8=0.0;
  for (i=1;i<=4;i++){
    u=c2*x[i];
    s8=s8+w[i]*((*f)(c1+u)+(*f)(c1-u));
  }
  s8 = c2*s8;
  s16=0.0;
  for (i=5;i<=12;i++){
    u=c2*x[i];
    s16=s16+w[i]*((*f)(c1+u)+(*f)(c1-u));
  }
  s16=c2*s16;
  if(abs(s16-s8) <= eps*(1.+abs(s16)) ) goto label5;
  bb=c1;
  if(1.+cconst*abs(c2) != 1.) goto label2;
  jbwgausst = 0.0;
  std::cout << "too high accuracy required, jbwgauss\n";
  return jbwgausst;
label5:
  jbwgausst=jbwgausst+s16;
  if(bb != b) goto label1;
  return jbwgausst;
}
