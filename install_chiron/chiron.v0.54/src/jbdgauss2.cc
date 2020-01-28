// jbdgauss2.cc  is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* adaption of CERNLIB dgauss fortran routine to C++ by J.~Bijnens,
  first done 22/12/2009. Changed subdivision strategy,
  better for high precision */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

double jbdgauss2(double (*f)(const double),const double a,const double b,
		const double eps)
  /*
     adaptive gaussian quadrature.

     jbdgauss2 is set equal to the approximate value of the integral of
     the function f over the interval (a,b), with accuracy parameter
     eps.

     eps is relative if the absolute value of the integral is larger
           than on one and absolute otherwise
     ******************************************************************
   */
{ const double x[13]={0.0
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
 double jbdgausst; 
 double aa,bb,c1,c2,s8,s16,cconst,u;
 int i,ndiv;
 /*  start. */
  jbdgausst=0.;
  if(b == a) return jbdgausst;
  cconst=0.005/(b-a);
  bb=a;
  ndiv = 1;
  /*  computational loop. */
label1: 
  aa=bb;
  if(ndiv<=2) {
    bb=b;
    ndiv = 1;}
  else{
    bb = aa+(b-aa)/pow(2,ndiv-1);
    ndiv = ndiv-1;
  }
label2:
  c1=0.5*(bb+aa);
  c2=0.5*(bb-aa);
  s8=0.0;
  for (i=1; i<=4;i++)
    { u=c2*x[i];
    s8=s8+w[i]*((*f)(c1+u)+(*f)(c1-u));
    }
  s8=c2*s8;
  s16=0.0;
  for (i=5;i<=12;i++)
    { u=c2*x[i];
    s16=s16+w[i]*((*f)(c1+u)+(*f)(c1-u));
    }
  s16=c2*s16;
  if( fabs(s16-s8)<= eps*(1.0+fabs(s16)) ) goto label5;
  ndiv += 1;
  bb=c1;
  if( 1.0+fabs(cconst*c2) != 1.0 ) goto label2;
  jbdgausst=0.0;
  std::cout << "jbdgauss2 too high accuracy required\n";
  return jbdgausst;
label5:
  jbdgausst=jbdgausst+s16;
  if(bb != b) goto label1;
  return jbdgausst;
}
