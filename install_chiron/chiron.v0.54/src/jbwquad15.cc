// jbwquad15.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// and the quadpack numbers which are public domain
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* uses a 15 point Gauss-Kronrod rule with the data for weights and points
from the public domain quadpack
*/

#include <cmath>
#include <iostream>

#include "jbnumlib.h"

dcomplex jbwquad15(dcomplex (*f)(const dcomplex),const dcomplex a,
		   const dcomplex b, const double eps){
  // uses a Gauss Kronrod rule with 15 points (21 = (2*3+1)+(2*3+1+1))
  // Gauss term for first estimate
  // Gauss-Kronrod with 15 for next estimate
  // points and weights from quadpack

  const double xgk[8]={ 0.991455371120812639206854697526329e0
		    ,0.949107912342758524526189684047851e0
		    ,0.864864423359769072789712788640926e0
		    ,0.741531185599394439863864773280788e0
		    ,0.586087235467691130294144838258730e0
		    ,0.405845151377397166906606412076961e0
		    ,0.207784955007898467600689403773245e0
		    ,0.000000000000000000000000000000000e0};
  const double wg[4]={ 0.129484966168869693270611432679082e0
		       ,0.279705391489276667901467771423780e0
		       ,0.381830050505118944950369775488975e0
		       ,0.417959183673469387755102040816327e0};
  const double wgk[8]={ 0.022935322010529224963732008058970e0
			  ,0.063092092629978553290700663189204e0
			  ,0.104790010322250183839876322541518e0
			  ,0.140653259715525918745189590510238e0
			  ,0.169004726639267902826583426598550e0
			  ,0.190350578064785409913256402421014e0
			  ,0.204432940075298892414161999234649e0
			  ,0.209482141084727828012999174891714e0};
  dcomplex result;
  dcomplex aa,bb,x,center,width,s1,s2,s0;
  int i,ndiv;
  // start
  result = dcomplex(0.);
  if (b == a) return result;
  double cconst = 0.005/abs(b-a);
  /* computational loop: does two things, divides constantly by half the size
 of the interval and when accuracy reached, starts again with the remaining
 part divided by 2^(n-1) when n was the number of subdivisions needed
 for the previous interval */
  bb=a;
  ndiv = 1;
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
  center=0.5*(bb+aa);
  width=0.5*(bb-aa);
  s0 = (*f)(center);
  s1 = wg[3]*s0;
  s2 = wgk[7]*s0;
  for (i=0;i<=2;i++){
    x=width*xgk[2*i+1];
    s0 = ((*f)(center+x)+(*f)(center-x));
    s1 = s1+ wg[i]*s0;
    s2 = s2+wgk[2*i+1]*s0;
  }
  for (i=0;i<=3;i++){
    x=width*xgk[2*i];
    s0 = ((*f)(center+x)+(*f)(center-x));
    s2 = s2+wgk[2*i]*s0;
  }
  s1 = width*s1;
  s2 = width*s2;
  if (abs(s2-s1) <= eps*(1.+abs(s2)) ) goto label5;
  ndiv = ndiv+1;
  bb=center;
  if(1.+cconst*abs(width) != 1.) goto label2;
  std::cout << "jbwquad15 too high accuracy required\n";
  return result;
label5:
  result = result+s2;
  if(bb != b) goto label1;
  return result;
}
