// jbdsing15.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* similar to dcauch in cernlib but with changed division strategy
   and 15 point gauss kronrod rule 5/1/2015 */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

/* integration around a singularity,
    a,b boundaries, s place of the singularity,
   checks the naive cancellation around the singularity */

double jbdsing15(double (*f)(double),const double a,const double b,
		const double s,const double eps){
 
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

  const double z1 = 1., hf = z1/2., cst = 5.*z1/1000.;
  double b0,c,aa,bb,c1,c2,c3,c4,s8,s16,u,s0;
  double h=0.;
  int ndiv;
  if( (s==a) || (s==b)){
      h=0.;
      std::cout << "singularity at border in jbdsing15" << std::endl;
      return h;}
  else {
    if( (s< fmin(a,b)) || (s > fmax(a,b))){
      h= jbdquad15(f,a,b,eps);
      return h;}
    else{
      if(2.*s<(a+b)) {
        h= jbdquad15(f,2.*s-a,b,eps);
        b0=s-a;}
      else{
        h= jbdquad15(f,a,2.*s-b,eps);
        b0=b-s;}
      c=cst/b0;
      bb=0.;
      ndiv = 1;      
    label1:
      aa=bb;
  if(ndiv<=2) {
    bb=b0;
    ndiv = 1;}
  else{
    bb = aa+(b0-aa)/pow(2,ndiv-1);
    ndiv = ndiv-1;
  }
    label2:
      c1=hf*(bb+aa);
      c2=hf*(bb-aa);
      c3=s+c1;
      c4=s-c1;
      s0=((*f)(c3)+(*f)(c4));
      s8 = wg[3]*s0;
      s16 = wgk[7]*s0;
      for (int i=0; i<=2;i++){
	u=c2*xgk[2*i+1];
	s0=((*f)(c3+u)+(*f)(c4-u))+((*f)(c3-u)+(*f)(c4+u));
	s8=s8+wg[i]*s0;
	s16=s16+wgk[2*i+1]*s0;
      }
      s8=c2*s8;
      for (int i=0; i<=3;i++){
	u=c2*xgk[2*i];
	s16=s16+wgk[2*i]*(((*f)(c3+u)+(*f)(c4-u))+((*f)(c3-u)+(*f)(c4+u)));
      }
      s16=c2*s16;
      if(fabs(s16-s8)<(eps*(1.+fabs(s16)))) goto label5;
      ndiv += 1;
      bb=c1;
      if(1.+fabs(c*c2) != 1) goto label2;
      std::cout << "jbdsing15 too high accuracy required" << std::endl;
      goto label9;
    label5:
      h=h+s16;
      if(bb !=  b0) goto label1;
    }}
 label9:
  return h;
}
  
