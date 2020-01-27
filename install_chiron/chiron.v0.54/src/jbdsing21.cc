// jbdsing21.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* similar to dcauch in cernlib but with changed division strategy
   and 21 point gauss kronrod rule 5/1/2015 */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

/* integration around a singularity, simple adaptation of the corresponding
   dcauch.for routine of cernlib, a,b boundaries, s place of the singularity,
   checks the naive cancellation around the singularity */

double jbdsing21(double (*f)(double),const double a,const double b,
		const double s,const double eps){
   const int n1 = 5;
  const double xgk[2*n1+1]={ 0.995657163025808080735527280689003e0
			    ,0.973906528517171720077964012084452e0
			    ,0.930157491355708226001207180059508e0
			    ,0.865063366688984510732096688423493e0
			    ,0.780817726586416897063717578345042e0
			    ,0.679409568299024406234327365114874e0
			    ,0.562757134668604683339000099272694e0
			    ,0.433395394129247190799265943165784e0
			    ,0.294392862701460198131126603103866e0
			    ,0.148874338981631210884826001129720e0
			    ,0.000000000000000000000000000000000e0};
  const double wg[n1]={ 0.066671344308688137593568809893332e0
		       ,0.149451349150580593145776339657697e0
		       ,0.219086362515982043995534934228163e0
		       ,0.269266719309996355091226921569469e0
		       ,0.295524224714752870173892994651338e0};
  const double wgk[2*n1+1]={ 0.011694638867371874278064396062192e0
			    ,0.032558162307964727478818972459390e0
			    ,0.054755896574351996031381300244580e0
			    ,0.075039674810919952767043140916190e0
			    ,0.093125454583697605535065465083366e0
			    ,0.109387158802297641899210590325805e0
			    ,0.123491976262065851077958109831074e0
			    ,0.134709217311473325928054001771707e0
			    ,0.142775938577060080797094273138717e0
			    ,0.147739104901338491374841515972068e0
			    ,0.149445554002916905664936468389821e0};

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
      h= jbdquad21(f,a,b,eps);
      return h;}
    else{
      if(2.*s<(a+b)) {
        h= jbdquad21(f,2.*s-a,b,eps);
        b0=s-a;}
      else{
        h= jbdquad21(f,a,2.*s-b,eps);
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
      s8 = 0;
      s16 = wgk[10]*s0;
      for (int i=0; i<5;i++){
	u=c2*xgk[2*i+1];
	s0=((*f)(c3+u)+(*f)(c4-u))+((*f)(c3-u)+(*f)(c4+u));
	s8=s8+wg[i]*s0;
	s16=s16+wgk[2*i+1]*s0;
	u=c2*xgk[2*i];
	s16=s16+wgk[2*i]*(((*f)(c3+u)+(*f)(c4-u))+((*f)(c3-u)+(*f)(c4+u)));
      }
      s8=c2*s8;
      s16=c2*s16;
      if(fabs(s16-s8)<(eps*(1.+fabs(s16)))) goto label5;
      ndiv += 1;
      bb=c1;
      if(1.+fabs(c*c2) != 1) goto label2;
      std::cout << "jbdsing21 too high accuracy required" << std::endl;
      goto label9;
    label5:
      h=h+s16;
      if(bb !=  b0) goto label1;
    }}
 label9:
  return h;
}
  
