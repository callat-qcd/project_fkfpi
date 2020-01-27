// jbdquad21.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// and the quadpack numbers which are public domain
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* J.~Bijnens, first programmed 21/12/2008
adaptive integration routine similar to jbdgauss2
subdivision strategy larger interval but not back to full interval
uses a 21 point Gauss-Kronrod rule with the data for weights and points
from the public domain quadpack
*/

#include <cmath>
#include <iostream>

#include "jbnumlib.h"

double jbdquad21(double (*f)(const double),const double a,const double b,
		const double eps){
  // uses a Gauss Kronrod rule with 21 points (21 = 2*5+ (2*5+1))
  // first Gauss tem for first estimate
  // Gauss-Kronrod with 21 for next estimate
  // points and weights taken from quadpack
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
  double result = 0.;
  double aa,bb,x,center,width,s1,s2,s0;
  int i,ndiv;
  // start
  if (b == a) return result;;
  double cconst = 0.005/fabs(b-a);
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
  s1 = 0.0;
  s2 = 0.0;
  for (i=0;i<n1;i++){
    x=width*xgk[2*i+1];
    s0 = ((*f)(center+x)+(*f)(center-x));
    s1 = s1+ wg[i]*s0;
    s2 = s2+wgk[2*i+1]*s0;
    x=width*xgk[2*i];
    s0 = ((*f)(center+x)+(*f)(center-x));
    s2 = s2+wgk[2*i]*s0;
  }
  x = width*xgk[2*n1];
  s2 = s2+wgk[2*n1]*(*f)(center+x);
  s1 = width*s1;
  s2 = width*s2;
  if (fabs(s2-s1) <= eps*(1.+fabs(s2)) ) goto label5;
  ndiv = ndiv+1;
  bb=center;
  if(1.+cconst*fabs(width) != 1.) goto label2;
  std::cout << "jbdquad21 too high accuracy required\n";
  return result;
label5:
  result= result+s2;
  if(bb != b) goto label1;
  return result;
}
