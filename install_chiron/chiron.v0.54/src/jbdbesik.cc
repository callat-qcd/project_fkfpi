// jbdbesik.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

/* naive port to C++ for cernlib dbesio, first done
   J. Bijnens 7 September 2012 */

#include <cmath>
#include <iostream>
#include "jbnumlib.h"

double jbdbesi0(double x){
  const double eps = 1e-15, z1 = 1., hf = z1/2., pi = M_PI,
    ce = 0.57721566490153286e0, pih = pi/2., rpih = 2./pi, rpi2 = 1/(2.*pi);
  const double ci[25]={
     +1.008279205458740032e0
    ,+0.008445122624920943e0
    ,+0.000172700630777567e0
    ,+0.000007247591099959e0
    ,+0.000000513587726878e0
    ,+0.000000056816965808e0
    ,+0.000000008513091223e0
    ,+0.000000001238425364e0
    ,+0.000000000029801672e0
    ,-0.000000000078956698e0
    ,-0.000000000033127128e0
    ,-0.000000000004497339e0
    ,+0.000000000001799790e0
    ,+0.000000000000965748e0
    ,+0.000000000000038604e0
    ,-0.000000000000104039e0
    ,-0.000000000000023950e0
    ,+0.000000000000009554e0
    ,+0.000000000000004443e0
    ,-0.000000000000000859e0
    ,-0.000000000000000709e0
    ,+0.000000000000000087e0
    ,+0.000000000000000112e0
    ,-0.000000000000000012e0
    ,-0.000000000000000018e0};

  double v=fabs(x);
  if (v < 8.){
    double y= pow(hf*v,2);
    double xl=2.;// (nu+2)
    double a0=1.;
    double a1=1.+2.*y/((xl+1.)*(xl-1.));
    double a2=1.+y*(4.+3.*y/((xl+2.)*xl))/((xl+3.)*(xl-1.));
    double b0=1.;
    double b1=1.-y/(xl+1.);
    double b2=1.-y*(1.-y/(2.*(xl+2.)))/(xl+3.);
    double w1=3.+xl;
    double v1=3.-xl;
    double v3=xl-1.;
    double v2=v3+v3;
    double c=0.;
    for (int n=3; n<=30;n++){
      double c0=c;
      double fn=double(n);
      w1=w1+2.;
      double w2=w1-1.;
      double w3=w2-1.;
      double w4=w3-1.;
      double w5=w4-1.;
      double w6=w5-1.;
      v1=v1+1.;
      v2=v2+1.;
      v3=v3+1.;
      double u1=fn*w4;
      double e=v3/(u1*w3);
      double u2=e*y;
      double f1=1.+y*v1/(u1*w1);
      double f2=(1.+y*v2/(v3*w2*w5))*u2;
      double f3=-y*y*u2/(w4*w5*w5*w6);
      double a=f1*a2+f2*a1+f3*a0;
      double b=f1*b2+f2*b1+f3*b0;
      c=a/b;
      if(fabs(c0-c) < (eps*fabs(c))) break;
      a0=a1;
      a1=a2;
      a2=a;
      b0=b1;
      b1=b2;
      b2=b;
    }
    //h = c;
    //       if(nu .eq. 1) h=hf*x*h
    //    if(lex) h=exp(-v)*h
    return c;
  }
  double r=1/v;
  double h=16.*r-1.;
  double alfa=h+h;
  double b1=0.;
  double b2=0.;
  double b0;
  for (int i=24; i>=0; i--){
    b0=ci[i]+alfa*b1-b2;
    b2=b1;
    b1=b0;
  }
  h=sqrt(rpi2*r)*(b0-h*b2);
  h = exp(v)*h;
  //     if(nu*x .lt. 0) h=-h
  //     if(.not.lex) h=exp(v)*h
  return h;
}


double jbdbesi1(double x){
  const double eps = 1e-15, z1 = 1., hf = z1/2., pi = M_PI,
    ce = 0.57721566490153286e0, pih = pi/2., rpih = 2./pi, rpi2 = 1/(2.*pi);
  const double ci[25]={
     +0.975800602326285926e0
    ,-0.024467442963276385e0
    ,-0.000277205360763829e0
    ,-0.000009732146728020e0
    ,-0.000000629724238640e0
    ,-0.000000065961142154e0
    ,-0.000000009613872919e0
    ,-0.000000001401140901e0
    ,-0.000000000047563167e0
    ,+0.000000000081530681e0
    ,+0.000000000035408148e0
    ,+0.000000000005102564e0
    ,-0.000000000001804409e0
    ,-0.000000000001023594e0
    ,-0.000000000000052678e0
    ,+0.000000000000107094e0
    ,+0.000000000000026120e0
    ,-0.000000000000009561e0
    ,-0.000000000000004713e0
    ,+0.000000000000000829e0
    ,+0.000000000000000743e0
    ,-0.000000000000000080e0
    ,-0.000000000000000117e0
    ,+0.000000000000000011e0
    ,+0.000000000000000019e0};

  double v=fabs(x);
  if (v < 8.){
    double y= pow(hf*v,2);
    double xl=3.;// (nu+2)
    double a0=1.;
    double a1=1.+2.*y/((xl+1.)*(xl-1.));
    double a2=1.+y*(4.+3.*y/((xl+2.)*xl))/((xl+3.)*(xl-1.));
    double b0=1.;
    double b1=1.-y/(xl+1.);
    double b2=1.-y*(1.-y/(2.*(xl+2.)))/(xl+3.);
    double w1=3.+xl;
    double v1=3.-xl;
    double v3=xl-1.;
    double v2=v3+v3;
    double c=0.;
    for (int n=3; n<=30;n++){
      double c0=c;
      double fn=double(n);
      w1=w1+2.;
      double w2=w1-1.;
      double w3=w2-1.;
      double w4=w3-1.;
      double w5=w4-1.;
      double w6=w5-1.;
      v1=v1+1.;
      v2=v2+1.;
      v3=v3+1.;
      double u1=fn*w4;
      double e=v3/(u1*w3);
      double u2=e*y;
      double f1=1.+y*v1/(u1*w1);
      double f2=(1.+y*v2/(v3*w2*w5))*u2;
      double f3=-y*y*u2/(w4*w5*w5*w6);
      double a=f1*a2+f2*a1+f3*a0;
      double b=f1*b2+f2*b1+f3*b0;
      c=a/b;
      if(fabs(c0-c) < (eps*fabs(c))) break;
      a0=a1;
      a1=a2;
      a2=a;
      b0=b1;
      b1=b2;
      b2=b;
    }
    double h = c;
    h = hf*x*h;
    //    if(lex) h=exp(-v)*h
    return h;
  }
  double r=1/v;
  double h=16.*r-1.;
  double alfa=h+h;
  double b1=0.;
  double b2=0.;
  double b0;
  for (int i=24; i>=0; i--){
    b0=ci[i]+alfa*b1-b2;
    b2=b1;
    b1=b0;
  }
  h=sqrt(rpi2*r)*(b0-h*b2);
  h = exp(v)*h;
  if(x <  0) h=-h;
  return h;
}

double jbdbesk0(double x){

  const double eps = 1e-15, z1 = 1., hf = z1/2., pi = M_PI,
    ce = 0.57721566490153286e0, pih = pi/2., rpih = 2./pi,
    rpi2 = 1/(2.*pi);
  const double ck[17] ={
     +0.988408174230825800e0
    ,-0.011310504646928281e0
    ,+0.000269532612762724e0
    ,-0.000011106685196665e0
    ,+0.000000632575108500e0
    ,-0.000000045047337641e0
    ,+0.000000003792996456e0
    ,-0.000000000364547179e0
    ,+0.000000000039043756e0
    ,-0.000000000004579936e0
    ,+0.000000000000580811e0
    ,-0.000000000000078832e0
    ,+0.000000000000011360e0
    ,-0.000000000000001727e0
    ,+0.000000000000000275e0
    ,-0.000000000000000046e0
    ,+0.000000000000000008e0};

  double h = 0.;
  if(x < 0){
    std::cout << "dbesk0 called with negative argument " << x << std::endl;
    return 0.;
  }
  if (x < 1.){
    double b=hf*x;
    double bk=-(log(b)+ce);
    double f=bk;
    double p=hf;
    double q=hf;
    double c=1.;
    double d=b*b;
    double bk1=p;
    double fn,rfn,g;
    for (int n=1; n<= 15; n++){
      fn=double(n);
      rfn=1./fn;
      p=p*rfn;
      q=q*rfn;
      f=(f+p+q)*rfn;
      c=c*d*rfn;
      g=c*(p-fn*f);
      h=c*f;
      bk=bk+h;
      bk1=bk1+g;
      if( (bk1*h+fabs(g)*bk) <= (eps*bk*bk1) ) break;
    }
    h=bk;
    return h;
  }
  if(x <= 5.){
    double xn=0.;
    double a=9.-xn;
    double b=25.-xn;
    double c=768.*x*x;
    double c0=48.*x;
    double a0=1.;
    double a1=(16.*x+7.+xn)/a;
    double a2=(c+c0*(xn+23.)+xn*(xn+62.)+129.)/(a*b);
    double b0=1.;
    double b1=(16.*x+9.-xn)/a;
    double b2=(c+c0*b)/(a*b)+1.;
    c=0.;
    double fn,fn1,fn2,fn3,fn4,fn5,ran,f1,f2,f3;
    for (int n=3; n<=30;n++){
      c0=c;
      fn=double(n);
      fn2=fn+fn;
      fn1=fn2-1.;
      fn3=fn1/(fn2-3);
      fn4=12.*fn*fn-(1.-xn);
      fn5=16.*fn1*x;
      ran=1./(pow(fn2+1.,2)-xn);
      f1=fn3*(fn4-20.*fn)+fn5;
      f2=28.*fn-fn4-8+fn5;
      f3=fn3*(pow(fn2-5.,2)-xn);
      a=(f1*a2+f2*a1+f3*a0)*ran;
      b=(f1*b2+f2*b1+f3*b0)*ran;
      c=a/b;
      if(fabs(c0-c) <= (eps*fabs(c))) break;
      a0=a1;
      a1=a2;
      a2=a;
      b0=b1;
      b1=b2;
      b2=b;
    }
    h=c/sqrt(rpih*x);
    h = exp(-x)*h;
    return h;
  }
  double r=1/x;
  h=10.*r-1.;
  double alfa=h+h;
  double b0,b1=0., b2=0.;
  for (int i=16; i >= 0; i--){
    b0=ck[i]+alfa*b1-b2;
    b2=b1;
    b1=b0;
  }
  h=sqrt(pih*r)*(b0-h*b2);
  h = exp(-x)*h;
  return h;
}

double jbdbesk1(double x){

  const double eps = 1e-15, z1 = 1., hf = z1/2., pi = M_PI,
    ce = 0.57721566490153286e0, pih = pi/2., rpih = 2./pi,
    rpi2 = 1/(2.*pi);
  const double ck[17] ={
       +1.035950858772358331e0
      ,+0.035465291243331114e0
      ,-0.000468475028166889e0
      ,+0.000016185063810053e0
      ,-0.000000845172048124e0
      ,+0.000000057132218103e0
      ,-0.000000004645554607e0
      ,+0.000000000435417339e0
      ,-0.000000000045757297e0
      ,+0.000000000005288133e0
      ,-0.000000000000662613e0
      ,+0.000000000000089048e0
      ,-0.000000000000012726e0
      ,+0.000000000000001921e0
      ,-0.000000000000000305e0
      ,+0.000000000000000050e0
      ,-0.000000000000000009e0};

  double h = 0.;
  if(x < 0){
    std::cout << "dbesk1 called with negative argument " << x << std::endl;
    return 0.;
  }
  if (x < 1.){
    double b=hf*x;
    double bk=-(log(b)+ce);
    double f=bk;
    double p=hf;
    double q=hf;
    double c=1.;
    double d=b*b;
    double bk1=p;
    double fn,rfn,g;
    for (int n=1; n<= 15; n++){
      fn=double(n);
      rfn=1./fn;
      p=p*rfn;
      q=q*rfn;
      f=(f+p+q)*rfn;
      c=c*d*rfn;
      g=c*(p-fn*f);
      h=c*f;
      bk=bk+h;
      bk1=bk1+g;
      if( (bk1*h+fabs(g)*bk) <= (eps*bk*bk1) ) break;
    }
    h=bk;
    h=bk1/b;
    return h;
  }
  if(x <= 5.){
    double xn=4.;
    double a=9.-xn;
    double b=25.-xn;
    double c=768.*x*x;
    double c0=48.*x;
    double a0=1.;
    double a1=(16.*x+7.+xn)/a;
    double a2=(c+c0*(xn+23.)+xn*(xn+62.)+129.)/(a*b);
    double b0=1.;
    double b1=(16.*x+9.-xn)/a;
    double b2=(c+c0*b)/(a*b)+1.;
    c=0.;
    double fn,fn1,fn2,fn3,fn4,fn5,ran,f1,f2,f3;
    for (int n=3; n<=30;n++){
      c0=c;
      fn=double(n);
      fn2=fn+fn;
      fn1=fn2-1.;
      fn3=fn1/(fn2-3);
      fn4=12.*fn*fn-(1.-xn);
      fn5=16.*fn1*x;
      ran=1./(pow(fn2+1.,2)-xn);
      f1=fn3*(fn4-20.*fn)+fn5;
      f2=28.*fn-fn4-8+fn5;
      f3=fn3*(pow(fn2-5.,2)-xn);
      a=(f1*a2+f2*a1+f3*a0)*ran;
      b=(f1*b2+f2*b1+f3*b0)*ran;
      c=a/b;
      if(fabs(c0-c) <= (eps*fabs(c))) break;
      a0=a1;
      a1=a2;
      a2=a;
      b0=b1;
      b1=b2;
      b2=b;
    }
    h=c/sqrt(rpih*x);
    h = exp(-x)*h;
    return h;
  }
  double r=1/x;
  h=10.*r-1.;
  double alfa=h+h;
  double b0,b1=0., b2=0.;
  for (int i=16; i >= 0; i--){
    b0=ck[i]+alfa*b1-b2;
    b2=b1;
    b1=b0;
  }
  h=sqrt(pih*r)*(b0-h*b2);
  h = exp(-x)*h;
  return h;
}

// recursion for the others
double jbdbesk2(const double x){
  return 2./x*jbdbesk1(x)+jbdbesk0(x);
}

double jbdbesk3(const double x){
  return (8./(x*x)+1.)*jbdbesk1(x)+4./x*jbdbesk0(x);
}

double jbdbesk4(const double x){
  return (48./(x*x*x)+8./x)*jbdbesk1(x)+(24./(x*x)+1.)*jbdbesk0(x);
}
