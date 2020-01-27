// sunsetintegrals.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains all the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]
//
// the function phi(x,y) defined in Davydychev Tausk Nucl. Phys. B397(1993)123
// except multiplied with m3**2 to make it fully symmetric in the masses.
// also contains the barred functions and the full functions
// plus a version with the imaginary part correctly. zhh,zhh1 and zhh21
// derivatives are hhd,hh1d and hh21d.
// above threshold, use zhhd,zhh1d,zhh21d
//
// sometimes you will get NaN even if the integral is really well defined
// this comes from zeros in the Kählen function lm


// allow to easily change integration routines
// real integral (CHIRON v0.51 possible jbdgauss,jbdgauss2,jbquad15,jbdquad21)
#ifndef DINTEGRAL
#define DINTEGRAL jbdgauss
#endif
// real integral with singularity (CHIRON v0.51 possible jbdcauch,jbdcauch2,
//                                                       jbdsing15,jbdsing21) 
#ifndef SINTEGRAL
#define SINTEGRAL jbdcauch
#endif

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include "sunsetintegrals.h"
#include "jbnumlib.h"

const double pi = M_PI;
const double pi2 = pi*pi;
const double pi16 = 1./(16.*pi2);
const double pi162 = pi16*pi16;

// this is the real precision pi162 etc included
// assumes m1sq,..,qsq etc are all of order 1 and no crazy other things
double precisionsunsetintegrals = 1e-10;

void setprecisionsunsetintegrals(const double eps){
  precisionsunsetintegrals = eps;
}

double getprecisionsunsetintegrals(void){
  return precisionsunsetintegrals;
}

// the sunsetintegrals, 0 = hbar, 1 = hbar_1, 2 = hbar_21
// version with k2 is hbar
// dispersive version is hbarp
// derivative with k2 is hbard
// dispersive derivative is hbarpd : notice that we use the dispersion
// relation for the derivative here. (otherwise it contains 1/(z-s)^2 in the
// integral)
// needs epsh set beforehand and pi162
// they do not check that psq < (m1+m2+m3)**2
// i = 0 hbar
// i = 1 hbar_1
// i = 2 hbar_21

struct hbarcom{
  double sigm1,psqh,xm12h,xm22h,xm32h;
  int ih;
};

struct hbarpcom{
  double zm1,psqhh,xm12hh,xm22hh,xm32hh;
  int ihh;
};
double dlam(const double x,const double y,const double z){
  return  pow(x-y-z,2)-4.*y*z;
}

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
hbarcom hbardat;

double hbar2(const double x){
  //const double stretch1 = 1.;  //  x near zero
  //const double stretch1p = 1.;  //  x near 1
  //const double stretch2 = 1.;  // sigma near infinity
  //const double stretch2p = 1.;  // sigma near threshold
  //double xxt = pow(x,stretch1);
  //double xx = 1.-pow(1.-xxt,stretch1p);
  //double yyt = pow(hbardat.sigm1,stretch2);
  //double yy = 1.-pow(1.-yyt,stretch2p);
  //double hbar2t = (pow(hbardat.sigm1,(stretch2-1.))*stretch2)
  //  *(pow(x,(stretch1-1.))*stretch1)*(pow(1.-yyt,stretch2p-1.)*stretch2p)
  //  *(pow(1.-xxt,stretch1p-1.)*stretch1p);
  // no strectching at all
  double xx = x;
  double yy = hbardat.sigm1;
  double hbar2t = 1.;

  double siglow = pow(sqrt(hbardat.xm22h)+sqrt(hbardat.xm32h),2);
  double sigma = siglow/yy;
  hbar2t = hbar2t*siglow/pow(yy,2);
  hbar2t = hbar2t*sqrt(dlam(1.,hbardat.xm22h/sigma,hbardat.xm32h/sigma));
  double z1 = hbardat.xm12h*(1.-xx)+sigma*xx;
  double z2 = xx*(1.-xx)*hbardat.psqh;
  double xk2 = pow(xx,hbardat.ih) * pi162*(log(1.-z2/z1)+z2/z1);
  return hbar2t*xk2;
}

double hbar1(const double y){
  hbardat.sigm1 = y;
  return DINTEGRAL(hbar2,0.,1.,precisionsunsetintegrals/5.);
}

double hbar(const double xm12,const double xm22,const double xm32
	    ,const double psq,const int i){
  hbardat.ih = i;
  hbardat.psqh = psq;
  hbardat.xm12h = xm12;
  hbardat.xm22h = xm22;
  hbardat.xm32h = xm32;
  return DINTEGRAL(hbar1,0.,1.,precisionsunsetintegrals);
}

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
hbarcom hbarddat;

double hbard2(const double x){
  const double stretch1 = 1.;  //  x near zero
  const double stretch2 = 1.;  // sigma near infinity
  double xx = pow(x,stretch1);
  double yy = pow(hbarddat.sigm1,stretch2);
  double hbar2t = (pow(hbarddat.sigm1,(stretch2-1.))*stretch2)
    *(pow(x,(stretch1-1.))*stretch1);
  double siglow = pow(sqrt(hbarddat.xm22h)+sqrt(hbarddat.xm32h),2);
  double sigma = siglow/yy;
  hbar2t = hbar2t*siglow/pow(yy,2);
  hbar2t = hbar2t*sqrt(dlam(1.,hbarddat.xm22h/sigma,hbarddat.xm32h/sigma));
  double z1 = hbarddat.xm12h*(1.-xx)+sigma*xx;
  double z2 = -xx*(1.-xx)*hbarddat.psqh;
  double xk2 = pow(xx,hbarddat.ih) * pi162*xx*(1.-xx)*z2/(z1*(z1+z2));
  return hbar2t*xk2;
}

double hbard1(const double y){
  hbarddat.sigm1 = y;
  return DINTEGRAL(hbard2,0.,1.,precisionsunsetintegrals/5.);
}

double hbard(const double xm12,const double xm22,const double xm32
	    ,const double psq,const int i){
  hbarddat.ih = i;
  hbarddat.psqh = psq;
  hbarddat.xm12h = xm12;
  hbarddat.xm22h = xm22;
  hbarddat.xm32h = xm32;
  return DINTEGRAL(hbard1,0.,1.,precisionsunsetintegrals);
}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// complex part also correct and usable above threshold

hbarpcom hbarpdat;

double zimhbar(const double x){
  const double stretch1 = 1.; // 
  // const double stretch2 = 1.; // z near infinity, not used here
  double xx = pow(x,stretch1);
  double hbarp2 = (pow(x,(stretch1-1.))*stretch1);
  // zlow = pow(sqrt(hbarpdat.xm12hh)+sqrt(hbarpdat.xm22hh)
  //	    +sqrt(hbarpdat.xm32hh),2);
  double zet = hbarpdat.psqhh;
  double xm1p = sqrt(hbarpdat.xm12hh);
  double p1 = sqrt(zet);
  double e1max = 1./(2.*p1)*(zet+hbarpdat.xm12hh
       -pow(sqrt(hbarpdat.xm22hh)+sqrt(hbarpdat.xm32hh),2));
  double e1 = xm1p+xx*(e1max-xm1p);
  hbarp2 = hbarp2*(e1max-xm1p);
  double s23 = zet+hbarpdat.xm12hh-2.*p1*e1;
  double e2maxmin = 
    sqrt(dlam(zet,hbarpdat.xm12hh,s23)*
	 dlam(s23,hbarpdat.xm22hh,hbarpdat.xm32hh))/
    (p1*s23);
  double xt;
  if (hbarpdat.ihh == 0){
    xt = 1.;}
  else{
    if(hbarpdat.ihh == 1){
      xt = e1/p1;}
    else{
      if (hbarpdat.ihh == 2){
	xt = (4.*e1*e1-hbarpdat.xm12hh)/(3.*zet);}
      else{
	if (hbarpdat.ihh == 3){
	  xt = e1*(4.*e1*e1-hbarpdat.xm12hh)/(3.*zet*p1);}
	else{
	  std::cout<<"stupid ihh in zimhbar"<<std::endl;
	  assert(0);}
      }
    }
  }
  return hbarp2/(-16.)/pow(2.*pi,3) *e2maxmin*xt;
}

double hbarp2(const double x){
  double stretch1 = 1.; //
  double stretch2 = 1.;  // z near infinity
  double xx = pow(x,stretch1);
  double yy = pow(hbarpdat.zm1,stretch2);
  double hbarp2p = (pow(hbarpdat.zm1,(stretch2-1.))*stretch2)
    *(pow(x,(stretch1-1.))*stretch1);
  double zlow = pow(sqrt(hbarpdat.xm12hh)+sqrt(hbarpdat.xm22hh)
		    +sqrt(hbarpdat.xm32hh),2);
  double zet = zlow/yy;
  hbarp2p = hbarp2p*zlow/pow(yy,2);
  hbarp2p = hbarp2p*pow(hbarpdat.psqhh,2)/pi/pow(zet,2)/(zet-hbarpdat.psqhh);
  double xm1p = sqrt(hbarpdat.xm12hh);
  double p1 = sqrt(zet);
  double e1max = 1./(2.*p1)*(zet+hbarpdat.xm12hh
      -pow(sqrt(hbarpdat.xm22hh)+sqrt(hbarpdat.xm32hh),2));
  double e1 = xm1p+xx*(e1max-xm1p);
  hbarp2p = hbarp2p*(e1max-xm1p);
  double s23 = zet+hbarpdat.xm12hh-2.*p1*e1;
  double e2maxmin = 
    sqrt(dlam(zet,hbarpdat.xm12hh,s23)*
	 dlam(s23,hbarpdat.xm22hh,hbarpdat.xm32hh))/
    (p1*s23);
  double xt;
  if (hbarpdat.ihh == 0){
    xt = 1.;}
  else{
    if (hbarpdat.ihh == 1){
      xt = e1/p1;}
    else{
      if (hbarpdat.ihh == 2){
        xt = (4.*e1*e1-hbarpdat.xm12hh)/(3.*zet);}
      else{
	if (hbarpdat.ihh == 3){
	  xt = e1*(2.*e1*e1-hbarpdat.xm12hh)/(zet*p1);}
	else{
        std::cout<<"stupid ihh in hbarp2"<<std::endl;
        assert(0);}
      }
    }
  }
  return hbarp2p/(-16.)/pow(2.*pi,3) *e2maxmin*xt;
}

double hbarp1(const double y){
  hbarpdat.zm1 = y;
  return DINTEGRAL(hbarp2,0.,1.,precisionsunsetintegrals/5.);
}

dcomplex zhbar(const double xm12,const double xm22,const double xm32,
	       const double psq,const int i){
//      if (psq <= pow(sqrt(xm12)+sqrt(xm22)+sqrt(xm32),2)){
//        return hbar(xm12,xm22,xm32,psq,i);}
// valid for stretch1 = 1
  hbarpdat.ihh = i;
  hbarpdat.psqhh = psq;
  hbarpdat.xm12hh = xm12;
  hbarpdat.xm22hh = xm22;
  hbarpdat.xm32hh = xm32;
  double zlow = pow(sqrt(xm12)+sqrt(xm22)+sqrt(xm32),2);
  double ysing = zlow/psq;
  dcomplex zhbar2 = SINTEGRAL(hbarp1,0.,1.,ysing,precisionsunsetintegrals);
  if (psq >= zlow){
    hbarpdat.zm1 = ysing;
    double zim = DINTEGRAL(zimhbar,0.,1.,precisionsunsetintegrals);
    zhbar2 = zhbar2+dcomplex(0.,zim);}
  return zhbar2;
}
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

hbarpcom hbarpddat;

double zimhbard(const double xxx){
  double stretch1 = 2.; //   e1 near m1
  double stretch1p = 2.;//   e1 near e1max
  //double stretch2 = 1.; // z near infinity, not used here
  double x = 1.-pow((1.-xxx),stretch1p);
  double xx = pow(x,stretch1);
  //double yy = pow(hbarpddat.zm1,stretch2);
  double hbarp2p = (pow(x,(stretch1-1.))*stretch1)
    *(pow((1.-xxx),(stretch1p-1.))*stretch1p);
  //double zlow = pow(sqrt(hbarpddat.xm12hh)+sqrt(hbarpddat.xm22hh)
  //		    +sqrt(hbarpddat.xm32hh),2);
  double zet = hbarpddat.psqhh;
  double xm1p = sqrt(hbarpddat.xm12hh);
  double p1 = sqrt(zet);
  double e1max = 1./(2.*p1)*(zet+hbarpddat.xm12hh
      -pow(sqrt(hbarpddat.xm22hh)+sqrt(hbarpddat.xm32hh),2));
  double e1 = xm1p+xx*(e1max-xm1p);
  hbarp2p = hbarp2p*(e1max-xm1p);
  double s23 = zet+hbarpddat.xm12hh-2.*p1*e1;
  double dm23ds = 1.-e1/p1;
  double xlam1 = sqrt(dlam(zet,hbarpddat.xm12hh,s23));
  double xlam2 = sqrt(dlam(s23,hbarpddat.xm22hh,hbarpddat.xm32hh));
  //double e2maxmin = xlam1*xlam2/(p1*s23);
  double xt;
  if (hbarpddat.ihh == 0){
    xt = dm23ds*(-xlam1*xlam2/(pow(s23,2)*p1)
     +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
       +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		 (s23*xlam1*xlam2*p1))
      +(-xlam1*xlam2/(s23*2.*pow(p1,3))
	+xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*p1));}
  else{
    if (hbarpddat.ihh == 1){
      xt = e1*(dm23ds*(-xlam1*xlam2/(pow(s23,2)*zet)
	     +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
	       +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		       (s23*xlam1*xlam2*zet))
	       +(-xlam1*xlam2/(s23*pow(zet,2))
		 +xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*zet)) );}
    else{
      if (hbarpddat.ihh == 2){
        xt = (4.*pow(e1,2)-hbarpddat.xm12hh)/3.*
	  (dm23ds*(-xlam1*xlam2/(pow(s23,2)*pow(p1,3))
		   +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
		     +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		   (s23*xlam1*xlam2*pow(p1,3)))
	   +(-xlam1*xlam2/(s23*2./3.*pow(p1,5))
	     +xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*pow(p1,3)) ) );}
      else{
	  std::cout<<"stupid ihh in zimhbard"<<std::endl;
	  assert(0);}
    }
  }
  return hbarp2p/(-16.)/pow(2.*pi,3) *xt;
}

double hbarpd2(const double xxx){
  double stretch1 = 2.; //   e1 near m1
  double stretch1p = 2.;//   e1 near e1max
  double stretch2 = 1.; // z near infinity
  double x = 1.-pow((1.-xxx),stretch1p);
  double xx = pow(x,stretch1);
  double yy = pow(hbarpddat.zm1,stretch2);
  double hbarp2p = (pow(hbarpddat.zm1,(stretch2-1.))*stretch2)
    *(pow(x,(stretch1-1.))*stretch1)
    *(pow((1.-xxx),(stretch1p-1.))*stretch1p);
  double zlow = pow(sqrt(hbarpddat.xm12hh)+sqrt(hbarpddat.xm22hh)
		    +sqrt(hbarpddat.xm32hh),2);
  double zet = zlow/yy;
  hbarp2p = hbarp2p*zlow/pow(yy,2);
  hbarp2p = hbarp2p*hbarpddat.psqhh/pi/zet/(zet-hbarpddat.psqhh);
  double xm1p = sqrt(hbarpddat.xm12hh);
  double p1 = sqrt(zet);
  double e1max = 1./(2.*p1)*(zet+hbarpddat.xm12hh
      -pow(sqrt(hbarpddat.xm22hh)+sqrt(hbarpddat.xm32hh),2));
  double e1 = xm1p+xx*(e1max-xm1p);
  hbarp2p = hbarp2p*(e1max-xm1p);
  double s23 = zet+hbarpddat.xm12hh-2.*p1*e1;
  double dm23ds = 1.-e1/p1;
  double xlam1 = sqrt(dlam(zet,hbarpddat.xm12hh,s23));
  double xlam2 = sqrt(dlam(s23,hbarpddat.xm22hh,hbarpddat.xm32hh));
  //double e2maxmin = xlam1*xlam2/(p1*s23);
  double xt;
  if (hbarpddat.ihh == 0){
    xt = dm23ds*(-xlam1*xlam2/(pow(s23,2)*p1)
     +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
       +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		 (s23*xlam1*xlam2*p1))
      +(-xlam1*xlam2/(s23*2.*pow(p1,3))
	+xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*p1));}
  else{
    if (hbarpddat.ihh == 1){
      xt = e1*(dm23ds*(-xlam1*xlam2/(pow(s23,2)*zet)
	     +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
	       +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		       (s23*xlam1*xlam2*zet))
	       +(-xlam1*xlam2/(s23*pow(zet,2))
		 +xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*zet)) );}
    else{
      if (hbarpddat.ihh == 2){
        xt = (4.*pow(e1,2)-hbarpddat.xm12hh)/3.*
	  (dm23ds*(-xlam1*xlam2/(pow(s23,2)*pow(p1,3))
		   +(pow(xlam2,2)*(s23-zet-hbarpddat.xm12hh)
		     +pow(xlam1,2)*(s23-hbarpddat.xm22hh-hbarpddat.xm32hh))/
		   (s23*xlam1*xlam2*pow(p1,3)))
	   +(-xlam1*xlam2/(s23*2./3.*pow(p1,5))
	     +xlam2*(zet-s23-hbarpddat.xm12hh)/(s23*xlam1*pow(p1,3)) ) );}
      else{
	  std::cout<<"stupid ihh in hbarpd2"<<std::endl;
	  assert(0);}
    }
  }
return hbarp2p/(-16.)/pow(2.*pi,3)*xt;
}

double hbarpd1(const double y){
  hbarpddat.zm1 = y;
  return DINTEGRAL(hbarpd2,0.,1.,precisionsunsetintegrals/5.);
}

dcomplex zhbard(const double xm12,const double xm22,const double xm32,
	       const double psq,const int i){
//      if (psq <= pow(sqrt(xm12)+sqrt(xm22)+sqrt(xm32),2)){
//        return hbar(xm12,xm22,xm32,psq,i);}
// valid for stretch1 = 1
  hbarpddat.ihh = i;
  hbarpddat.psqhh = psq;
  hbarpddat.xm12hh = xm12;
  hbarpddat.xm22hh = xm22;
  hbarpddat.xm32hh = xm32;
  double zlow = pow(sqrt(xm12)+sqrt(xm22)+sqrt(xm32),2);
  double ysing = zlow/psq;
  dcomplex zhbar2 = SINTEGRAL(hbarpd1,0.,1.,ysing,precisionsunsetintegrals);
  if (psq >= zlow){
    hbarpddat.zm1 = ysing;
    double zim = DINTEGRAL(zimhbard,0.,1.,precisionsunsetintegrals);
    zhbar2 = zhbar2+dcomplex(0.,zim);}
  return zhbar2;
}




//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// clausen's function in terms of li2 and with the definition.
double  claus2(const double x){
  const double dli2one =  pi2/6.;//1.64493406684823;
  dcomplex expt=exp(dcomplex(0.,x));
  return real(dcomplex(0.,-1.)*(jbdli2(expt)-dli2one)
	      +dcomplex(0.,-pi/2.*x+x*x/4.));
}

dcomplex zphi(const double xm1,const double xm2,const double xm3){
  double xm12 = xm1*xm1;
  double xm22 = xm2*xm2;
  double xm32 = xm3*xm3;
  double xlam2 = pow(xm12-xm22-xm32,2)-4.*xm22*xm32;
  double xx,yy,xlam,xmass2;
  if (xlam2 < 0.){
    //  no separate cases needed, is fully symmetric
    xx = xm12/xm32;
    yy = xm22/xm32;
    xlam = sqrt(-xlam2)/xm32;
    xmass2 = xm32;
    double x1 = (-1.+xx+yy)/(2.*sqrt(xx*yy));
    double x2 = (1.+xx-yy)/(2.*sqrt(xx));
    double x3 = (1.-xx+yy)/(2.*sqrt(yy));
    return 2./xlam/xmass2*(claus2(2.*acos(x1))+claus2(2.*acos(x2))
			   +claus2(2.*acos(x3)));}
  else{
    if ((xm1+xm2)<= xm3){
      xx = xm12/xm32;
      yy = xm22/xm32;
      xlam = sqrt(xlam2)/xm32;
      xmass2 = xm32;}
    else{
      if ((xm1+xm3)<=xm2){
	xx = xm12/xm22;
	yy = xm32/xm22;
	xlam = sqrt(xlam2)/xm22;
	xmass2 = xm22;}
      else
	if ((xm2+xm3)<=xm1){
          xx = xm22/xm12;
          yy = xm32/xm12;
          xlam = sqrt(xlam2)/xm12;
          xmass2 = xm12;}
        else{
          std::cout << "crazy mass combination in phi" << std::endl;
          assert(0);}//makes program stop
    }
    dcomplex zx1 = 0.5*(1.+xx-yy-xlam);
    dcomplex zx2 = 0.5*(1.-xx+yy-xlam);     
    return 1./xlam*(2.*log(zx1)*log(zx2)-2.*jbdli2(zx1)-2.*jbdli2(zx2)
		    -log(xx)*log(yy)+pi2/3.)/xmass2;
  }
};

double psi(const double xm1,const double xm2,const double xm3){
  double xm12 = xm1*xm1;
  double xm22 = xm2*xm2;
  double xm32 = xm3*xm3;
  double xlam2 = pow(xm12-xm22-xm32,2)-4.*xm22*xm32;
  return  -xlam2*real(zphi(xm1,xm2,xm3));
}

double h0(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  //double dm = m1sq-m2sq-m3sq;
  double hh0 = -0.5*psi(m1,m2,m3)
    + m1sq*(pi2/12.+1.5-ln1+0.5*(-ln2*ln3+ln1*ln4))
    + m2sq*(pi2/12.+1.5-ln2+0.5*(-ln1*ln3+ln2*ln4))
    + m3sq*(pi2/12.+1.5-ln3+0.5*(-ln1*ln2+ln3*ln4));
  return pi162*hh0;
}


double h0p(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = pow(m1sq-m2sq-m3sq,2)-4.*m2sq*m3sq;
  double hh0p = m1sq*m2sq*m3sq/pow(lm,2)*psi(m1,m2,m3)+1./8.
    + m1sq*dm*ln1/2./lm
    + m2sq/2./lm*(m2sq-m1sq-m3sq)*ln2
    + m3sq/2./lm*(m3sq-m1sq-m2sq)*ln3;
  return pi162*hh0p;
}

double h10(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  return (0.25*(-1.+m1sq*dm/lm)*psi(m1,m2,m3)
      +m1sq*(3./8.-0.5*ln1)
      +m2sq*(pi2/24.+9./16.-0.25*ln2+0.25*(ln2*ln4-ln1*ln3))
	  +m3sq*(pi2/24.+9./16.-0.25*ln3+0.25*(ln3*ln4-ln1*ln2)))
      *pi162;
}

double h10p(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  return (pow(m1sq,2)*m2sq*m3sq*dm/pow(lm,3)*psi(m1,m2,m3)
     -m1sq*dm/(6.*lm)+pow(m1sq,2)*ln1/(6.*lm)*(1.+12.*m2sq*m3sq/lm)
     +m2sq*ln2/(6.*lm)*
           (m2sq-m3sq-2.*m1sq-6.*m1sq*m3sq/lm*(dm+2.*m2sq))
      +m3sq*ln3/(6.*lm)*
           (m3sq-m2sq-2.*m1sq-6.*m1sq*m2sq/lm*(dm+2.*m3sq))
	  + 7./72.)*pi162;
}

double h210(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  return (1./(6.*lm*lm)*
       (-lm*lm+lm*m1sq*dm+2.*pow(m1sq,2)*m2sq*m3sq)*psi(m1,m2,m3)
     + m1sq*(17./72.-ln1/3.+m1sq*dm*ln1/(6.*lm))
     + m2sq*(pi2/36.+19./54.-ln2/9.
            -m1sq/(6.*lm)*(dm+2.*m3sq)*ln2+1./6.*(ln2*ln4-ln1*ln3))
      + m3sq*(pi2/36.+19./54.-ln3/9.
	      -m1sq/(6.*lm)*(dm+2.*m2sq)*ln3+1./6.*(ln3*ln4-ln1*ln2)))*pi162;
}

double h210p(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  return ( 17./288.-m1sq/(24.*lm)*(dm+2.*m1sq)
	   -5.*pow(m1sq,2)*m2sq*m3sq/(6.*lm*lm)
	   +pow(m1sq,3)*m2sq*m3sq/pow(lm,3)*(1.+5.*m2sq*m3sq/lm)*psi(m1,m2,m3)
	   +pow(m1sq,3)*dm/(12.*lm*lm)*ln1*(1.+30.*m2sq*m3sq/lm)
	   +pow(m1sq,2)*m2sq*pow(m3sq,2)/pow(lm,3) *ln3*5./2.*(m3sq-m1sq-m2sq)
	   +m3sq/(12.*lm)*ln3*(m3sq-m2sq)
	   +m1sq*m3sq/(12.*pow(lm,2))*ln3*(-3.*pow(m1sq,2)-12.*m1sq*m2sq
		    +5.*m1sq*m3sq+3.*pow(m2sq,2)-m2sq*m3sq-2.*pow(m3,4))
	   +pow(m1sq,2)*m3sq*pow(m2sq,2)/pow(lm,3)*ln2*5./2.*(m2sq-m1sq-m3sq)
	   +m2sq/(12.*lm)*ln2*(m2sq-m3sq)
	   +m1sq*m2sq/(12.*pow(lm,2))*ln2*(-3.*pow(m1sq,2)-12.*m1sq*m3sq
	       +5.*m1sq*m2sq+3.*pow(m3sq,2)-m2sq*m3sq-2.*pow(m2,4))
	   )*pi162;
}


double h310(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  double hh31 = 0.;
   hh31 +=
       + pow(lm,-2) * (  - 1./12.*ln3*pow(m1sq,2)*pow(m3sq,3) - 1./12.*
         ln3*pow(m1sq,2)*m2sq*pow(m3sq,2) + 1./6.*ln3*pow(m1sq,2)*pow(
         m2sq,2)*m3sq + 1./6.*ln3*pow(m1sq,3)*pow(m3sq,2) - 1./12.*ln3*
         pow(m1sq,3)*m2sq*m3sq - 1./12.*ln3*pow(m1sq,4)*m3sq + 1./12.*
         ln2*pow(m1sq,2)*pow(m3sq,3) + 1./12.*ln2*pow(m1sq,2)*m2sq*pow(
         m3sq,2) - 1./6.*ln2*pow(m1sq,2)*pow(m2sq,2)*m3sq - 1./6.*ln2*
         pow(m1sq,3)*pow(m3sq,2) - 5./12.*ln2*pow(m1sq,3)*m2sq*m3sq + 1.
         /12.*ln2*pow(m1sq,4)*m3sq + 1./2.*ln1*pow(m1sq,3)*m2sq*m3sq );

      hh31 +=  + pow(lm,-1) * ( 1./24.*pow(m1sq,2)*m3sq + 1./24.*pow(
         m1sq,2)*m2sq - 1./24.*pow(m1sq,3) + 1./6.*ln3*m1sq*pow(m3sq,2)
          - 1./6.*ln3*m1sq*m2sq*m3sq - 1./8.*ln3*pow(m1sq,2)*m3sq - 1./
         6.*ln2*m1sq*pow(m3sq,2) + 1./6.*ln2*m1sq*m2sq*m3sq + 1./4.*ln2
         *pow(m1sq,2)*m3sq + 1./8.*ln2*pow(m1sq,2)*m2sq - 1./6.*ln2*
         pow(m1sq,3) - 1./8.*ln1*pow(m1sq,2)*m3sq - 1./8.*ln1*pow(
         m1sq,2)*m2sq + 1./6.*ln1*pow(m1sq,3) );

      hh31 +=  + psi(m1,m2,m3)*(pow(lm,-3) * 
				   (  - 1./4.*pow(m1sq,3)*m2sq*pow(
         m3sq,2) - 1./4.*pow(m1sq,3)*pow(m2sq,2)*m3sq + 1./4.*pow(
         m1sq,4)*m2sq*m3sq )

                    +pow(lm,-2) * ( 1./4.*pow(m1sq,2)*m2sq*m3sq )

                     +pow(lm,-1) * (  - 1./8.*m1sq*m3sq - 1./8.*m1sq*
                                    m2sq + 1./8.*pow(m1sq,2) )

		   + (  - 1./8. ));

      hh31 +=  + 33./128.*m3sq + 33./128.*m2sq + 55./288.*m1sq - 1./16.
         *ln3*m3sq + 1./8.*pow(ln3,2)*m3sq - 1./16.*ln2*m2sq + 1./6.*
         ln2*m1sq + 1./8.*ln2*ln3*m3sq + 1./8.*ln2*ln3*m2sq + 1./8.*
         pow(ln2,2)*m2sq - 1./4.*ln1*m1sq + 1./8.*ln1*ln3*m3sq - 1./8.*
         ln1*ln3*m2sq - 1./8.*ln1*ln2*m3sq + 1./8.*ln1*ln2*m2sq + 1./48.
         *pow(pi,2)*m3sq + 1./48.*pow(pi,2)*m2sq;

  return pi162*hh31;
}

double h310p(const double m1sq, const double m2sq,const double m3sq,
	  const double xmu2){
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  double hh31 = 0.;

   hh31 =
       + pow(lm,-4) * (  - 2./5.*ln3*pow(m1sq,3)*pow(m3sq,5) - 3./2.*
         ln3*pow(m1sq,3)*m2sq*pow(m3sq,4) - 17./10.*ln3*pow(m1sq,3)* 
         pow(m2sq,2)*pow(m3sq,3) + 5./2.*ln3*pow(m1sq,3)*pow(m2sq,3)*
         pow(m3sq,2) + 11./10.*ln3*pow(m1sq,3)*pow(m2sq,4)*m3sq + 4./5.
         *ln3*pow(m1sq,4)*pow(m3sq,4) + 27./5.*ln3*pow(m1sq,4)*m2sq*   
         pow(m3sq,3) - 23./5.*ln3*pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)  
          - 11./5.*ln3*pow(m1sq,4)*pow(m2sq,3)*m3sq - 2./5.*ln3*pow(   
         m1sq,5)*pow(m3sq,3) - 23./10.*ln3*pow(m1sq,5)*m2sq*pow(m3sq,2)
          + 11./10.*ln3*pow(m1sq,5)*pow(m2sq,2)*m3sq + 2./5.*ln2*pow(  
         m1sq,3)*pow(m3sq,5) + 3./2.*ln2*pow(m1sq,3)*m2sq*pow(m3sq,4)  
          + 17./10.*ln2*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3) - 5./2.*ln2
         *pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2) - 11./10.*ln2*pow(m1sq,3)
         *pow(m2sq,4)*m3sq - 4./5.*ln2*pow(m1sq,4)*pow(m3sq,4) - 27./5.
         *ln2*pow(m1sq,4)*m2sq*pow(m3sq,3) - 47./5.*ln2*pow(m1sq,4)*   
         pow(m2sq,2)*pow(m3sq,2) + 11./5.*ln2*pow(m1sq,4)*pow(m2sq,3)*
         m3sq );

      hh31 +=  + pow(lm,-4) * ( 2./5.*ln2*pow(m1sq,5)*pow(m3sq,3) + 23./
         10.*ln2*pow(m1sq,5)*m2sq*pow(m3sq,2) - 11./10.*ln2*pow(m1sq,5)
         *pow(m2sq,2)*m3sq + 14*ln1*pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)
          );

      hh31 +=  + pow(lm,-3) * ( 1./15.*pow(m1sq,3)*pow(m3sq,3) + 11./10.
         *pow(m1sq,3)*m2sq*pow(m3sq,2) + 11./10.*pow(m1sq,3)*pow(
         m2sq,2)*m3sq + 1./15.*pow(m1sq,3)*pow(m2sq,3) - 2./15.*pow(
         m1sq,4)*pow(m3sq,2) - 43./30.*pow(m1sq,4)*m2sq*m3sq - 2./15.*
         pow(m1sq,4)*pow(m2sq,2) + 1./15.*pow(m1sq,5)*m3sq + 1./15.*
         pow(m1sq,5)*m2sq + 1./5.*ln3*pow(m1sq,2)*pow(m3sq,4) + 23./60.
         *ln3*pow(m1sq,2)*m2sq*pow(m3sq,3) - 1./5.*ln3*pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,2) - 23./60.*ln3*pow(m1sq,2)*pow(m2sq,3)*m3sq
          - 1./5.*ln3*pow(m1sq,3)*pow(m3sq,3) - 4./15.*ln3*pow(m1sq,3)*
         m2sq*pow(m3sq,2) + 7./15.*ln3*pow(m1sq,3)*pow(m2sq,2)*m3sq + 3.
         /5.*ln3*pow(m1sq,4)*pow(m3sq,2) - 59./60.*ln3*pow(m1sq,4)*m2sq
         *m3sq - 1./5.*ln3*pow(m1sq,5)*m3sq - 1./5.*ln2*pow(m1sq,2)*
         pow(m3sq,4) - 23./60.*ln2*pow(m1sq,2)*m2sq*pow(m3sq,3) + 1./5.
         *ln2*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2) + 23./60.*ln2*pow(
         m1sq,2)*pow(m2sq,3)*m3sq + 1./5.*ln2*pow(m1sq,3)*pow(m3sq,3)
          + 4./15.*ln2*pow(m1sq,3)*m2sq*pow(m3sq,2) );

      hh31 +=  + pow(lm,-3) * (  - 7./15.*ln2*pow(m1sq,3)*pow(m2sq,2)*
         m3sq - 3./5.*ln2*pow(m1sq,4)*pow(m3sq,2) - 131./60.*ln2*pow(
         m1sq,4)*m2sq*m3sq + 1./5.*ln2*pow(m1sq,5)*m3sq + 19./6.*ln1*
         pow(m1sq,4)*m2sq*m3sq );

      hh31 +=  + pow(lm,-2) * (  - 1./30.*pow(m1sq,2)*pow(m3sq,2) - 1./
         6.*pow(m1sq,2)*m2sq*m3sq - 1./30.*pow(m1sq,2)*pow(m2sq,2) + 1./
         20.*pow(m1sq,3)*m3sq + 1./20.*pow(m1sq,3)*m2sq - 1./12.*pow(
         m1sq,4) - 1./10.*ln3*m1sq*pow(m3sq,3) - 1./30.*ln3*m1sq*m2sq*
         pow(m3sq,2) + 2./15.*ln3*m1sq*pow(m2sq,2)*m3sq + 1./20.*ln3*
         pow(m1sq,2)*pow(m3sq,2) - 1./20.*ln3*pow(m1sq,2)*m2sq*m3sq + 1.
         /10.*ln2*m1sq*pow(m3sq,3) + 1./30.*ln2*m1sq*m2sq*pow(m3sq,2)
          - 2./15.*ln2*m1sq*pow(m2sq,2)*m3sq - 1./20.*ln2*pow(m1sq,2)*
         pow(m3sq,2) + 1./20.*ln2*pow(m1sq,2)*m2sq*m3sq - 1./20.*ln2*
         pow(m1sq,4) + 1./20.*ln1*pow(m1sq,4) );

      hh31 +=  + pow(lm,-1) * ( 1./60.*m1sq*m3sq + 1./60.*m1sq*m2sq - 1.
         /120.*pow(m1sq,2) + 1./20.*ln3*pow(m3sq,2) - 1./20.*ln3*m2sq*
         m3sq - 1./20.*ln2*pow(m3sq,2) + 1./20.*ln2*m2sq*m3sq );

      hh31 +=  + psi(m1,m2,m3)*pow(lm,-5) * (  - 7*pow(m1sq,4)*pow(m2sq,2)*pow(
         m3sq,3) - 7*pow(m1sq,4)*pow(m2sq,3)*pow(m3sq,2) + 7*pow(
         m1sq,5)*pow(m2sq,2)*pow(m3sq,2) );

      hh31 +=  + psi(m1,m2,m3)*pow(lm,-4) * (  - pow(m1sq,4)*m2sq*pow(m3sq,2) -
         pow(m1sq,4)*pow(m2sq,2)*m3sq + pow(m1sq,5)*m2sq*m3sq );

      hh31 +=  + 31./800. + 1./20.*ln2;



  return pi162*hh31;
}

double hh(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh"<<std::endl;
    return 0.;}
  return h0(m1sq,m2sq,m3sq,xmu2)+qsq*h0p(m1sq,m2sq,m3sq,xmu2)
    +hbar(m1sq,m2sq,m3sq,qsq,0);
}

double hh1(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh1"<<std::endl;
    return 0.;}
  return h10(m1sq,m2sq,m3sq,xmu2)+qsq*h10p(m1sq,m2sq,m3sq,xmu2)
    +hbar(m1sq,m2sq,m3sq,qsq,1);
}

double hh21(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh21"<<std::endl;
    return 0.;}
  return h210(m1sq,m2sq,m3sq,xmu2)+qsq*h210p(m1sq,m2sq,m3sq,xmu2)
    +hbar(m1sq,m2sq,m3sq,qsq,2);
}

double hh31(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh31"<<std::endl;
    return 0.;}
  return h310(m1sq,m2sq,m3sq,xmu2)+qsq*h310p(m1sq,m2sq,m3sq,xmu2)
    +hbar(m1sq,m2sq,m3sq,qsq,3);
}

dcomplex zhh(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h0(m1sq,m2sq,m3sq,xmu2)+qsq*h0p(m1sq,m2sq,m3sq,xmu2)
      +zhbar(m1sq,m2sq,m3sq,qsq,0);
}

dcomplex zhh1(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h10(m1sq,m2sq,m3sq,xmu2)+qsq*h10p(m1sq,m2sq,m3sq,xmu2)
      +zhbar(m1sq,m2sq,m3sq,qsq,1);
}

dcomplex zhh21(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h210(m1sq,m2sq,m3sq,xmu2)+qsq*h210p(m1sq,m2sq,m3sq,xmu2)
      +zhbar(m1sq,m2sq,m3sq,qsq,2);
}

dcomplex zhh31(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h310(m1sq,m2sq,m3sq,xmu2)+qsq*h310p(m1sq,m2sq,m3sq,xmu2)
      +zhbar(m1sq,m2sq,m3sq,qsq,3);
}


double hhd(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hhd"<<std::endl;
    return 0.;}
  return h0p(m1sq,m2sq,m3sq,xmu2)
    +hbard(m1sq,m2sq,m3sq,qsq,0);
}

double hh1d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh1d"<<std::endl;
    return 0.;}
  return h10p(m1sq,m2sq,m3sq,xmu2)
    +hbard(m1sq,m2sq,m3sq,qsq,1);
}

double hh21d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh21d"<<std::endl;
    return 0.;}
  return h210p(m1sq,m2sq,m3sq,xmu2)
    +hbard(m1sq,m2sq,m3sq,qsq,2);
}

dcomplex zhhd(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h0p(m1sq,m2sq,m3sq,xmu2)
      +zhbard(m1sq,m2sq,m3sq,qsq,0);
}

dcomplex zhh1d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h10p(m1sq,m2sq,m3sq,xmu2)
      +zhbard(m1sq,m2sq,m3sq,qsq,1);
}

dcomplex zhh21d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double xmu2){
  return h210p(m1sq,m2sq,m3sq,xmu2)
      +zhbard(m1sq,m2sq,m3sq,qsq,2);
}
