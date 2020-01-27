// quenchedsunsetintegrals.cc is part of the CHIRON ChPT
// at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains all the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]
// has the extra ones needed for the partially quenched case with the
// extra first integer index indicating the power of the propagators added
// Derived in:
//
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 73} (2006) 074509
//  [hep-lat/0602003].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 72} (2005) 074502
//  [hep-lat/0506004].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 71} (2005) 094502
//  [hep-lat/0501014].
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 70} (2004) 111503
//  [hep-lat/0406017].


// so far only below threshold /////////////////////////////////////////
//
// the function phi(x,y) defined in Davydychev Tausk Nucl. Phys. B397(1993)123
// except multiplied with m3**2 to make it fully symmetric in the masses.
// also contains the barred functions and the full functions
// derivatives are hhd,hh1d and hh21d.
// only valid below threshold
//
// the cases with the Kählen function lm=0 are implemented separately
// however no expansions in the neighbourhood are done
// sometimes ones also looses quite some precision due to large cancellations



#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>

#include "quenchedsunsetintegrals.h"
#include "jbnumlib.h"

#ifndef DINTEGRAL
#define DINTEGRAL jbdgauss
#endif
// real integral with singularity (CHIRON v0.51 possible jbdcauch,jbdcauch2,
//                                                       jbdsing15,jbdsing21) 
//#ifndef SINTEGRAL
//#define SINTEGRAL jbdcauch
//#endif

const double pi = M_PI;
const double pi2 = pi*pi;
const double pi16 = 1./(16.*pi2);
const double pi162 = pi16*pi16;

// this is the real precision pi162 etc included
// assumes m1sq,..,qsq etc are all of order 1 and no crazy other things
double precisionquenchedsunsetintegrals = 1e-10;

void setprecisionquenchedsunsetintegrals(const double eps){
  precisionquenchedsunsetintegrals = eps;
}

double getprecisionquenchedsunsetintegrals(void){
  return precisionquenchedsunsetintegrals;
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

struct hbarquenchedcom{
  double sigm1,psqh,xm12h,xm22h,xm32h;
  int ih;
  int iprop;
};

// for above threshold
struct hbarpquenchedcom{
  double zm1,psqhh,xm12hh,xm22hh,xm32hh;
  int ihh;
  int iprop;
};
double dlamq(const double x,const double y,const double z){
  return  pow(x-y-z,2)-4.*y*z;
}

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
hbarquenchedcom hbardatq;


double hbarquenched2(const double x){
  const double stretch1 = 1.;  //  x near zero
  const double stretch2 = 1.;  // sigma near infinity
  double xx = pow(x,stretch1);
  double yy = pow(hbardatq.sigm1,stretch2);
  double hbar2t = (pow(hbardatq.sigm1,(stretch2-1.))*stretch2)
    *(pow(x,(stretch1-1.))*stretch1);
  double xm2 = sqrt(hbardatq.xm22h);
  double xm3 = sqrt(hbardatq.xm32h);
  double siglow = pow(xm2+xm3,2);
  double sigmat = 1./yy;
  double sigma = siglow*sigmat;
  hbar2t = hbar2t*sigmat*sigmat;
  double z1,z2,xk2,xlam,dlamm2,dlamm3,dk2dsig,dk2dm1,dk2dm1sig,dlam23,dk2dsig2;
  switch(hbardatq.iprop){
  case 1:
    hbar2t = hbar2t*sqrt(dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy));
    z1 = hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = pow(xx,hbardatq.ih) * pi162*(log(1.-z2/z1)+z2/z1);
    return hbar2t*xk2;
  case 2:
    hbar2t = hbar2t*sqrt(dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy));
    z1 = hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = pow(xx,hbardatq.ih) * pi162*(1.-xx)*z2*z2/(z1*z1*(z1-z2));
    return hbar2t*xk2;
  case 3:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
	(1.-1.*yy+xm3/xm2) -4*hbardatq.xm32h*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = log(1.-z2/z1)+z2/z1;
    dk2dsig = xx*z2*z2/(z1*z1*(z1-z2));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(xk2*dlamm2/(2.*sqrt(xlam))
	      +sqrt(xlam)*(1.+xm3/xm2)*sigmat*dk2dsig );
    return hbar2t;
  case 4:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm3 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbardatq.xm22h*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = log(1.-z2/z1)+z2/z1;
    dk2dsig = xx*z2*z2/(z1*z1*(z1-z2));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(xk2*dlamm3/(2.*sqrt(xlam))
	      +sqrt(xlam)*(1.+xm2/xm3)*sigmat*dk2dsig );
    return hbar2t;
  case 5:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
      (1.-yy+xm3/xm2) -4*hbardatq.xm32h*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    dk2dm1 =(1.-xx)*z2*z2/(z1*z1*(z1-z2));
    dk2dm1sig = xx*(1.-xx)*((-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,2)));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(dk2dm1*dlamm2/(2.*sqrt(xlam))
	  +sqrt(xlam)*(1.+xm3/xm2)*sigmat*dk2dm1sig );
    return hbar2t;
  case 6:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm3 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbardatq.xm22h*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    dk2dm1 =(1.-xx)*z2*z2/(z1*z1*(z1-z2));
    dk2dm1sig = xx*(1.-xx)*((-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,2)));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(dk2dm1*dlamm3/(2.*sqrt(xlam))
	 +sqrt(xlam)*(1.+xm2/xm3)*sigmat*dk2dm1sig );
    return hbar2t;
  case 7:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
           (1.-yy+xm3/xm2) -4*hbardatq.xm32h*yy*yy;
    dlamm3 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbardatq.xm22h*yy*yy;
    dlam23 = (siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)/(xm2*xm3)
      +2.*(1.-yy+xm3/xm2)*(1.-yy+xm2/xm3)-4.*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = log(1.-z2/z1)+z2/z1;
    dk2dsig = xx*z2*z2/(z1*z1*(z1-z2));
    dk2dsig2= xx*xx*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,2));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(
       -xk2*dlamm2*dlamm3/(4.*xlam*sqrt(xlam))
       +xk2*dlam23/(2.*sqrt(xlam))
       +dlamm3*(1.+xm3/xm2)*sigmat*dk2dsig/(2.*sqrt(xlam))
       +dlamm2*(1.+xm2/xm3)*sigmat*dk2dsig/(2.*sqrt(xlam))
       +sqrt(xlam)*sigmat/(2.*xm2*xm3)*dk2dsig
       +sqrt(xlam)*siglow/(xm2*xm3)*pow(sigmat,2)*dk2dsig2
					      );
    return hbar2t;
  case 8:
    xlam = dlamq(siglow,hbardatq.xm22h*yy,hbardatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
           (1.-yy+xm3/xm2) -4*hbardatq.xm32h*yy*yy;
    dlamm3 = 2.*(siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbardatq.xm22h*yy*yy;
    dlam23 = (siglow-hbardatq.xm22h*yy-hbardatq.xm32h*yy)/(xm2*xm3)
      +2.*(1.-1.*yy+xm3/xm2)*(1.-1.*yy+xm2/xm3)
      -4.*yy*yy;
    z1 =  hbardatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbardatq.psqh;
    xk2 = (1.-xx)*z2*z2/(z1*z1*(z1-z2));
    dk2dsig = xx*(1.-xx)*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,2));
    dk2dsig2=xx*xx*(1.-xx)*(12.*z1*z1-16.*z1*z2+6.*z2*z2)*z2*z2
        /(pow(z1,4)*pow(z1-z2,3));
    hbar2t = hbar2t*pow(xx,hbardatq.ih)*pi162*(
       -xk2*dlamm2*dlamm3/(4.*xlam*sqrt(xlam))
       +xk2*dlam23/(2.*sqrt(xlam))
       +dlamm3*(1.+xm3/xm2)*sigmat*dk2dsig/(2.*sqrt(xlam))
       +dlamm2*(1.+xm2/xm3)*sigmat*dk2dsig/(2.*sqrt(xlam))
       +sqrt(xlam)*sigmat/(2.*xm2*xm3)*dk2dsig
       +sqrt(xlam)*siglow/(xm2*xm3)*pow(sigmat,2)*dk2dsig2
						);
    return hbar2t;
  default:
    std::cout << "wrong iprop = "<<hbardatq.iprop<<" in hbarquenched2\n";
    return 0.;
  }
  return 0.;
}

double hbarquenched1(const double y){
  hbardatq.sigm1 = y;
  return DINTEGRAL(hbarquenched2,0.,1.,precisionquenchedsunsetintegrals/5.);
}

double hbar(const int iprop, const double xm12,const double xm22,
	    const double xm32, const double psq,const int i){
  hbardatq.ih = i;
  hbardatq.iprop = iprop;
  hbardatq.psqh = psq;
  hbardatq.xm12h = xm12;
  hbardatq.xm22h = xm22;
  hbardatq.xm32h = xm32;
  return DINTEGRAL(hbarquenched1,0.,1.,precisionquenchedsunsetintegrals);
}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
hbarquenchedcom hbarddatq;

double hbard2q(const double x){
  const double stretch1 = 2.;  //  x near zero
  const double stretch2 = 1.;  // sigma near infinity
  double xx = pow(x,stretch1);
  double yy = pow(hbarddatq.sigm1,stretch2);
  double hbard2q = (pow(hbarddatq.sigm1,(stretch2-1.))*stretch2)
    *(pow(x,(stretch1-1.))*stretch1);
  double xm2 = sqrt(hbarddatq.xm22h);
  double xm3 = sqrt(hbarddatq.xm32h);
  double siglow = pow(xm2+xm3,2);
  double sigma = siglow/yy;
  hbard2q = hbard2q/pow(yy,2);
  double z1,z2,xk2,xlam,dlamm2,dlamm3,dk2dsig,dk2dm1,dk2dm1sig,dlam23,dk2dsig2;
  switch(hbarddatq.iprop){
  case 1:
    hbard2q = hbard2q*sqrt(dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy));
    z1 = hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    xk2 = pow(xx,hbarddatq.ih) *xx*(1.-xx)* pi162*(-z2)/(z1*(z1-z2));
    return hbard2q*xk2;
  case 2:
    hbard2q = hbard2q*sqrt(dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy));
    z1 = hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    xk2 = pow(xx,hbarddatq.ih) * pi162*xx*pow(1.-xx,2)
      *(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2)); 
    return hbard2q*xk2;
  case 3:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm3/xm2) -4*hbarddatq.xm32h*yy*yy;
    z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    xk2 = -z2/(z1*(z1-z2));
    dk2dsig = xx*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
    hbard2q = hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*
        (xk2*dlamm2/(2.*sqrt(xlam))
         +sqrt(xlam)*(1.+xm3/xm2)/yy*dk2dsig );
    return hbard2q;
  case 4:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
      dlamm3 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
	(1.-yy+xm2/xm3) -4*hbarddatq.xm22h*yy*yy;
      z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
      z2 = xx*(1.-xx)*hbarddatq.psqh;
      xk2 = -z2/(z1*(z1-z2));
      dk2dsig = xx*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
      hbard2q =hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*
       (xk2*dlamm3/(2.*sqrt(xlam))
	+sqrt(xlam)*(1.+xm2/xm3)/yy*dk2dsig );
    return hbard2q;
  case 5:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm3/xm2) -4*hbarddatq.xm32h*yy*yy;
    z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    dk2dm1 =(1.-xx)*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
    dk2dm1sig = xx*(1.-xx)*(
           (-6.*z1+6.*z2)*z2/(pow(z1,3)*pow(z1-z2,2))
	   +2.*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,3)) );
    hbard2q = hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*
      (dk2dm1*dlamm2/(2.*sqrt(xlam))
       +sqrt(xlam)*(1.+xm3/xm2)/yy*dk2dm1sig );
    return hbard2q;
  case 6:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
    dlamm3 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbarddatq.xm22h*yy*yy;
    z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    dk2dm1 = (1.-xx)*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
    dk2dm1sig =  xx*(1.-xx)*(
	     (-6.*z1+6.*z2)*z2/(pow(z1,3)*pow(z1-z2,2))
	     +2.*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,3)) );
    hbard2q = hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*
      (dk2dm1*dlamm3/(2.*sqrt(xlam))
       +sqrt(xlam)*(1.+xm2/xm3)/yy*dk2dm1sig );
    return hbard2q;
  case 7:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm3/xm2) -4*hbarddatq.xm32h*yy*yy;
    dlamm3 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbarddatq.xm22h*yy*yy;
    dlam23 = (siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)/(xm2*xm3)
      +2.*(1.-yy+xm3/xm2)*(1.-yy+xm2/xm3)
      -4.*yy*yy;
    z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    xk2 = -z2/(z1*(z1-z2));
    dk2dsig = xx*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
    dk2dsig2= xx*xx*(
      (-6.*z1+6.*z2)*z2/(pow(z1,3)*pow(z1-z2,2))
      +2.*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,3)) );
    hbard2q = hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*(
       -xk2*dlamm2*dlamm3/(4.*xlam*sqrt(xlam))
       +xk2*dlam23/(2.*sqrt(xlam))
       +dlamm3*(1.+xm3/xm2)/yy*dk2dsig/(2.*sqrt(xlam))
       +dlamm2*(1.+xm2/xm3)/yy*dk2dsig/(2.*sqrt(xlam))
       +sqrt(xlam)/yy/(2.*xm2*xm3)*dk2dsig
       +sqrt(xlam)*siglow/(xm2*xm3)/(yy*yy)*dk2dsig2
						   );
    return hbard2q;
  case 8:
    xlam = dlamq(siglow,hbarddatq.xm22h*yy,hbarddatq.xm32h*yy);
    dlamm2 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm3/xm2) -4*hbarddatq.xm32h*yy*yy;
    dlamm3 = 2.*(siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)*
      (1.-yy+xm2/xm3) -4*hbarddatq.xm22h*yy*yy;
    dlam23 = (siglow-hbarddatq.xm22h*yy-hbarddatq.xm32h*yy)/(xm2*xm3)
      +2.*(1.-yy+xm3/xm2)*(1.-yy+xm2/xm3)
      -4.*yy*yy;
    z1 =  hbarddatq.xm12h*(1.-xx)+sigma*xx;
    z2 = xx*(1.-xx)*hbarddatq.psqh;
    xk2 =  (1.-xx)*(2.*z1-z2)*z2/(z1*z1*pow(z1-z2,2));
    dk2dsig = xx*(1.-xx)*(
	  (-6.*z1+6.*z2)*z2/(pow(z1,3)*pow(z1-z2,2))
	  +2.*(-3.*z1+2.*z2)*z2*z2/(pow(z1,3)*pow(z1-z2,3)));
    dk2dsig2=xx*xx*(1.-xx)*(
    (24.*z1*z1-48.*z1*z2+24.*z2*z2)*z2/(pow(z1,4)*pow(z1-z2,3))
    +3.*(12.*z1*z1-16.*z1*z2+6.*z2*z2)*z2*z2/(pow(z1*(z1-z2),4))
			    );
    hbard2q = hbard2q*pow(xx,hbarddatq.ih)*pi162*xx*(1.-xx)*(
       -xk2*dlamm2*dlamm3/(4.*xlam*sqrt(xlam))
       +xk2*dlam23/(2.*sqrt(xlam))
       +dlamm3*(1.+xm3/xm2)/yy*dk2dsig/(2.*sqrt(xlam))
       +dlamm2*(1.+xm2/xm3)/yy*dk2dsig/(2.*sqrt(xlam))
       +sqrt(xlam)/yy/(2.*xm2*xm3)*dk2dsig
       +sqrt(xlam)*siglow/(xm2*xm3)/(yy*yy)*dk2dsig2
						 );
    return hbard2q;
  default:
    std::cout << "wrong iprop = "<<hbardatq.iprop<<" in hbarquenched2\n";
    return 0.;
  }
  return 0.;
}

double hbard1q(const double y){
  hbarddatq.sigm1 = y;
  return DINTEGRAL(hbard2q,0.,1.,precisionquenchedsunsetintegrals/5.);
}

double hbard(const int iprop, const double xm12,const double xm22,
	     const double xm32, const double psq,const int i){
  hbarddatq.ih = i;
  hbarddatq.iprop = iprop;
  hbarddatq.psqh = psq;
  hbarddatq.xm12h = xm12;
  hbarddatq.xm22h = xm22;
  hbarddatq.xm32h = xm32;
  return DINTEGRAL(hbard1q,0.,1.,precisionquenchedsunsetintegrals);
}
//xxxxxxxxxxxxx lm=0 cases xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// these functions ASSUME lm = 0
double h0sing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h0x1,h0x2,h0x3,h0x4,h0x5,h0x6,h0x7,h0x8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0x1 =
       + ln3 * (  - pow(mm3,2) );

      h0x1 +=  + pow(ln3,2) * ( 1./2.*pow(mm3,2) );

      h0x1 +=  + ln2 * (  - pow(mm2,2) );

      h0x1 +=  + ln2*ln3 * (  - mm2*mm3 );

      h0x1 +=  + pow(ln2,2) * ( 1./2.*pow(mm2,2) );

      h0x1 +=  + ln1 * (  - pow(mm3,2) - 2*mm2*mm3 - pow(mm2,2) );

      h0x1 +=  + ln1*ln3 * ( pow(mm3,2) + mm2*mm3 );

      h0x1 +=  + ln1*ln2 * ( mm2*mm3 + pow(mm2,2) );

      h0x1 +=  + pow(ln1,2) * ( 1./2.*pow(mm3,2) + mm2*mm3 + 1./2.*pow(
         mm2,2) );

      h0x1 +=  + 3*pow(mm3,2) + 1./6.*pow(mm3,2)*pi2 + 3*mm2*mm3 + 1./6.
         *mm2*mm3*pi2 + 3*pow(mm2,2) + 1./6.*pow(mm2,2)*pi2;
      return h0x1*pi162;
    case 2:
      h0x2 =
       + ln3 * (  - mm3 );

      h0x2 +=  + ln2 * (  - mm2 );

      h0x2 +=  + ln2*ln3 * (  - 1./2.*mm3 - 1./2.*mm2 );

      h0x2 +=  + ln1 * ( 2*mm3 + 2*mm2 );

      h0x2 +=  + ln1*ln3 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x2 +=  + ln1*ln2 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x2 +=  + pow(ln1,2) * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x2 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm2 + 1./12.*mm2*
         pi2;
      return h0x2*pi162/mm1;
    case 3:
      h0x3 =
       + ln3 * ( pow(mm2,-1)*pow(mm3,2) + mm3 );

      h0x3 +=  + ln2 * ( 2*mm3 + 2*mm2 );

      h0x3 +=  + ln2*ln3 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x3 +=  + pow(ln2,2) * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x3 +=  + ln1 * (  - pow(mm2,-1)*pow(mm3,2) - 2*mm3 - mm2 );

      h0x3 +=  + ln1*ln3 * (  - 1./2.*mm3 - 1./2.*mm2 );

      h0x3 +=  + ln1*ln2 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x3 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm2 + 1./12.*mm2*
         pi2;
      return h0x3*pi162/mm1;
    case 4:
      h0x4 =
       + ln3 * ( 2*mm3 + 2*mm2 );

      h0x4 +=  + pow(ln3,2) * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x4 +=  + ln2 * ( mm2 + pow(mm2,2)*pow(mm3,-1) );

      h0x4 +=  + ln2*ln3 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x4 +=  + ln1 * (  - mm3 - 2*mm2 - pow(mm2,2)*pow(mm3,-1) );

      h0x4 +=  + ln1*ln3 * ( 1./2.*mm3 + 1./2.*mm2 );

      h0x4 +=  + ln1*ln2 * (  - 1./2.*mm3 - 1./2.*mm2 );

      h0x4 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm2 + 1./12.*mm2*
         pi2;
      return h0x4*pi162/mm1;
    case 5:
      h0x5 =
       + ln3 * (  - 1./6.*pow(mm2,-2)*pow(mm3,4) - 1./3.*pow(mm2,-1)*
         pow(mm3,3) - 1./6.*pow(mm3,2) );

      h0x5 +=  + ln2 * ( 1./2.*pow(mm3,2) + 4./3.*mm2*mm3 + 7./6.*pow(
         mm2,2) + 1./3.*pow(mm2,3)*pow(mm3,-1) );

      h0x5 +=  + ln1 * ( 1./6.*pow(mm2,-2)*pow(mm3,4) + 1./3.*pow(
         mm2,-1)*pow(mm3,3) - 1./3.*pow(mm3,2) - 4./3.*mm2*mm3 - 7./6.*
         pow(mm2,2) - 1./3.*pow(mm2,3)*pow(mm3,-1) );

      h0x5 +=  - 1./3.*pow(mm2,-1)*pow(mm3,3) - pow(mm3,2) - mm2*mm3 - 
         1./3.*pow(mm2,2);
      return h0x5*pi162/pow(mm1,4);
    case 6:
   h0x6 =
       + ln3 * ( 1./3.*pow(mm2,-1)*pow(mm3,3) + 7./6.*pow(mm3,2) + 4./3.
         *mm2*mm3 + 1./2.*pow(mm2,2) );

      h0x6 +=  + ln2 * (  - 1./6.*pow(mm2,2) - 1./3.*pow(mm2,3)*pow(
         mm3,-1) - 1./6.*pow(mm2,4)*pow(mm3,-2) );

      h0x6 +=  + ln1 * (  - 1./3.*pow(mm2,-1)*pow(mm3,3) - 7./6.*pow(
         mm3,2) - 4./3.*mm2*mm3 - 1./3.*pow(mm2,2) + 1./3.*pow(mm2,3)*
         pow(mm3,-1) + 1./6.*pow(mm2,4)*pow(mm3,-2) );

      h0x6 +=  - 1./3.*pow(mm3,2) - mm2*mm3 - pow(mm2,2) - 1./3.*pow(
         mm2,3)*pow(mm3,-1);
      return h0x6*pi162/pow(mm1,4);
    case 7:
   h0x7 =
       + ln3 * ( 1./6.*pow(mm2,-2)*pow(mm3,4) + pow(mm2,-1)*pow(mm3,3)
          + 2*pow(mm3,2) + 5./3.*mm2*mm3 + 1./2.*pow(mm2,2) );

      h0x7 +=  + ln2 * ( 1./2.*pow(mm3,2) + 5./3.*mm2*mm3 + 2*pow(
         mm2,2) + pow(mm2,3)*pow(mm3,-1) + 1./6.*pow(mm2,4)*pow(mm3,-2)
          );

      h0x7 +=  + ln1 * (  - 1./6.*pow(mm2,-2)*pow(mm3,4) - pow(mm2,-1)*
         pow(mm3,3) - 5./2.*pow(mm3,2) - 10./3.*mm2*mm3 - 5./2.*pow(
         mm2,2) - pow(mm2,3)*pow(mm3,-1) - 1./6.*pow(mm2,4)*pow(mm3,-2)
          );

      h0x7 +=  + 1./3.*pow(mm2,-1)*pow(mm3,3) + 4./3.*pow(mm3,2) + 2*
         mm2*mm3 + 4./3.*pow(mm2,2) + 1./3.*pow(mm2,3)*pow(mm3,-1);
      return h0x7*pi162/pow(mm1,4);
    case 8:
      h0x8 =
       + ln3 * (  - 1./15.*pow(mm2,-3)*pow(mm3,5) - 11./30.*pow(mm2,-2)
         *pow(mm3,4) - 13./15.*pow(mm2,-1)*pow(mm3,3) - 16./15.*pow(
         mm3,2) - 2./3.*mm2*mm3 - 1./6.*pow(mm2,2) );

      h0x8 +=  + ln2 * (  - 1./6.*pow(mm3,2) - 2./3.*mm2*mm3 - 16./15.*
         pow(mm2,2) - 13./15.*pow(mm2,3)*pow(mm3,-1) - 11./30.*pow(
         mm2,4)*pow(mm3,-2) - 1./15.*pow(mm2,5)*pow(mm3,-3) );

      h0x8 +=  + ln1 * ( 1./15.*pow(mm2,-3)*pow(mm3,5) + 11./30.*pow(
         mm2,-2)*pow(mm3,4) + 13./15.*pow(mm2,-1)*pow(mm3,3) + 37./30.*
         pow(mm3,2) + 4./3.*mm2*mm3 + 37./30.*pow(mm2,2) + 13./15.*pow(
         mm2,3)*pow(mm3,-1) + 11./30.*pow(mm2,4)*pow(mm3,-2) + 1./15.*
         pow(mm2,5)*pow(mm3,-3) );

      h0x8 +=  - 2./15.*pow(mm2,-2)*pow(mm3,4) - 2./3.*pow(mm2,-1)*pow(
         mm3,3) - 22./15.*pow(mm3,2) - 28./15.*mm2*mm3 - 22./15.*pow(
         mm2,2) - 2./3.*pow(mm2,3)*pow(mm3,-1) - 2./15.*pow(mm2,4)*pow(
         mm3,-2);
      return h0x8*pi162/pow(mm1,6);
    default:
      std::cout << "h0sing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0x1 =
       + ln3 * (  - pow(mm3,2) );

      h0x1 +=  + pow(ln3,2) * ( 1./2.*pow(mm3,2) );

      h0x1 +=  + ln2 * (  - pow(mm3,2) - 2*mm1*mm3 - pow(mm1,2) );

      h0x1 +=  + ln2*ln3 * ( pow(mm3,2) + mm1*mm3 );

      h0x1 +=  + pow(ln2,2) * ( 1./2.*pow(mm3,2) + mm1*mm3 + 1./2.*pow(
         mm1,2) );

      h0x1 +=  + ln1 * (  - pow(mm1,2) );

      h0x1 +=  + ln1*ln3 * (  - mm1*mm3 );

      h0x1 +=  + ln1*ln2 * ( mm1*mm3 + pow(mm1,2) );

      h0x1 +=  + pow(ln1,2) * ( 1./2.*pow(mm1,2) );

      h0x1 +=  + 3*pow(mm3,2) + 1./6.*pow(mm3,2)*pi2 + 3*mm1*mm3 + 1./6.
         *mm1*mm3*pi2 + 3*pow(mm1,2) + 1./6.*pow(mm1,2)*pi2;
      return h0x1*pi162;
    case 2:
      h0x2 =
       + ln3 * ( pow(mm1,-1)*pow(mm3,2) + mm3 );

      h0x2 +=  + ln2 * (  - pow(mm1,-1)*pow(mm3,2) - 2*mm3 - mm1 );

      h0x2 +=  + ln2*ln3 * (  - 1./2.*mm3 - 1./2.*mm1 );

      h0x2 +=  + ln1 * ( 2*mm3 + 2*mm1 );

      h0x2 +=  + ln1*ln3 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x2 +=  + ln1*ln2 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x2 +=  + pow(ln1,2) * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x2 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x2*pi162/mm2;
    case 3:
      h0x3 =
       + ln3 * (  - mm3 );

      h0x3 +=  + ln2 * ( 2*mm3 + 2*mm1 );

      h0x3 +=  + ln2*ln3 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x3 +=  + pow(ln2,2) * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x3 +=  + ln1 * (  - mm1 );

      h0x3 +=  + ln1*ln3 * (  - 1./2.*mm3 - 1./2.*mm1 );

      h0x3 +=  + ln1*ln2 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x3 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x3*pi162/mm2;
    case 4:
      h0x4 =
       + ln3 * ( 2*mm3 + 2*mm1 );

      h0x4 +=  + pow(ln3,2) * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x4 +=  + ln2 * (  - mm3 - 2*mm1 - pow(mm1,2)*pow(mm3,-1) );

      h0x4 +=  + ln2*ln3 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x4 +=  + ln1 * ( mm1 + pow(mm1,2)*pow(mm3,-1) );

      h0x4 +=  + ln1*ln3 * ( 1./2.*mm3 + 1./2.*mm1 );

      h0x4 +=  + ln1*ln2 * (  - 1./2.*mm3 - 1./2.*mm1 );

      h0x4 +=  + 1./2.*mm3 + 1./12.*mm3*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x4*pi162/mm2;
    case 5:
      h0x5 =
       + ln3 * (  - 1./6.*pow(mm1,-2)*pow(mm3,4) - 1./3.*pow(mm1,-1)*
         pow(mm3,3) - 1./6.*pow(mm3,2) );

      h0x5 +=  + ln2 * ( 1./6.*pow(mm1,-2)*pow(mm3,4) + 1./3.*pow(
         mm1,-1)*pow(mm3,3) - 1./3.*pow(mm3,2) - 4./3.*mm1*mm3 - 7./6.*
         pow(mm1,2) - 1./3.*pow(mm1,3)*pow(mm3,-1) );

      h0x5 +=  + ln1 * ( 1./2.*pow(mm3,2) + 4./3.*mm1*mm3 + 7./6.*pow(
         mm1,2) + 1./3.*pow(mm1,3)*pow(mm3,-1) );

      h0x5 +=  - 1./3.*pow(mm1,-1)*pow(mm3,3) - pow(mm3,2) - mm1*mm3 - 
         1./3.*pow(mm1,2);
      return h0x5*pi162/pow(mm2,4);
    case 6:
      h0x6 =
       + ln3 * ( 1./6.*pow(mm1,-2)*pow(mm3,4) + pow(mm1,-1)*pow(mm3,3)
          + 2*pow(mm3,2) + 5./3.*mm1*mm3 + 1./2.*pow(mm1,2) );

      h0x6 +=  + ln2 * (  - 1./6.*pow(mm1,-2)*pow(mm3,4) - pow(mm1,-1)*
         pow(mm3,3) - 5./2.*pow(mm3,2) - 10./3.*mm1*mm3 - 5./2.*pow(
         mm1,2) - pow(mm1,3)*pow(mm3,-1) - 1./6.*pow(mm1,4)*pow(mm3,-2)
          );

      h0x6 +=  + ln1 * ( 1./2.*pow(mm3,2) + 5./3.*mm1*mm3 + 2*pow(
         mm1,2) + pow(mm1,3)*pow(mm3,-1) + 1./6.*pow(mm1,4)*pow(mm3,-2)
          );

      h0x6 +=  + 1./3.*pow(mm1,-1)*pow(mm3,3) + 4./3.*pow(mm3,2) + 2*
         mm1*mm3 + 4./3.*pow(mm1,2) + 1./3.*pow(mm1,3)*pow(mm3,-1);
      return h0x6*pi162/pow(mm2,4);
    case 7:
      h0x7 =
       + ln3 * ( 1./3.*pow(mm1,-1)*pow(mm3,3) + 7./6.*pow(mm3,2) + 4./3.
         *mm1*mm3 + 1./2.*pow(mm1,2) );

      h0x7 +=  + ln2 * (  - 1./3.*pow(mm1,-1)*pow(mm3,3) - 7./6.*pow(
         mm3,2) - 4./3.*mm1*mm3 - 1./3.*pow(mm1,2) + 1./3.*pow(mm1,3)*
         pow(mm3,-1) + 1./6.*pow(mm1,4)*pow(mm3,-2) );

      h0x7 +=  + ln1 * (  - 1./6.*pow(mm1,2) - 1./3.*pow(mm1,3)*pow(
         mm3,-1) - 1./6.*pow(mm1,4)*pow(mm3,-2) );

      h0x7 +=  - 1./3.*pow(mm3,2) - mm1*mm3 - pow(mm1,2) - 1./3.*pow(
         mm1,3)*pow(mm3,-1);
      return h0x7*pi162/pow(mm2,4);
    case 8:
      h0x8 =
       + ln3 * (  - 1./15.*pow(mm1,-3)*pow(mm3,5) - 11./30.*pow(mm1,-2)
         *pow(mm3,4) - 13./15.*pow(mm1,-1)*pow(mm3,3) - 16./15.*pow(
         mm3,2) - 2./3.*mm1*mm3 - 1./6.*pow(mm1,2) );

      h0x8 +=  + ln2 * ( 1./15.*pow(mm1,-3)*pow(mm3,5) + 11./30.*pow(
         mm1,-2)*pow(mm3,4) + 13./15.*pow(mm1,-1)*pow(mm3,3) + 37./30.*
         pow(mm3,2) + 4./3.*mm1*mm3 + 37./30.*pow(mm1,2) + 13./15.*pow(
         mm1,3)*pow(mm3,-1) + 11./30.*pow(mm1,4)*pow(mm3,-2) + 1./15.*
         pow(mm1,5)*pow(mm3,-3) );

      h0x8 +=  + ln1 * (  - 1./6.*pow(mm3,2) - 2./3.*mm1*mm3 - 16./15.*
         pow(mm1,2) - 13./15.*pow(mm1,3)*pow(mm3,-1) - 11./30.*pow(
         mm1,4)*pow(mm3,-2) - 1./15.*pow(mm1,5)*pow(mm3,-3) );

      h0x8 +=  - 2./15.*pow(mm1,-2)*pow(mm3,4) - 2./3.*pow(mm1,-1)*pow(
         mm3,3) - 22./15.*pow(mm3,2) - 28./15.*mm1*mm3 - 22./15.*pow(
         mm1,2) - 2./3.*pow(mm1,3)*pow(mm3,-1) - 2./15.*pow(mm1,4)*pow(
         mm3,-2);
      return h0x8*pi162/pow(mm2,6);
    default:
      std::cout << "h0sing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0x1 =
       + ln3 * (  - pow(mm2,2) - 2*mm1*mm2 - pow(mm1,2) );

      h0x1 +=  + pow(ln3,2) * ( 1./2.*pow(mm2,2) + mm1*mm2 + 1./2.*pow(
         mm1,2) );

      h0x1 +=  + ln2 * (  - pow(mm2,2) );

      h0x1 +=  + ln2*ln3 * ( pow(mm2,2) + mm1*mm2 );

      h0x1 +=  + pow(ln2,2) * ( 1./2.*pow(mm2,2) );

      h0x1 +=  + ln1 * (  - pow(mm1,2) );

      h0x1 +=  + ln1*ln3 * ( mm1*mm2 + pow(mm1,2) );

      h0x1 +=  + ln1*ln2 * (  - mm1*mm2 );

      h0x1 +=  + pow(ln1,2) * ( 1./2.*pow(mm1,2) );

      h0x1 +=  + 3*pow(mm2,2) + 1./6.*pow(mm2,2)*pi2 + 3*mm1*mm2 + 1./6.
         *mm1*mm2*pi2 + 3*pow(mm1,2) + 1./6.*pow(mm1,2)*pi2;
      return h0x1*pi162;
    case 2:
      h0x2 =
       + ln3 * (  - pow(mm1,-1)*pow(mm2,2) - 2*mm2 - mm1 );

      h0x2 +=  + ln2 * ( pow(mm1,-1)*pow(mm2,2) + mm2 );

      h0x2 +=  + ln2*ln3 * (  - 1./2.*mm2 - 1./2.*mm1 );

      h0x2 +=  + ln1 * ( 2*mm2 + 2*mm1 );

      h0x2 +=  + ln1*ln3 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x2 +=  + ln1*ln2 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x2 +=  + pow(ln1,2) * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x2 +=  + 1./2.*mm2 + 1./12.*mm2*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x2*pi162/mm3;
    case 3:
      h0x3 =
       + ln3 * (  - mm2 - 2*mm1 - pow(mm1,2)*pow(mm2,-1) );

      h0x3 +=  + ln2 * ( 2*mm2 + 2*mm1 );

      h0x3 +=  + ln2*ln3 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x3 +=  + pow(ln2,2) * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x3 +=  + ln1 * ( mm1 + pow(mm1,2)*pow(mm2,-1) );

      h0x3 +=  + ln1*ln3 * (  - 1./2.*mm2 - 1./2.*mm1 );

      h0x3 +=  + ln1*ln2 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x3 +=  + 1./2.*mm2 + 1./12.*mm2*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x3*pi162/mm3;
    case 4:
      h0x4 =
       + ln3 * ( 2*mm2 + 2*mm1 );

      h0x4 +=  + pow(ln3,2) * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x4 +=  + ln2 * (  - mm2 );

      h0x4 +=  + ln2*ln3 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x4 +=  + ln1 * (  - mm1 );

      h0x4 +=  + ln1*ln3 * ( 1./2.*mm2 + 1./2.*mm1 );

      h0x4 +=  + ln1*ln2 * (  - 1./2.*mm2 - 1./2.*mm1 );

      h0x4 +=  + 1./2.*mm2 + 1./12.*mm2*pi2 + 1./2.*mm1 + 1./12.*mm1*
         pi2;
      return h0x4*pi162/mm3;
    case 5:
      h0x5 =
       + ln3 * (  - 1./6.*pow(mm1,-2)*pow(mm2,4) - pow(mm1,-1)*pow(
         mm2,3) - 5./2.*pow(mm2,2) - 10./3.*mm1*mm2 - 5./2.*pow(mm1,2)
          - pow(mm1,3)*pow(mm2,-1) - 1./6.*pow(mm1,4)*pow(mm2,-2) );

      h0x5 +=  + ln2 * ( 1./6.*pow(mm1,-2)*pow(mm2,4) + pow(mm1,-1)*
         pow(mm2,3) + 2*pow(mm2,2) + 5./3.*mm1*mm2 + 1./2.*pow(mm1,2) )
         ;

      h0x5 +=  + ln1 * ( 1./2.*pow(mm2,2) + 5./3.*mm1*mm2 + 2*pow(
         mm1,2) + pow(mm1,3)*pow(mm2,-1) + 1./6.*pow(mm1,4)*pow(mm2,-2)
          );

      h0x5 +=  + 1./3.*pow(mm1,-1)*pow(mm2,3) + 4./3.*pow(mm2,2) + 2*
         mm1*mm2 + 4./3.*pow(mm1,2) + 1./3.*pow(mm1,3)*pow(mm2,-1);
      return h0x5*pi162/pow(mm3,4);
    case 6:
      h0x6 =
       + ln3 * ( 1./6.*pow(mm1,-2)*pow(mm2,4) + 1./3.*pow(mm1,-1)*pow(
         mm2,3) - 1./3.*pow(mm2,2) - 4./3.*mm1*mm2 - 7./6.*pow(mm1,2)
          - 1./3.*pow(mm1,3)*pow(mm2,-1) );

      h0x6 +=  + ln2 * (  - 1./6.*pow(mm1,-2)*pow(mm2,4) - 1./3.*pow(
         mm1,-1)*pow(mm2,3) - 1./6.*pow(mm2,2) );

      h0x6 +=  + ln1 * ( 1./2.*pow(mm2,2) + 4./3.*mm1*mm2 + 7./6.*pow(
         mm1,2) + 1./3.*pow(mm1,3)*pow(mm2,-1) );

      h0x6 +=  - 1./3.*pow(mm1,-1)*pow(mm2,3) - pow(mm2,2) - mm1*mm2 - 
         1./3.*pow(mm1,2);
      return h0x6*pi162/pow(mm3,4);
    case 7:
      h0x7 =
       + ln3 * (  - 1./3.*pow(mm1,-1)*pow(mm2,3) - 7./6.*pow(mm2,2) - 4.
         /3.*mm1*mm2 - 1./3.*pow(mm1,2) + 1./3.*pow(mm1,3)*pow(mm2,-1)
          + 1./6.*pow(mm1,4)*pow(mm2,-2) );

      h0x7 +=  + ln2 * ( 1./3.*pow(mm1,-1)*pow(mm2,3) + 7./6.*pow(
         mm2,2) + 4./3.*mm1*mm2 + 1./2.*pow(mm1,2) );

      h0x7 +=  + ln1 * (  - 1./6.*pow(mm1,2) - 1./3.*pow(mm1,3)*pow(
         mm2,-1) - 1./6.*pow(mm1,4)*pow(mm2,-2) );

      h0x7 +=  - 1./3.*pow(mm2,2) - mm1*mm2 - pow(mm1,2) - 1./3.*pow(
         mm1,3)*pow(mm2,-1);
      return h0x7*pi162/pow(mm3,4);
    case 8:
      h0x8 =
       + ln3 * ( 1./15.*pow(mm1,-3)*pow(mm2,5) + 11./30.*pow(mm1,-2)*
         pow(mm2,4) + 13./15.*pow(mm1,-1)*pow(mm2,3) + 37./30.*pow(
         mm2,2) + 4./3.*mm1*mm2 + 37./30.*pow(mm1,2) + 13./15.*pow(
         mm1,3)*pow(mm2,-1) + 11./30.*pow(mm1,4)*pow(mm2,-2) + 1./15.*
         pow(mm1,5)*pow(mm2,-3) );

      h0x8 +=  + ln2 * (  - 1./15.*pow(mm1,-3)*pow(mm2,5) - 11./30.*
         pow(mm1,-2)*pow(mm2,4) - 13./15.*pow(mm1,-1)*pow(mm2,3) - 16./
         15.*pow(mm2,2) - 2./3.*mm1*mm2 - 1./6.*pow(mm1,2) );

      h0x8 +=  + ln1 * (  - 1./6.*pow(mm2,2) - 2./3.*mm1*mm2 - 16./15.*
         pow(mm1,2) - 13./15.*pow(mm1,3)*pow(mm2,-1) - 11./30.*pow(
         mm1,4)*pow(mm2,-2) - 1./15.*pow(mm1,5)*pow(mm2,-3) );

      h0x8 +=  - 2./15.*pow(mm1,-2)*pow(mm2,4) - 2./3.*pow(mm1,-1)*pow(
         mm2,3) - 22./15.*pow(mm2,2) - 28./15.*mm1*mm2 - 22./15.*pow(
         mm1,2) - 2./3.*pow(mm1,3)*pow(mm2,-1) - 2./15.*pow(mm1,4)*pow(
         mm2,-2);
      return h0x8*pi162/pow(mm3,6);
    default:
      std::cout << "h0sing wrong iprop ="<<iprop<<'\n';
    }
  }
  std::cout << "funny masses in h0sing\n";
  return 0.;
}

double h0psing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h0px1,h0px2,h0px3,h0px4,h0px5,h0px6,h0px7,h0px8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0px1 =
       + ln3 * (  - 1./6.*pow(mm2,-1)*pow(mm3,2) );

      h0px1 +=  + ln2 * (  - 1./6.*pow(mm2,2)*pow(mm3,-1) );

      h0px1 +=  + ln1 * ( 1./6.*pow(mm2,-1)*pow(mm3,2) + 1./2.*mm3 + 1./
         2.*mm2 + 1./6.*pow(mm2,2)*pow(mm3,-1) );

      h0px1 +=  + 7./24.*mm3 + 7./24.*mm2;
      return h0px1*pi162/mm1;
    case 2:
      h0px2 =
       + ln3 * ( 1./20.*pow(mm2,-2)*pow(mm3,3) + 1./12.*pow(mm2,-1)*
         pow(mm3,2) );

      h0px2 +=  + ln2 * ( 1./12.*pow(mm2,2)*pow(mm3,-1) + 1./20.*pow(
         mm2,3)*pow(mm3,-2) );

      h0px2 +=  + ln1 * (  - 1./20.*pow(mm2,-2)*pow(mm3,3) - 1./12.*
         pow(mm2,-1)*pow(mm3,2) - 1./12.*pow(mm2,2)*pow(mm3,-1) - 1./20.
         *pow(mm2,3)*pow(mm3,-2) );

      h0px2 +=  + 1./10.*pow(mm2,-1)*pow(mm3,2) + 11./30.*mm3 + 11./30.
         *mm2 + 1./10.*pow(mm2,2)*pow(mm3,-1);
      return h0px2*pi162/pow(mm1,3);
    case 3:
      h0px3 =
       + ln3 * ( 1./30.*pow(mm2,-3)*pow(mm3,4) + 7./60.*pow(mm2,-2)*
         pow(mm3,3) + 1./12.*pow(mm2,-1)*pow(mm3,2) );

      h0px3 +=  + ln2 * (  - 1./6.*mm3 - 1./3.*mm2 - 13./60.*pow(mm2,2)
         *pow(mm3,-1) - 1./20.*pow(mm2,3)*pow(mm3,-2) );

      h0px3 +=  + ln1 * (  - 1./30.*pow(mm2,-3)*pow(mm3,4) - 7./60.*
         pow(mm2,-2)*pow(mm3,3) - 1./12.*pow(mm2,-1)*pow(mm3,2) + 1./6.
         *mm3 + 1./3.*mm2 + 13./60.*pow(mm2,2)*pow(mm3,-1) + 1./20.*
         pow(mm2,3)*pow(mm3,-2) );

      h0px3 +=  + 1./15.*pow(mm2,-2)*pow(mm3,3) + 1./5.*pow(mm2,-1)*
         pow(mm3,2) + 1./10.*mm3 - 2./15.*mm2 - 1./10.*pow(mm2,2)*pow(
         mm3,-1);
      return h0px3*pi162/pow(mm1,3);
    case 4:
      h0px4 =
       + ln3 * (  - 1./20.*pow(mm2,-2)*pow(mm3,3) - 13./60.*pow(mm2,-1)
         *pow(mm3,2) - 1./3.*mm3 - 1./6.*mm2 );

      h0px4 +=  + ln2 * ( 1./12.*pow(mm2,2)*pow(mm3,-1) + 7./60.*pow(
         mm2,3)*pow(mm3,-2) + 1./30.*pow(mm2,4)*pow(mm3,-3) );

      h0px4 +=  + ln1 * ( 1./20.*pow(mm2,-2)*pow(mm3,3) + 13./60.*pow(
         mm2,-1)*pow(mm3,2) + 1./3.*mm3 + 1./6.*mm2 - 1./12.*pow(mm2,2)
         *pow(mm3,-1) - 7./60.*pow(mm2,3)*pow(mm3,-2) - 1./30.*pow(
         mm2,4)*pow(mm3,-3) );

      h0px4 +=  - 1./10.*pow(mm2,-1)*pow(mm3,2) - 2./15.*mm3 + 1./10.*
         mm2 + 1./5.*pow(mm2,2)*pow(mm3,-1) + 1./15.*pow(mm2,3)*pow(
         mm3,-2);
      return h0px4*pi162/pow(mm1,3);
    case 5:
      h0px5 =
       + ln3 * (  - 3./140.*pow(mm2,-4)*pow(mm3,5) - 37./420.*pow(
         mm2,-3)*pow(mm3,4) - 2./15.*pow(mm2,-2)*pow(mm3,3) - 1./15.*
         pow(mm2,-1)*pow(mm3,2) );

      h0px5 +=  + ln2 * ( 1./12.*mm3 + 13./60.*mm2 + 7./30.*pow(mm2,2)*
         pow(mm3,-1) + 9./70.*pow(mm2,3)*pow(mm3,-2) + 1./35.*pow(
         mm2,4)*pow(mm3,-3) );

      h0px5 +=  + ln1 * ( 3./140.*pow(mm2,-4)*pow(mm3,5) + 37./420.*
         pow(mm2,-3)*pow(mm3,4) + 2./15.*pow(mm2,-2)*pow(mm3,3) + 1./15.
         *pow(mm2,-1)*pow(mm3,2) - 1./12.*mm3 - 13./60.*mm2 - 7./30.*
         pow(mm2,2)*pow(mm3,-1) - 9./70.*pow(mm2,3)*pow(mm3,-2) - 1./35.
         *pow(mm2,4)*pow(mm3,-3) );

      h0px5 +=  - 3./70.*pow(mm2,-3)*pow(mm3,4) - 13./84.*pow(mm2,-2)*
         pow(mm3,3) - 27./140.*pow(mm2,-1)*pow(mm3,2) + 3./140.*mm3 + 
         23./84.*mm2 + 8./35.*pow(mm2,2)*pow(mm3,-1) + 2./35.*pow(
         mm2,3)*pow(mm3,-2);
      return h0px5*pi162/pow(mm1,5);
    case 6:
      h0px6 =
       + ln3 * ( 1./35.*pow(mm2,-3)*pow(mm3,4) + 9./70.*pow(mm2,-2)*
         pow(mm3,3) + 7./30.*pow(mm2,-1)*pow(mm3,2) + 13./60.*mm3 + 1./
         12.*mm2 );

      h0px6 +=  + ln2 * (  - 1./15.*pow(mm2,2)*pow(mm3,-1) - 2./15.*
         pow(mm2,3)*pow(mm3,-2) - 37./420.*pow(mm2,4)*pow(mm3,-3) - 3./
         140.*pow(mm2,5)*pow(mm3,-4) );

      h0px6 +=  + ln1 * (  - 1./35.*pow(mm2,-3)*pow(mm3,4) - 9./70.*
         pow(mm2,-2)*pow(mm3,3) - 7./30.*pow(mm2,-1)*pow(mm3,2) - 13./
         60.*mm3 - 1./12.*mm2 + 1./15.*pow(mm2,2)*pow(mm3,-1) + 2./15.*
         pow(mm2,3)*pow(mm3,-2) + 37./420.*pow(mm2,4)*pow(mm3,-3) + 3./
         140.*pow(mm2,5)*pow(mm3,-4) );

      h0px6 +=  + 2./35.*pow(mm2,-2)*pow(mm3,3) + 8./35.*pow(mm2,-1)*
         pow(mm3,2) + 23./84.*mm3 + 3./140.*mm2 - 27./140.*pow(mm2,2)*
         pow(mm3,-1) - 13./84.*pow(mm2,3)*pow(mm3,-2) - 3./70.*pow(
         mm2,4)*pow(mm3,-3);
      return h0px6*pi162/pow(mm1,5);
    case 7:
      h0px7 =
       + ln3 * ( 3./140.*pow(mm2,-4)*pow(mm3,5) + 53./420.*pow(mm2,-3)*
         pow(mm3,4) + 32./105.*pow(mm2,-2)*pow(mm3,3) + 2./5.*pow(
         mm2,-1)*pow(mm3,2) + 17./60.*mm3 + 1./12.*mm2 );

      h0px7 +=  + ln2 * ( 1./12.*mm3 + 17./60.*mm2 + 2./5.*pow(mm2,2)*
         pow(mm3,-1) + 32./105.*pow(mm2,3)*pow(mm3,-2) + 53./420.*pow(
         mm2,4)*pow(mm3,-3) + 3./140.*pow(mm2,5)*pow(mm3,-4) );

      h0px7 +=  + ln1 * (  - 3./140.*pow(mm2,-4)*pow(mm3,5) - 53./420.*
         pow(mm2,-3)*pow(mm3,4) - 32./105.*pow(mm2,-2)*pow(mm3,3) - 2./
         5.*pow(mm2,-1)*pow(mm3,2) - 11./30.*mm3 - 11./30.*mm2 - 2./5.*
         pow(mm2,2)*pow(mm3,-1) - 32./105.*pow(mm2,3)*pow(mm3,-2) - 53./
         420.*pow(mm2,4)*pow(mm3,-3) - 3./140.*pow(mm2,5)*pow(mm3,-4) )
         ;

      h0px7 +=  + 3./70.*pow(mm2,-3)*pow(mm3,4) + 97./420.*pow(mm2,-2)*
         pow(mm3,3) + 209./420.*pow(mm2,-1)*pow(mm3,2) + 67./105.*mm3
          + 67./105.*mm2 + 209./420.*pow(mm2,2)*pow(mm3,-1) + 97./420.*
         pow(mm2,3)*pow(mm3,-2) + 3./70.*pow(mm2,4)*pow(mm3,-3);
      return h0px7*pi162/pow(mm1,5);
    case 8:
      h0px8 =
       + ln3 * (  - 2./105.*pow(mm2,-5)*pow(mm3,6) - 13./105.*pow(
         mm2,-4)*pow(mm3,5) - 12./35.*pow(mm2,-3)*pow(mm3,4) - 11./21.*
         pow(mm2,-2)*pow(mm3,3) - 17./35.*pow(mm2,-1)*pow(mm3,2) - 4./
         15.*mm3 - 1./15.*mm2 );

      h0px8 +=  + ln2 * (  - 1./15.*mm3 - 4./15.*mm2 - 17./35.*pow(
         mm2,2)*pow(mm3,-1) - 11./21.*pow(mm2,3)*pow(mm3,-2) - 12./35.*
         pow(mm2,4)*pow(mm3,-3) - 13./105.*pow(mm2,5)*pow(mm3,-4) - 2./
         105.*pow(mm2,6)*pow(mm3,-5) );

      h0px8 +=  + ln1 * ( 2./105.*pow(mm2,-5)*pow(mm3,6) + 13./105.*
         pow(mm2,-4)*pow(mm3,5) + 12./35.*pow(mm2,-3)*pow(mm3,4) + 11./
         21.*pow(mm2,-2)*pow(mm3,3) + 17./35.*pow(mm2,-1)*pow(mm3,2) + 
         1./3.*mm3 + 1./3.*mm2 + 17./35.*pow(mm2,2)*pow(mm3,-1) + 11./
         21.*pow(mm2,3)*pow(mm3,-2) + 12./35.*pow(mm2,4)*pow(mm3,-3) + 
         13./105.*pow(mm2,5)*pow(mm3,-4) + 2./105.*pow(mm2,6)*pow(
         mm3,-5) );

      h0px8 +=  - 4./105.*pow(mm2,-4)*pow(mm3,5) - 8./35.*pow(mm2,-3)*
         pow(mm3,4) - 181./315.*pow(mm2,-2)*pow(mm3,3) - 7./9.*pow(
         mm2,-1)*pow(mm3,2) - 226./315.*mm3 - 226./315.*mm2 - 7./9.*
         pow(mm2,2)*pow(mm3,-1) - 181./315.*pow(mm2,3)*pow(mm3,-2) - 8./
         35.*pow(mm2,4)*pow(mm3,-3) - 4./105.*pow(mm2,5)*pow(mm3,-4);
      return h0px8*pi162/pow(mm1,7);

    default:
      std::cout << "h0psing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0px1 =
       + ln3 * (  - 1./6.*pow(mm1,-1)*pow(mm3,2) );

      h0px1 +=  + ln2 * ( 1./6.*pow(mm1,-1)*pow(mm3,2) + 1./2.*mm3 + 1./
         2.*mm1 + 1./6.*pow(mm1,2)*pow(mm3,-1) );

      h0px1 +=  + ln1 * (  - 1./6.*pow(mm1,2)*pow(mm3,-1) );

      h0px1 +=  + 7./24.*mm3 + 7./24.*mm1;
      return h0px1*pi162/mm2;
    case 2:
      h0px2 =
       + ln3 * ( 1./30.*pow(mm1,-3)*pow(mm3,4) + 7./60.*pow(mm1,-2)*
         pow(mm3,3) + 1./12.*pow(mm1,-1)*pow(mm3,2) );

      h0px2 +=  + ln2 * (  - 1./30.*pow(mm1,-3)*pow(mm3,4) - 7./60.*
         pow(mm1,-2)*pow(mm3,3) - 1./12.*pow(mm1,-1)*pow(mm3,2) + 1./6.
         *mm3 + 1./3.*mm1 + 13./60.*pow(mm1,2)*pow(mm3,-1) + 1./20.*
         pow(mm1,3)*pow(mm3,-2) );

      h0px2 +=  + ln1 * (  - 1./6.*mm3 - 1./3.*mm1 - 13./60.*pow(mm1,2)
         *pow(mm3,-1) - 1./20.*pow(mm1,3)*pow(mm3,-2) );

      h0px2 +=  + 1./15.*pow(mm1,-2)*pow(mm3,3) + 1./5.*pow(mm1,-1)*
         pow(mm3,2) + 1./10.*mm3 - 2./15.*mm1 - 1./10.*pow(mm1,2)*pow(
         mm3,-1);
      return h0px2*pi162/pow(mm2,3);
    case 3:
      h0px3 =
       + ln3 * ( 1./20.*pow(mm1,-2)*pow(mm3,3) + 1./12.*pow(mm1,-1)*
         pow(mm3,2) );

      h0px3 +=  + ln2 * (  - 1./20.*pow(mm1,-2)*pow(mm3,3) - 1./12.*
         pow(mm1,-1)*pow(mm3,2) - 1./12.*pow(mm1,2)*pow(mm3,-1) - 1./20.
         *pow(mm1,3)*pow(mm3,-2) );

      h0px3 +=  + ln1 * ( 1./12.*pow(mm1,2)*pow(mm3,-1) + 1./20.*pow(
         mm1,3)*pow(mm3,-2) );

      h0px3 +=  + 1./10.*pow(mm1,-1)*pow(mm3,2) + 11./30.*mm3 + 11./30.
         *mm1 + 1./10.*pow(mm1,2)*pow(mm3,-1);
      return h0px3*pi162/pow(mm2,3);
    case 4:
      h0px4 =
       + ln3 * (  - 1./20.*pow(mm1,-2)*pow(mm3,3) - 13./60.*pow(mm1,-1)
         *pow(mm3,2) - 1./3.*mm3 - 1./6.*mm1 );

      h0px4 +=  + ln2 * ( 1./20.*pow(mm1,-2)*pow(mm3,3) + 13./60.*pow(
         mm1,-1)*pow(mm3,2) + 1./3.*mm3 + 1./6.*mm1 - 1./12.*pow(mm1,2)
         *pow(mm3,-1) - 7./60.*pow(mm1,3)*pow(mm3,-2) - 1./30.*pow(
         mm1,4)*pow(mm3,-3) );

      h0px4 +=  + ln1 * ( 1./12.*pow(mm1,2)*pow(mm3,-1) + 7./60.*pow(
         mm1,3)*pow(mm3,-2) + 1./30.*pow(mm1,4)*pow(mm3,-3) );

      h0px4 +=  - 1./10.*pow(mm1,-1)*pow(mm3,2) - 2./15.*mm3 + 1./10.*
         mm1 + 1./5.*pow(mm1,2)*pow(mm3,-1) + 1./15.*pow(mm1,3)*pow(
         mm3,-2);
      return h0px4*pi162/pow(mm2,3);
    case 5:
      h0px5 =
       + ln3 * (  - 3./140.*pow(mm1,-4)*pow(mm3,5) - 37./420.*pow(
         mm1,-3)*pow(mm3,4) - 2./15.*pow(mm1,-2)*pow(mm3,3) - 1./15.*
         pow(mm1,-1)*pow(mm3,2) );

      h0px5 +=  + ln2 * ( 3./140.*pow(mm1,-4)*pow(mm3,5) + 37./420.*
         pow(mm1,-3)*pow(mm3,4) + 2./15.*pow(mm1,-2)*pow(mm3,3) + 1./15.
         *pow(mm1,-1)*pow(mm3,2) - 1./12.*mm3 - 13./60.*mm1 - 7./30.*
         pow(mm1,2)*pow(mm3,-1) - 9./70.*pow(mm1,3)*pow(mm3,-2) - 1./35.
         *pow(mm1,4)*pow(mm3,-3) );

      h0px5 +=  + ln1 * ( 1./12.*mm3 + 13./60.*mm1 + 7./30.*pow(mm1,2)*
         pow(mm3,-1) + 9./70.*pow(mm1,3)*pow(mm3,-2) + 1./35.*pow(
         mm1,4)*pow(mm3,-3) );

      h0px5 +=  - 3./70.*pow(mm1,-3)*pow(mm3,4) - 13./84.*pow(mm1,-2)*
         pow(mm3,3) - 27./140.*pow(mm1,-1)*pow(mm3,2) + 3./140.*mm3 + 
         23./84.*mm1 + 8./35.*pow(mm1,2)*pow(mm3,-1) + 2./35.*pow(
         mm1,3)*pow(mm3,-2);
      return h0px5*pi162/pow(mm2,5);
    case 6:
      h0px6 =
       + ln3 * ( 3./140.*pow(mm1,-4)*pow(mm3,5) + 53./420.*pow(mm1,-3)*
         pow(mm3,4) + 32./105.*pow(mm1,-2)*pow(mm3,3) + 2./5.*pow(
         mm1,-1)*pow(mm3,2) + 17./60.*mm3 + 1./12.*mm1 );

      h0px6 +=  + ln2 * (  - 3./140.*pow(mm1,-4)*pow(mm3,5) - 53./420.*
         pow(mm1,-3)*pow(mm3,4) - 32./105.*pow(mm1,-2)*pow(mm3,3) - 2./
         5.*pow(mm1,-1)*pow(mm3,2) - 11./30.*mm3 - 11./30.*mm1 - 2./5.*
         pow(mm1,2)*pow(mm3,-1) - 32./105.*pow(mm1,3)*pow(mm3,-2) - 53./
         420.*pow(mm1,4)*pow(mm3,-3) - 3./140.*pow(mm1,5)*pow(mm3,-4) )
         ;

      h0px6 +=  + ln1 * ( 1./12.*mm3 + 17./60.*mm1 + 2./5.*pow(mm1,2)*
         pow(mm3,-1) + 32./105.*pow(mm1,3)*pow(mm3,-2) + 53./420.*pow(
         mm1,4)*pow(mm3,-3) + 3./140.*pow(mm1,5)*pow(mm3,-4) );

      h0px6 +=  + 3./70.*pow(mm1,-3)*pow(mm3,4) + 97./420.*pow(mm1,-2)*
         pow(mm3,3) + 209./420.*pow(mm1,-1)*pow(mm3,2) + 67./105.*mm3
          + 67./105.*mm1 + 209./420.*pow(mm1,2)*pow(mm3,-1) + 97./420.*
         pow(mm1,3)*pow(mm3,-2) + 3./70.*pow(mm1,4)*pow(mm3,-3);
      return h0px6*pi162/pow(mm2,5);
    case 7:
      h0px7 =
       + ln3 * ( 1./35.*pow(mm1,-3)*pow(mm3,4) + 9./70.*pow(mm1,-2)*
         pow(mm3,3) + 7./30.*pow(mm1,-1)*pow(mm3,2) + 13./60.*mm3 + 1./
         12.*mm1 );

      h0px7 +=  + ln2 * (  - 1./35.*pow(mm1,-3)*pow(mm3,4) - 9./70.*
         pow(mm1,-2)*pow(mm3,3) - 7./30.*pow(mm1,-1)*pow(mm3,2) - 13./
         60.*mm3 - 1./12.*mm1 + 1./15.*pow(mm1,2)*pow(mm3,-1) + 2./15.*
         pow(mm1,3)*pow(mm3,-2) + 37./420.*pow(mm1,4)*pow(mm3,-3) + 3./
         140.*pow(mm1,5)*pow(mm3,-4) );

      h0px7 +=  + ln1 * (  - 1./15.*pow(mm1,2)*pow(mm3,-1) - 2./15.*
         pow(mm1,3)*pow(mm3,-2) - 37./420.*pow(mm1,4)*pow(mm3,-3) - 3./
         140.*pow(mm1,5)*pow(mm3,-4) );

      h0px7 +=  + 2./35.*pow(mm1,-2)*pow(mm3,3) + 8./35.*pow(mm1,-1)*
         pow(mm3,2) + 23./84.*mm3 + 3./140.*mm1 - 27./140.*pow(mm1,2)*
         pow(mm3,-1) - 13./84.*pow(mm1,3)*pow(mm3,-2) - 3./70.*pow(
         mm1,4)*pow(mm3,-3);
      return h0px7*pi162/pow(mm2,5);
    case 8:
      h0px8 =
       + ln3 * (  - 2./105.*pow(mm1,-5)*pow(mm3,6) - 13./105.*pow(
         mm1,-4)*pow(mm3,5) - 12./35.*pow(mm1,-3)*pow(mm3,4) - 11./21.*
         pow(mm1,-2)*pow(mm3,3) - 17./35.*pow(mm1,-1)*pow(mm3,2) - 4./
         15.*mm3 - 1./15.*mm1 );

      h0px8 +=  + ln2 * ( 2./105.*pow(mm1,-5)*pow(mm3,6) + 13./105.*
         pow(mm1,-4)*pow(mm3,5) + 12./35.*pow(mm1,-3)*pow(mm3,4) + 11./
         21.*pow(mm1,-2)*pow(mm3,3) + 17./35.*pow(mm1,-1)*pow(mm3,2) + 
         1./3.*mm3 + 1./3.*mm1 + 17./35.*pow(mm1,2)*pow(mm3,-1) + 11./
         21.*pow(mm1,3)*pow(mm3,-2) + 12./35.*pow(mm1,4)*pow(mm3,-3) + 
         13./105.*pow(mm1,5)*pow(mm3,-4) + 2./105.*pow(mm1,6)*pow(
         mm3,-5) );

      h0px8 +=  + ln1 * (  - 1./15.*mm3 - 4./15.*mm1 - 17./35.*pow(
         mm1,2)*pow(mm3,-1) - 11./21.*pow(mm1,3)*pow(mm3,-2) - 12./35.*
         pow(mm1,4)*pow(mm3,-3) - 13./105.*pow(mm1,5)*pow(mm3,-4) - 2./
         105.*pow(mm1,6)*pow(mm3,-5) );

      h0px8 +=  - 4./105.*pow(mm1,-4)*pow(mm3,5) - 8./35.*pow(mm1,-3)*
         pow(mm3,4) - 181./315.*pow(mm1,-2)*pow(mm3,3) - 7./9.*pow(
         mm1,-1)*pow(mm3,2) - 226./315.*mm3 - 226./315.*mm1 - 7./9.*
         pow(mm1,2)*pow(mm3,-1) - 181./315.*pow(mm1,3)*pow(mm3,-2) - 8./
         35.*pow(mm1,4)*pow(mm3,-3) - 4./105.*pow(mm1,5)*pow(mm3,-4);
      return h0px8*pi162/pow(mm2,7);
    default:
      std::cout << "h0psing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h0px1 =
       + ln3 * ( 1./6.*pow(mm1,-1)*pow(mm2,2) + 1./2.*mm2 + 1./2.*mm1
          + 1./6.*pow(mm1,2)*pow(mm2,-1) );

      h0px1 +=  + ln2 * (  - 1./6.*pow(mm1,-1)*pow(mm2,2) );

      h0px1 +=  + ln1 * (  - 1./6.*pow(mm1,2)*pow(mm2,-1) );

      h0px1 +=  + 7./24.*mm2 + 7./24.*mm1;
      return h0px1*pi162/mm3;
    case 2:
      h0px2 =
       + ln3 * (  - 1./30.*pow(mm1,-3)*pow(mm2,4) - 7./60.*pow(mm1,-2)*
         pow(mm2,3) - 1./12.*pow(mm1,-1)*pow(mm2,2) + 1./6.*mm2 + 1./3.
         *mm1 + 13./60.*pow(mm1,2)*pow(mm2,-1) + 1./20.*pow(mm1,3)*pow(
         mm2,-2) );

      h0px2 +=  + ln2 * ( 1./30.*pow(mm1,-3)*pow(mm2,4) + 7./60.*pow(
         mm1,-2)*pow(mm2,3) + 1./12.*pow(mm1,-1)*pow(mm2,2) );

      h0px2 +=  + ln1 * (  - 1./6.*mm2 - 1./3.*mm1 - 13./60.*pow(mm1,2)
         *pow(mm2,-1) - 1./20.*pow(mm1,3)*pow(mm2,-2) );

      h0px2 +=  + 1./15.*pow(mm1,-2)*pow(mm2,3) + 1./5.*pow(mm1,-1)*
         pow(mm2,2) + 1./10.*mm2 - 2./15.*mm1 - 1./10.*pow(mm1,2)*pow(
         mm2,-1);
      return h0px2*pi162/pow(mm3,3);
    case 3:
      h0px3 =
       + ln3 * ( 1./20.*pow(mm1,-2)*pow(mm2,3) + 13./60.*pow(mm1,-1)*
         pow(mm2,2) + 1./3.*mm2 + 1./6.*mm1 - 1./12.*pow(mm1,2)*pow(
         mm2,-1) - 7./60.*pow(mm1,3)*pow(mm2,-2) - 1./30.*pow(mm1,4)*
         pow(mm2,-3) );

      h0px3 +=  + ln2 * (  - 1./20.*pow(mm1,-2)*pow(mm2,3) - 13./60.*
         pow(mm1,-1)*pow(mm2,2) - 1./3.*mm2 - 1./6.*mm1 );

      h0px3 +=  + ln1 * ( 1./12.*pow(mm1,2)*pow(mm2,-1) + 7./60.*pow(
         mm1,3)*pow(mm2,-2) + 1./30.*pow(mm1,4)*pow(mm2,-3) );

      h0px3 +=  - 1./10.*pow(mm1,-1)*pow(mm2,2) - 2./15.*mm2 + 1./10.*
         mm1 + 1./5.*pow(mm1,2)*pow(mm2,-1) + 1./15.*pow(mm1,3)*pow(
         mm2,-2);
      return h0px3*pi162/pow(mm3,3);
    case 4:
      h0px4 =
       + ln3 * (  - 1./20.*pow(mm1,-2)*pow(mm2,3) - 1./12.*pow(mm1,-1)*
         pow(mm2,2) - 1./12.*pow(mm1,2)*pow(mm2,-1) - 1./20.*pow(mm1,3)
         *pow(mm2,-2) );

      h0px4 +=  + ln2 * ( 1./20.*pow(mm1,-2)*pow(mm2,3) + 1./12.*pow(
         mm1,-1)*pow(mm2,2) );

      h0px4 +=  + ln1 * ( 1./12.*pow(mm1,2)*pow(mm2,-1) + 1./20.*pow(
         mm1,3)*pow(mm2,-2) );

      h0px4 +=  + 1./10.*pow(mm1,-1)*pow(mm2,2) + 11./30.*mm2 + 11./30.
         *mm1 + 1./10.*pow(mm1,2)*pow(mm2,-1);
     return h0px4*pi162/pow(mm3,3);
    case 5:
      h0px5 =
       + ln3 * (  - 3./140.*pow(mm1,-4)*pow(mm2,5) - 53./420.*pow(
         mm1,-3)*pow(mm2,4) - 32./105.*pow(mm1,-2)*pow(mm2,3) - 2./5.*
         pow(mm1,-1)*pow(mm2,2) - 11./30.*mm2 - 11./30.*mm1 - 2./5.*
         pow(mm1,2)*pow(mm2,-1) - 32./105.*pow(mm1,3)*pow(mm2,-2) - 53./
         420.*pow(mm1,4)*pow(mm2,-3) - 3./140.*pow(mm1,5)*pow(mm2,-4) )
         ;

      h0px5 +=  + ln2 * ( 3./140.*pow(mm1,-4)*pow(mm2,5) + 53./420.*
         pow(mm1,-3)*pow(mm2,4) + 32./105.*pow(mm1,-2)*pow(mm2,3) + 2./
         5.*pow(mm1,-1)*pow(mm2,2) + 17./60.*mm2 + 1./12.*mm1 );

      h0px5 +=  + ln1 * ( 1./12.*mm2 + 17./60.*mm1 + 2./5.*pow(mm1,2)*
         pow(mm2,-1) + 32./105.*pow(mm1,3)*pow(mm2,-2) + 53./420.*pow(
         mm1,4)*pow(mm2,-3) + 3./140.*pow(mm1,5)*pow(mm2,-4) );

      h0px5 +=  + 3./70.*pow(mm1,-3)*pow(mm2,4) + 97./420.*pow(mm1,-2)*
         pow(mm2,3) + 209./420.*pow(mm1,-1)*pow(mm2,2) + 67./105.*mm2
          + 67./105.*mm1 + 209./420.*pow(mm1,2)*pow(mm2,-1) + 97./420.*
         pow(mm1,3)*pow(mm2,-2) + 3./70.*pow(mm1,4)*pow(mm2,-3);
     return h0px5*pi162/pow(mm3,5);
    case 6:
      h0px6 =
       + ln3 * ( 3./140.*pow(mm1,-4)*pow(mm2,5) + 37./420.*pow(mm1,-3)*
         pow(mm2,4) + 2./15.*pow(mm1,-2)*pow(mm2,3) + 1./15.*pow(
         mm1,-1)*pow(mm2,2) - 1./12.*mm2 - 13./60.*mm1 - 7./30.*pow(
         mm1,2)*pow(mm2,-1) - 9./70.*pow(mm1,3)*pow(mm2,-2) - 1./35.*
         pow(mm1,4)*pow(mm2,-3) );

      h0px6 +=  + ln2 * (  - 3./140.*pow(mm1,-4)*pow(mm2,5) - 37./420.*
         pow(mm1,-3)*pow(mm2,4) - 2./15.*pow(mm1,-2)*pow(mm2,3) - 1./15.
         *pow(mm1,-1)*pow(mm2,2) );

      h0px6 +=  + ln1 * ( 1./12.*mm2 + 13./60.*mm1 + 7./30.*pow(mm1,2)*
         pow(mm2,-1) + 9./70.*pow(mm1,3)*pow(mm2,-2) + 1./35.*pow(
         mm1,4)*pow(mm2,-3) );

      h0px6 +=  - 3./70.*pow(mm1,-3)*pow(mm2,4) - 13./84.*pow(mm1,-2)*
         pow(mm2,3) - 27./140.*pow(mm1,-1)*pow(mm2,2) + 3./140.*mm2 + 
         23./84.*mm1 + 8./35.*pow(mm1,2)*pow(mm2,-1) + 2./35.*pow(
         mm1,3)*pow(mm2,-2);
     return h0px6*pi162/pow(mm3,5);
    case 7:
      h0px7 =
       + ln3 * (  - 1./35.*pow(mm1,-3)*pow(mm2,4) - 9./70.*pow(mm1,-2)*
         pow(mm2,3) - 7./30.*pow(mm1,-1)*pow(mm2,2) - 13./60.*mm2 - 1./
         12.*mm1 + 1./15.*pow(mm1,2)*pow(mm2,-1) + 2./15.*pow(mm1,3)*
         pow(mm2,-2) + 37./420.*pow(mm1,4)*pow(mm2,-3) + 3./140.*pow(
         mm1,5)*pow(mm2,-4) );

      h0px7 +=  + ln2 * ( 1./35.*pow(mm1,-3)*pow(mm2,4) + 9./70.*pow(
         mm1,-2)*pow(mm2,3) + 7./30.*pow(mm1,-1)*pow(mm2,2) + 13./60.*
         mm2 + 1./12.*mm1 );

      h0px7 +=  + ln1 * (  - 1./15.*pow(mm1,2)*pow(mm2,-1) - 2./15.*
         pow(mm1,3)*pow(mm2,-2) - 37./420.*pow(mm1,4)*pow(mm2,-3) - 3./
         140.*pow(mm1,5)*pow(mm2,-4) );

      h0px7 +=  + 2./35.*pow(mm1,-2)*pow(mm2,3) + 8./35.*pow(mm1,-1)*
         pow(mm2,2) + 23./84.*mm2 + 3./140.*mm1 - 27./140.*pow(mm1,2)*
         pow(mm2,-1) - 13./84.*pow(mm1,3)*pow(mm2,-2) - 3./70.*pow(
         mm1,4)*pow(mm2,-3);
     return h0px7*pi162/pow(mm3,5);
    case 8:
      h0px8 =
       + ln3 * ( 2./105.*pow(mm1,-5)*pow(mm2,6) + 13./105.*pow(mm1,-4)*
         pow(mm2,5) + 12./35.*pow(mm1,-3)*pow(mm2,4) + 11./21.*pow(
         mm1,-2)*pow(mm2,3) + 17./35.*pow(mm1,-1)*pow(mm2,2) + 1./3.*
         mm2 + 1./3.*mm1 + 17./35.*pow(mm1,2)*pow(mm2,-1) + 11./21.*
         pow(mm1,3)*pow(mm2,-2) + 12./35.*pow(mm1,4)*pow(mm2,-3) + 13./
         105.*pow(mm1,5)*pow(mm2,-4) + 2./105.*pow(mm1,6)*pow(mm2,-5) )
         ;

      h0px8 +=  + ln2 * (  - 2./105.*pow(mm1,-5)*pow(mm2,6) - 13./105.*
         pow(mm1,-4)*pow(mm2,5) - 12./35.*pow(mm1,-3)*pow(mm2,4) - 11./
         21.*pow(mm1,-2)*pow(mm2,3) - 17./35.*pow(mm1,-1)*pow(mm2,2) - 
         4./15.*mm2 - 1./15.*mm1 );

      h0px8 +=  + ln1 * (  - 1./15.*mm2 - 4./15.*mm1 - 17./35.*pow(
         mm1,2)*pow(mm2,-1) - 11./21.*pow(mm1,3)*pow(mm2,-2) - 12./35.*
         pow(mm1,4)*pow(mm2,-3) - 13./105.*pow(mm1,5)*pow(mm2,-4) - 2./
         105.*pow(mm1,6)*pow(mm2,-5) );

      h0px8 +=  - 4./105.*pow(mm1,-4)*pow(mm2,5) - 8./35.*pow(mm1,-3)*
         pow(mm2,4) - 181./315.*pow(mm1,-2)*pow(mm2,3) - 7./9.*pow(
         mm1,-1)*pow(mm2,2) - 226./315.*mm2 - 226./315.*mm1 - 7./9.*
         pow(mm1,2)*pow(mm2,-1) - 181./315.*pow(mm1,3)*pow(mm2,-2) - 8./
         35.*pow(mm1,4)*pow(mm2,-3) - 4./105.*pow(mm1,5)*pow(mm2,-4);

      return h0px8*pi162/pow(mm3,7);
    default:
      std::cout << "h0psing, wrong iprop = "<<iprop<<'\n';
      return 0.;
    }
  }
  std::cout << "funny masses in h0psing\n";
  return 0.;
  }

double h10sing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h10x1,h10x2,h10x3,h10x4,h10x5,h10x6,h10x7,h10x8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10x1 =
       + ln3 * ( 1./4.*pow(mm3,3) + 3./4.*mm2*pow(mm3,2) + 1./2.*pow(
         mm2,2)*mm3 );

      h10x1 +=  + pow(ln3,2) * ( 1./4.*pow(mm3,3) + 1./4.*mm2*pow(
         mm3,2) );

      h10x1 +=  + ln2 * ( 1./2.*mm2*pow(mm3,2) + 3./4.*pow(mm2,2)*mm3
          + 1./4.*pow(mm2,3) );

      h10x1 +=  + ln2*ln3 * ( 1./4.*pow(mm3,3) + 1./4.*mm2*pow(mm3,2)
          + 1./4.*pow(mm2,2)*mm3 + 1./4.*pow(mm2,3) );

      h10x1 +=  + pow(ln2,2) * ( 1./4.*pow(mm2,2)*mm3 + 1./4.*pow(
         mm2,3) );

      h10x1 +=  + ln1 * (  - pow(mm3,3) - 3*mm2*pow(mm3,2) - 3*pow(
         mm2,2)*mm3 - pow(mm2,3) );

      h10x1 +=  + ln1*ln3 * ( 1./4.*pow(mm3,3) + 1./4.*mm2*pow(mm3,2)
          - 1./4.*pow(mm2,2)*mm3 - 1./4.*pow(mm2,3) );

      h10x1 +=  + ln1*ln2 * (  - 1./4.*pow(mm3,3) - 1./4.*mm2*pow(
         mm3,2) + 1./4.*pow(mm2,2)*mm3 + 1./4.*pow(mm2,3) );

      h10x1 +=  + 15./16.*pow(mm3,3) + 1./24.*pow(mm3,3)*pi2 + 27./16.*
         mm2*pow(mm3,2) + 1./24.*mm2*pow(mm3,2)*pi2 + 27./16.*pow(
         mm2,2)*mm3 + 1./24.*pow(mm2,2)*mm3*pi2 + 15./16.*pow(mm2,3) + 
         1./24.*pow(mm2,3)*pi2;
      return h10x1*pi162/mm1;
    case 2:
      h10x2 =
       + ln3 * ( 1./6.*pow(mm2,-1)*pow(mm3,5) + 1./2.*pow(mm3,4) + 1./2.
         *mm2*pow(mm3,3) + 1./6.*pow(mm2,2)*pow(mm3,2) );

      h10x2 +=  + ln2 * ( 1./6.*pow(mm2,2)*pow(mm3,2) + 1./2.*pow(
         mm2,3)*mm3 + 1./2.*pow(mm2,4) + 1./6.*pow(mm2,5)*pow(mm3,-1) )
         ;

      h10x2 +=  + ln1 * (  - 1./6.*pow(mm2,-1)*pow(mm3,5) - pow(mm3,4)
          - 5./2.*mm2*pow(mm3,3) - 10./3.*pow(mm2,2)*pow(mm3,2) - 5./2.
         *pow(mm2,3)*mm3 - pow(mm2,4) - 1./6.*pow(mm2,5)*pow(mm3,-1) );

      h10x2 +=  - 7./24.*pow(mm3,4) - 7./6.*mm2*pow(mm3,3) - 7./4.*pow(
         mm2,2)*pow(mm3,2) - 7./6.*pow(mm2,3)*mm3 - 7./24.*pow(mm2,4);
      return h10x2*pi162/pow(mm1,4);
    case 3:
      h10x3 =
       + ln3 * ( 1./12.*pow(mm2,-2)*pow(mm3,6) + 5./6.*pow(mm2,-1)*pow(
         mm3,5) + 5./2.*pow(mm3,4) + 10./3.*mm2*pow(mm3,3) + 25./12.*
         pow(mm2,2)*pow(mm3,2) + 1./2.*pow(mm2,3)*mm3 );

      h10x3 +=  + ln2 * ( pow(mm3,4) + 23./6.*mm2*pow(mm3,3) + 16./3.*
         pow(mm2,2)*pow(mm3,2) + 3*pow(mm2,3)*mm3 + 1./3.*pow(mm2,4) - 
         1./6.*pow(mm2,5)*pow(mm3,-1) );

      h10x3 +=  + ln2*ln3 * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./2.
         *pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) );

      h10x3 +=  + pow(ln2,2) * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./
         2.*pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) )
         ;

      h10x3 +=  + ln1 * (  - 1./12.*pow(mm2,-2)*pow(mm3,6) - 5./6.*pow(
         mm2,-1)*pow(mm3,5) - 11./4.*pow(mm3,4) - 25./6.*mm2*pow(mm3,3)
          - 35./12.*pow(mm2,2)*pow(mm3,2) - 1./2.*pow(mm2,3)*mm3 + 5./
         12.*pow(mm2,4) + 1./6.*pow(mm2,5)*pow(mm3,-1) );

      h10x3 +=  + ln1*ln3 * (  - 1./4.*pow(mm3,4) - mm2*pow(mm3,3) - 3./
         2.*pow(mm2,2)*pow(mm3,2) - pow(mm2,3)*mm3 - 1./4.*pow(mm2,4) )
         ;

      h10x3 +=  + ln1*ln2 * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./2.
         *pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) );

      h10x3 +=  + 1./6.*pow(mm2,-1)*pow(mm3,5) + 55./48.*pow(mm3,4) + 1.
         /24.*pow(mm3,4)*pi2 + 35./12.*mm2*pow(mm3,3) + 1./6.*mm2*pow(
         mm3,3)*pi2 + 85./24.*pow(mm2,2)*pow(mm3,2) + 1./4.*pow(mm2,2)*
         pow(mm3,2)*pi2 + 25./12.*pow(mm2,3)*mm3 + 1./6.*pow(mm2,3)*mm3
         *pi2 + 23./48.*pow(mm2,4) + 1./24.*pow(mm2,4)*pi2;
      return h10x3*pi162/pow(mm1,4);
    case 4:
      h10x4 =
       + ln3 * (  - 1./6.*pow(mm2,-1)*pow(mm3,5) + 1./3.*pow(mm3,4) + 3
         *mm2*pow(mm3,3) + 16./3.*pow(mm2,2)*pow(mm3,2) + 23./6.*pow(
         mm2,3)*mm3 + pow(mm2,4) );

      h10x4 +=  + pow(ln3,2) * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./
         2.*pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) )
         ;

      h10x4 +=  + ln2 * ( 1./2.*mm2*pow(mm3,3) + 25./12.*pow(mm2,2)*
         pow(mm3,2) + 10./3.*pow(mm2,3)*mm3 + 5./2.*pow(mm2,4) + 5./6.*
         pow(mm2,5)*pow(mm3,-1) + 1./12.*pow(mm2,6)*pow(mm3,-2) );

      h10x4 +=  + ln2*ln3 * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./2.
         *pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) );

      h10x4 +=  + ln1 * ( 1./6.*pow(mm2,-1)*pow(mm3,5) + 5./12.*pow(
         mm3,4) - 1./2.*mm2*pow(mm3,3) - 35./12.*pow(mm2,2)*pow(mm3,2)
          - 25./6.*pow(mm2,3)*mm3 - 11./4.*pow(mm2,4) - 5./6.*pow(
         mm2,5)*pow(mm3,-1) - 1./12.*pow(mm2,6)*pow(mm3,-2) );

      h10x4 +=  + ln1*ln3 * ( 1./4.*pow(mm3,4) + mm2*pow(mm3,3) + 3./2.
         *pow(mm2,2)*pow(mm3,2) + pow(mm2,3)*mm3 + 1./4.*pow(mm2,4) );

      h10x4 +=  + ln1*ln2 * (  - 1./4.*pow(mm3,4) - mm2*pow(mm3,3) - 3./
         2.*pow(mm2,2)*pow(mm3,2) - pow(mm2,3)*mm3 - 1./4.*pow(mm2,4) )
         ;

      h10x4 +=  + 23./48.*pow(mm3,4) + 1./24.*pow(mm3,4)*pi2 + 25./12.*
         mm2*pow(mm3,3) + 1./6.*mm2*pow(mm3,3)*pi2 + 85./24.*pow(mm2,2)
         *pow(mm3,2) + 1./4.*pow(mm2,2)*pow(mm3,2)*pi2 + 35./12.*pow(
         mm2,3)*mm3 + 1./6.*pow(mm2,3)*mm3*pi2 + 55./48.*pow(mm2,4) + 1.
         /24.*pow(mm2,4)*pi2 + 1./6.*pow(mm2,5)*pow(mm3,-1);
      return h10x4*pi162/pow(mm1,4);
    case 5:
      h10x5 =
       + ln3 * (  - 1./30.*pow(mm2,-3)*pow(mm3,7) - 13./60.*pow(mm2,-2)
         *pow(mm3,6) - 8./15.*pow(mm2,-1)*pow(mm3,5) - 19./30.*pow(
         mm3,4) - 11./30.*mm2*pow(mm3,3) - 1./12.*pow(mm2,2)*pow(mm3,2)
          );

      h10x5 +=  + ln2 * ( 1./6.*pow(mm3,4) + 5./6.*mm2*pow(mm3,3) + 103.
         /60.*pow(mm2,2)*pow(mm3,2) + 28./15.*pow(mm2,3)*mm3 + 17./15.*
         pow(mm2,4) + 11./30.*pow(mm2,5)*pow(mm3,-1) + 1./20.*pow(
         mm2,6)*pow(mm3,-2) );

      h10x5 +=  + ln1 * ( 1./30.*pow(mm2,-3)*pow(mm3,7) + 13./60.*pow(
         mm2,-2)*pow(mm3,6) + 8./15.*pow(mm2,-1)*pow(mm3,5) + 7./15.*
         pow(mm3,4) - 7./15.*mm2*pow(mm3,3) - 49./30.*pow(mm2,2)*pow(
         mm3,2) - 28./15.*pow(mm2,3)*mm3 - 17./15.*pow(mm2,4) - 11./30.
         *pow(mm2,5)*pow(mm3,-1) - 1./20.*pow(mm2,6)*pow(mm3,-2) );

      h10x5 +=  - 1./15.*pow(mm2,-2)*pow(mm3,6) - 2./5.*pow(mm2,-1)*
         pow(mm3,5) - 9./10.*pow(mm3,4) - 5./6.*mm2*pow(mm3,3) + 3./5.*
         pow(mm2,3)*mm3 + 13./30.*pow(mm2,4) + 1./10.*pow(mm2,5)*pow(
         mm3,-1);
      return h10x5*pi162/pow(mm1,6);
    case 6:
      h10x6 =
       + ln3 * ( 1./20.*pow(mm2,-2)*pow(mm3,6) + 11./30.*pow(mm2,-1)*
         pow(mm3,5) + 17./15.*pow(mm3,4) + 28./15.*mm2*pow(mm3,3) + 103.
         /60.*pow(mm2,2)*pow(mm3,2) + 5./6.*pow(mm2,3)*mm3 + 1./6.*pow(
         mm2,4) );

      h10x6 +=  + ln2 * (  - 1./12.*pow(mm2,2)*pow(mm3,2) - 11./30.*
         pow(mm2,3)*mm3 - 19./30.*pow(mm2,4) - 8./15.*pow(mm2,5)*pow(
         mm3,-1) - 13./60.*pow(mm2,6)*pow(mm3,-2) - 1./30.*pow(mm2,7)*
         pow(mm3,-3) );

      h10x6 +=  + ln1 * (  - 1./20.*pow(mm2,-2)*pow(mm3,6) - 11./30.*
         pow(mm2,-1)*pow(mm3,5) - 17./15.*pow(mm3,4) - 28./15.*mm2*pow(
         mm3,3) - 49./30.*pow(mm2,2)*pow(mm3,2) - 7./15.*pow(mm2,3)*mm3
          + 7./15.*pow(mm2,4) + 8./15.*pow(mm2,5)*pow(mm3,-1) + 13./60.
         *pow(mm2,6)*pow(mm3,-2) + 1./30.*pow(mm2,7)*pow(mm3,-3) );

      h10x6 +=  + 1./10.*pow(mm2,-1)*pow(mm3,5) + 13./30.*pow(mm3,4) + 
         3./5.*mm2*pow(mm3,3) - 5./6.*pow(mm2,3)*mm3 - 9./10.*pow(
         mm2,4) - 2./5.*pow(mm2,5)*pow(mm3,-1) - 1./15.*pow(mm2,6)*pow(
         mm3,-2);
      return h10x6*pi162/pow(mm1,6);
    case 7:
      h10x7 =
       + ln3 * ( 1./30.*pow(mm2,-3)*pow(mm3,7) + 1./3.*pow(mm2,-2)*pow(
         mm3,6) + 3./2.*pow(mm2,-1)*pow(mm3,5) + 11./3.*pow(mm3,4) + 31.
         /6.*mm2*pow(mm3,3) + 21./5.*pow(mm2,2)*pow(mm3,2) + 11./6.*
         pow(mm2,3)*mm3 + 1./3.*pow(mm2,4) );

      h10x7 +=  + ln2 * ( 1./3.*pow(mm3,4) + 11./6.*mm2*pow(mm3,3) + 21.
         /5.*pow(mm2,2)*pow(mm3,2) + 31./6.*pow(mm2,3)*mm3 + 11./3.*
         pow(mm2,4) + 3./2.*pow(mm2,5)*pow(mm3,-1) + 1./3.*pow(mm2,6)*
         pow(mm3,-2) + 1./30.*pow(mm2,7)*pow(mm3,-3) );

      h10x7 +=  + ln1 * (  - 1./30.*pow(mm2,-3)*pow(mm3,7) - 1./3.*pow(
         mm2,-2)*pow(mm3,6) - 3./2.*pow(mm2,-1)*pow(mm3,5) - 4*pow(
         mm3,4) - 7*mm2*pow(mm3,3) - 42./5.*pow(mm2,2)*pow(mm3,2) - 7*
         pow(mm2,3)*mm3 - 4*pow(mm2,4) - 3./2.*pow(mm2,5)*pow(mm3,-1)
          - 1./3.*pow(mm2,6)*pow(mm3,-2) - 1./30.*pow(mm2,7)*pow(
         mm3,-3) );

      h10x7 +=  + 1./15.*pow(mm2,-2)*pow(mm3,6) + 19./30.*pow(mm2,-1)*
         pow(mm3,5) + 37./15.*pow(mm3,4) + 157./30.*mm2*pow(mm3,3) + 20.
         /3.*pow(mm2,2)*pow(mm3,2) + 157./30.*pow(mm2,3)*mm3 + 37./15.*
         pow(mm2,4) + 19./30.*pow(mm2,5)*pow(mm3,-1) + 1./15.*pow(
         mm2,6)*pow(mm3,-2);
      return h10x7*pi162/pow(mm1,6);
    case 8:
      h10x8 =
       + ln3 * (  - 3./140.*pow(mm2,-4)*pow(mm3,8) - 4./21.*pow(mm2,-3)
         *pow(mm3,7) - 157./210.*pow(mm2,-2)*pow(mm3,6) - 12./7.*pow(
         mm2,-1)*pow(mm3,5) - 53./21.*pow(mm3,4) - 256./105.*mm2*pow(
         mm3,3) - 3./2.*pow(mm2,2)*pow(mm3,2) - 8./15.*pow(mm2,3)*mm3
          - 1./12.*pow(mm2,4) );

      h10x8 +=  + ln2 * (  - 1./12.*pow(mm3,4) - 8./15.*mm2*pow(mm3,3)
          - 3./2.*pow(mm2,2)*pow(mm3,2) - 256./105.*pow(mm2,3)*mm3 - 53.
         /21.*pow(mm2,4) - 12./7.*pow(mm2,5)*pow(mm3,-1) - 157./210.*
         pow(mm2,6)*pow(mm3,-2) - 4./21.*pow(mm2,7)*pow(mm3,-3) - 3./
         140.*pow(mm2,8)*pow(mm3,-4) );

      h10x8 +=  + ln1 * ( 3./140.*pow(mm2,-4)*pow(mm3,8) + 4./21.*pow(
         mm2,-3)*pow(mm3,7) + 157./210.*pow(mm2,-2)*pow(mm3,6) + 12./7.
         *pow(mm2,-1)*pow(mm3,5) + 73./28.*pow(mm3,4) + 104./35.*mm2*
         pow(mm3,3) + 3*pow(mm2,2)*pow(mm3,2) + 104./35.*pow(mm2,3)*mm3
          + 73./28.*pow(mm2,4) + 12./7.*pow(mm2,5)*pow(mm3,-1) + 157./
         210.*pow(mm2,6)*pow(mm3,-2) + 4./21.*pow(mm2,7)*pow(mm3,-3) + 
         3./140.*pow(mm2,8)*pow(mm3,-4) );

      h10x8 +=  - 3./70.*pow(mm2,-3)*pow(mm3,7) - 151./420.*pow(mm2,-2)
         *pow(mm3,6) - 277./210.*pow(mm2,-1)*pow(mm3,5) - 43./15.*pow(
         mm3,4) - 449./105.*mm2*pow(mm3,3) - 1013./210.*pow(mm2,2)*pow(
         mm3,2) - 449./105.*pow(mm2,3)*mm3 - 43./15.*pow(mm2,4) - 277./
         210.*pow(mm2,5)*pow(mm3,-1) - 151./420.*pow(mm2,6)*pow(mm3,-2)
          - 3./70.*pow(mm2,7)*pow(mm3,-3);
      return h10x8*pi162/pow(mm1,8);
    default:
      std::cout << "h10sing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10x1 =
       + ln3 * (  - 1./4.*pow(mm3,3) - 3./4.*mm1*pow(mm3,2) - 1./2.*
         pow(mm1,2)*mm3 );

      h10x1 +=  + pow(ln3,2) * ( 1./4.*pow(mm3,3) + 1./4.*mm1*pow(
         mm3,2) );

      h10x1 +=  + ln2 * (  - 1./4.*pow(mm3,3) - 1./4.*mm1*pow(mm3,2) + 
         1./4.*pow(mm1,2)*mm3 + 1./4.*pow(mm1,3) );

      h10x1 +=  + ln2*ln3 * ( 1./2.*pow(mm3,3) + mm1*pow(mm3,2) + 3./4.
         *pow(mm1,2)*mm3 + 1./4.*pow(mm1,3) );

      h10x1 +=  + pow(ln2,2) * ( 1./4.*pow(mm3,3) + 3./4.*mm1*pow(
         mm3,2) + 3./4.*pow(mm1,2)*mm3 + 1./4.*pow(mm1,3) );

      h10x1 +=  + ln1 * (  - pow(mm1,2)*mm3 - pow(mm1,3) );

      h10x1 +=  + ln1*ln3 * (  - 1./2.*mm1*pow(mm3,2) - 3./4.*pow(
         mm1,2)*mm3 - 1./4.*pow(mm1,3) );

      h10x1 +=  + ln1*ln2 * ( 1./2.*mm1*pow(mm3,2) + 3./4.*pow(mm1,2)*
         mm3 + 1./4.*pow(mm1,3) );

      h10x1 +=  + 9./8.*pow(mm3,3) + 1./12.*pow(mm3,3)*pi2 + 9./4.*mm1*
         pow(mm3,2) + 1./6.*mm1*pow(mm3,2)*pi2 + 33./16.*pow(mm1,2)*mm3
          + 1./8.*pow(mm1,2)*mm3*pi2 + 15./16.*pow(mm1,3) + 1./24.*pow(
         mm1,3)*pi2;
      return h10x1*pi162/mm2;
    case 2:
      h10x2 =
       + ln3 * ( 1./6.*pow(mm1,-1)*pow(mm3,5) + 1./2.*pow(mm3,4) + 1./2.
         *mm1*pow(mm3,3) + 1./6.*pow(mm1,2)*pow(mm3,2) );

      h10x2 +=  + ln2 * (  - 1./6.*pow(mm1,-1)*pow(mm3,5) - pow(mm3,4)
          - 5./2.*mm1*pow(mm3,3) - 10./3.*pow(mm1,2)*pow(mm3,2) - 5./2.
         *pow(mm1,3)*mm3 - pow(mm1,4) - 1./6.*pow(mm1,5)*pow(mm3,-1) );

      h10x2 +=  + ln1 * ( 1./6.*pow(mm1,2)*pow(mm3,2) + 1./2.*pow(
         mm1,3)*mm3 + 1./2.*pow(mm1,4) + 1./6.*pow(mm1,5)*pow(mm3,-1) )
         ;

      h10x2 +=  - 7./24.*pow(mm3,4) - 7./6.*mm1*pow(mm3,3) - 7./4.*pow(
         mm1,2)*pow(mm3,2) - 7./6.*pow(mm1,3)*mm3 - 7./24.*pow(mm1,4);
      return h10x2*pi162/pow(mm2,4);
    case 3:
      h10x3 =
       + ln3 * (  - 5./12.*pow(mm3,4) - 4./3.*mm1*pow(mm3,3) - 17./12.*
         pow(mm1,2)*pow(mm3,2) - 1./2.*pow(mm1,3)*mm3 );

      h10x3 +=  + ln2 * ( 7./6.*pow(mm3,4) + 29./6.*mm1*pow(mm3,3) + 23.
         /3.*pow(mm1,2)*pow(mm3,2) + 17./3.*pow(mm1,3)*mm3 + 11./6.*
         pow(mm1,4) + 1./6.*pow(mm1,5)*pow(mm3,-1) );

      h10x3 +=  + ln2*ln3 * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./2.
         *pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) );

      h10x3 +=  + pow(ln2,2) * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./
         2.*pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) )
         ;

      h10x3 +=  + ln1 * (  - 1./2.*mm1*pow(mm3,3) - 7./4.*pow(mm1,2)*
         pow(mm3,2) - 13./6.*pow(mm1,3)*mm3 - 13./12.*pow(mm1,4) - 1./6.
         *pow(mm1,5)*pow(mm3,-1) );

      h10x3 +=  + ln1*ln3 * (  - 1./4.*pow(mm3,4) - mm1*pow(mm3,3) - 3./
         2.*pow(mm1,2)*pow(mm3,2) - pow(mm1,3)*mm3 - 1./4.*pow(mm1,4) )
         ;

      h10x3 +=  + ln1*ln2 * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./2.
         *pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) );

      h10x3 +=  + 5./16.*pow(mm3,4) + 1./24.*pow(mm3,4)*pi2 + 17./12.*
         mm1*pow(mm3,3) + 1./6.*mm1*pow(mm3,3)*pi2 + 19./8.*pow(mm1,2)*
         pow(mm3,2) + 1./4.*pow(mm1,2)*pow(mm3,2)*pi2 + 7./4.*pow(
         mm1,3)*mm3 + 1./6.*pow(mm1,3)*mm3*pi2 + 23./48.*pow(mm1,4) + 1.
         /24.*pow(mm1,4)*pi2;
      return h10x3*pi162/pow(mm2,4);
    case 4:
      h10x4 =
       + ln3 * ( 7./6.*pow(mm3,4) + 9./2.*mm1*pow(mm3,3) + 13./2.*pow(
         mm1,2)*pow(mm3,2) + 25./6.*pow(mm1,3)*mm3 + pow(mm1,4) );

      h10x4 +=  + pow(ln3,2) * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./
         2.*pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) )
         ;

      h10x4 +=  + ln2 * (  - 5./12.*pow(mm3,4) - 2*mm1*pow(mm3,3) - 15./
         4.*pow(mm1,2)*pow(mm3,2) - 10./3.*pow(mm1,3)*mm3 - 5./4.*pow(
         mm1,4) + 1./12.*pow(mm1,6)*pow(mm3,-2) );

      h10x4 +=  + ln2*ln3 * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./2.
         *pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) );

      h10x4 +=  + ln1 * ( 1./2.*mm1*pow(mm3,3) + 7./4.*pow(mm1,2)*pow(
         mm3,2) + 13./6.*pow(mm1,3)*mm3 + pow(mm1,4) - 1./12.*pow(
         mm1,6)*pow(mm3,-2) );

      h10x4 +=  + ln1*ln3 * ( 1./4.*pow(mm3,4) + mm1*pow(mm3,3) + 3./2.
         *pow(mm1,2)*pow(mm3,2) + pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) );

      h10x4 +=  + ln1*ln2 * (  - 1./4.*pow(mm3,4) - mm1*pow(mm3,3) - 3./
         2.*pow(mm1,2)*pow(mm3,2) - pow(mm1,3)*mm3 - 1./4.*pow(mm1,4) )
         ;

      h10x4 +=  + 5./16.*pow(mm3,4) + 1./24.*pow(mm3,4)*pi2 + 13./12.*
         mm1*pow(mm3,3) + 1./6.*mm1*pow(mm3,3)*pi2 + 29./24.*pow(mm1,2)
         *pow(mm3,2) + 1./4.*pow(mm1,2)*pow(mm3,2)*pi2 + 1./4.*pow(
         mm1,3)*mm3 + 1./6.*pow(mm1,3)*mm3*pi2 - 17./48.*pow(mm1,4) + 1.
         /24.*pow(mm1,4)*pi2 - 1./6.*pow(mm1,5)*pow(mm3,-1);
      return h10x4*pi162/pow(mm2,4);
    case 5:
      h10x5 =
       + ln3 * (  - 1./20.*pow(mm1,-2)*pow(mm3,6) - 7./30.*pow(mm1,-1)*
         pow(mm3,5) - 2./5.*pow(mm3,4) - 3./10.*mm1*pow(mm3,3) - 1./12.
         *pow(mm1,2)*pow(mm3,2) );

      h10x5 +=  + ln2 * ( 1./20.*pow(mm1,-2)*pow(mm3,6) + 7./30.*pow(
         mm1,-1)*pow(mm3,5) + 2./5.*pow(mm3,4) + 3./10.*mm1*pow(mm3,3)
          + 1./6.*pow(mm1,2)*pow(mm3,2) + 3./10.*pow(mm1,3)*mm3 + 2./5.
         *pow(mm1,4) + 7./30.*pow(mm1,5)*pow(mm3,-1) + 1./20.*pow(
         mm1,6)*pow(mm3,-2) );

      h10x5 +=  + ln1 * (  - 1./12.*pow(mm1,2)*pow(mm3,2) - 3./10.*pow(
         mm1,3)*mm3 - 2./5.*pow(mm1,4) - 7./30.*pow(mm1,5)*pow(mm3,-1)
          - 1./20.*pow(mm1,6)*pow(mm3,-2) );

      h10x5 +=  - 1./10.*pow(mm1,-1)*pow(mm3,5) - 2./3.*pow(mm3,4) - 53.
         /30.*mm1*pow(mm3,3) - 12./5.*pow(mm1,2)*pow(mm3,2) - 53./30.*
         pow(mm1,3)*mm3 - 2./3.*pow(mm1,4) - 1./10.*pow(mm1,5)*pow(
         mm3,-1);
      return h10x5*pi162/pow(mm2,6);
    case 6:
      h10x6 =
       + ln3 * ( 1./20.*pow(mm1,-2)*pow(mm3,6) + 11./30.*pow(mm1,-1)*
         pow(mm3,5) + 17./15.*pow(mm3,4) + 28./15.*mm1*pow(mm3,3) + 103.
         /60.*pow(mm1,2)*pow(mm3,2) + 5./6.*pow(mm1,3)*mm3 + 1./6.*pow(
         mm1,4) );

      h10x6 +=  + ln2 * (  - 1./20.*pow(mm1,-2)*pow(mm3,6) - 11./30.*
         pow(mm1,-1)*pow(mm3,5) - 17./15.*pow(mm3,4) - 28./15.*mm1*pow(
         mm3,3) - 49./30.*pow(mm1,2)*pow(mm3,2) - 7./15.*pow(mm1,3)*mm3
          + 7./15.*pow(mm1,4) + 8./15.*pow(mm1,5)*pow(mm3,-1) + 13./60.
         *pow(mm1,6)*pow(mm3,-2) + 1./30.*pow(mm1,7)*pow(mm3,-3) );

      h10x6 +=  + ln1 * (  - 1./12.*pow(mm1,2)*pow(mm3,2) - 11./30.*
         pow(mm1,3)*mm3 - 19./30.*pow(mm1,4) - 8./15.*pow(mm1,5)*pow(
         mm3,-1) - 13./60.*pow(mm1,6)*pow(mm3,-2) - 1./30.*pow(mm1,7)*
         pow(mm3,-3) );

      h10x6 +=  + 1./10.*pow(mm1,-1)*pow(mm3,5) + 13./30.*pow(mm3,4) + 
         3./5.*mm1*pow(mm3,3) - 5./6.*pow(mm1,3)*mm3 - 9./10.*pow(
         mm1,4) - 2./5.*pow(mm1,5)*pow(mm3,-1) - 1./15.*pow(mm1,6)*pow(
         mm3,-2);
      return h10x6*pi162/pow(mm2,6);
    case 7:
      h10x7 =
       + ln3 * ( 1./5.*pow(mm1,-1)*pow(mm3,5) + 11./10.*pow(mm3,4) + 73.
         /30.*mm1*pow(mm3,3) + 27./10.*pow(mm1,2)*pow(mm3,2) + 3./2.*
         pow(mm1,3)*mm3 + 1./3.*pow(mm1,4) );

      h10x7 +=  + ln2 * (  - 1./5.*pow(mm1,-1)*pow(mm3,5) - 11./10.*
         pow(mm3,4) - 73./30.*mm1*pow(mm3,3) - 27./10.*pow(mm1,2)*pow(
         mm3,2) - 3./2.*pow(mm1,3)*mm3 - 11./30.*pow(mm1,4) - 1./10.*
         pow(mm1,5)*pow(mm3,-1) - 1./10.*pow(mm1,6)*pow(mm3,-2) - 1./30.
         *pow(mm1,7)*pow(mm3,-3) );

      h10x7 +=  + ln1 * ( 1./30.*pow(mm1,4) + 1./10.*pow(mm1,5)*pow(
         mm3,-1) + 1./10.*pow(mm1,6)*pow(mm3,-2) + 1./30.*pow(mm1,7)*
         pow(mm3,-3) );

      h10x7 +=  - 1./10.*pow(mm3,4) - 1./2.*mm1*pow(mm3,3) - 14./15.*
         pow(mm1,2)*pow(mm3,2) - 11./15.*pow(mm1,3)*mm3 - 1./10.*pow(
         mm1,4) + 1./6.*pow(mm1,5)*pow(mm3,-1) + 1./15.*pow(mm1,6)*pow(
         mm3,-2);
      return h10x7*pi162/pow(mm2,6);
    case 8:
      h10x8 =
       + ln3 * (  - 1./35.*pow(mm1,-3)*pow(mm3,7) - 3./14.*pow(mm1,-2)*
         pow(mm3,6) - 74./105.*pow(mm1,-1)*pow(mm3,5) - 559./420.*pow(
         mm3,4) - 164./105.*mm1*pow(mm3,3) - 17./15.*pow(mm1,2)*pow(
         mm3,2) - 7./15.*pow(mm1,3)*mm3 - 1./12.*pow(mm1,4) );

      h10x8 +=  + ln2 * ( 1./35.*pow(mm1,-3)*pow(mm3,7) + 3./14.*pow(
         mm1,-2)*pow(mm3,6) + 74./105.*pow(mm1,-1)*pow(mm3,5) + 559./
         420.*pow(mm3,4) + 164./105.*mm1*pow(mm3,3) + 16./15.*pow(
         mm1,2)*pow(mm3,2) + 2./15.*pow(mm1,3)*mm3 - 127./210.*pow(
         mm1,4) - 79./105.*pow(mm1,5)*pow(mm3,-1) - 97./210.*pow(mm1,6)
         *pow(mm3,-2) - 16./105.*pow(mm1,7)*pow(mm3,-3) - 3./140.*pow(
         mm1,8)*pow(mm3,-4) );

      h10x8 +=  + ln1 * ( 1./15.*pow(mm1,2)*pow(mm3,2) + 1./3.*pow(
         mm1,3)*mm3 + 289./420.*pow(mm1,4) + 79./105.*pow(mm1,5)*pow(
         mm3,-1) + 97./210.*pow(mm1,6)*pow(mm3,-2) + 16./105.*pow(
         mm1,7)*pow(mm3,-3) + 3./140.*pow(mm1,8)*pow(mm3,-4) );

      h10x8 +=  - 2./35.*pow(mm1,-2)*pow(mm3,6) - 2./5.*pow(mm1,-1)*
         pow(mm3,5) - 95./84.*pow(mm3,4) - 111./70.*mm1*pow(mm3,3) - 
         129./140.*pow(mm1,2)*pow(mm3,2) + 83./210.*pow(mm1,3)*mm3 + 
         149./140.*pow(mm1,4) + 11./14.*pow(mm1,5)*pow(mm3,-1) + 17./60.
         *pow(mm1,6)*pow(mm3,-2) + 3./70.*pow(mm1,7)*pow(mm3,-3);
      return h10x8*pi162/pow(mm2,8);
    default:
      std::cout << "h10sing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10x1 =
       + ln3 * (  - 1./4.*pow(mm2,3) - 1./4.*mm1*pow(mm2,2) + 1./4.*
         pow(mm1,2)*mm2 + 1./4.*pow(mm1,3) );

      h10x1 +=  + pow(ln3,2) * ( 1./4.*pow(mm2,3) + 3./4.*mm1*pow(
         mm2,2) + 3./4.*pow(mm1,2)*mm2 + 1./4.*pow(mm1,3) );

      h10x1 +=  + ln2 * (  - 1./4.*pow(mm2,3) - 3./4.*mm1*pow(mm2,2) - 
         1./2.*pow(mm1,2)*mm2 );

      h10x1 +=  + ln2*ln3 * ( 1./2.*pow(mm2,3) + mm1*pow(mm2,2) + 3./4.
         *pow(mm1,2)*mm2 + 1./4.*pow(mm1,3) );

      h10x1 +=  + pow(ln2,2) * ( 1./4.*pow(mm2,3) + 1./4.*mm1*pow(
         mm2,2) );

      h10x1 +=  + ln1 * (  - pow(mm1,2)*mm2 - pow(mm1,3) );

      h10x1 +=  + ln1*ln3 * ( 1./2.*mm1*pow(mm2,2) + 3./4.*pow(mm1,2)*
         mm2 + 1./4.*pow(mm1,3) );

      h10x1 +=  + ln1*ln2 * (  - 1./2.*mm1*pow(mm2,2) - 3./4.*pow(
         mm1,2)*mm2 - 1./4.*pow(mm1,3) );

      h10x1 +=  + 9./8.*pow(mm2,3) + 1./12.*pow(mm2,3)*pi2 + 9./4.*mm1*
         pow(mm2,2) + 1./6.*mm1*pow(mm2,2)*pi2 + 33./16.*pow(mm1,2)*mm2
          + 1./8.*pow(mm1,2)*mm2*pi2 + 15./16.*pow(mm1,3) + 1./24.*pow(
         mm1,3)*pi2;
      return h10x1*pi162/mm3;
    case 2:
      h10x2 =
       + ln3 * (  - 1./6.*pow(mm1,-1)*pow(mm2,5) - pow(mm2,4) - 5./2.*
         mm1*pow(mm2,3) - 10./3.*pow(mm1,2)*pow(mm2,2) - 5./2.*pow(
         mm1,3)*mm2 - pow(mm1,4) - 1./6.*pow(mm1,5)*pow(mm2,-1) );

      h10x2 +=  + ln2 * ( 1./6.*pow(mm1,-1)*pow(mm2,5) + 1./2.*pow(
         mm2,4) + 1./2.*mm1*pow(mm2,3) + 1./6.*pow(mm1,2)*pow(mm2,2) );

      h10x2 +=  + ln1 * ( 1./6.*pow(mm1,2)*pow(mm2,2) + 1./2.*pow(
         mm1,3)*mm2 + 1./2.*pow(mm1,4) + 1./6.*pow(mm1,5)*pow(mm2,-1) )
         ;

      h10x2 +=  - 7./24.*pow(mm2,4) - 7./6.*mm1*pow(mm2,3) - 7./4.*pow(
         mm1,2)*pow(mm2,2) - 7./6.*pow(mm1,3)*mm2 - 7./24.*pow(mm1,4);
      return h10x2*pi162/pow(mm3,4);
    case 3:
      h10x3 =
       + ln3 * (  - 5./12.*pow(mm2,4) - 2*mm1*pow(mm2,3) - 15./4.*pow(
         mm1,2)*pow(mm2,2) - 10./3.*pow(mm1,3)*mm2 - 5./4.*pow(mm1,4)
          + 1./12.*pow(mm1,6)*pow(mm2,-2) );

      h10x3 +=  + ln2 * ( 7./6.*pow(mm2,4) + 9./2.*mm1*pow(mm2,3) + 13./
         2.*pow(mm1,2)*pow(mm2,2) + 25./6.*pow(mm1,3)*mm2 + pow(mm1,4)
          );

      h10x3 +=  + ln2*ln3 * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./2.
         *pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) );

      h10x3 +=  + pow(ln2,2) * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./
         2.*pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) )
         ;

      h10x3 +=  + ln1 * ( 1./2.*mm1*pow(mm2,3) + 7./4.*pow(mm1,2)*pow(
         mm2,2) + 13./6.*pow(mm1,3)*mm2 + pow(mm1,4) - 1./12.*pow(
         mm1,6)*pow(mm2,-2) );

      h10x3 +=  + ln1*ln3 * (  - 1./4.*pow(mm2,4) - mm1*pow(mm2,3) - 3./
         2.*pow(mm1,2)*pow(mm2,2) - pow(mm1,3)*mm2 - 1./4.*pow(mm1,4) )
         ;

      h10x3 +=  + ln1*ln2 * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./2.
         *pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) );

      h10x3 +=  + 5./16.*pow(mm2,4) + 1./24.*pow(mm2,4)*pi2 + 13./12.*
         mm1*pow(mm2,3) + 1./6.*mm1*pow(mm2,3)*pi2 + 29./24.*pow(mm1,2)
         *pow(mm2,2) + 1./4.*pow(mm1,2)*pow(mm2,2)*pi2 + 1./4.*pow(
         mm1,3)*mm2 + 1./6.*pow(mm1,3)*mm2*pi2 - 17./48.*pow(mm1,4) + 1.
         /24.*pow(mm1,4)*pi2 - 1./6.*pow(mm1,5)*pow(mm2,-1);
      return h10x3*pi162/pow(mm3,4);
    case 4:
      h10x4 =
       + ln3 * ( 7./6.*pow(mm2,4) + 29./6.*mm1*pow(mm2,3) + 23./3.*pow(
         mm1,2)*pow(mm2,2) + 17./3.*pow(mm1,3)*mm2 + 11./6.*pow(mm1,4)
          + 1./6.*pow(mm1,5)*pow(mm2,-1) );

      h10x4 +=  + pow(ln3,2) * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./
         2.*pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) )
         ;

      h10x4 +=  + ln2 * (  - 5./12.*pow(mm2,4) - 4./3.*mm1*pow(mm2,3)
          - 17./12.*pow(mm1,2)*pow(mm2,2) - 1./2.*pow(mm1,3)*mm2 );

      h10x4 +=  + ln2*ln3 * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./2.
         *pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) );

      h10x4 +=  + ln1 * (  - 1./2.*mm1*pow(mm2,3) - 7./4.*pow(mm1,2)*
         pow(mm2,2) - 13./6.*pow(mm1,3)*mm2 - 13./12.*pow(mm1,4) - 1./6.
         *pow(mm1,5)*pow(mm2,-1) );

      h10x4 +=  + ln1*ln3 * ( 1./4.*pow(mm2,4) + mm1*pow(mm2,3) + 3./2.
         *pow(mm1,2)*pow(mm2,2) + pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) );

      h10x4 +=  + ln1*ln2 * (  - 1./4.*pow(mm2,4) - mm1*pow(mm2,3) - 3./
         2.*pow(mm1,2)*pow(mm2,2) - pow(mm1,3)*mm2 - 1./4.*pow(mm1,4) )
         ;

      h10x4 +=  + 5./16.*pow(mm2,4) + 1./24.*pow(mm2,4)*pi2 + 17./12.*
         mm1*pow(mm2,3) + 1./6.*mm1*pow(mm2,3)*pi2 + 19./8.*pow(mm1,2)*
         pow(mm2,2) + 1./4.*pow(mm1,2)*pow(mm2,2)*pi2 + 7./4.*pow(
         mm1,3)*mm2 + 1./6.*pow(mm1,3)*mm2*pi2 + 23./48.*pow(mm1,4) + 1.
         /24.*pow(mm1,4)*pi2;
      return h10x4*pi162/pow(mm3,4);
    case 5:
      h10x5 =
       + ln3 * (  - 1./20.*pow(mm1,-2)*pow(mm2,6) - 11./30.*pow(mm1,-1)
         *pow(mm2,5) - 17./15.*pow(mm2,4) - 28./15.*mm1*pow(mm2,3) - 49.
         /30.*pow(mm1,2)*pow(mm2,2) - 7./15.*pow(mm1,3)*mm2 + 7./15.*
         pow(mm1,4) + 8./15.*pow(mm1,5)*pow(mm2,-1) + 13./60.*pow(
         mm1,6)*pow(mm2,-2) + 1./30.*pow(mm1,7)*pow(mm2,-3) );

      h10x5 +=  + ln2 * ( 1./20.*pow(mm1,-2)*pow(mm2,6) + 11./30.*pow(
         mm1,-1)*pow(mm2,5) + 17./15.*pow(mm2,4) + 28./15.*mm1*pow(
         mm2,3) + 103./60.*pow(mm1,2)*pow(mm2,2) + 5./6.*pow(mm1,3)*mm2
          + 1./6.*pow(mm1,4) );

      h10x5 +=  + ln1 * (  - 1./12.*pow(mm1,2)*pow(mm2,2) - 11./30.*
         pow(mm1,3)*mm2 - 19./30.*pow(mm1,4) - 8./15.*pow(mm1,5)*pow(
         mm2,-1) - 13./60.*pow(mm1,6)*pow(mm2,-2) - 1./30.*pow(mm1,7)*
         pow(mm2,-3) );

      h10x5 +=  + 1./10.*pow(mm1,-1)*pow(mm2,5) + 13./30.*pow(mm2,4) + 
         3./5.*mm1*pow(mm2,3) - 5./6.*pow(mm1,3)*mm2 - 9./10.*pow(
         mm1,4) - 2./5.*pow(mm1,5)*pow(mm2,-1) - 1./15.*pow(mm1,6)*pow(
         mm2,-2);
      return h10x5*pi162/pow(mm3,6);
    case 6:
      h10x6 =
       + ln3 * ( 1./20.*pow(mm1,-2)*pow(mm2,6) + 7./30.*pow(mm1,-1)*
         pow(mm2,5) + 2./5.*pow(mm2,4) + 3./10.*mm1*pow(mm2,3) + 1./6.*
         pow(mm1,2)*pow(mm2,2) + 3./10.*pow(mm1,3)*mm2 + 2./5.*pow(
         mm1,4) + 7./30.*pow(mm1,5)*pow(mm2,-1) + 1./20.*pow(mm1,6)*
         pow(mm2,-2) );

      h10x6 +=  + ln2 * (  - 1./20.*pow(mm1,-2)*pow(mm2,6) - 7./30.*
         pow(mm1,-1)*pow(mm2,5) - 2./5.*pow(mm2,4) - 3./10.*mm1*pow(
         mm2,3) - 1./12.*pow(mm1,2)*pow(mm2,2) );

      h10x6 +=  + ln1 * (  - 1./12.*pow(mm1,2)*pow(mm2,2) - 3./10.*pow(
         mm1,3)*mm2 - 2./5.*pow(mm1,4) - 7./30.*pow(mm1,5)*pow(mm2,-1)
          - 1./20.*pow(mm1,6)*pow(mm2,-2) );

      h10x6 +=  - 1./10.*pow(mm1,-1)*pow(mm2,5) - 2./3.*pow(mm2,4) - 53.
         /30.*mm1*pow(mm2,3) - 12./5.*pow(mm1,2)*pow(mm2,2) - 53./30.*
         pow(mm1,3)*mm2 - 2./3.*pow(mm1,4) - 1./10.*pow(mm1,5)*pow(
         mm2,-1);
      return h10x6*pi162/pow(mm3,6);
    case 7:
      h10x7 =
       + ln3 * (  - 1./5.*pow(mm1,-1)*pow(mm2,5) - 11./10.*pow(mm2,4)
          - 73./30.*mm1*pow(mm2,3) - 27./10.*pow(mm1,2)*pow(mm2,2) - 3./
         2.*pow(mm1,3)*mm2 - 11./30.*pow(mm1,4) - 1./10.*pow(mm1,5)*
         pow(mm2,-1) - 1./10.*pow(mm1,6)*pow(mm2,-2) - 1./30.*pow(
         mm1,7)*pow(mm2,-3) );

      h10x7 +=  + ln2 * ( 1./5.*pow(mm1,-1)*pow(mm2,5) + 11./10.*pow(
         mm2,4) + 73./30.*mm1*pow(mm2,3) + 27./10.*pow(mm1,2)*pow(
         mm2,2) + 3./2.*pow(mm1,3)*mm2 + 1./3.*pow(mm1,4) );

      h10x7 +=  + ln1 * ( 1./30.*pow(mm1,4) + 1./10.*pow(mm1,5)*pow(
         mm2,-1) + 1./10.*pow(mm1,6)*pow(mm2,-2) + 1./30.*pow(mm1,7)*
         pow(mm2,-3) );

      h10x7 +=  - 1./10.*pow(mm2,4) - 1./2.*mm1*pow(mm2,3) - 14./15.*
         pow(mm1,2)*pow(mm2,2) - 11./15.*pow(mm1,3)*mm2 - 1./10.*pow(
         mm1,4) + 1./6.*pow(mm1,5)*pow(mm2,-1) + 1./15.*pow(mm1,6)*pow(
         mm2,-2);
      return h10x7*pi162/pow(mm3,6);
    case 8:
      h10x8 =
       + ln3 * ( 1./35.*pow(mm1,-3)*pow(mm2,7) + 3./14.*pow(mm1,-2)*
         pow(mm2,6) + 74./105.*pow(mm1,-1)*pow(mm2,5) + 559./420.*pow(
         mm2,4) + 164./105.*mm1*pow(mm2,3) + 16./15.*pow(mm1,2)*pow(
         mm2,2) + 2./15.*pow(mm1,3)*mm2 - 127./210.*pow(mm1,4) - 79./
         105.*pow(mm1,5)*pow(mm2,-1) - 97./210.*pow(mm1,6)*pow(mm2,-2)
          - 16./105.*pow(mm1,7)*pow(mm2,-3) - 3./140.*pow(mm1,8)*pow(
         mm2,-4) );

      h10x8 +=  + ln2 * (  - 1./35.*pow(mm1,-3)*pow(mm2,7) - 3./14.*
         pow(mm1,-2)*pow(mm2,6) - 74./105.*pow(mm1,-1)*pow(mm2,5) - 559.
         /420.*pow(mm2,4) - 164./105.*mm1*pow(mm2,3) - 17./15.*pow(
         mm1,2)*pow(mm2,2) - 7./15.*pow(mm1,3)*mm2 - 1./12.*pow(mm1,4)
          );

      h10x8 +=  + ln1 * ( 1./15.*pow(mm1,2)*pow(mm2,2) + 1./3.*pow(
         mm1,3)*mm2 + 289./420.*pow(mm1,4) + 79./105.*pow(mm1,5)*pow(
         mm2,-1) + 97./210.*pow(mm1,6)*pow(mm2,-2) + 16./105.*pow(
         mm1,7)*pow(mm2,-3) + 3./140.*pow(mm1,8)*pow(mm2,-4) );

      h10x8 +=  - 2./35.*pow(mm1,-2)*pow(mm2,6) - 2./5.*pow(mm1,-1)*
         pow(mm2,5) - 95./84.*pow(mm2,4) - 111./70.*mm1*pow(mm2,3) - 
         129./140.*pow(mm1,2)*pow(mm2,2) + 83./210.*pow(mm1,3)*mm2 + 
         149./140.*pow(mm1,4) + 11./14.*pow(mm1,5)*pow(mm2,-1) + 17./60.
         *pow(mm1,6)*pow(mm2,-2) + 3./70.*pow(mm1,7)*pow(mm2,-3);
      return h10x8*pi162/pow(mm3,8);
    default:
      std::cout << "h10sing, wrong iprop = "<<iprop<<'\n';
      return 0.;
     }
  }
  std::cout << "funny masses in h10sing\n";
  return 0.;
}

double h10psing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h10px1,h10px2,h10px3,h10px4,h10px5,h10px6,h10px7,h10px8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10px1 =
       + ln3 * (  - 1./60.*pow(mm2,-2)*pow(mm3,5) - 7./60.*pow(mm2,-1)*
         pow(mm3,4) - 11./60.*pow(mm3,3) - 1./12.*mm2*pow(mm3,2) );

      h10px1 +=  + ln2 * (  - 1./12.*pow(mm2,2)*mm3 - 11./60.*pow(
         mm2,3) - 7./60.*pow(mm2,4)*pow(mm3,-1) - 1./60.*pow(mm2,5)*
         pow(mm3,-2) );

      h10px1 +=  + ln1 * ( 1./60.*pow(mm2,-2)*pow(mm3,5) + 7./60.*pow(
         mm2,-1)*pow(mm3,4) + 7./20.*pow(mm3,3) + 7./12.*mm2*pow(mm3,2)
          + 7./12.*pow(mm2,2)*mm3 + 7./20.*pow(mm2,3) + 7./60.*pow(
         mm2,4)*pow(mm3,-1) + 1./60.*pow(mm2,5)*pow(mm3,-2) );

      h10px1 +=  - 1./30.*pow(mm2,-1)*pow(mm3,4) - 13./360.*pow(mm3,3)
          + 7./120.*mm2*pow(mm3,2) + 7./120.*pow(mm2,2)*mm3 - 13./360.*
         pow(mm2,3) - 1./30.*pow(mm2,4)*pow(mm3,-1);
      return h10px1*pi162/pow(mm1,3);
    case 2:
      h10px2 =
       + ln3 * ( 1./105.*pow(mm2,-3)*pow(mm3,6) + 11./210.*pow(mm2,-2)*
         pow(mm3,5) + 23./210.*pow(mm2,-1)*pow(mm3,4) + 1./10.*pow(
         mm3,3) + 1./30.*mm2*pow(mm3,2) );

      h10px2 +=  + ln2 * ( 1./30.*pow(mm2,2)*mm3 + 1./10.*pow(mm2,3) + 
         23./210.*pow(mm2,4)*pow(mm3,-1) + 11./210.*pow(mm2,5)*pow(
         mm3,-2) + 1./105.*pow(mm2,6)*pow(mm3,-3) );

      h10px2 +=  + ln1 * (  - 1./105.*pow(mm2,-3)*pow(mm3,6) - 11./210.
         *pow(mm2,-2)*pow(mm3,5) - 23./210.*pow(mm2,-1)*pow(mm3,4) - 1./
         10.*pow(mm3,3) - 1./30.*mm2*pow(mm3,2) - 1./30.*pow(mm2,2)*mm3
          - 1./10.*pow(mm2,3) - 23./210.*pow(mm2,4)*pow(mm3,-1) - 11./
         210.*pow(mm2,5)*pow(mm3,-2) - 1./105.*pow(mm2,6)*pow(mm3,-3) )
         ;

      h10px2 +=  + 2./105.*pow(mm2,-2)*pow(mm3,5) + 2./21.*pow(mm2,-1)*
         pow(mm3,4) + 8./35.*pow(mm3,3) + 12./35.*mm2*pow(mm3,2) + 12./
         35.*pow(mm2,2)*mm3 + 8./35.*pow(mm2,3) + 2./21.*pow(mm2,4)*
         pow(mm3,-1) + 2./105.*pow(mm2,5)*pow(mm3,-2);
      return h10px2*pi162/pow(mm1,5);
    case 3:
      h10px3 =
       + ln3 * ( 1./140.*pow(mm2,-4)*pow(mm3,7) + 23./420.*pow(mm2,-3)*
         pow(mm3,6) + 6./35.*pow(mm2,-2)*pow(mm3,5) + 9./35.*pow(
         mm2,-1)*pow(mm3,4) + 11./60.*pow(mm3,3) + 1./20.*mm2*pow(
         mm3,2) );

      h10px3 +=  + ln2 * (  - 1./12.*pow(mm3,3) - 7./20.*mm2*pow(mm3,2)
          - 3./5.*pow(mm2,2)*mm3 - 19./35.*pow(mm2,3) - 39./140.*pow(
         mm2,4)*pow(mm3,-1) - 11./140.*pow(mm2,5)*pow(mm3,-2) - 1./105.
         *pow(mm2,6)*pow(mm3,-3) );

      h10px3 +=  + ln1 * (  - 1./140.*pow(mm2,-4)*pow(mm3,7) - 23./420.
         *pow(mm2,-3)*pow(mm3,6) - 6./35.*pow(mm2,-2)*pow(mm3,5) - 9./
         35.*pow(mm2,-1)*pow(mm3,4) - 1./10.*pow(mm3,3) + 3./10.*mm2*
         pow(mm3,2) + 3./5.*pow(mm2,2)*mm3 + 19./35.*pow(mm2,3) + 39./
         140.*pow(mm2,4)*pow(mm3,-1) + 11./140.*pow(mm2,5)*pow(mm3,-2)
          + 1./105.*pow(mm2,6)*pow(mm3,-3) );

      h10px3 +=  + 1./70.*pow(mm2,-3)*pow(mm3,6) + 43./420.*pow(mm2,-2)
         *pow(mm3,5) + 41./140.*pow(mm2,-1)*pow(mm3,4) + 38./105.*pow(
         mm3,3) + 1./21.*mm2*pow(mm3,2) - 149./420.*pow(mm2,2)*mm3 - 
         157./420.*pow(mm2,3) - 31./210.*pow(mm2,4)*pow(mm3,-1) - 2./
         105.*pow(mm2,5)*pow(mm3,-2);
      return h10px3*pi162/pow(mm1,5);
    case 4:
      h10px4 =
       + ln3 * (  - 1./105.*pow(mm2,-3)*pow(mm3,6) - 11./140.*pow(
         mm2,-2)*pow(mm3,5) - 39./140.*pow(mm2,-1)*pow(mm3,4) - 19./35.
         *pow(mm3,3) - 3./5.*mm2*pow(mm3,2) - 7./20.*pow(mm2,2)*mm3 - 1.
         /12.*pow(mm2,3) );

      h10px4 +=  + ln2 * ( 1./20.*pow(mm2,2)*mm3 + 11./60.*pow(mm2,3)
          + 9./35.*pow(mm2,4)*pow(mm3,-1) + 6./35.*pow(mm2,5)*pow(
         mm3,-2) + 23./420.*pow(mm2,6)*pow(mm3,-3) + 1./140.*pow(mm2,7)
         *pow(mm3,-4) );

      h10px4 +=  + ln1 * ( 1./105.*pow(mm2,-3)*pow(mm3,6) + 11./140.*
         pow(mm2,-2)*pow(mm3,5) + 39./140.*pow(mm2,-1)*pow(mm3,4) + 19./
         35.*pow(mm3,3) + 3./5.*mm2*pow(mm3,2) + 3./10.*pow(mm2,2)*mm3
          - 1./10.*pow(mm2,3) - 9./35.*pow(mm2,4)*pow(mm3,-1) - 6./35.*
         pow(mm2,5)*pow(mm3,-2) - 23./420.*pow(mm2,6)*pow(mm3,-3) - 1./
         140.*pow(mm2,7)*pow(mm3,-4) );

      h10px4 +=  - 2./105.*pow(mm2,-2)*pow(mm3,5) - 31./210.*pow(
         mm2,-1)*pow(mm3,4) - 157./420.*pow(mm3,3) - 149./420.*mm2*pow(
         mm3,2) + 1./21.*pow(mm2,2)*mm3 + 38./105.*pow(mm2,3) + 41./140.
         *pow(mm2,4)*pow(mm3,-1) + 43./420.*pow(mm2,5)*pow(mm3,-2) + 1./
         70.*pow(mm2,6)*pow(mm3,-3);
      return h10px4*pi162/pow(mm1,5);
    case 5:
      h10px5 =
       + ln3 * (  - 2./315.*pow(mm2,-5)*pow(mm3,8) - 1./21.*pow(mm2,-4)
         *pow(mm3,7) - 16./105.*pow(mm2,-3)*pow(mm3,6) - 169./630.*pow(
         mm2,-2)*pow(mm3,5) - 19./70.*pow(mm2,-1)*pow(mm3,4) - 31./210.
         *pow(mm3,3) - 1./30.*mm2*pow(mm3,2) );

      h10px5 +=  + ln2 * ( 1./30.*pow(mm3,3) + 1./6.*mm2*pow(mm3,2) + 
         13./35.*pow(mm2,2)*mm3 + 17./35.*pow(mm2,3) + 127./315.*pow(
         mm2,4)*pow(mm3,-1) + 22./105.*pow(mm2,5)*pow(mm3,-2) + 13./210.
         *pow(mm2,6)*pow(mm3,-3) + 1./126.*pow(mm2,7)*pow(mm3,-4) );

      h10px5 +=  + ln1 * ( 2./315.*pow(mm2,-5)*pow(mm3,8) + 1./21.*pow(
         mm2,-4)*pow(mm3,7) + 16./105.*pow(mm2,-3)*pow(mm3,6) + 169./
         630.*pow(mm2,-2)*pow(mm3,5) + 19./70.*pow(mm2,-1)*pow(mm3,4)
          + 4./35.*pow(mm3,3) - 2./15.*mm2*pow(mm3,2) - 13./35.*pow(
         mm2,2)*mm3 - 17./35.*pow(mm2,3) - 127./315.*pow(mm2,4)*pow(
         mm3,-1) - 22./105.*pow(mm2,5)*pow(mm3,-2) - 13./210.*pow(
         mm2,6)*pow(mm3,-3) - 1./126.*pow(mm2,7)*pow(mm3,-4) );

      h10px5 +=  - 4./315.*pow(mm2,-4)*pow(mm3,7) - 4./45.*pow(mm2,-3)*
         pow(mm3,6) - 247./945.*pow(mm2,-2)*pow(mm3,5) - 26./63.*pow(
         mm2,-1)*pow(mm3,4) - 197./630.*pow(mm3,3) + 191./1890.*mm2*
         pow(mm3,2) + 167./315.*pow(mm2,2)*mm3 + 191./315.*pow(mm2,3)
          + 137./378.*pow(mm2,4)*pow(mm3,-1) + 73./630.*pow(mm2,5)*pow(
         mm3,-2) + 1./63.*pow(mm2,6)*pow(mm3,-3);
      return h10px5*pi162/pow(mm1,7);
    case 6:
      h10px6 =
       + ln3 * ( 1./126.*pow(mm2,-4)*pow(mm3,7) + 13./210.*pow(mm2,-3)*
         pow(mm3,6) + 22./105.*pow(mm2,-2)*pow(mm3,5) + 127./315.*pow(
         mm2,-1)*pow(mm3,4) + 17./35.*pow(mm3,3) + 13./35.*mm2*pow(
         mm3,2) + 1./6.*pow(mm2,2)*mm3 + 1./30.*pow(mm2,3) );

      h10px6 +=  + ln2 * (  - 1./30.*pow(mm2,2)*mm3 - 31./210.*pow(
         mm2,3) - 19./70.*pow(mm2,4)*pow(mm3,-1) - 169./630.*pow(mm2,5)
         *pow(mm3,-2) - 16./105.*pow(mm2,6)*pow(mm3,-3) - 1./21.*pow(
         mm2,7)*pow(mm3,-4) - 2./315.*pow(mm2,8)*pow(mm3,-5) );

      h10px6 +=  + ln1 * (  - 1./126.*pow(mm2,-4)*pow(mm3,7) - 13./210.
         *pow(mm2,-3)*pow(mm3,6) - 22./105.*pow(mm2,-2)*pow(mm3,5) - 
         127./315.*pow(mm2,-1)*pow(mm3,4) - 17./35.*pow(mm3,3) - 13./35.
         *mm2*pow(mm3,2) - 2./15.*pow(mm2,2)*mm3 + 4./35.*pow(mm2,3) + 
         19./70.*pow(mm2,4)*pow(mm3,-1) + 169./630.*pow(mm2,5)*pow(
         mm3,-2) + 16./105.*pow(mm2,6)*pow(mm3,-3) + 1./21.*pow(mm2,7)*
         pow(mm3,-4) + 2./315.*pow(mm2,8)*pow(mm3,-5) );

      h10px6 +=  + 1./63.*pow(mm2,-3)*pow(mm3,6) + 73./630.*pow(mm2,-2)
         *pow(mm3,5) + 137./378.*pow(mm2,-1)*pow(mm3,4) + 191./315.*
         pow(mm3,3) + 167./315.*mm2*pow(mm3,2) + 191./1890.*pow(mm2,2)*
         mm3 - 197./630.*pow(mm2,3) - 26./63.*pow(mm2,4)*pow(mm3,-1) - 
         247./945.*pow(mm2,5)*pow(mm3,-2) - 4./45.*pow(mm2,6)*pow(
         mm3,-3) - 4./315.*pow(mm2,7)*pow(mm3,-4);
      return h10px6*pi162/pow(mm1,7);
    case 7:
      h10px7 =
       + ln3 * ( 2./315.*pow(mm2,-5)*pow(mm3,8) + 11./180.*pow(mm2,-4)*
         pow(mm3,7) + 109./420.*pow(mm2,-3)*pow(mm3,6) + 803./1260.*
         pow(mm2,-2)*pow(mm3,5) + 253./252.*pow(mm2,-1)*pow(mm3,4) + 21.
         /20.*pow(mm3,3) + 299./420.*mm2*pow(mm3,2) + 17./60.*pow(
         mm2,2)*mm3 + 1./20.*pow(mm2,3) );

      h10px7 +=  + ln2 * ( 1./20.*pow(mm3,3) + 17./60.*mm2*pow(mm3,2)
          + 299./420.*pow(mm2,2)*mm3 + 21./20.*pow(mm2,3) + 253./252.*
         pow(mm2,4)*pow(mm3,-1) + 803./1260.*pow(mm2,5)*pow(mm3,-2) + 
         109./420.*pow(mm2,6)*pow(mm3,-3) + 11./180.*pow(mm2,7)*pow(
         mm3,-4) + 2./315.*pow(mm2,8)*pow(mm3,-5) );

      h10px7 +=  + ln1 * (  - 2./315.*pow(mm2,-5)*pow(mm3,8) - 11./180.
         *pow(mm2,-4)*pow(mm3,7) - 109./420.*pow(mm2,-3)*pow(mm3,6) - 
         803./1260.*pow(mm2,-2)*pow(mm3,5) - 253./252.*pow(mm2,-1)*pow(
         mm3,4) - 11./10.*pow(mm3,3) - 209./210.*mm2*pow(mm3,2) - 209./
         210.*pow(mm2,2)*mm3 - 11./10.*pow(mm2,3) - 253./252.*pow(
         mm2,4)*pow(mm3,-1) - 803./1260.*pow(mm2,5)*pow(mm3,-2) - 109./
         420.*pow(mm2,6)*pow(mm3,-3) - 11./180.*pow(mm2,7)*pow(mm3,-4)
          - 2./315.*pow(mm2,8)*pow(mm3,-5) );

      h10px7 +=  + 4./315.*pow(mm2,-4)*pow(mm3,7) + 73./630.*pow(
         mm2,-3)*pow(mm3,6) + 1747./3780.*pow(mm2,-2)*pow(mm3,5) + 3979.
         /3780.*pow(mm2,-1)*pow(mm3,4) + 1979./1260.*pow(mm3,3) + 6731./
         3780.*mm2*pow(mm3,2) + 6731./3780.*pow(mm2,2)*mm3 + 1979./1260.
         *pow(mm2,3) + 3979./3780.*pow(mm2,4)*pow(mm3,-1) + 1747./3780.
         *pow(mm2,5)*pow(mm3,-2) + 73./630.*pow(mm2,6)*pow(mm3,-3) + 4./
         315.*pow(mm2,7)*pow(mm3,-4);
      return h10px7*pi162/pow(mm1,7);
    case 8:
      h10px8 =
       + ln3 * (  - 5./693.*pow(mm2,-6)*pow(mm3,9) - 27./385.*pow(
         mm2,-5)*pow(mm3,8) - 47./154.*pow(mm2,-4)*pow(mm3,7) - 155./
         198.*pow(mm2,-3)*pow(mm3,6) - 1009./770.*pow(mm2,-2)*pow(
         mm3,5) - 3./2.*pow(mm2,-1)*pow(mm3,4) - 751./630.*pow(mm3,3)
          - 9./14.*mm2*pow(mm3,2) - 3./14.*pow(mm2,2)*mm3 - 1./30.*pow(
         mm2,3) );

      h10px8 +=  + ln2 * (  - 1./30.*pow(mm3,3) - 3./14.*mm2*pow(mm3,2)
          - 9./14.*pow(mm2,2)*mm3 - 751./630.*pow(mm2,3) - 3./2.*pow(
         mm2,4)*pow(mm3,-1) - 1009./770.*pow(mm2,5)*pow(mm3,-2) - 155./
         198.*pow(mm2,6)*pow(mm3,-3) - 47./154.*pow(mm2,7)*pow(mm3,-4)
          - 27./385.*pow(mm2,8)*pow(mm3,-5) - 5./693.*pow(mm2,9)*pow(
         mm3,-6) );

      h10px8 +=  + ln1 * ( 5./693.*pow(mm2,-6)*pow(mm3,9) + 27./385.*
         pow(mm2,-5)*pow(mm3,8) + 47./154.*pow(mm2,-4)*pow(mm3,7) + 155.
         /198.*pow(mm2,-3)*pow(mm3,6) + 1009./770.*pow(mm2,-2)*pow(
         mm3,5) + 3./2.*pow(mm2,-1)*pow(mm3,4) + 386./315.*pow(mm3,3)
          + 6./7.*mm2*pow(mm3,2) + 6./7.*pow(mm2,2)*mm3 + 386./315.*
         pow(mm2,3) + 3./2.*pow(mm2,4)*pow(mm3,-1) + 1009./770.*pow(
         mm2,5)*pow(mm3,-2) + 155./198.*pow(mm2,6)*pow(mm3,-3) + 47./
         154.*pow(mm2,7)*pow(mm3,-4) + 27./385.*pow(mm2,8)*pow(mm3,-5)
          + 5./693.*pow(mm2,9)*pow(mm3,-6) );

      h10px8 +=  - 10./693.*pow(mm2,-5)*pow(mm3,8) - 461./3465.*pow(
         mm2,-4)*pow(mm3,7) - 5666./10395.*pow(mm2,-3)*pow(mm3,6) - 
         4517./3465.*pow(mm2,-2)*pow(mm3,5) - 6962./3465.*pow(mm2,-1)*
         pow(mm3,4) - 1061./495.*pow(mm3,3) - 2171./1155.*mm2*pow(
         mm3,2) - 2171./1155.*pow(mm2,2)*mm3 - 1061./495.*pow(mm2,3) - 
         6962./3465.*pow(mm2,4)*pow(mm3,-1) - 4517./3465.*pow(mm2,5)*
         pow(mm3,-2) - 5666./10395.*pow(mm2,6)*pow(mm3,-3) - 461./3465.
         *pow(mm2,7)*pow(mm3,-4) - 10./693.*pow(mm2,8)*pow(mm3,-5);
      return h10px8*pi162/pow(mm1,9);
    default:
      std::cout << "h10psing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10px1 =
       + ln3 * (  - 1./15.*pow(mm1,-1)*pow(mm3,4) - 3./20.*pow(mm3,3)
          - 1./12.*mm1*pow(mm3,2) );

      h10px1 +=  + ln2 * ( 1./15.*pow(mm1,-1)*pow(mm3,4) + 19./60.*pow(
         mm3,3) + 7./12.*mm1*pow(mm3,2) + 1./2.*pow(mm1,2)*mm3 + 1./6.*
         pow(mm1,3) - 1./60.*pow(mm1,4)*pow(mm3,-1) - 1./60.*pow(mm1,5)
         *pow(mm3,-2) );

      h10px1 +=  + ln1 * ( 1./60.*pow(mm1,4)*pow(mm3,-1) + 1./60.*pow(
         mm1,5)*pow(mm3,-2) );

      h10px1 +=  + 47./360.*pow(mm3,3) + 47./120.*mm1*pow(mm3,2) + 17./
         40.*pow(mm1,2)*mm3 + 71./360.*pow(mm1,3) + 1./30.*pow(mm1,4)*
         pow(mm3,-1);
      return h10px1*pi162/pow(mm2,3);
    case 2:
      h10px2 =
       + ln3 * ( 1./105.*pow(mm1,-3)*pow(mm3,6) + 11./210.*pow(mm1,-2)*
         pow(mm3,5) + 23./210.*pow(mm1,-1)*pow(mm3,4) + 1./10.*pow(
         mm3,3) + 1./30.*mm1*pow(mm3,2) );

      h10px2 +=  + ln2 * (  - 1./105.*pow(mm1,-3)*pow(mm3,6) - 11./210.
         *pow(mm1,-2)*pow(mm3,5) - 23./210.*pow(mm1,-1)*pow(mm3,4) - 1./
         10.*pow(mm3,3) - 1./30.*mm1*pow(mm3,2) - 1./30.*pow(mm1,2)*mm3
          - 1./10.*pow(mm1,3) - 23./210.*pow(mm1,4)*pow(mm3,-1) - 11./
         210.*pow(mm1,5)*pow(mm3,-2) - 1./105.*pow(mm1,6)*pow(mm3,-3) )
         ;

      h10px2 +=  + ln1 * ( 1./30.*pow(mm1,2)*mm3 + 1./10.*pow(mm1,3) + 
         23./210.*pow(mm1,4)*pow(mm3,-1) + 11./210.*pow(mm1,5)*pow(
         mm3,-2) + 1./105.*pow(mm1,6)*pow(mm3,-3) );

      h10px2 +=  + 2./105.*pow(mm1,-2)*pow(mm3,5) + 2./21.*pow(mm1,-1)*
         pow(mm3,4) + 8./35.*pow(mm3,3) + 12./35.*mm1*pow(mm3,2) + 12./
         35.*pow(mm1,2)*mm3 + 8./35.*pow(mm1,3) + 2./21.*pow(mm1,4)*
         pow(mm3,-1) + 2./105.*pow(mm1,5)*pow(mm3,-2);
      return h10px2*pi162/pow(mm2,5);
    case 3:
      h10px3 =
       + ln3 * ( 1./42.*pow(mm1,-2)*pow(mm3,5) + 19./210.*pow(mm1,-1)*
         pow(mm3,4) + 7./60.*pow(mm3,3) + 1./20.*mm1*pow(mm3,2) );

      h10px3 +=  + ln2 * (  - 1./42.*pow(mm1,-2)*pow(mm3,5) - 19./210.*
         pow(mm1,-1)*pow(mm3,4) - 7./60.*pow(mm3,3) - 1./20.*mm1*pow(
         mm3,2) + 1./60.*pow(mm1,4)*pow(mm3,-1) + 11./420.*pow(mm1,5)*
         pow(mm3,-2) + 1./105.*pow(mm1,6)*pow(mm3,-3) );

      h10px3 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm3,-1) - 11./420.*
         pow(mm1,5)*pow(mm3,-2) - 1./105.*pow(mm1,6)*pow(mm3,-3) );

      h10px3 +=  + 1./21.*pow(mm1,-1)*pow(mm3,4) + 101./420.*pow(mm3,3)
          + 13./28.*mm1*pow(mm3,2) + 11./28.*pow(mm1,2)*mm3 + 41./420.*
         pow(mm1,3) - 3./70.*pow(mm1,4)*pow(mm3,-1) - 2./105.*pow(
         mm1,5)*pow(mm3,-2);
      return h10px3*pi162/pow(mm2,5);
    case 4:
      h10px4 =
       + ln3 * (  - 1./42.*pow(mm1,-2)*pow(mm3,5) - 31./210.*pow(
         mm1,-1)*pow(mm3,4) - 157./420.*pow(mm3,3) - 29./60.*mm1*pow(
         mm3,2) - 19./60.*pow(mm1,2)*mm3 - 1./12.*pow(mm1,3) );

      h10px4 +=  + ln2 * ( 1./42.*pow(mm1,-2)*pow(mm3,5) + 31./210.*
         pow(mm1,-1)*pow(mm3,4) + 157./420.*pow(mm3,3) + 29./60.*mm1*
         pow(mm3,2) + 19./60.*pow(mm1,2)*mm3 + 1./12.*pow(mm1,3) + 1./
         60.*pow(mm1,4)*pow(mm3,-1) + 17./420.*pow(mm1,5)*pow(mm3,-2)
          + 13./420.*pow(mm1,6)*pow(mm3,-3) + 1./140.*pow(mm1,7)*pow(
         mm3,-4) );

      h10px4 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm3,-1) - 17./420.*
         pow(mm1,5)*pow(mm3,-2) - 13./420.*pow(mm1,6)*pow(mm3,-3) - 1./
         140.*pow(mm1,7)*pow(mm3,-4) );

      h10px4 +=  - 1./21.*pow(mm1,-1)*pow(mm3,4) - 79./420.*pow(mm3,3)
          - 107./420.*mm1*pow(mm3,2) - 13./105.*pow(mm1,2)*mm3 - 1./42.
         *pow(mm1,3) - 23./420.*pow(mm1,4)*pow(mm3,-1) - 23./420.*pow(
         mm1,5)*pow(mm3,-2) - 1./70.*pow(mm1,6)*pow(mm3,-3);
      return h10px4*pi162/pow(mm2,5);
    case 5:
      h10px5 =
       + ln3 * (  - 1./126.*pow(mm1,-4)*pow(mm3,7) - 31./630.*pow(
         mm1,-3)*pow(mm3,6) - 8./63.*pow(mm1,-2)*pow(mm3,5) - 6./35.*
         pow(mm1,-1)*pow(mm3,4) - 5./42.*pow(mm3,3) - 1./30.*mm1*pow(
         mm3,2) );

      h10px5 +=  + ln2 * ( 1./126.*pow(mm1,-4)*pow(mm3,7) + 31./630.*
         pow(mm1,-3)*pow(mm3,6) + 8./63.*pow(mm1,-2)*pow(mm3,5) + 6./35.
         *pow(mm1,-1)*pow(mm3,4) + 5./42.*pow(mm3,3) + 1./30.*mm1*pow(
         mm3,2) + 1./30.*pow(mm1,2)*mm3 + 5./42.*pow(mm1,3) + 6./35.*
         pow(mm1,4)*pow(mm3,-1) + 8./63.*pow(mm1,5)*pow(mm3,-2) + 31./
         630.*pow(mm1,6)*pow(mm3,-3) + 1./126.*pow(mm1,7)*pow(mm3,-4) )
         ;

      h10px5 +=  + ln1 * (  - 1./30.*pow(mm1,2)*mm3 - 5./42.*pow(mm1,3)
          - 6./35.*pow(mm1,4)*pow(mm3,-1) - 8./63.*pow(mm1,5)*pow(
         mm3,-2) - 31./630.*pow(mm1,6)*pow(mm3,-3) - 1./126.*pow(mm1,7)
         *pow(mm3,-4) );

      h10px5 +=  - 1./63.*pow(mm1,-3)*pow(mm3,6) - 19./210.*pow(mm1,-2)
         *pow(mm3,5) - 397./1890.*pow(mm1,-1)*pow(mm3,4) - 103./378.*
         pow(mm3,3) - 487./1890.*mm1*pow(mm3,2) - 487./1890.*pow(mm1,2)
         *mm3 - 103./378.*pow(mm1,3) - 397./1890.*pow(mm1,4)*pow(
         mm3,-1) - 19./210.*pow(mm1,5)*pow(mm3,-2) - 1./63.*pow(mm1,6)*
         pow(mm3,-3);
      return h10px5*pi162/pow(mm2,7);
    case 6:
      h10px6 =
       + ln3 * ( 1./126.*pow(mm1,-4)*pow(mm3,7) + 13./210.*pow(mm1,-3)*
         pow(mm3,6) + 22./105.*pow(mm1,-2)*pow(mm3,5) + 127./315.*pow(
         mm1,-1)*pow(mm3,4) + 17./35.*pow(mm3,3) + 13./35.*mm1*pow(
         mm3,2) + 1./6.*pow(mm1,2)*mm3 + 1./30.*pow(mm1,3) );

      h10px6 +=  + ln2 * (  - 1./126.*pow(mm1,-4)*pow(mm3,7) - 13./210.
         *pow(mm1,-3)*pow(mm3,6) - 22./105.*pow(mm1,-2)*pow(mm3,5) - 
         127./315.*pow(mm1,-1)*pow(mm3,4) - 17./35.*pow(mm3,3) - 13./35.
         *mm1*pow(mm3,2) - 2./15.*pow(mm1,2)*mm3 + 4./35.*pow(mm1,3) + 
         19./70.*pow(mm1,4)*pow(mm3,-1) + 169./630.*pow(mm1,5)*pow(
         mm3,-2) + 16./105.*pow(mm1,6)*pow(mm3,-3) + 1./21.*pow(mm1,7)*
         pow(mm3,-4) + 2./315.*pow(mm1,8)*pow(mm3,-5) );

      h10px6 +=  + ln1 * (  - 1./30.*pow(mm1,2)*mm3 - 31./210.*pow(
         mm1,3) - 19./70.*pow(mm1,4)*pow(mm3,-1) - 169./630.*pow(mm1,5)
         *pow(mm3,-2) - 16./105.*pow(mm1,6)*pow(mm3,-3) - 1./21.*pow(
         mm1,7)*pow(mm3,-4) - 2./315.*pow(mm1,8)*pow(mm3,-5) );

      h10px6 +=  + 1./63.*pow(mm1,-3)*pow(mm3,6) + 73./630.*pow(mm1,-2)
         *pow(mm3,5) + 137./378.*pow(mm1,-1)*pow(mm3,4) + 191./315.*
         pow(mm3,3) + 167./315.*mm1*pow(mm3,2) + 191./1890.*pow(mm1,2)*
         mm3 - 197./630.*pow(mm1,3) - 26./63.*pow(mm1,4)*pow(mm3,-1) - 
         247./945.*pow(mm1,5)*pow(mm3,-2) - 4./45.*pow(mm1,6)*pow(
         mm3,-3) - 4./315.*pow(mm1,7)*pow(mm3,-4);
      return h10px6*pi162/pow(mm2,7);
    case 7:
      h10px7 =
       + ln3 * ( 1./63.*pow(mm1,-3)*pow(mm3,6) + 13./126.*pow(mm1,-2)*
         pow(mm3,5) + 181./630.*pow(mm1,-1)*pow(mm3,4) + 187./420.*pow(
         mm3,3) + 173./420.*mm1*pow(mm3,2) + 13./60.*pow(mm1,2)*mm3 + 1.
         /20.*pow(mm1,3) );

      h10px7 +=  + ln2 * (  - 1./63.*pow(mm1,-3)*pow(mm3,6) - 13./126.*
         pow(mm1,-2)*pow(mm3,5) - 181./630.*pow(mm1,-1)*pow(mm3,4) - 
         187./420.*pow(mm3,3) - 173./420.*mm1*pow(mm3,2) - 13./60.*pow(
         mm1,2)*mm3 - 1./20.*pow(mm1,3) - 3./140.*pow(mm1,4)*pow(
         mm3,-1) - 9./140.*pow(mm1,5)*pow(mm3,-2) - 89./1260.*pow(
         mm1,6)*pow(mm3,-3) - 43./1260.*pow(mm1,7)*pow(mm3,-4) - 2./315.
         *pow(mm1,8)*pow(mm3,-5) );

      h10px7 +=  + ln1 * ( 3./140.*pow(mm1,4)*pow(mm3,-1) + 9./140.*
         pow(mm1,5)*pow(mm3,-2) + 89./1260.*pow(mm1,6)*pow(mm3,-3) + 43.
         /1260.*pow(mm1,7)*pow(mm3,-4) + 2./315.*pow(mm1,8)*pow(mm3,-5)
          );

      h10px7 +=  + 2./63.*pow(mm1,-2)*pow(mm3,5) + 4./21.*pow(mm1,-1)*
         pow(mm3,4) + 1717./3780.*pow(mm3,3) + 397./756.*mm1*pow(mm3,2)
          + 53./189.*pow(mm1,2)*mm3 + 25./378.*pow(mm1,3) + 293./3780.*
         pow(mm1,4)*pow(mm3,-1) + 421./3780.*pow(mm1,5)*pow(mm3,-2) + 
         13./210.*pow(mm1,6)*pow(mm3,-3) + 4./315.*pow(mm1,7)*pow(
         mm3,-4);
      return h10px7*pi162/pow(mm2,7);
    case 8:
      h10px8 =
       + ln3 * (  - 2./231.*pow(mm1,-5)*pow(mm3,8) - 17./231.*pow(
         mm1,-4)*pow(mm3,7) - 137./495.*pow(mm1,-3)*pow(mm3,6) - 93./
         154.*pow(mm1,-2)*pow(mm3,5) - 59./70.*pow(mm1,-1)*pow(mm3,4)
          - 247./315.*pow(mm3,3) - 17./35.*mm1*pow(mm3,2) - 13./70.*
         pow(mm1,2)*mm3 - 1./30.*pow(mm1,3) );

      h10px8 +=  + ln2 * ( 2./231.*pow(mm1,-5)*pow(mm3,8) + 17./231.*
         pow(mm1,-4)*pow(mm3,7) + 137./495.*pow(mm1,-3)*pow(mm3,6) + 93.
         /154.*pow(mm1,-2)*pow(mm3,5) + 59./70.*pow(mm1,-1)*pow(mm3,4)
          + 247./315.*pow(mm3,3) + 17./35.*mm1*pow(mm3,2) + 1./7.*pow(
         mm1,2)*mm3 - 19./105.*pow(mm1,3) - 7./15.*pow(mm1,4)*pow(
         mm3,-1) - 61./105.*pow(mm1,5)*pow(mm3,-2) - 3119./6930.*pow(
         mm1,6)*pow(mm3,-3) - 167./770.*pow(mm1,7)*pow(mm3,-4) - 23./
         385.*pow(mm1,8)*pow(mm3,-5) - 5./693.*pow(mm1,9)*pow(mm3,-6) )
         ;

      h10px8 +=  + ln1 * ( 3./70.*pow(mm1,2)*mm3 + 3./14.*pow(mm1,3) + 
         7./15.*pow(mm1,4)*pow(mm3,-1) + 61./105.*pow(mm1,5)*pow(
         mm3,-2) + 3119./6930.*pow(mm1,6)*pow(mm3,-3) + 167./770.*pow(
         mm1,7)*pow(mm3,-4) + 23./385.*pow(mm1,8)*pow(mm3,-5) + 5./693.
         *pow(mm1,9)*pow(mm3,-6) );

      h10px8 +=  - 4./231.*pow(mm1,-4)*pow(mm3,7) - 32./231.*pow(
         mm1,-3)*pow(mm3,6) - 17./35.*pow(mm1,-2)*pow(mm3,5) - 161./165.
         *pow(mm1,-1)*pow(mm3,4) - 1798./1485.*pow(mm3,3) - 437./495.*
         mm1*pow(mm3,2) - 83./495.*pow(mm1,2)*mm3 + 5266./10395.*pow(
         mm1,3) + 958./1155.*pow(mm1,4)*pow(mm3,-1) + 277./385.*pow(
         mm1,5)*pow(mm3,-2) + 358./945.*pow(mm1,6)*pow(mm3,-3) + 389./
         3465.*pow(mm1,7)*pow(mm3,-4) + 10./693.*pow(mm1,8)*pow(mm3,-5)
         ;
      return h10px8*pi162/pow(mm2,9);
    default:
      std::cout << "h10psing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h10px1 =
       + ln3 * ( 1./15.*pow(mm1,-1)*pow(mm2,4) + 19./60.*pow(mm2,3) + 7.
         /12.*mm1*pow(mm2,2) + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3)
          - 1./60.*pow(mm1,4)*pow(mm2,-1) - 1./60.*pow(mm1,5)*pow(
         mm2,-2) );

      h10px1 +=  + ln2 * (  - 1./15.*pow(mm1,-1)*pow(mm2,4) - 3./20.*
         pow(mm2,3) - 1./12.*mm1*pow(mm2,2) );

      h10px1 +=  + ln1 * ( 1./60.*pow(mm1,4)*pow(mm2,-1) + 1./60.*pow(
         mm1,5)*pow(mm2,-2) );

      h10px1 +=  + 47./360.*pow(mm2,3) + 47./120.*mm1*pow(mm2,2) + 17./
         40.*pow(mm1,2)*mm2 + 71./360.*pow(mm1,3) + 1./30.*pow(mm1,4)*
         pow(mm2,-1);
      return h10px1*pi162/pow(mm3,3);
    case 2:
      h10px2 =
       + ln3 * (  - 1./105.*pow(mm1,-3)*pow(mm2,6) - 11./210.*pow(
         mm1,-2)*pow(mm2,5) - 23./210.*pow(mm1,-1)*pow(mm2,4) - 1./10.*
         pow(mm2,3) - 1./30.*mm1*pow(mm2,2) - 1./30.*pow(mm1,2)*mm2 - 1.
         /10.*pow(mm1,3) - 23./210.*pow(mm1,4)*pow(mm2,-1) - 11./210.*
         pow(mm1,5)*pow(mm2,-2) - 1./105.*pow(mm1,6)*pow(mm2,-3) );

      h10px2 +=  + ln2 * ( 1./105.*pow(mm1,-3)*pow(mm2,6) + 11./210.*
         pow(mm1,-2)*pow(mm2,5) + 23./210.*pow(mm1,-1)*pow(mm2,4) + 1./
         10.*pow(mm2,3) + 1./30.*mm1*pow(mm2,2) );

      h10px2 +=  + ln1 * ( 1./30.*pow(mm1,2)*mm2 + 1./10.*pow(mm1,3) + 
         23./210.*pow(mm1,4)*pow(mm2,-1) + 11./210.*pow(mm1,5)*pow(
         mm2,-2) + 1./105.*pow(mm1,6)*pow(mm2,-3) );

      h10px2 +=  + 2./105.*pow(mm1,-2)*pow(mm2,5) + 2./21.*pow(mm1,-1)*
         pow(mm2,4) + 8./35.*pow(mm2,3) + 12./35.*mm1*pow(mm2,2) + 12./
         35.*pow(mm1,2)*mm2 + 8./35.*pow(mm1,3) + 2./21.*pow(mm1,4)*
         pow(mm2,-1) + 2./105.*pow(mm1,5)*pow(mm2,-2);
      return h10px2*pi162/pow(mm3,5);
    case 3:
      h10px3 =
       + ln3 * ( 1./42.*pow(mm1,-2)*pow(mm2,5) + 31./210.*pow(mm1,-1)*
         pow(mm2,4) + 157./420.*pow(mm2,3) + 29./60.*mm1*pow(mm2,2) + 
         19./60.*pow(mm1,2)*mm2 + 1./12.*pow(mm1,3) + 1./60.*pow(mm1,4)
         *pow(mm2,-1) + 17./420.*pow(mm1,5)*pow(mm2,-2) + 13./420.*pow(
         mm1,6)*pow(mm2,-3) + 1./140.*pow(mm1,7)*pow(mm2,-4) );

      h10px3 +=  + ln2 * (  - 1./42.*pow(mm1,-2)*pow(mm2,5) - 31./210.*
         pow(mm1,-1)*pow(mm2,4) - 157./420.*pow(mm2,3) - 29./60.*mm1*
         pow(mm2,2) - 19./60.*pow(mm1,2)*mm2 - 1./12.*pow(mm1,3) );

      h10px3 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm2,-1) - 17./420.*
         pow(mm1,5)*pow(mm2,-2) - 13./420.*pow(mm1,6)*pow(mm2,-3) - 1./
         140.*pow(mm1,7)*pow(mm2,-4) );

      h10px3 +=  - 1./21.*pow(mm1,-1)*pow(mm2,4) - 79./420.*pow(mm2,3)
          - 107./420.*mm1*pow(mm2,2) - 13./105.*pow(mm1,2)*mm2 - 1./42.
         *pow(mm1,3) - 23./420.*pow(mm1,4)*pow(mm2,-1) - 23./420.*pow(
         mm1,5)*pow(mm2,-2) - 1./70.*pow(mm1,6)*pow(mm2,-3);
      return h10px3*pi162/pow(mm3,5);
    case 4:
      h10px4 =
       + ln3 * (  - 1./42.*pow(mm1,-2)*pow(mm2,5) - 19./210.*pow(
         mm1,-1)*pow(mm2,4) - 7./60.*pow(mm2,3) - 1./20.*mm1*pow(mm2,2)
          + 1./60.*pow(mm1,4)*pow(mm2,-1) + 11./420.*pow(mm1,5)*pow(
         mm2,-2) + 1./105.*pow(mm1,6)*pow(mm2,-3) );

      h10px4 +=  + ln2 * ( 1./42.*pow(mm1,-2)*pow(mm2,5) + 19./210.*
         pow(mm1,-1)*pow(mm2,4) + 7./60.*pow(mm2,3) + 1./20.*mm1*pow(
         mm2,2) );

      h10px4 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm2,-1) - 11./420.*
         pow(mm1,5)*pow(mm2,-2) - 1./105.*pow(mm1,6)*pow(mm2,-3) );

      h10px4 +=  + 1./21.*pow(mm1,-1)*pow(mm2,4) + 101./420.*pow(mm2,3)
          + 13./28.*mm1*pow(mm2,2) + 11./28.*pow(mm1,2)*mm2 + 41./420.*
         pow(mm1,3) - 3./70.*pow(mm1,4)*pow(mm2,-1) - 2./105.*pow(
         mm1,5)*pow(mm2,-2);
      return h10px4*pi162/pow(mm3,5);
    case 5:
      h10px5 =
       + ln3 * (  - 1./126.*pow(mm1,-4)*pow(mm2,7) - 13./210.*pow(
         mm1,-3)*pow(mm2,6) - 22./105.*pow(mm1,-2)*pow(mm2,5) - 127./
         315.*pow(mm1,-1)*pow(mm2,4) - 17./35.*pow(mm2,3) - 13./35.*mm1
         *pow(mm2,2) - 2./15.*pow(mm1,2)*mm2 + 4./35.*pow(mm1,3) + 19./
         70.*pow(mm1,4)*pow(mm2,-1) + 169./630.*pow(mm1,5)*pow(mm2,-2)
          + 16./105.*pow(mm1,6)*pow(mm2,-3) + 1./21.*pow(mm1,7)*pow(
         mm2,-4) + 2./315.*pow(mm1,8)*pow(mm2,-5) );

      h10px5 +=  + ln2 * ( 1./126.*pow(mm1,-4)*pow(mm2,7) + 13./210.*
         pow(mm1,-3)*pow(mm2,6) + 22./105.*pow(mm1,-2)*pow(mm2,5) + 127.
         /315.*pow(mm1,-1)*pow(mm2,4) + 17./35.*pow(mm2,3) + 13./35.*
         mm1*pow(mm2,2) + 1./6.*pow(mm1,2)*mm2 + 1./30.*pow(mm1,3) );

      h10px5 +=  + ln1 * (  - 1./30.*pow(mm1,2)*mm2 - 31./210.*pow(
         mm1,3) - 19./70.*pow(mm1,4)*pow(mm2,-1) - 169./630.*pow(mm1,5)
         *pow(mm2,-2) - 16./105.*pow(mm1,6)*pow(mm2,-3) - 1./21.*pow(
         mm1,7)*pow(mm2,-4) - 2./315.*pow(mm1,8)*pow(mm2,-5) );

      h10px5 +=  + 1./63.*pow(mm1,-3)*pow(mm2,6) + 73./630.*pow(mm1,-2)
         *pow(mm2,5) + 137./378.*pow(mm1,-1)*pow(mm2,4) + 191./315.*
         pow(mm2,3) + 167./315.*mm1*pow(mm2,2) + 191./1890.*pow(mm1,2)*
         mm2 - 197./630.*pow(mm1,3) - 26./63.*pow(mm1,4)*pow(mm2,-1) - 
         247./945.*pow(mm1,5)*pow(mm2,-2) - 4./45.*pow(mm1,6)*pow(
         mm2,-3) - 4./315.*pow(mm1,7)*pow(mm2,-4);
      return h10px5*pi162/pow(mm3,7);
    case 6:
      h10px6 =
       + ln3 * ( 1./126.*pow(mm1,-4)*pow(mm2,7) + 31./630.*pow(mm1,-3)*
         pow(mm2,6) + 8./63.*pow(mm1,-2)*pow(mm2,5) + 6./35.*pow(
         mm1,-1)*pow(mm2,4) + 5./42.*pow(mm2,3) + 1./30.*mm1*pow(mm2,2)
          + 1./30.*pow(mm1,2)*mm2 + 5./42.*pow(mm1,3) + 6./35.*pow(
         mm1,4)*pow(mm2,-1) + 8./63.*pow(mm1,5)*pow(mm2,-2) + 31./630.*
         pow(mm1,6)*pow(mm2,-3) + 1./126.*pow(mm1,7)*pow(mm2,-4) );

      h10px6 +=  + ln2 * (  - 1./126.*pow(mm1,-4)*pow(mm2,7) - 31./630.
         *pow(mm1,-3)*pow(mm2,6) - 8./63.*pow(mm1,-2)*pow(mm2,5) - 6./
         35.*pow(mm1,-1)*pow(mm2,4) - 5./42.*pow(mm2,3) - 1./30.*mm1*
         pow(mm2,2) );

      h10px6 +=  + ln1 * (  - 1./30.*pow(mm1,2)*mm2 - 5./42.*pow(mm1,3)
          - 6./35.*pow(mm1,4)*pow(mm2,-1) - 8./63.*pow(mm1,5)*pow(
         mm2,-2) - 31./630.*pow(mm1,6)*pow(mm2,-3) - 1./126.*pow(mm1,7)
         *pow(mm2,-4) );

      h10px6 +=  - 1./63.*pow(mm1,-3)*pow(mm2,6) - 19./210.*pow(mm1,-2)
         *pow(mm2,5) - 397./1890.*pow(mm1,-1)*pow(mm2,4) - 103./378.*
         pow(mm2,3) - 487./1890.*mm1*pow(mm2,2) - 487./1890.*pow(mm1,2)
         *mm2 - 103./378.*pow(mm1,3) - 397./1890.*pow(mm1,4)*pow(
         mm2,-1) - 19./210.*pow(mm1,5)*pow(mm2,-2) - 1./63.*pow(mm1,6)*
         pow(mm2,-3);
      return h10px6*pi162/pow(mm3,7);
    case 7:
      h10px7 =
       + ln3 * (  - 1./63.*pow(mm1,-3)*pow(mm2,6) - 13./126.*pow(
         mm1,-2)*pow(mm2,5) - 181./630.*pow(mm1,-1)*pow(mm2,4) - 187./
         420.*pow(mm2,3) - 173./420.*mm1*pow(mm2,2) - 13./60.*pow(
         mm1,2)*mm2 - 1./20.*pow(mm1,3) - 3./140.*pow(mm1,4)*pow(
         mm2,-1) - 9./140.*pow(mm1,5)*pow(mm2,-2) - 89./1260.*pow(
         mm1,6)*pow(mm2,-3) - 43./1260.*pow(mm1,7)*pow(mm2,-4) - 2./315.
         *pow(mm1,8)*pow(mm2,-5) );

      h10px7 +=  + ln2 * ( 1./63.*pow(mm1,-3)*pow(mm2,6) + 13./126.*
         pow(mm1,-2)*pow(mm2,5) + 181./630.*pow(mm1,-1)*pow(mm2,4) + 
         187./420.*pow(mm2,3) + 173./420.*mm1*pow(mm2,2) + 13./60.*pow(
         mm1,2)*mm2 + 1./20.*pow(mm1,3) );

      h10px7 +=  + ln1 * ( 3./140.*pow(mm1,4)*pow(mm2,-1) + 9./140.*
         pow(mm1,5)*pow(mm2,-2) + 89./1260.*pow(mm1,6)*pow(mm2,-3) + 43.
         /1260.*pow(mm1,7)*pow(mm2,-4) + 2./315.*pow(mm1,8)*pow(mm2,-5)
          );

      h10px7 +=  + 2./63.*pow(mm1,-2)*pow(mm2,5) + 4./21.*pow(mm1,-1)*
         pow(mm2,4) + 1717./3780.*pow(mm2,3) + 397./756.*mm1*pow(mm2,2)
          + 53./189.*pow(mm1,2)*mm2 + 25./378.*pow(mm1,3) + 293./3780.*
         pow(mm1,4)*pow(mm2,-1) + 421./3780.*pow(mm1,5)*pow(mm2,-2) + 
         13./210.*pow(mm1,6)*pow(mm2,-3) + 4./315.*pow(mm1,7)*pow(
         mm2,-4);
      return h10px7*pi162/pow(mm3,7);
    case 8:
      h10px8 =
       + ln3 * ( 2./231.*pow(mm1,-5)*pow(mm2,8) + 17./231.*pow(mm1,-4)*
         pow(mm2,7) + 137./495.*pow(mm1,-3)*pow(mm2,6) + 93./154.*pow(
         mm1,-2)*pow(mm2,5) + 59./70.*pow(mm1,-1)*pow(mm2,4) + 247./315.
         *pow(mm2,3) + 17./35.*mm1*pow(mm2,2) + 1./7.*pow(mm1,2)*mm2 - 
         19./105.*pow(mm1,3) - 7./15.*pow(mm1,4)*pow(mm2,-1) - 61./105.
         *pow(mm1,5)*pow(mm2,-2) - 3119./6930.*pow(mm1,6)*pow(mm2,-3)
          - 167./770.*pow(mm1,7)*pow(mm2,-4) - 23./385.*pow(mm1,8)*pow(
         mm2,-5) - 5./693.*pow(mm1,9)*pow(mm2,-6) );

      h10px8 +=  + ln2 * (  - 2./231.*pow(mm1,-5)*pow(mm2,8) - 17./231.
         *pow(mm1,-4)*pow(mm2,7) - 137./495.*pow(mm1,-3)*pow(mm2,6) - 
         93./154.*pow(mm1,-2)*pow(mm2,5) - 59./70.*pow(mm1,-1)*pow(
         mm2,4) - 247./315.*pow(mm2,3) - 17./35.*mm1*pow(mm2,2) - 13./
         70.*pow(mm1,2)*mm2 - 1./30.*pow(mm1,3) );

      h10px8 +=  + ln1 * ( 3./70.*pow(mm1,2)*mm2 + 3./14.*pow(mm1,3) + 
         7./15.*pow(mm1,4)*pow(mm2,-1) + 61./105.*pow(mm1,5)*pow(
         mm2,-2) + 3119./6930.*pow(mm1,6)*pow(mm2,-3) + 167./770.*pow(
         mm1,7)*pow(mm2,-4) + 23./385.*pow(mm1,8)*pow(mm2,-5) + 5./693.
         *pow(mm1,9)*pow(mm2,-6) );

      h10px8 +=  - 4./231.*pow(mm1,-4)*pow(mm2,7) - 32./231.*pow(
         mm1,-3)*pow(mm2,6) - 17./35.*pow(mm1,-2)*pow(mm2,5) - 161./165.
         *pow(mm1,-1)*pow(mm2,4) - 1798./1485.*pow(mm2,3) - 437./495.*
         mm1*pow(mm2,2) - 83./495.*pow(mm1,2)*mm2 + 5266./10395.*pow(
         mm1,3) + 958./1155.*pow(mm1,4)*pow(mm2,-1) + 277./385.*pow(
         mm1,5)*pow(mm2,-2) + 358./945.*pow(mm1,6)*pow(mm2,-3) + 389./
         3465.*pow(mm1,7)*pow(mm2,-4) + 10./693.*pow(mm1,8)*pow(mm2,-5)
         ;
      return h10px8*pi162/pow(mm3,9);
    default:
      std::cout << "h10psing, wrong iprop = "<<iprop<<'\n';
      return 0.;
     }
  }
  std::cout << "funny masses in h10psing\n";
  return 0.;
}

double h210sing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h210x1,h210x2,h210x3,h210x4,h210x5,h210x6,h210x7,h210x8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210x1 =
       + ln3 * (  - 1./18.*pow(mm2,-1)*pow(mm3,4) + 1./9.*pow(mm3,3) + 
         1./2.*mm2*pow(mm3,2) + 1./3.*pow(mm2,2)*mm3 );

      h210x1 +=  + pow(ln3,2) * ( 1./6.*pow(mm3,3) + 1./6.*mm2*pow(
         mm3,2) );

      h210x1 +=  + ln2 * ( 1./3.*mm2*pow(mm3,2) + 1./2.*pow(mm2,2)*mm3
          + 1./9.*pow(mm2,3) - 1./18.*pow(mm2,4)*pow(mm3,-1) );

      h210x1 +=  + ln2*ln3 * ( 1./6.*pow(mm3,3) + 1./6.*mm2*pow(mm3,2)
          + 1./6.*pow(mm2,2)*mm3 + 1./6.*pow(mm2,3) );

      h210x1 +=  + pow(ln2,2) * ( 1./6.*pow(mm2,2)*mm3 + 1./6.*pow(
         mm2,3) );

      h210x1 +=  + ln1 * ( 1./18.*pow(mm2,-1)*pow(mm3,4) - 7./18.*pow(
         mm3,3) - 13./9.*mm2*pow(mm3,2) - 13./9.*pow(mm2,2)*mm3 - 7./18.
         *pow(mm2,3) + 1./18.*pow(mm2,4)*pow(mm3,-1) );

      h210x1 +=  + ln1*ln3 * ( 1./6.*pow(mm3,3) + 1./6.*mm2*pow(mm3,2)
          - 1./6.*pow(mm2,2)*mm3 - 1./6.*pow(mm2,3) );

      h210x1 +=  + ln1*ln2 * (  - 1./6.*pow(mm3,3) - 1./6.*mm2*pow(
         mm3,2) + 1./6.*pow(mm2,2)*mm3 + 1./6.*pow(mm2,3) );

      h210x1 +=  + 139./216.*pow(mm3,3) + 1./36.*pow(mm3,3)*pi2 + 265./
         216.*mm2*pow(mm3,2) + 1./36.*mm2*pow(mm3,2)*pi2 + 265./216.*
         pow(mm2,2)*mm3 + 1./36.*pow(mm2,2)*mm3*pi2 + 139./216.*pow(
         mm2,3) + 1./36.*pow(mm2,3)*pi2;
      return h210x1*pi162/mm1;
    case 2:
      h210x2 =
       + ln3 * ( 1./60.*pow(mm2,-2)*pow(mm3,5) + 7./60.*pow(mm2,-1)*
         pow(mm3,4) + 11./60.*pow(mm3,3) + 1./12.*mm2*pow(mm3,2) );

      h210x2 +=  + ln2 * ( 1./12.*pow(mm2,2)*mm3 + 11./60.*pow(mm2,3)
          + 7./60.*pow(mm2,4)*pow(mm3,-1) + 1./60.*pow(mm2,5)*pow(
         mm3,-2) );

      h210x2 +=  + ln1 * (  - 1./60.*pow(mm2,-2)*pow(mm3,5) - 7./60.*
         pow(mm2,-1)*pow(mm3,4) - 7./20.*pow(mm3,3) - 7./12.*mm2*pow(
         mm3,2) - 7./12.*pow(mm2,2)*mm3 - 7./20.*pow(mm2,3) - 7./60.*
         pow(mm2,4)*pow(mm3,-1) - 1./60.*pow(mm2,5)*pow(mm3,-2) );

      h210x2 +=  + 1./30.*pow(mm2,-1)*pow(mm3,4) + 13./360.*pow(mm3,3)
          - 7./120.*mm2*pow(mm3,2) - 7./120.*pow(mm2,2)*mm3 + 13./360.*
         pow(mm2,3) + 1./30.*pow(mm2,4)*pow(mm3,-1);
      return h210x2*pi162/pow(mm1,3);
    case 3:
      h210x3 =
       + ln3 * ( 1./90.*pow(mm2,-3)*pow(mm3,6) + 7./60.*pow(mm2,-2)*
         pow(mm3,5) + 37./60.*pow(mm2,-1)*pow(mm3,4) + 227./180.*pow(
         mm3,3) + 13./12.*mm2*pow(mm3,2) + 1./3.*pow(mm2,2)*mm3 );

      h210x3 +=  + ln2 * ( 2./3.*pow(mm3,3) + 11./6.*mm2*pow(mm3,2) + 
         89./60.*pow(mm2,2)*mm3 + 7./60.*pow(mm2,3) - 13./60.*pow(
         mm2,4)*pow(mm3,-1) - 1./60.*pow(mm2,5)*pow(mm3,-2) );

      h210x3 +=  + ln2*ln3 * ( 1./6.*pow(mm3,3) + 1./2.*mm2*pow(mm3,2)
          + 1./2.*pow(mm2,2)*mm3 + 1./6.*pow(mm2,3) );

      h210x3 +=  + pow(ln2,2) * ( 1./6.*pow(mm3,3) + 1./2.*mm2*pow(
         mm3,2) + 1./2.*pow(mm2,2)*mm3 + 1./6.*pow(mm2,3) );

      h210x3 +=  + ln1 * (  - 1./90.*pow(mm2,-3)*pow(mm3,6) - 7./60.*
         pow(mm2,-2)*pow(mm3,5) - 37./60.*pow(mm2,-1)*pow(mm3,4) - 247./
         180.*pow(mm3,3) - 5./4.*mm2*pow(mm3,2) - 3./20.*pow(mm2,2)*mm3
          + 79./180.*pow(mm2,3) + 13./60.*pow(mm2,4)*pow(mm3,-1) + 1./
         60.*pow(mm2,5)*pow(mm3,-2) );

      h210x3 +=  + ln1*ln3 * (  - 1./6.*pow(mm3,3) - 1./2.*mm2*pow(
         mm3,2) - 1./2.*pow(mm2,2)*mm3 - 1./6.*pow(mm2,3) );

      h210x3 +=  + ln1*ln2 * ( 1./6.*pow(mm3,3) + 1./2.*mm2*pow(mm3,2)
          + 1./2.*pow(mm2,2)*mm3 + 1./6.*pow(mm2,3) );

      h210x3 +=  + 1./45.*pow(mm2,-2)*pow(mm3,5) + 2./9.*pow(mm2,-1)*
         pow(mm3,4) + 118./135.*pow(mm3,3) + 1./36.*pow(mm3,3)*pi2 + 
         133./90.*mm2*pow(mm3,2) + 1./12.*mm2*pow(mm3,2)*pi2 + 97./90.*
         pow(mm2,2)*mm3 + 1./12.*pow(mm2,2)*mm3*pi2 + 13./54.*pow(
         mm2,3) + 1./36.*pow(mm2,3)*pi2 - 1./30.*pow(mm2,4)*pow(mm3,-1)
         ;
      return h210x3*pi162/pow(mm1,3);
    case 4:
      h210x4 =
       + ln3 * (  - 1./60.*pow(mm2,-2)*pow(mm3,7) - 1./4.*pow(mm2,-1)*
         pow(mm3,6) - 1./3.*pow(mm3,5) + 3./2.*mm2*pow(mm3,4) + 59./12.
         *pow(mm2,2)*pow(mm3,3) + 349./60.*pow(mm2,3)*pow(mm3,2) + 19./
         6.*pow(mm2,4)*mm3 + 2./3.*pow(mm2,5) );

      h210x4 +=  + pow(ln3,2) * ( 1./6.*pow(mm3,5) + 5./6.*mm2*pow(
         mm3,4) + 5./3.*pow(mm2,2)*pow(mm3,3) + 5./3.*pow(mm2,3)*pow(
         mm3,2) + 5./6.*pow(mm2,4)*mm3 + 1./6.*pow(mm2,5) );

      h210x4 +=  + ln2 * ( 1./3.*mm2*pow(mm3,4) + 7./4.*pow(mm2,2)*pow(
         mm3,3) + 677./180.*pow(mm2,3)*pow(mm3,2) + 38./9.*pow(mm2,4)*
         mm3 + 47./18.*pow(mm2,5) + 31./36.*pow(mm2,6)*pow(mm3,-1) + 5./
         36.*pow(mm2,7)*pow(mm3,-2) + 1./90.*pow(mm2,8)*pow(mm3,-3) );

      h210x4 +=  + ln2*ln3 * ( 1./6.*pow(mm3,5) + 5./6.*mm2*pow(mm3,4)
          + 5./3.*pow(mm2,2)*pow(mm3,3) + 5./3.*pow(mm2,3)*pow(mm3,2)
          + 5./6.*pow(mm2,4)*mm3 + 1./6.*pow(mm2,5) );

      h210x4 +=  + ln1 * ( 1./60.*pow(mm2,-2)*pow(mm3,7) + 1./4.*pow(
         mm2,-1)*pow(mm3,6) + 8./9.*pow(mm3,5) + 17./18.*mm2*pow(mm3,4)
          - 10./9.*pow(mm2,2)*pow(mm3,3) - 181./45.*pow(mm2,3)*pow(
         mm3,2) - 83./18.*pow(mm2,4)*mm3 - 49./18.*pow(mm2,5) - 31./36.
         *pow(mm2,6)*pow(mm3,-1) - 5./36.*pow(mm2,7)*pow(mm3,-2) - 1./
         90.*pow(mm2,8)*pow(mm3,-3) );

      h210x4 +=  + ln1*ln3 * ( 1./6.*pow(mm3,5) + 5./6.*mm2*pow(mm3,4)
          + 5./3.*pow(mm2,2)*pow(mm3,3) + 5./3.*pow(mm2,3)*pow(mm3,2)
          + 5./6.*pow(mm2,4)*mm3 + 1./6.*pow(mm2,5) );

      h210x4 +=  + ln1*ln2 * (  - 1./6.*pow(mm3,5) - 5./6.*mm2*pow(
         mm3,4) - 5./3.*pow(mm2,2)*pow(mm3,3) - 5./3.*pow(mm2,3)*pow(
         mm3,2) - 5./6.*pow(mm2,4)*mm3 - 1./6.*pow(mm2,5) );

      h210x4 +=  - 1./30.*pow(mm2,-1)*pow(mm3,6) + 47./270.*pow(mm3,5)
          + 1./36.*pow(mm3,5)*pi2 + 206./135.*mm2*pow(mm3,4) + 5./36.*
         mm2*pow(mm3,4)*pi2 + 523./135.*pow(mm2,2)*pow(mm3,3) + 5./18.*
         pow(mm2,2)*pow(mm3,3)*pi2 + 265./54.*pow(mm2,3)*pow(mm3,2) + 5.
         /18.*pow(mm2,3)*pow(mm3,2)*pi2 + 931./270.*pow(mm2,4)*mm3 + 5./
         36.*pow(mm2,4)*mm3*pi2 + 181./135.*pow(mm2,5) + 1./36.*pow(
         mm2,5)*pi2 + 4./15.*pow(mm2,6)*pow(mm3,-1) + 1./45.*pow(mm2,7)
         *pow(mm3,-2);
      return h210x4*pi162/pow(mm1,5);
    case 5:
      h210x5 =
       + ln3 * (  - 1./140.*pow(mm2,-4)*pow(mm3,8) - 13./210.*pow(
         mm2,-3)*pow(mm3,7) - 19./84.*pow(mm2,-2)*pow(mm3,6) - 3./7.*
         pow(mm2,-1)*pow(mm3,5) - 37./84.*pow(mm3,4) - 7./30.*mm2*pow(
         mm3,3) - 1./20.*pow(mm2,2)*pow(mm3,2) );

      h210x5 +=  + ln2 * ( 1./12.*pow(mm3,4) + 13./30.*mm2*pow(mm3,3)
          + 19./20.*pow(mm2,2)*pow(mm3,2) + 8./7.*pow(mm2,3)*mm3 + 23./
         28.*pow(mm2,4) + 5./14.*pow(mm2,5)*pow(mm3,-1) + 37./420.*pow(
         mm2,6)*pow(mm3,-2) + 1./105.*pow(mm2,7)*pow(mm3,-3) );

      h210x5 +=  + ln1 * ( 1./140.*pow(mm2,-4)*pow(mm3,8) + 13./210.*
         pow(mm2,-3)*pow(mm3,7) + 19./84.*pow(mm2,-2)*pow(mm3,6) + 3./7.
         *pow(mm2,-1)*pow(mm3,5) + 5./14.*pow(mm3,4) - 1./5.*mm2*pow(
         mm3,3) - 9./10.*pow(mm2,2)*pow(mm3,2) - 8./7.*pow(mm2,3)*mm3
          - 23./28.*pow(mm2,4) - 5./14.*pow(mm2,5)*pow(mm3,-1) - 37./
         420.*pow(mm2,6)*pow(mm3,-2) - 1./105.*pow(mm2,7)*pow(mm3,-3) )
         ;

      h210x5 +=  - 1./70.*pow(mm2,-3)*pow(mm3,7) - 7./60.*pow(mm2,-2)*
         pow(mm3,6) - 83./210.*pow(mm2,-1)*pow(mm3,5) - 55./84.*pow(
         mm3,4) - 43./105.*mm2*pow(mm3,3) + 43./140.*pow(mm2,2)*pow(
         mm3,2) + 51./70.*pow(mm2,3)*mm3 + 73./140.*pow(mm2,4) + 1./6.*
         pow(mm2,5)*pow(mm3,-1) + 2./105.*pow(mm2,6)*pow(mm3,-2);
      return h210x5*pi162/pow(mm1,6);
    case 6:
      h210x6 =
       + ln3 * ( 1./105.*pow(mm2,-3)*pow(mm3,7) + 37./420.*pow(mm2,-2)*
         pow(mm3,6) + 5./14.*pow(mm2,-1)*pow(mm3,5) + 23./28.*pow(
         mm3,4) + 8./7.*mm2*pow(mm3,3) + 19./20.*pow(mm2,2)*pow(mm3,2)
          + 13./30.*pow(mm2,3)*mm3 + 1./12.*pow(mm2,4) );

      h210x6 +=  + ln2 * (  - 1./20.*pow(mm2,2)*pow(mm3,2) - 7./30.*
         pow(mm2,3)*mm3 - 37./84.*pow(mm2,4) - 3./7.*pow(mm2,5)*pow(
         mm3,-1) - 19./84.*pow(mm2,6)*pow(mm3,-2) - 13./210.*pow(mm2,7)
         *pow(mm3,-3) - 1./140.*pow(mm2,8)*pow(mm3,-4) );

      h210x6 +=  + ln1 * (  - 1./105.*pow(mm2,-3)*pow(mm3,7) - 37./420.
         *pow(mm2,-2)*pow(mm3,6) - 5./14.*pow(mm2,-1)*pow(mm3,5) - 23./
         28.*pow(mm3,4) - 8./7.*mm2*pow(mm3,3) - 9./10.*pow(mm2,2)*pow(
         mm3,2) - 1./5.*pow(mm2,3)*mm3 + 5./14.*pow(mm2,4) + 3./7.*pow(
         mm2,5)*pow(mm3,-1) + 19./84.*pow(mm2,6)*pow(mm3,-2) + 13./210.
         *pow(mm2,7)*pow(mm3,-3) + 1./140.*pow(mm2,8)*pow(mm3,-4) );

      h210x6 +=  + 2./105.*pow(mm2,-2)*pow(mm3,6) + 1./6.*pow(mm2,-1)*
         pow(mm3,5) + 73./140.*pow(mm3,4) + 51./70.*mm2*pow(mm3,3) + 43.
         /140.*pow(mm2,2)*pow(mm3,2) - 43./105.*pow(mm2,3)*mm3 - 55./84.
         *pow(mm2,4) - 83./210.*pow(mm2,5)*pow(mm3,-1) - 7./60.*pow(
         mm2,6)*pow(mm3,-2) - 1./70.*pow(mm2,7)*pow(mm3,-3);
      return h210x6*pi162/pow(mm1,6);
    case 7:
      h210x7 =
       + ln3 * ( 1./140.*pow(mm2,-4)*pow(mm3,8) + 3./35.*pow(mm2,-3)*
         pow(mm3,7) + 33./70.*pow(mm2,-2)*pow(mm3,6) + 11./7.*pow(
         mm2,-1)*pow(mm3,5) + 23./7.*pow(mm3,4) + 149./35.*mm2*pow(
         mm3,3) + 33./10.*pow(mm2,2)*pow(mm3,2) + 7./5.*pow(mm2,3)*mm3
          + 1./4.*pow(mm2,4) );

      h210x7 +=  + ln2 * ( 1./4.*pow(mm3,4) + 7./5.*mm2*pow(mm3,3) + 33.
         /10.*pow(mm2,2)*pow(mm3,2) + 149./35.*pow(mm2,3)*mm3 + 23./7.*
         pow(mm2,4) + 11./7.*pow(mm2,5)*pow(mm3,-1) + 33./70.*pow(
         mm2,6)*pow(mm3,-2) + 3./35.*pow(mm2,7)*pow(mm3,-3) + 1./140.*
         pow(mm2,8)*pow(mm3,-4) );

      h210x7 +=  + ln1 * (  - 1./140.*pow(mm2,-4)*pow(mm3,8) - 3./35.*
         pow(mm2,-3)*pow(mm3,7) - 33./70.*pow(mm2,-2)*pow(mm3,6) - 11./
         7.*pow(mm2,-1)*pow(mm3,5) - 99./28.*pow(mm3,4) - 198./35.*mm2*
         pow(mm3,3) - 33./5.*pow(mm2,2)*pow(mm3,2) - 198./35.*pow(
         mm2,3)*mm3 - 99./28.*pow(mm2,4) - 11./7.*pow(mm2,5)*pow(
         mm3,-1) - 33./70.*pow(mm2,6)*pow(mm3,-2) - 3./35.*pow(mm2,7)*
         pow(mm3,-3) - 1./140.*pow(mm2,8)*pow(mm3,-4) );

      h210x7 +=  + 1./70.*pow(mm2,-3)*pow(mm3,7) + 23./140.*pow(mm2,-2)
         *pow(mm3,6) + 181./210.*pow(mm2,-1)*pow(mm3,5) + 13./5.*pow(
         mm3,4) + 172./35.*mm2*pow(mm3,3) + 1271./210.*pow(mm2,2)*pow(
         mm3,2) + 172./35.*pow(mm2,3)*mm3 + 13./5.*pow(mm2,4) + 181./
         210.*pow(mm2,5)*pow(mm3,-1) + 23./140.*pow(mm2,6)*pow(mm3,-2)
          + 1./70.*pow(mm2,7)*pow(mm3,-3);
      return h210x7*pi162/pow(mm1,6);
    case 8:
      h210x8 =
       + ln3 * (  - 2./315.*pow(mm2,-5)*pow(mm3,9) - 17./252.*pow(
         mm2,-4)*pow(mm3,8) - 101./315.*pow(mm2,-3)*pow(mm3,7) - 113./
         126.*pow(mm2,-2)*pow(mm3,6) - 517./315.*pow(mm2,-1)*pow(mm3,5)
          - 647./315.*pow(mm3,4) - 37./21.*mm2*pow(mm3,3) - 209./210.*
         pow(mm2,2)*pow(mm3,2) - 1./3.*pow(mm2,3)*mm3 - 1./20.*pow(
         mm2,4) );

      h210x8 +=  + ln2 * (  - 1./20.*pow(mm3,4) - 1./3.*mm2*pow(mm3,3)
          - 209./210.*pow(mm2,2)*pow(mm3,2) - 37./21.*pow(mm2,3)*mm3 - 
         647./315.*pow(mm2,4) - 517./315.*pow(mm2,5)*pow(mm3,-1) - 113./
         126.*pow(mm2,6)*pow(mm3,-2) - 101./315.*pow(mm2,7)*pow(mm3,-3)
          - 17./252.*pow(mm2,8)*pow(mm3,-4) - 2./315.*pow(mm2,9)*pow(
         mm3,-5) );

      h210x8 +=  + ln1 * ( 2./315.*pow(mm2,-5)*pow(mm3,9) + 17./252.*
         pow(mm2,-4)*pow(mm3,8) + 101./315.*pow(mm2,-3)*pow(mm3,7) + 
         113./126.*pow(mm2,-2)*pow(mm3,6) + 517./315.*pow(mm2,-1)*pow(
         mm3,5) + 2651./1260.*pow(mm3,4) + 44./21.*mm2*pow(mm3,3) + 209.
         /105.*pow(mm2,2)*pow(mm3,2) + 44./21.*pow(mm2,3)*mm3 + 2651./
         1260.*pow(mm2,4) + 517./315.*pow(mm2,5)*pow(mm3,-1) + 113./126.
         *pow(mm2,6)*pow(mm3,-2) + 101./315.*pow(mm2,7)*pow(mm3,-3) + 
         17./252.*pow(mm2,8)*pow(mm3,-4) + 2./315.*pow(mm2,9)*pow(
         mm3,-5) );

      h210x8 +=  - 4./315.*pow(mm2,-4)*pow(mm3,8) - 9./70.*pow(mm2,-3)*
         pow(mm3,7) - 437./756.*pow(mm2,-2)*pow(mm3,6) - 409./270.*pow(
         mm2,-1)*pow(mm3,5) - 2479./945.*pow(mm3,4) - 3167./945.*mm2*
         pow(mm3,3) - 6731./1890.*pow(mm2,2)*pow(mm3,2) - 3167./945.*
         pow(mm2,3)*mm3 - 2479./945.*pow(mm2,4) - 409./270.*pow(mm2,5)*
         pow(mm3,-1) - 437./756.*pow(mm2,6)*pow(mm3,-2) - 9./70.*pow(
         mm2,7)*pow(mm3,-3) - 4./315.*pow(mm2,8)*pow(mm3,-4);
      return h210x8*pi162/pow(mm1,8);
    default:
      std::cout << "h210sing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210x1 =
       + ln3 * (  - 1./9.*pow(mm3,3) - 1./2.*mm1*pow(mm3,2) - 1./3.*
         pow(mm1,2)*mm3 );

      h210x1 +=  + pow(ln3,2) * ( 1./6.*pow(mm3,3) + 1./6.*mm1*pow(
         mm3,2) );

      h210x1 +=  + ln2 * (  - 1./9.*pow(mm3,3) + 1./18.*mm1*pow(mm3,2)
          + 1./2.*pow(mm1,2)*mm3 + 7./18.*pow(mm1,3) + 1./18.*pow(
         mm1,4)*pow(mm3,-1) );

      h210x1 +=  + ln2*ln3 * ( 1./3.*pow(mm3,3) + 2./3.*mm1*pow(mm3,2)
          + 1./2.*pow(mm1,2)*mm3 + 1./6.*pow(mm1,3) );

      h210x1 +=  + pow(ln2,2) * ( 1./6.*pow(mm3,3) + 1./2.*mm1*pow(
         mm3,2) + 1./2.*pow(mm1,2)*mm3 + 1./6.*pow(mm1,3) );

      h210x1 +=  + ln1 * (  - 2./3.*pow(mm1,2)*mm3 - 2./3.*pow(mm1,3)
          - 1./18.*pow(mm1,4)*pow(mm3,-1) );

      h210x1 +=  + ln1*ln3 * (  - 1./3.*mm1*pow(mm3,2) - 1./2.*pow(
         mm1,2)*mm3 - 1./6.*pow(mm1,3) );

      h210x1 +=  + ln1*ln2 * ( 1./3.*mm1*pow(mm3,2) + 1./2.*pow(mm1,2)*
         mm3 + 1./6.*pow(mm1,3) );

      h210x1 +=  + 19./27.*pow(mm3,3) + 1./18.*pow(mm3,3)*pi2 + 38./27.
         *mm1*pow(mm3,2) + 1./9.*mm1*pow(mm3,2)*pi2 + 97./72.*pow(
         mm1,2)*mm3 + 1./12.*pow(mm1,2)*mm3*pi2 + 139./216.*pow(mm1,3)
          + 1./36.*pow(mm1,3)*pi2;
      return h210x1*pi162/mm2;
    case 2:
      h210x2 =
       + ln3 * ( 1./15.*pow(mm1,-1)*pow(mm3,4) + 3./20.*pow(mm3,3) + 1./
         12.*mm1*pow(mm3,2) );

      h210x2 +=  + ln2 * (  - 1./15.*pow(mm1,-1)*pow(mm3,4) - 19./60.*
         pow(mm3,3) - 7./12.*mm1*pow(mm3,2) - 1./2.*pow(mm1,2)*mm3 - 1./
         6.*pow(mm1,3) + 1./60.*pow(mm1,4)*pow(mm3,-1) + 1./60.*pow(
         mm1,5)*pow(mm3,-2) );

      h210x2 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm3,-1) - 1./60.*
         pow(mm1,5)*pow(mm3,-2) );

      h210x2 +=  - 47./360.*pow(mm3,3) - 47./120.*mm1*pow(mm3,2) - 17./
         40.*pow(mm1,2)*mm3 - 71./360.*pow(mm1,3) - 1./30.*pow(mm1,4)*
         pow(mm3,-1);
      return h210x2*pi162/pow(mm2,3);
    case 3:
      h210x3 =
       + ln3 * (  - 47./180.*pow(mm3,3) - 5./8.*mm1*pow(mm3,2) + 1./24.
         *mm1*pow(mm3,3)*pow(mm2,-1) - 7./24.*pow(mm1,2)*mm3 + 1./24.*
         pow(mm1,3) - 1./12.*pow(mm1,3)*mm3*pow(mm2,-1) - 1./24.*pow(
         mm1,4)*pow(mm3,-1) + 1./24.*pow(mm1,5)*pow(mm3,-1)*pow(
         mm2,-1) );

      h210x3 +=  + ln2 * ( 49./60.*pow(mm3,3) + 31./12.*mm1*pow(mm3,2)
          + 17./6.*pow(mm1,2)*mm3 + 7./6.*pow(mm1,3) + 1./12.*pow(
         mm1,4)*pow(mm3,-1) - 1./60.*pow(mm1,5)*pow(mm3,-2) );

      h210x3 +=  + ln2*ln3 * ( 1./6.*pow(mm3,3) + 1./2.*mm1*pow(mm3,2)
          + 1./2.*pow(mm1,2)*mm3 + 1./6.*pow(mm1,3) );

      h210x3 +=  + pow(ln2,2) * ( 1./6.*pow(mm3,3) + 1./2.*mm1*pow(
         mm3,2) + 1./2.*pow(mm1,2)*mm3 + 1./6.*pow(mm1,3) );

      h210x3 +=  + ln1 * (  - 7./24.*mm1*pow(mm3,2) - 1./24.*mm1*pow(
         mm3,3)*pow(mm2,-1) - 7./8.*pow(mm1,2)*mm3 - 47./72.*pow(
         mm1,3) + 1./12.*pow(mm1,3)*mm3*pow(mm2,-1) - 1./24.*pow(
         mm1,4)*pow(mm3,-1) + 1./60.*pow(mm1,5)*pow(mm3,-2) - 1./24.*
         pow(mm1,5)*pow(mm3,-1)*pow(mm2,-1) );

      h210x3 +=  + ln1*ln3 * (  - 1./6.*pow(mm3,3) - 1./2.*mm1*pow(
         mm3,2) - 1./2.*pow(mm1,2)*mm3 - 1./6.*pow(mm1,3) );

      h210x3 +=  + ln1*ln2 * ( 1./6.*pow(mm3,3) + 1./2.*mm1*pow(mm3,2)
          + 1./2.*pow(mm1,2)*mm3 + 1./6.*pow(mm1,3) );

      h210x3 +=  + 13./54.*pow(mm3,3) + 1./36.*pow(mm3,3)*pi2 + 13./15.
         *mm1*pow(mm3,2) + 1./12.*mm1*pow(mm3,2)*pi2 + 16./15.*pow(
         mm1,2)*mm3 + 1./12.*pow(mm1,2)*mm3*pi2 + 64./135.*pow(mm1,3)
          + 1./36.*pow(mm1,3)*pi2 + 1./30.*pow(mm1,4)*pow(mm3,-1);
      return h210x3*pi162/pow(mm2,3);
    case 4:
      h210x4 =
       + ln3 * ( 49./60.*pow(mm3,5) + 79./20.*mm1*pow(mm3,4) + 457./60.
         *pow(mm1,2)*pow(mm3,3) + 439./60.*pow(mm1,3)*pow(mm3,2) + 7./2.
         *pow(mm1,4)*mm3 + 2./3.*pow(mm1,5) );

      h210x4 +=  + pow(ln3,2) * ( 1./6.*pow(mm3,5) + 5./6.*mm1*pow(
         mm3,4) + 5./3.*pow(mm1,2)*pow(mm3,3) + 5./3.*pow(mm1,3)*pow(
         mm3,2) + 5./6.*pow(mm1,4)*mm3 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln2 * (  - 47./180.*pow(mm3,5) - 271./180.*mm1*pow(
         mm3,4) - 641./180.*pow(mm1,2)*pow(mm3,3) - 787./180.*pow(
         mm1,3)*pow(mm3,2) - 103./36.*pow(mm1,4)*mm3 - 157./180.*pow(
         mm1,5) - 11./180.*pow(mm1,6)*pow(mm3,-1) - 1./180.*pow(mm1,7)*
         pow(mm3,-2) - 1./90.*pow(mm1,8)*pow(mm3,-3) );

      h210x4 +=  + ln2*ln3 * ( 1./6.*pow(mm3,5) + 5./6.*mm1*pow(mm3,4)
          + 5./3.*pow(mm1,2)*pow(mm3,3) + 5./3.*pow(mm1,3)*pow(mm3,2)
          + 5./6.*pow(mm1,4)*mm3 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln1 * ( 1./3.*mm1*pow(mm3,4) + 3./2.*pow(mm1,2)*pow(
         mm3,3) + 47./18.*pow(mm1,3)*pow(mm3,2) + 77./36.*pow(mm1,4)*
         mm3 + 137./180.*pow(mm1,5) + 11./180.*pow(mm1,6)*pow(mm3,-1)
          + 1./180.*pow(mm1,7)*pow(mm3,-2) + 1./90.*pow(mm1,8)*pow(
         mm3,-3) );

      h210x4 +=  + ln1*ln3 * ( 1./6.*pow(mm3,5) + 5./6.*mm1*pow(mm3,4)
          + 5./3.*pow(mm1,2)*pow(mm3,3) + 5./3.*pow(mm1,3)*pow(mm3,2)
          + 5./6.*pow(mm1,4)*mm3 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln1*ln2 * (  - 1./6.*pow(mm3,5) - 5./6.*mm1*pow(
         mm3,4) - 5./3.*pow(mm1,2)*pow(mm3,3) - 5./3.*pow(mm1,3)*pow(
         mm3,2) - 5./6.*pow(mm1,4)*mm3 - 1./6.*pow(mm1,5) );

      h210x4 +=  + 13./54.*pow(mm3,5) + 1./36.*pow(mm3,5)*pi2 + 143./
         135.*mm1*pow(mm3,4) + 5./36.*mm1*pow(mm3,4)*pi2 + 47./27.*pow(
         mm1,2)*pow(mm3,3) + 5./18.*pow(mm1,2)*pow(mm3,3)*pi2 + 163./
         135.*pow(mm1,3)*pow(mm3,2) + 5./18.*pow(mm1,3)*pow(mm3,2)*pi2
          + 49./270.*pow(mm1,4)*mm3 + 5./36.*pow(mm1,4)*mm3*pi2 - 17./
         135.*pow(mm1,5) + 1./36.*pow(mm1,5)*pi2 + 1./45.*pow(mm1,7)*
         pow(mm3,-2);
      return h210x4*pi162/pow(mm2,5);
    case 5:
      h210x5 =
       + ln3 * (  - 1./42.*pow(mm1,-2)*pow(mm3,6) - 4./35.*pow(mm1,-1)*
         pow(mm3,5) - 29./140.*pow(mm3,4) - 1./6.*mm1*pow(mm3,3) - 1./
         20.*pow(mm1,2)*pow(mm3,2) );

      h210x5 +=  + ln2 * ( 1./42.*pow(mm1,-2)*pow(mm3,6) + 4./35.*pow(
         mm1,-1)*pow(mm3,5) + 29./140.*pow(mm3,4) + 1./6.*mm1*pow(
         mm3,3) + 1./20.*pow(mm1,2)*pow(mm3,2) - 1./60.*pow(mm1,4) - 3./
         70.*pow(mm1,5)*pow(mm3,-1) - 1./28.*pow(mm1,6)*pow(mm3,-2) - 1.
         /105.*pow(mm1,7)*pow(mm3,-3) );

      h210x5 +=  + ln1 * ( 1./60.*pow(mm1,4) + 3./70.*pow(mm1,5)*pow(
         mm3,-1) + 1./28.*pow(mm1,6)*pow(mm3,-2) + 1./105.*pow(mm1,7)*
         pow(mm3,-3) );

      h210x5 +=  - 1./21.*pow(mm1,-1)*pow(mm3,5) - 121./420.*pow(mm3,4)
          - 74./105.*mm1*pow(mm3,3) - 6./7.*pow(mm1,2)*pow(mm3,2) - 103.
         /210.*pow(mm1,3)*mm3 - 23./420.*pow(mm1,4) + 13./210.*pow(
         mm1,5)*pow(mm3,-1) + 2./105.*pow(mm1,6)*pow(mm3,-2);
      return h210x5*pi162/pow(mm2,6);
    case 6:
      h210x6 =
       + ln3 * ( 1./42.*pow(mm1,-2)*pow(mm3,6) + 6./35.*pow(mm1,-1)*
         pow(mm3,5) + 73./140.*pow(mm3,4) + 6./7.*mm1*pow(mm3,3) + 4./5.
         *pow(mm1,2)*pow(mm3,2) + 2./5.*pow(mm1,3)*mm3 + 1./12.*pow(
         mm1,4) );

      h210x6 +=  + ln2 * (  - 1./42.*pow(mm1,-2)*pow(mm3,6) - 6./35.*
         pow(mm1,-1)*pow(mm3,5) - 73./140.*pow(mm3,4) - 6./7.*mm1*pow(
         mm3,3) - 4./5.*pow(mm1,2)*pow(mm3,2) - 2./5.*pow(mm1,3)*mm3 - 
         1./10.*pow(mm1,4) - 2./35.*pow(mm1,5)*pow(mm3,-1) - 1./14.*
         pow(mm1,6)*pow(mm3,-2) - 4./105.*pow(mm1,7)*pow(mm3,-3) - 1./
         140.*pow(mm1,8)*pow(mm3,-4) );

      h210x6 +=  + ln1 * ( 1./60.*pow(mm1,4) + 2./35.*pow(mm1,5)*pow(
         mm3,-1) + 1./14.*pow(mm1,6)*pow(mm3,-2) + 4./105.*pow(mm1,7)*
         pow(mm3,-3) + 1./140.*pow(mm1,8)*pow(mm3,-4) );

      h210x6 +=  + 1./21.*pow(mm1,-1)*pow(mm3,5) + 33./140.*pow(mm3,4)
          + 31./70.*mm1*pow(mm3,3) + 53./140.*pow(mm1,2)*pow(mm3,2) + 
         31./210.*pow(mm1,3)*mm3 + 11./140.*pow(mm1,4) + 23./210.*pow(
         mm1,5)*pow(mm3,-1) + 29./420.*pow(mm1,6)*pow(mm3,-2) + 1./70.*
         pow(mm1,7)*pow(mm3,-3);
      return h210x6*pi162/pow(mm2,6);
    case 7:
      h210x7 =
       + ln3 * ( 1./7.*pow(mm1,-1)*pow(mm3,5) + 11./14.*pow(mm3,4) + 61.
         /35.*mm1*pow(mm3,3) + 39./20.*pow(mm1,2)*pow(mm3,2) + 11./10.*
         pow(mm1,3)*mm3 + 1./4.*pow(mm1,4) );

      h210x7 +=  + ln2 * (  - 1./7.*pow(mm1,-1)*pow(mm3,5) - 11./14.*
         pow(mm3,4) - 61./35.*mm1*pow(mm3,3) - 39./20.*pow(mm1,2)*pow(
         mm3,2) - 11./10.*pow(mm1,3)*mm3 - 1./4.*pow(mm1,4) + 1./140.*
         pow(mm1,6)*pow(mm3,-2) + 1./70.*pow(mm1,7)*pow(mm3,-3) + 1./
         140.*pow(mm1,8)*pow(mm3,-4) );

      h210x7 +=  + ln1 * (  - 1./140.*pow(mm1,6)*pow(mm3,-2) - 1./70.*
         pow(mm1,7)*pow(mm3,-3) - 1./140.*pow(mm1,8)*pow(mm3,-4) );

      h210x7 +=  - 1./21.*pow(mm3,4) - 5./21.*mm1*pow(mm3,3) - 191./420.
         *pow(mm1,2)*pow(mm3,2) - 41./105.*pow(mm1,3)*mm3 - 13./105.*
         pow(mm1,4) - 1./210.*pow(mm1,5)*pow(mm3,-1) - 3./140.*pow(
         mm1,6)*pow(mm3,-2) - 1./70.*pow(mm1,7)*pow(mm3,-3);
      return h210x7*pi162/pow(mm2,6);
    case 8:
      h210x8 =
       + ln3 * (  - 1./63.*pow(mm1,-3)*pow(mm3,7) - 5./42.*pow(mm1,-2)*
         pow(mm3,6) - 41./105.*pow(mm1,-1)*pow(mm3,5) - 923./1260.*pow(
         mm3,4) - 6./7.*mm1*pow(mm3,3) - 22./35.*pow(mm1,2)*pow(mm3,2)
          - 4./15.*pow(mm1,3)*mm3 - 1./20.*pow(mm1,4) );

      h210x8 +=  + ln2 * ( 1./63.*pow(mm1,-3)*pow(mm3,7) + 5./42.*pow(
         mm1,-2)*pow(mm3,6) + 41./105.*pow(mm1,-1)*pow(mm3,5) + 923./
         1260.*pow(mm3,4) + 6./7.*mm1*pow(mm3,3) + 22./35.*pow(mm1,2)*
         pow(mm3,2) + 4./15.*pow(mm1,3)*mm3 + 1./14.*pow(mm1,4) + 3./35.
         *pow(mm1,5)*pow(mm3,-1) + 17./126.*pow(mm1,6)*pow(mm3,-2) + 11.
         /105.*pow(mm1,7)*pow(mm3,-3) + 17./420.*pow(mm1,8)*pow(mm3,-4)
          + 2./315.*pow(mm1,9)*pow(mm3,-5) );

      h210x8 +=  + ln1 * (  - 3./140.*pow(mm1,4) - 3./35.*pow(mm1,5)*
         pow(mm3,-1) - 17./126.*pow(mm1,6)*pow(mm3,-2) - 11./105.*pow(
         mm1,7)*pow(mm3,-3) - 17./420.*pow(mm1,8)*pow(mm3,-4) - 2./315.
         *pow(mm1,9)*pow(mm3,-5) );

      h210x8 +=  - 2./63.*pow(mm1,-2)*pow(mm3,6) - 2./9.*pow(mm1,-1)*
         pow(mm3,5) - 2437./3780.*pow(mm3,4) - 617./630.*mm1*pow(mm3,3)
          - 29./36.*pow(mm1,2)*pow(mm3,2) - 131./378.*pow(mm1,3)*mm3 - 
         181./1260.*pow(mm1,4) - 17./90.*pow(mm1,5)*pow(mm3,-1) - 131./
         756.*pow(mm1,6)*pow(mm3,-2) - 47./630.*pow(mm1,7)*pow(mm3,-3)
          - 4./315.*pow(mm1,8)*pow(mm3,-4);
      return h210x8*pi162/pow(mm2,8);
    default:
      std::cout << "h210sing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210x1 =
       + ln3 * (  - 1./9.*pow(mm2,3) + 1./18.*mm1*pow(mm2,2) + 1./2.*
         pow(mm1,2)*mm2 + 7./18.*pow(mm1,3) + 1./18.*pow(mm1,4)*pow(
         mm2,-1) );

      h210x1 +=  + pow(ln3,2) * ( 1./6.*pow(mm2,3) + 1./2.*mm1*pow(
         mm2,2) + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3) );

      h210x1 +=  + ln2 * (  - 1./9.*pow(mm2,3) - 1./2.*mm1*pow(mm2,2)
          - 1./3.*pow(mm1,2)*mm2 );

      h210x1 +=  + ln2*ln3 * ( 1./3.*pow(mm2,3) + 2./3.*mm1*pow(mm2,2)
          + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3) );

      h210x1 +=  + pow(ln2,2) * ( 1./6.*pow(mm2,3) + 1./6.*mm1*pow(
         mm2,2) );

      h210x1 +=  + ln1 * (  - 2./3.*pow(mm1,2)*mm2 - 2./3.*pow(mm1,3)
          - 1./18.*pow(mm1,4)*pow(mm2,-1) );

      h210x1 +=  + ln1*ln3 * ( 1./3.*mm1*pow(mm2,2) + 1./2.*pow(mm1,2)*
         mm2 + 1./6.*pow(mm1,3) );

      h210x1 +=  + ln1*ln2 * (  - 1./3.*mm1*pow(mm2,2) - 1./2.*pow(
         mm1,2)*mm2 - 1./6.*pow(mm1,3) );

      h210x1 +=  + 19./27.*pow(mm2,3) + 1./18.*pow(mm2,3)*pi2 + 38./27.
         *mm1*pow(mm2,2) + 1./9.*mm1*pow(mm2,2)*pi2 + 97./72.*pow(
         mm1,2)*mm2 + 1./12.*pow(mm1,2)*mm2*pi2 + 139./216.*pow(mm1,3)
          + 1./36.*pow(mm1,3)*pi2;
      return h210x1*pi162/mm3;
    case 2:
      h210x2 =
       + ln3 * (  - 1./15.*pow(mm1,-1)*pow(mm2,4) - 19./60.*pow(mm2,3)
          - 7./12.*mm1*pow(mm2,2) - 1./2.*pow(mm1,2)*mm2 - 1./6.*pow(
         mm1,3) + 1./60.*pow(mm1,4)*pow(mm2,-1) + 1./60.*pow(mm1,5)*
         pow(mm2,-2) );

      h210x2 +=  + ln2 * ( 1./15.*pow(mm1,-1)*pow(mm2,4) + 3./20.*pow(
         mm2,3) + 1./12.*mm1*pow(mm2,2) );

      h210x2 +=  + ln1 * (  - 1./60.*pow(mm1,4)*pow(mm2,-1) - 1./60.*
         pow(mm1,5)*pow(mm2,-2) );

      h210x2 +=  - 47./360.*pow(mm2,3) - 47./120.*mm1*pow(mm2,2) - 17./
         40.*pow(mm1,2)*mm2 - 71./360.*pow(mm1,3) - 1./30.*pow(mm1,4)*
         pow(mm2,-1);
      return h210x2*pi162/pow(mm3,3);
    case 3:
      h210x3 =
       + ln3 * (  - 47./180.*pow(mm2,3) - 59./60.*mm1*pow(mm2,2) - 4./3.
         *pow(mm1,2)*mm2 - 13./18.*pow(mm1,3) - 1./12.*pow(mm1,4)*pow(
         mm2,-1) + 1./60.*pow(mm1,5)*pow(mm2,-2) - 1./90.*pow(mm1,6)*
         pow(mm2,-3) );

      h210x3 +=  + ln2 * ( 49./60.*pow(mm2,3) + 139./60.*mm1*pow(mm2,2)
          + 13./6.*pow(mm1,2)*mm2 + 2./3.*pow(mm1,3) );

      h210x3 +=  + ln2*ln3 * ( 1./6.*pow(mm2,3) + 1./2.*mm1*pow(mm2,2)
          + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3) );

      h210x3 +=  + pow(ln2,2) * ( 1./6.*pow(mm2,3) + 1./2.*mm1*pow(
         mm2,2) + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3) );

      h210x3 +=  + ln1 * ( 1./3.*mm1*pow(mm2,2) + 5./6.*pow(mm1,2)*mm2
          + 11./18.*pow(mm1,3) + 1./12.*pow(mm1,4)*pow(mm2,-1) - 1./60.
         *pow(mm1,5)*pow(mm2,-2) + 1./90.*pow(mm1,6)*pow(mm2,-3) );

      h210x3 +=  + ln1*ln3 * (  - 1./6.*pow(mm2,3) - 1./2.*mm1*pow(
         mm2,2) - 1./2.*pow(mm1,2)*mm2 - 1./6.*pow(mm1,3) );

      h210x3 +=  + ln1*ln2 * ( 1./6.*pow(mm2,3) + 1./2.*mm1*pow(mm2,2)
          + 1./2.*pow(mm1,2)*mm2 + 1./6.*pow(mm1,3) );

      h210x3 +=  + 13./54.*pow(mm2,3) + 1./36.*pow(mm2,3)*pi2 + 26./45.
         *mm1*pow(mm2,2) + 1./12.*mm1*pow(mm2,2)*pi2 + 31./90.*pow(
         mm1,2)*mm2 + 1./12.*pow(mm1,2)*mm2*pi2 - 8./135.*pow(mm1,3) + 
         1./36.*pow(mm1,3)*pi2 - 2./45.*pow(mm1,4)*pow(mm2,-1) + 1./45.
         *pow(mm1,5)*pow(mm2,-2);
      return h210x3*pi162/pow(mm3,3);
    case 4:
      h210x4 =
       + ln3 * ( 49./60.*pow(mm2,5) + 253./60.*mm1*pow(mm2,4) + 529./60.
         *pow(mm1,2)*pow(mm2,3) + 113./12.*pow(mm1,3)*pow(mm2,2) + 21./
         4.*pow(mm1,4)*mm2 + 79./60.*pow(mm1,5) + 1./20.*pow(mm1,6)*
         pow(mm2,-1) - 1./60.*pow(mm1,7)*pow(mm2,-2) );

      h210x4 +=  + pow(ln3,2) * ( 1./6.*pow(mm2,5) + 5./6.*mm1*pow(
         mm2,4) + 5./3.*pow(mm1,2)*pow(mm2,3) + 5./3.*pow(mm1,3)*pow(
         mm2,2) + 5./6.*pow(mm1,4)*mm2 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln2 * (  - 47./180.*pow(mm2,5) - 199./180.*mm1*pow(
         mm2,4) - 317./180.*pow(mm1,2)*pow(mm2,3) - 5./4.*pow(mm1,3)*
         pow(mm2,2) - 1./3.*pow(mm1,4)*mm2 );

      h210x4 +=  + ln2*ln3 * ( 1./6.*pow(mm2,5) + 5./6.*mm1*pow(mm2,4)
          + 5./3.*pow(mm1,2)*pow(mm2,3) + 5./3.*pow(mm1,3)*pow(mm2,2)
          + 5./6.*pow(mm1,4)*mm2 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln1 * (  - 1./3.*mm1*pow(mm2,4) - 3./2.*pow(mm1,2)*
         pow(mm2,3) - 47./18.*pow(mm1,3)*pow(mm2,2) - 77./36.*pow(
         mm1,4)*mm2 - 137./180.*pow(mm1,5) - 1./20.*pow(mm1,6)*pow(
         mm2,-1) + 1./60.*pow(mm1,7)*pow(mm2,-2) );

      h210x4 +=  + ln1*ln3 * ( 1./6.*pow(mm2,5) + 5./6.*mm1*pow(mm2,4)
          + 5./3.*pow(mm1,2)*pow(mm2,3) + 5./3.*pow(mm1,3)*pow(mm2,2)
          + 5./6.*pow(mm1,4)*mm2 + 1./6.*pow(mm1,5) );

      h210x4 +=  + ln1*ln2 * (  - 1./6.*pow(mm2,5) - 5./6.*mm1*pow(
         mm2,4) - 5./3.*pow(mm1,2)*pow(mm2,3) - 5./3.*pow(mm1,3)*pow(
         mm2,2) - 5./6.*pow(mm1,4)*mm2 - 1./6.*pow(mm1,5) );

      h210x4 +=  + 13./54.*pow(mm2,5) + 1./36.*pow(mm2,5)*pi2 + 182./
         135.*mm1*pow(mm2,4) + 5./36.*mm1*pow(mm2,4)*pi2 + 821./270.*
         pow(mm1,2)*pow(mm2,3) + 5./18.*pow(mm1,2)*pow(mm2,3)*pi2 + 469.
         /135.*pow(mm1,3)*pow(mm2,2) + 5./18.*pow(mm1,3)*pow(mm2,2)*pi2
          + 553./270.*pow(mm1,4)*mm2 + 5./36.*pow(mm1,4)*mm2*pi2 + 73./
         135.*pow(mm1,5) + 1./36.*pow(mm1,5)*pi2 + 1./30.*pow(mm1,6)*
         pow(mm2,-1);
      return h210x4*pi162/pow(mm3,5);
    case 5:
      h210x5 =
       + ln3 * (  - 1./42.*pow(mm1,-2)*pow(mm2,6) - 6./35.*pow(mm1,-1)*
         pow(mm2,5) - 73./140.*pow(mm2,4) - 6./7.*mm1*pow(mm2,3) - 4./5.
         *pow(mm1,2)*pow(mm2,2) - 2./5.*pow(mm1,3)*mm2 - 1./10.*pow(
         mm1,4) - 2./35.*pow(mm1,5)*pow(mm2,-1) - 1./14.*pow(mm1,6)*
         pow(mm2,-2) - 4./105.*pow(mm1,7)*pow(mm2,-3) - 1./140.*pow(
         mm1,8)*pow(mm2,-4) );

      h210x5 +=  + ln2 * ( 1./42.*pow(mm1,-2)*pow(mm2,6) + 6./35.*pow(
         mm1,-1)*pow(mm2,5) + 73./140.*pow(mm2,4) + 6./7.*mm1*pow(
         mm2,3) + 4./5.*pow(mm1,2)*pow(mm2,2) + 2./5.*pow(mm1,3)*mm2 + 
         1./12.*pow(mm1,4) );

      h210x5 +=  + ln1 * ( 1./60.*pow(mm1,4) + 2./35.*pow(mm1,5)*pow(
         mm2,-1) + 1./14.*pow(mm1,6)*pow(mm2,-2) + 4./105.*pow(mm1,7)*
         pow(mm2,-3) + 1./140.*pow(mm1,8)*pow(mm2,-4) );

      h210x5 +=  + 1./21.*pow(mm1,-1)*pow(mm2,5) + 33./140.*pow(mm2,4)
          + 31./70.*mm1*pow(mm2,3) + 53./140.*pow(mm1,2)*pow(mm2,2) + 
         31./210.*pow(mm1,3)*mm2 + 11./140.*pow(mm1,4) + 23./210.*pow(
         mm1,5)*pow(mm2,-1) + 29./420.*pow(mm1,6)*pow(mm2,-2) + 1./70.*
         pow(mm1,7)*pow(mm2,-3);
      return h210x5*pi162/pow(mm3,6);
    case 6:
      h210x6 =
       + ln3 * ( 1./42.*pow(mm1,-2)*pow(mm2,6) + 4./35.*pow(mm1,-1)*
         pow(mm2,5) + 29./140.*pow(mm2,4) + 1./6.*mm1*pow(mm2,3) + 1./
         20.*pow(mm1,2)*pow(mm2,2) - 1./60.*pow(mm1,4) - 3./70.*pow(
         mm1,5)*pow(mm2,-1) - 1./28.*pow(mm1,6)*pow(mm2,-2) - 1./105.*
         pow(mm1,7)*pow(mm2,-3) );

      h210x6 +=  + ln2 * (  - 1./42.*pow(mm1,-2)*pow(mm2,6) - 4./35.*
         pow(mm1,-1)*pow(mm2,5) - 29./140.*pow(mm2,4) - 1./6.*mm1*pow(
         mm2,3) - 1./20.*pow(mm1,2)*pow(mm2,2) );

      h210x6 +=  + ln1 * ( 1./60.*pow(mm1,4) + 3./70.*pow(mm1,5)*pow(
         mm2,-1) + 1./28.*pow(mm1,6)*pow(mm2,-2) + 1./105.*pow(mm1,7)*
         pow(mm2,-3) );

      h210x6 +=  - 1./21.*pow(mm1,-1)*pow(mm2,5) - 121./420.*pow(mm2,4)
          - 74./105.*mm1*pow(mm2,3) - 6./7.*pow(mm1,2)*pow(mm2,2) - 103.
         /210.*pow(mm1,3)*mm2 - 23./420.*pow(mm1,4) + 13./210.*pow(
         mm1,5)*pow(mm2,-1) + 2./105.*pow(mm1,6)*pow(mm2,-2);
      return h210x6*pi162/pow(mm3,6);
    case 7:
      h210x7 =
       + ln3 * (  - 1./7.*pow(mm1,-1)*pow(mm2,5) - 11./14.*pow(mm2,4)
          - 61./35.*mm1*pow(mm2,3) - 39./20.*pow(mm1,2)*pow(mm2,2) - 11.
         /10.*pow(mm1,3)*mm2 - 1./4.*pow(mm1,4) + 1./140.*pow(mm1,6)*
         pow(mm2,-2) + 1./70.*pow(mm1,7)*pow(mm2,-3) + 1./140.*pow(
         mm1,8)*pow(mm2,-4) );

      h210x7 +=  + ln2 * ( 1./7.*pow(mm1,-1)*pow(mm2,5) + 11./14.*pow(
         mm2,4) + 61./35.*mm1*pow(mm2,3) + 39./20.*pow(mm1,2)*pow(
         mm2,2) + 11./10.*pow(mm1,3)*mm2 + 1./4.*pow(mm1,4) );

      h210x7 +=  + ln1 * (  - 1./140.*pow(mm1,6)*pow(mm2,-2) - 1./70.*
         pow(mm1,7)*pow(mm2,-3) - 1./140.*pow(mm1,8)*pow(mm2,-4) );

      h210x7 +=  - 1./21.*pow(mm2,4) - 5./21.*mm1*pow(mm2,3) - 191./420.
         *pow(mm1,2)*pow(mm2,2) - 41./105.*pow(mm1,3)*mm2 - 13./105.*
         pow(mm1,4) - 1./210.*pow(mm1,5)*pow(mm2,-1) - 3./140.*pow(
         mm1,6)*pow(mm2,-2) - 1./70.*pow(mm1,7)*pow(mm2,-3);
      return h210x7*pi162/pow(mm3,6);
    case 8:
      h210x8 =
       + ln3 * ( 1./63.*pow(mm1,-3)*pow(mm2,7) + 5./42.*pow(mm1,-2)*
         pow(mm2,6) + 41./105.*pow(mm1,-1)*pow(mm2,5) + 923./1260.*pow(
         mm2,4) + 6./7.*mm1*pow(mm2,3) + 22./35.*pow(mm1,2)*pow(mm2,2)
          + 4./15.*pow(mm1,3)*mm2 + 1./14.*pow(mm1,4) + 3./35.*pow(
         mm1,5)*pow(mm2,-1) + 17./126.*pow(mm1,6)*pow(mm2,-2) + 11./105.
         *pow(mm1,7)*pow(mm2,-3) + 17./420.*pow(mm1,8)*pow(mm2,-4) + 2./
         315.*pow(mm1,9)*pow(mm2,-5) );

      h210x8 +=  + ln2 * (  - 1./63.*pow(mm1,-3)*pow(mm2,7) - 5./42.*
         pow(mm1,-2)*pow(mm2,6) - 41./105.*pow(mm1,-1)*pow(mm2,5) - 923.
         /1260.*pow(mm2,4) - 6./7.*mm1*pow(mm2,3) - 22./35.*pow(mm1,2)*
         pow(mm2,2) - 4./15.*pow(mm1,3)*mm2 - 1./20.*pow(mm1,4) );

      h210x8 +=  + ln1 * (  - 3./140.*pow(mm1,4) - 3./35.*pow(mm1,5)*
         pow(mm2,-1) - 17./126.*pow(mm1,6)*pow(mm2,-2) - 11./105.*pow(
         mm1,7)*pow(mm2,-3) - 17./420.*pow(mm1,8)*pow(mm2,-4) - 2./315.
         *pow(mm1,9)*pow(mm2,-5) );

      h210x8 +=  - 2./63.*pow(mm1,-2)*pow(mm2,6) - 2./9.*pow(mm1,-1)*
         pow(mm2,5) - 2437./3780.*pow(mm2,4) - 617./630.*mm1*pow(mm2,3)
          - 29./36.*pow(mm1,2)*pow(mm2,2) - 131./378.*pow(mm1,3)*mm2 - 
         181./1260.*pow(mm1,4) - 17./90.*pow(mm1,5)*pow(mm2,-1) - 131./
         756.*pow(mm1,6)*pow(mm2,-2) - 47./630.*pow(mm1,7)*pow(mm2,-3)
          - 4./315.*pow(mm1,8)*pow(mm2,-4);
      return h210x8*pi162/pow(mm3,8);
    default:
      std::cout << "h210sing, wrong iprop = "<<iprop<<'\n';
      return 0.;
     }
  }
  std::cout << "funny masses in h210sing\n";
  return 0.;
}

double h210psing(const int iprop, const double m1sq, const double m2sq,
	      const double m3sq, const double xmu2){
  double mm1 = sqrt(m1sq);
  double mm2 = sqrt(m2sq);
  double mm3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double h210px1,h210px2,h210px3,h210px4,h210px5,h210px6,h210px7,h210px8;
  // three cases here m1+m2 = m3 (3), m1+m3 = m2 (2), m2 + m3 = m1 (1)
  // case 1
  if (fabs(mm2+mm3-mm1) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210px1 =
       + ln3 * (  - 1./420.*pow(mm2,-3)*pow(mm3,6) - 3./140.*pow(
         mm2,-2)*pow(mm3,5) - 3./35.*pow(mm2,-1)*pow(mm3,4) - 7./60.*
         pow(mm3,3) - 1./20.*mm2*pow(mm3,2) );

      h210px1 +=  + ln2 * (  - 1./20.*pow(mm2,2)*mm3 - 7./60.*pow(
         mm2,3) - 3./35.*pow(mm2,4)*pow(mm3,-1) - 3./140.*pow(mm2,5)*
         pow(mm3,-2) - 1./420.*pow(mm2,6)*pow(mm3,-3) );

      h210px1 +=  + ln1 * ( 1./420.*pow(mm2,-3)*pow(mm3,6) + 3./140.*
         pow(mm2,-2)*pow(mm3,5) + 3./35.*pow(mm2,-1)*pow(mm3,4) + 1./5.
         *pow(mm3,3) + 3./10.*mm2*pow(mm3,2) + 3./10.*pow(mm2,2)*mm3 + 
         1./5.*pow(mm2,3) + 3./35.*pow(mm2,4)*pow(mm3,-1) + 3./140.*
         pow(mm2,5)*pow(mm3,-2) + 1./420.*pow(mm2,6)*pow(mm3,-3) );

      h210px1 +=  - 1./210.*pow(mm2,-2)*pow(mm3,5) - 17./420.*pow(
         mm2,-1)*pow(mm3,4) - 653./10080.*pow(mm3,3) - 17./672.*mm2*
         pow(mm3,2) - 17./672.*pow(mm2,2)*mm3 - 653./10080.*pow(mm2,3)
          - 17./420.*pow(mm2,4)*pow(mm3,-1) - 1./210.*pow(mm2,5)*pow(
         mm3,-2);
      return h210px1*pi162/pow(mm1,3);
    case 2:
      h210px2 =
       + ln3 * ( 1./504.*pow(mm2,-4)*pow(mm3,7) + 37./2520.*pow(mm2,-3)
         *pow(mm3,6) + 113./2520.*pow(mm2,-2)*pow(mm3,5) + 59./840.*
         pow(mm2,-1)*pow(mm3,4) + 23./420.*pow(mm3,3) + 1./60.*mm2*pow(
         mm3,2) );

      h210px2 +=  + ln2 * ( 1./60.*pow(mm2,2)*mm3 + 23./420.*pow(mm2,3)
          + 59./840.*pow(mm2,4)*pow(mm3,-1) + 113./2520.*pow(mm2,5)*
         pow(mm3,-2) + 37./2520.*pow(mm2,6)*pow(mm3,-3) + 1./504.*pow(
         mm2,7)*pow(mm3,-4) );

      h210px2 +=  + ln1 * (  - 1./504.*pow(mm2,-4)*pow(mm3,7) - 37./
         2520.*pow(mm2,-3)*pow(mm3,6) - 113./2520.*pow(mm2,-2)*pow(
         mm3,5) - 59./840.*pow(mm2,-1)*pow(mm3,4) - 23./420.*pow(mm3,3)
          - 1./60.*mm2*pow(mm3,2) - 1./60.*pow(mm2,2)*mm3 - 23./420.*
         pow(mm2,3) - 59./840.*pow(mm2,4)*pow(mm3,-1) - 113./2520.*pow(
         mm2,5)*pow(mm3,-2) - 37./2520.*pow(mm2,6)*pow(mm3,-3) - 1./504.
         *pow(mm2,7)*pow(mm3,-4) );

      h210px2 +=  + 1./252.*pow(mm2,-3)*pow(mm3,6) + 23./840.*pow(
         mm2,-2)*pow(mm3,5) + 577./7560.*pow(mm2,-1)*pow(mm3,4) + 947./
         7560.*pow(mm3,3) + 227./1512.*mm2*pow(mm3,2) + 227./1512.*pow(
         mm2,2)*mm3 + 947./7560.*pow(mm2,3) + 577./7560.*pow(mm2,4)*
         pow(mm3,-1) + 23./840.*pow(mm2,5)*pow(mm3,-2) + 1./252.*pow(
         mm2,6)*pow(mm3,-3);
      return h210px2*pi162/pow(mm1,5);
    case 3:
      h210px3 =
       + ln3 * ( 1./630.*pow(mm2,-5)*pow(mm3,8) + 13./840.*pow(mm2,-4)*
         pow(mm3,7) + 11./168.*pow(mm2,-3)*pow(mm3,6) + 11./72.*pow(
         mm2,-2)*pow(mm3,5) + 11./56.*pow(mm2,-1)*pow(mm3,4) + 9./70.*
         pow(mm3,3) + 1./30.*mm2*pow(mm3,2) );

      h210px3 +=  + ln2 * (  - 1./20.*pow(mm3,3) - 13./60.*mm2*pow(
         mm3,2) - 11./28.*pow(mm2,2)*mm3 - 11./28.*pow(mm2,3) - 121./
         504.*pow(mm2,4)*pow(mm3,-1) - 11./120.*pow(mm2,5)*pow(mm3,-2)
          - 17./840.*pow(mm2,6)*pow(mm3,-3) - 1./504.*pow(mm2,7)*pow(
         mm3,-4) );

      h210px3 +=  + ln1 * (  - 1./630.*pow(mm2,-5)*pow(mm3,8) - 13./840.
         *pow(mm2,-4)*pow(mm3,7) - 11./168.*pow(mm2,-3)*pow(mm3,6) - 11.
         /72.*pow(mm2,-2)*pow(mm3,5) - 11./56.*pow(mm2,-1)*pow(mm3,4)
          - 11./140.*pow(mm3,3) + 11./60.*mm2*pow(mm3,2) + 11./28.*pow(
         mm2,2)*mm3 + 11./28.*pow(mm2,3) + 121./504.*pow(mm2,4)*pow(
         mm3,-1) + 11./120.*pow(mm2,5)*pow(mm3,-2) + 17./840.*pow(
         mm2,6)*pow(mm3,-3) + 1./504.*pow(mm2,7)*pow(mm3,-4) );

      h210px3 +=  + 1./315.*pow(mm2,-4)*pow(mm3,7) + 37./1260.*pow(
         mm2,-3)*pow(mm3,6) + 881./7560.*pow(mm2,-2)*pow(mm3,5) + 629./
         2520.*pow(mm2,-1)*pow(mm3,4) + 653./2520.*pow(mm3,3) - 11./
         7560.*mm2*pow(mm3,2) - 781./2520.*pow(mm2,2)*mm3 - 853./2520.*
         pow(mm2,3) - 1243./7560.*pow(mm2,4)*pow(mm3,-1) - 97./2520.*
         pow(mm2,5)*pow(mm3,-2) - 1./252.*pow(mm2,6)*pow(mm3,-3);
      return h210px3*pi162/pow(mm1,5);
    case 4:
      h210px4 =
       + ln3 * (  - 1./504.*pow(mm2,-4)*pow(mm3,7) - 17./840.*pow(
         mm2,-3)*pow(mm3,6) - 11./120.*pow(mm2,-2)*pow(mm3,5) - 121./
         504.*pow(mm2,-1)*pow(mm3,4) - 11./28.*pow(mm3,3) - 11./28.*mm2
         *pow(mm3,2) - 13./60.*pow(mm2,2)*mm3 - 1./20.*pow(mm2,3) );

      h210px4 +=  + ln2 * ( 1./30.*pow(mm2,2)*mm3 + 9./70.*pow(mm2,3)
          + 11./56.*pow(mm2,4)*pow(mm3,-1) + 11./72.*pow(mm2,5)*pow(
         mm3,-2) + 11./168.*pow(mm2,6)*pow(mm3,-3) + 13./840.*pow(
         mm2,7)*pow(mm3,-4) + 1./630.*pow(mm2,8)*pow(mm3,-5) );

      h210px4 +=  + ln1 * ( 1./504.*pow(mm2,-4)*pow(mm3,7) + 17./840.*
         pow(mm2,-3)*pow(mm3,6) + 11./120.*pow(mm2,-2)*pow(mm3,5) + 121.
         /504.*pow(mm2,-1)*pow(mm3,4) + 11./28.*pow(mm3,3) + 11./28.*
         mm2*pow(mm3,2) + 11./60.*pow(mm2,2)*mm3 - 11./140.*pow(mm2,3)
          - 11./56.*pow(mm2,4)*pow(mm3,-1) - 11./72.*pow(mm2,5)*pow(
         mm3,-2) - 11./168.*pow(mm2,6)*pow(mm3,-3) - 13./840.*pow(
         mm2,7)*pow(mm3,-4) - 1./630.*pow(mm2,8)*pow(mm3,-5) );

      h210px4 +=  - 1./252.*pow(mm2,-3)*pow(mm3,6) - 97./2520.*pow(
         mm2,-2)*pow(mm3,5) - 1243./7560.*pow(mm2,-1)*pow(mm3,4) - 853./
         2520.*pow(mm3,3) - 781./2520.*mm2*pow(mm3,2) - 11./7560.*pow(
         mm2,2)*mm3 + 653./2520.*pow(mm2,3) + 629./2520.*pow(mm2,4)*
         pow(mm3,-1) + 881./7560.*pow(mm2,5)*pow(mm3,-2) + 37./1260.*
         pow(mm2,6)*pow(mm3,-3) + 1./315.*pow(mm2,7)*pow(mm3,-4);
      return h210px4*pi162/pow(mm1,5);
    case 5:
      h210px5 =
       + ln3 * (  - 5./2772.*pow(mm2,-6)*pow(mm3,9) - 229./13860.*pow(
         mm2,-5)*pow(mm3,8) - 611./9240.*pow(mm2,-4)*pow(mm3,7) - 835./
         5544.*pow(mm2,-3)*pow(mm3,6) - 107./504.*pow(mm2,-2)*pow(
         mm3,5) - 31./168.*pow(mm2,-1)*pow(mm3,4) - 19./210.*pow(mm3,3)
          - 2./105.*mm2*pow(mm3,2) );

      h210px5 +=  + ln2 * ( 1./60.*pow(mm3,3) + 37./420.*mm2*pow(mm3,2)
          + 3./14.*pow(mm2,2)*mm3 + 20./63.*pow(mm2,3) + 157./504.*pow(
         mm2,4)*pow(mm3,-1) + 1879./9240.*pow(mm2,5)*pow(mm3,-2) + 2347.
         /27720.*pow(mm2,6)*pow(mm3,-3) + 113./5544.*pow(mm2,7)*pow(
         mm3,-4) + 1./462.*pow(mm2,8)*pow(mm3,-5) );

      h210px5 +=  + ln1 * ( 5./2772.*pow(mm2,-6)*pow(mm3,9) + 229./
         13860.*pow(mm2,-5)*pow(mm3,8) + 611./9240.*pow(mm2,-4)*pow(
         mm3,7) + 835./5544.*pow(mm2,-3)*pow(mm3,6) + 107./504.*pow(
         mm2,-2)*pow(mm3,5) + 31./168.*pow(mm2,-1)*pow(mm3,4) + 31./420.
         *pow(mm3,3) - 29./420.*mm2*pow(mm3,2) - 3./14.*pow(mm2,2)*mm3
          - 20./63.*pow(mm2,3) - 157./504.*pow(mm2,4)*pow(mm3,-1) - 
         1879./9240.*pow(mm2,5)*pow(mm3,-2) - 2347./27720.*pow(mm2,6)*
         pow(mm3,-3) - 113./5544.*pow(mm2,7)*pow(mm3,-4) - 1./462.*pow(
         mm2,8)*pow(mm3,-5) );

      h210px5 +=  - 5./1386.*pow(mm2,-5)*pow(mm3,8) - 433./13860.*pow(
         mm2,-4)*pow(mm3,7) - 221./1890.*pow(mm2,-3)*pow(mm3,6) - 2549./
         10395.*pow(mm2,-2)*pow(mm3,5) - 1076./3465.*pow(mm2,-1)*pow(
         mm3,4) - 17033./83160.*pow(mm3,3) + 5587./83160.*mm2*pow(
         mm3,2) + 136./385.*pow(mm2,2)*mm3 + 18889./41580.*pow(mm2,3)
          + 27821./83160.*pow(mm2,4)*pow(mm3,-1) + 379./2520.*pow(
         mm2,5)*pow(mm3,-2) + 107./2772.*pow(mm2,6)*pow(mm3,-3) + 1./
         231.*pow(mm2,7)*pow(mm3,-4);
      return h210px5*pi162/pow(mm1,7);
    case 6:
      h210px6 =
       + ln3 * ( 1./462.*pow(mm2,-5)*pow(mm3,8) + 113./5544.*pow(
         mm2,-4)*pow(mm3,7) + 2347./27720.*pow(mm2,-3)*pow(mm3,6) + 
         1879./9240.*pow(mm2,-2)*pow(mm3,5) + 157./504.*pow(mm2,-1)*
         pow(mm3,4) + 20./63.*pow(mm3,3) + 3./14.*mm2*pow(mm3,2) + 37./
         420.*pow(mm2,2)*mm3 + 1./60.*pow(mm2,3) );

      h210px6 +=  + ln2 * (  - 2./105.*pow(mm2,2)*mm3 - 19./210.*pow(
         mm2,3) - 31./168.*pow(mm2,4)*pow(mm3,-1) - 107./504.*pow(
         mm2,5)*pow(mm3,-2) - 835./5544.*pow(mm2,6)*pow(mm3,-3) - 611./
         9240.*pow(mm2,7)*pow(mm3,-4) - 229./13860.*pow(mm2,8)*pow(
         mm3,-5) - 5./2772.*pow(mm2,9)*pow(mm3,-6) );

      h210px6 +=  + ln1 * (  - 1./462.*pow(mm2,-5)*pow(mm3,8) - 113./
         5544.*pow(mm2,-4)*pow(mm3,7) - 2347./27720.*pow(mm2,-3)*pow(
         mm3,6) - 1879./9240.*pow(mm2,-2)*pow(mm3,5) - 157./504.*pow(
         mm2,-1)*pow(mm3,4) - 20./63.*pow(mm3,3) - 3./14.*mm2*pow(
         mm3,2) - 29./420.*pow(mm2,2)*mm3 + 31./420.*pow(mm2,3) + 31./
         168.*pow(mm2,4)*pow(mm3,-1) + 107./504.*pow(mm2,5)*pow(mm3,-2)
          + 835./5544.*pow(mm2,6)*pow(mm3,-3) + 611./9240.*pow(mm2,7)*
         pow(mm3,-4) + 229./13860.*pow(mm2,8)*pow(mm3,-5) + 5./2772.*
         pow(mm2,9)*pow(mm3,-6) );

      h210px6 +=  + 1./231.*pow(mm2,-4)*pow(mm3,7) + 107./2772.*pow(
         mm2,-3)*pow(mm3,6) + 379./2520.*pow(mm2,-2)*pow(mm3,5) + 27821.
         /83160.*pow(mm2,-1)*pow(mm3,4) + 18889./41580.*pow(mm3,3) + 
         136./385.*mm2*pow(mm3,2) + 5587./83160.*pow(mm2,2)*mm3 - 17033.
         /83160.*pow(mm2,3) - 1076./3465.*pow(mm2,4)*pow(mm3,-1) - 2549.
         /10395.*pow(mm2,5)*pow(mm3,-2) - 221./1890.*pow(mm2,6)*pow(
         mm3,-3) - 433./13860.*pow(mm2,7)*pow(mm3,-4) - 5./1386.*pow(
         mm2,8)*pow(mm3,-5);
      return h210px6*pi162/pow(mm1,7);
    case 7:
      h210px7 =
       + ln3 * ( 5./2772.*pow(mm2,-6)*pow(mm3,9) + 41./1980.*pow(
         mm2,-5)*pow(mm3,8) + 1481./13860.*pow(mm2,-4)*pow(mm3,7) + 
         4511./13860.*pow(mm2,-3)*pow(mm3,6) + 8957./13860.*pow(mm2,-2)
         *pow(mm3,5) + 221./252.*pow(mm2,-1)*pow(mm3,4) + 1037./1260.*
         pow(mm3,3) + 31./60.*mm2*pow(mm3,2) + 41./210.*pow(mm2,2)*mm3
          + 1./30.*pow(mm2,3) );

      h210px7 +=  + ln2 * ( 1./30.*pow(mm3,3) + 41./210.*mm2*pow(mm3,2)
          + 31./60.*pow(mm2,2)*mm3 + 1037./1260.*pow(mm2,3) + 221./252.
         *pow(mm2,4)*pow(mm3,-1) + 8957./13860.*pow(mm2,5)*pow(mm3,-2)
          + 4511./13860.*pow(mm2,6)*pow(mm3,-3) + 1481./13860.*pow(
         mm2,7)*pow(mm3,-4) + 41./1980.*pow(mm2,8)*pow(mm3,-5) + 5./
         2772.*pow(mm2,9)*pow(mm3,-6) );

      h210px7 +=  + ln1 * (  - 5./2772.*pow(mm2,-6)*pow(mm3,9) - 41./
         1980.*pow(mm2,-5)*pow(mm3,8) - 1481./13860.*pow(mm2,-4)*pow(
         mm3,7) - 4511./13860.*pow(mm2,-3)*pow(mm3,6) - 8957./13860.*
         pow(mm2,-2)*pow(mm3,5) - 221./252.*pow(mm2,-1)*pow(mm3,4) - 
         1079./1260.*pow(mm3,3) - 299./420.*mm2*pow(mm3,2) - 299./420.*
         pow(mm2,2)*mm3 - 1079./1260.*pow(mm2,3) - 221./252.*pow(mm2,4)
         *pow(mm3,-1) - 8957./13860.*pow(mm2,5)*pow(mm3,-2) - 4511./
         13860.*pow(mm2,6)*pow(mm3,-3) - 1481./13860.*pow(mm2,7)*pow(
         mm3,-4) - 41./1980.*pow(mm2,8)*pow(mm3,-5) - 5./2772.*pow(
         mm2,9)*pow(mm3,-6) );

      h210px7 +=  + 5./1386.*pow(mm2,-5)*pow(mm3,8) + 61./1540.*pow(
         mm2,-4)*pow(mm3,7) + 1615./8316.*pow(mm2,-3)*pow(mm3,6) + 6617.
         /11880.*pow(mm2,-2)*pow(mm3,5) + 85541./83160.*pow(mm2,-1)*
         pow(mm3,4) + 36623./27720.*pow(mm3,3) + 113119./83160.*mm2*
         pow(mm3,2) + 113119./83160.*pow(mm2,2)*mm3 + 36623./27720.*
         pow(mm2,3) + 85541./83160.*pow(mm2,4)*pow(mm3,-1) + 6617./
         11880.*pow(mm2,5)*pow(mm3,-2) + 1615./8316.*pow(mm2,6)*pow(
         mm3,-3) + 61./1540.*pow(mm2,7)*pow(mm3,-4) + 5./1386.*pow(
         mm2,8)*pow(mm3,-5);
      return h210px7*pi162/pow(mm1,7);
    case 8 :
      h210px8 =
       + ln3 * (  - 5./2002.*pow(mm2,-7)*pow(mm3,10) - 505./18018.*pow(
         mm2,-6)*pow(mm3,9) - 712./5005.*pow(mm2,-5)*pow(mm3,8) - 28./
         65.*pow(mm2,-4)*pow(mm3,7) - 38924./45045.*pow(mm2,-3)*pow(
         mm3,6) - 464./385.*pow(mm2,-2)*pow(mm3,5) - 92./77.*pow(
         mm2,-1)*pow(mm3,4) - 38./45.*pow(mm3,3) - 29./70.*mm2*pow(
         mm3,2) - 9./70.*pow(mm2,2)*mm3 - 2./105.*pow(mm2,3) );

      h210px8 +=  + ln2 * (  - 2./105.*pow(mm3,3) - 9./70.*mm2*pow(
         mm3,2) - 29./70.*pow(mm2,2)*mm3 - 38./45.*pow(mm2,3) - 92./77.
         *pow(mm2,4)*pow(mm3,-1) - 464./385.*pow(mm2,5)*pow(mm3,-2) - 
         38924./45045.*pow(mm2,6)*pow(mm3,-3) - 28./65.*pow(mm2,7)*pow(
         mm3,-4) - 712./5005.*pow(mm2,8)*pow(mm3,-5) - 505./18018.*pow(
         mm2,9)*pow(mm3,-6) - 5./2002.*pow(mm2,10)*pow(mm3,-7) );

      h210px8 +=  + ln1 * ( 5./2002.*pow(mm2,-7)*pow(mm3,10) + 505./
         18018.*pow(mm2,-6)*pow(mm3,9) + 712./5005.*pow(mm2,-5)*pow(
         mm3,8) + 28./65.*pow(mm2,-4)*pow(mm3,7) + 38924./45045.*pow(
         mm2,-3)*pow(mm3,6) + 464./385.*pow(mm2,-2)*pow(mm3,5) + 92./77.
         *pow(mm2,-1)*pow(mm3,4) + 272./315.*pow(mm3,3) + 19./35.*mm2*
         pow(mm3,2) + 19./35.*pow(mm2,2)*mm3 + 272./315.*pow(mm2,3) + 
         92./77.*pow(mm2,4)*pow(mm3,-1) + 464./385.*pow(mm2,5)*pow(
         mm3,-2) + 38924./45045.*pow(mm2,6)*pow(mm3,-3) + 28./65.*pow(
         mm2,7)*pow(mm3,-4) + 712./5005.*pow(mm2,8)*pow(mm3,-5) + 505./
         18018.*pow(mm2,9)*pow(mm3,-6) + 5./2002.*pow(mm2,10)*pow(
         mm3,-7) );

      h210px8 +=  - 5./1001.*pow(mm2,-6)*pow(mm3,9) - 965./18018.*pow(
         mm2,-5)*pow(mm3,8) - 1789./6930.*pow(mm2,-4)*pow(mm3,7) - 
         79645./108108.*pow(mm2,-3)*pow(mm3,6) - 248519./180180.*pow(
         mm2,-2)*pow(mm3,5) - 159703./90090.*pow(mm2,-1)*pow(mm3,4) - 
         146947./90090.*pow(mm3,3) - 38219./30030.*mm2*pow(mm3,2) - 
         38219./30030.*pow(mm2,2)*mm3 - 146947./90090.*pow(mm2,3) - 
         159703./90090.*pow(mm2,4)*pow(mm3,-1) - 248519./180180.*pow(
         mm2,5)*pow(mm3,-2) - 79645./108108.*pow(mm2,6)*pow(mm3,-3) - 
         1789./6930.*pow(mm2,7)*pow(mm3,-4) - 965./18018.*pow(mm2,8)*
         pow(mm3,-5) - 5./1001.*pow(mm2,9)*pow(mm3,-6);
      return h210px8*pi162/pow(mm1,9);
    default:
      std::cout << "h210psing wrong iprop ="<<iprop<<'\n';
    }
  }
  // case 2
  if (fabs(mm1+mm3-mm2) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210px1 =
       + ln3 * (  - 1./28.*pow(mm1,-1)*pow(mm3,4) - 1./12.*pow(mm3,3)
          - 1./20.*mm1*pow(mm3,2) );

      h210px1 +=  + ln2 * ( 1./28.*pow(mm1,-1)*pow(mm3,4) + 1./6.*pow(
         mm3,3) + 3./10.*mm1*pow(mm3,2) + 1./4.*pow(mm1,2)*mm3 + 1./12.
         *pow(mm1,3) + 1./420.*pow(mm1,6)*pow(mm3,-3) );

      h210px1 +=  + ln1 * (  - 1./420.*pow(mm1,6)*pow(mm3,-3) );

      h210px1 +=  + 143./2016.*pow(mm3,3) + 143./672.*mm1*pow(mm3,2) + 
         739./3360.*pow(mm1,2)*mm3 + 859./10080.*pow(mm1,3) + 1./420.*
         pow(mm1,4)*pow(mm3,-1) - 1./210.*pow(mm1,5)*pow(mm3,-2);
      return h210px1*pi162/pow(mm2,3);
    case 2:
      h210px2 =
       + ln3 * ( 1./252.*pow(mm1,-3)*pow(mm3,6) + 11./504.*pow(mm1,-2)*
         pow(mm3,5) + 13./280.*pow(mm1,-1)*pow(mm3,4) + 19./420.*pow(
         mm3,3) + 1./60.*mm1*pow(mm3,2) );

      h210px2 +=  + ln2 * (  - 1./252.*pow(mm1,-3)*pow(mm3,6) - 11./504.
         *pow(mm1,-2)*pow(mm3,5) - 13./280.*pow(mm1,-1)*pow(mm3,4) - 19.
         /420.*pow(mm3,3) - 1./60.*mm1*pow(mm3,2) + 1./140.*pow(mm1,4)*
         pow(mm3,-1) + 1./70.*pow(mm1,5)*pow(mm3,-2) + 23./2520.*pow(
         mm1,6)*pow(mm3,-3) + 1./504.*pow(mm1,7)*pow(mm3,-4) );

      h210px2 +=  + ln1 * (  - 1./140.*pow(mm1,4)*pow(mm3,-1) - 1./70.*
         pow(mm1,5)*pow(mm3,-2) - 23./2520.*pow(mm1,6)*pow(mm3,-3) - 1./
         504.*pow(mm1,7)*pow(mm3,-4) );

      h210px2 +=  + 1./126.*pow(mm1,-2)*pow(mm3,5) + 5./126.*pow(
         mm1,-1)*pow(mm3,4) + 331./3780.*pow(mm3,3) + 106./945.*mm1*
         pow(mm3,2) + 29./360.*pow(mm1,2)*mm3 + 17./1080.*pow(mm1,3) - 
         157./7560.*pow(mm1,4)*pow(mm3,-1) - 41./2520.*pow(mm1,5)*pow(
         mm3,-2) - 1./252.*pow(mm1,6)*pow(mm3,-3);
      return h210px2*pi162/pow(mm2,5);
    case 3:
      h210px3 =
       + ln3 * ( 1./72.*pow(mm1,-2)*pow(mm3,5) + 3./56.*pow(mm1,-1)*
         pow(mm3,4) + 1./14.*pow(mm3,3) + 1./30.*mm1*pow(mm3,2) );

      h210px3 +=  + ln2 * (  - 1./72.*pow(mm1,-2)*pow(mm3,5) - 3./56.*
         pow(mm1,-1)*pow(mm3,4) - 1./14.*pow(mm3,3) - 1./30.*mm1*pow(
         mm3,2) - 1./280.*pow(mm1,6)*pow(mm3,-3) - 1./504.*pow(mm1,7)*
         pow(mm3,-4) );

      h210px3 +=  + ln1 * ( 1./280.*pow(mm1,6)*pow(mm3,-3) + 1./504.*
         pow(mm1,7)*pow(mm3,-4) );

      h210px3 +=  + 1./36.*pow(mm1,-1)*pow(mm3,4) + 17./126.*pow(mm3,3)
          + 95./378.*mm1*pow(mm3,2) + 533./2520.*pow(mm1,2)*mm3 + 23./
         360.*pow(mm1,3) - 17./7560.*pow(mm1,4)*pow(mm3,-1) + 13./2520.
         *pow(mm1,5)*pow(mm3,-2) + 1./252.*pow(mm1,6)*pow(mm3,-3);
      return h210px3*pi162/pow(mm2,5);
    case 4:
      h210px4 =
       + ln3 * (  - 1./72.*pow(mm1,-2)*pow(mm3,5) - 43./504.*pow(
         mm1,-1)*pow(mm3,4) - 3./14.*pow(mm3,3) - 29./105.*mm1*pow(
         mm3,2) - 11./60.*pow(mm1,2)*mm3 - 1./20.*pow(mm1,3) );

      h210px4 +=  + ln2 * ( 1./72.*pow(mm1,-2)*pow(mm3,5) + 43./504.*
         pow(mm1,-1)*pow(mm3,4) + 3./14.*pow(mm3,3) + 29./105.*mm1*pow(
         mm3,2) + 11./60.*pow(mm1,2)*mm3 + 1./20.*pow(mm1,3) - 1./280.*
         pow(mm1,6)*pow(mm3,-3) - 13./2520.*pow(mm1,7)*pow(mm3,-4) - 1./
         630.*pow(mm1,8)*pow(mm3,-5) );

      h210px4 +=  + ln1 * ( 1./280.*pow(mm1,6)*pow(mm3,-3) + 13./2520.*
         pow(mm1,7)*pow(mm3,-4) + 1./630.*pow(mm1,8)*pow(mm3,-5) );

      h210px4 +=  - 1./36.*pow(mm1,-1)*pow(mm3,4) - 29./252.*pow(mm3,3)
          - 65./378.*mm1*pow(mm3,2) - 781./7560.*pow(mm1,2)*mm3 - 43./
         2520.*pow(mm1,3) - 1./1080.*pow(mm1,4)*pow(mm3,-1) + 23./7560.
         *pow(mm1,5)*pow(mm3,-2) + 11./1260.*pow(mm1,6)*pow(mm3,-3) + 1.
         /315.*pow(mm1,7)*pow(mm3,-4);
      return h210px4*pi162/pow(mm2,5);
    case 5:
      h210px5 =
       + ln3 * (  - 1./264.*pow(mm1,-4)*pow(mm3,7) - 131./5544.*pow(
         mm1,-3)*pow(mm3,6) - 31./504.*pow(mm1,-2)*pow(mm3,5) - 71./840.
         *pow(mm1,-1)*pow(mm3,4) - 13./210.*pow(mm3,3) - 2./105.*mm1*
         pow(mm3,2) );

      h210px5 +=  + ln2 * ( 1./264.*pow(mm1,-4)*pow(mm3,7) + 131./5544.
         *pow(mm1,-3)*pow(mm3,6) + 31./504.*pow(mm1,-2)*pow(mm3,5) + 71.
         /840.*pow(mm1,-1)*pow(mm3,4) + 13./210.*pow(mm3,3) + 2./105.*
         mm1*pow(mm3,2) - 3./280.*pow(mm1,4)*pow(mm3,-1) - 23./840.*
         pow(mm1,5)*pow(mm3,-2) - 67./2520.*pow(mm1,6)*pow(mm3,-3) - 67.
         /5544.*pow(mm1,7)*pow(mm3,-4) - 1./462.*pow(mm1,8)*pow(mm3,-5)
          );

      h210px5 +=  + ln1 * ( 3./280.*pow(mm1,4)*pow(mm3,-1) + 23./840.*
         pow(mm1,5)*pow(mm3,-2) + 67./2520.*pow(mm1,6)*pow(mm3,-3) + 67.
         /5544.*pow(mm1,7)*pow(mm3,-4) + 1./462.*pow(mm1,8)*pow(mm3,-5)
          );

      h210px5 +=  - 1./132.*pow(mm1,-3)*pow(mm3,6) - 241./5544.*pow(
         mm1,-2)*pow(mm3,5) - 565./5544.*pow(mm1,-1)*pow(mm3,4) - 10673.
         /83160.*pow(mm3,3) - 1211./11880.*mm1*pow(mm3,2) - 221./3960.*
         pow(mm1,2)*mm3 - 491./83160.*pow(mm1,3) + 2923./83160.*pow(
         mm1,4)*pow(mm3,-1) + 131./3080.*pow(mm1,5)*pow(mm3,-2) + 61./
         2772.*pow(mm1,6)*pow(mm3,-3) + 1./231.*pow(mm1,7)*pow(mm3,-4);
      return h210px5*pi162/pow(mm2,7);
    case 6:
      h210px6 =
       + ln3 * ( 1./264.*pow(mm1,-4)*pow(mm3,7) + 163./5544.*pow(
         mm1,-3)*pow(mm3,6) + 61./616.*pow(mm1,-2)*pow(mm3,5) + 53./280.
         *pow(mm1,-1)*pow(mm3,4) + 71./315.*pow(mm3,3) + 6./35.*mm1*
         pow(mm3,2) + 11./140.*pow(mm1,2)*mm3 + 1./60.*pow(mm1,3) );

      h210px6 +=  + ln2 * (  - 1./264.*pow(mm1,-4)*pow(mm3,7) - 163./
         5544.*pow(mm1,-3)*pow(mm3,6) - 61./616.*pow(mm1,-2)*pow(mm3,5)
          - 53./280.*pow(mm1,-1)*pow(mm3,4) - 71./315.*pow(mm3,3) - 6./
         35.*mm1*pow(mm3,2) - 11./140.*pow(mm1,2)*mm3 - 1./60.*pow(
         mm1,3) - 3./280.*pow(mm1,4)*pow(mm3,-1) - 31./840.*pow(mm1,5)*
         pow(mm3,-2) - 127./2520.*pow(mm1,6)*pow(mm3,-3) - 107./3080.*
         pow(mm1,7)*pow(mm3,-4) - 19./1540.*pow(mm1,8)*pow(mm3,-5) - 5./
         2772.*pow(mm1,9)*pow(mm3,-6) );

      h210px6 +=  + ln1 * ( 3./280.*pow(mm1,4)*pow(mm3,-1) + 31./840.*
         pow(mm1,5)*pow(mm3,-2) + 127./2520.*pow(mm1,6)*pow(mm3,-3) + 
         107./3080.*pow(mm1,7)*pow(mm3,-4) + 19./1540.*pow(mm1,8)*pow(
         mm3,-5) + 5./2772.*pow(mm1,9)*pow(mm3,-6) );

      h210px6 +=  + 1./132.*pow(mm1,-3)*pow(mm3,6) + 305./5544.*pow(
         mm1,-2)*pow(mm3,5) + 949./5544.*pow(mm1,-1)*pow(mm3,4) + 12071.
         /41580.*pow(mm3,3) + 1279./4620.*mm1*pow(mm3,2) + 3709./27720.
         *pow(mm1,2)*mm3 + 505./16632.*pow(mm1,3) + 19./462.*pow(mm1,4)
         *pow(mm3,-1) + 113./1540.*pow(mm1,5)*pow(mm3,-2) + 1213./20790.
         *pow(mm1,6)*pow(mm3,-3) + 317./13860.*pow(mm1,7)*pow(mm3,-4)
          + 5./1386.*pow(mm1,8)*pow(mm3,-5);
      return h210px6*pi162/pow(mm2,7);
    case 7:
      h210px7 =
       + ln3 * ( 1./99.*pow(mm1,-3)*pow(mm3,6) + 13./198.*pow(mm1,-2)*
         pow(mm3,5) + 23./126.*pow(mm1,-1)*pow(mm3,4) + 71./252.*pow(
         mm3,3) + 109./420.*mm1*pow(mm3,2) + 29./210.*pow(mm1,2)*mm3 + 
         1./30.*pow(mm1,3) );

      h210px7 +=  + ln2 * (  - 1./99.*pow(mm1,-3)*pow(mm3,6) - 13./198.
         *pow(mm1,-2)*pow(mm3,5) - 23./126.*pow(mm1,-1)*pow(mm3,4) - 71.
         /252.*pow(mm3,3) - 109./420.*mm1*pow(mm3,2) - 29./210.*pow(
         mm1,2)*mm3 - 1./30.*pow(mm1,3) + 2./315.*pow(mm1,6)*pow(
         mm3,-3) + 4./315.*pow(mm1,7)*pow(mm3,-4) + 113./13860.*pow(
         mm1,8)*pow(mm3,-5) + 5./2772.*pow(mm1,9)*pow(mm3,-6) );

      h210px7 +=  + ln1 * (  - 2./315.*pow(mm1,6)*pow(mm3,-3) - 4./315.
         *pow(mm1,7)*pow(mm3,-4) - 113./13860.*pow(mm1,8)*pow(mm3,-5)
          - 5./2772.*pow(mm1,9)*pow(mm3,-6) );

      h210px7 +=  + 2./99.*pow(mm1,-2)*pow(mm3,5) + 4./33.*pow(mm1,-1)*
         pow(mm3,4) + 4861./16632.*pow(mm3,3) + 5825./16632.*mm1*pow(
         mm3,2) + 8417./41580.*pow(mm1,2)*mm3 + 1733./41580.*pow(mm1,3)
          + 103./83160.*pow(mm1,4)*pow(mm3,-1) - 377./83160.*pow(mm1,5)
         *pow(mm3,-2) - 767./41580.*pow(mm1,6)*pow(mm3,-3) - 67./4620.*
         pow(mm1,7)*pow(mm3,-4) - 5./1386.*pow(mm1,8)*pow(mm3,-5);
      return h210px7*pi162/pow(mm2,7);
    case 8:
      h210px8 =
       + ln3 * (  - 2./429.*pow(mm1,-5)*pow(mm3,8) - 17./429.*pow(
         mm1,-4)*pow(mm3,7) - 122./819.*pow(mm1,-3)*pow(mm3,6) - 25./77.
         *pow(mm1,-2)*pow(mm3,5) - 174./385.*pow(mm1,-1)*pow(mm3,4) - 
         263./630.*pow(mm3,3) - 9./35.*mm1*pow(mm3,2) - 1./10.*pow(
         mm1,2)*mm3 - 2./105.*pow(mm1,3) );

      h210px8 +=  + ln2 * ( 2./429.*pow(mm1,-5)*pow(mm3,8) + 17./429.*
         pow(mm1,-4)*pow(mm3,7) + 122./819.*pow(mm1,-3)*pow(mm3,6) + 25.
         /77.*pow(mm1,-2)*pow(mm3,5) + 174./385.*pow(mm1,-1)*pow(mm3,4)
          + 263./630.*pow(mm3,3) + 9./35.*mm1*pow(mm3,2) + 1./10.*pow(
         mm1,2)*mm3 + 2./105.*pow(mm1,3) + 2./105.*pow(mm1,4)*pow(
         mm3,-1) + 8./105.*pow(mm1,5)*pow(mm3,-2) + 446./3465.*pow(
         mm1,6)*pow(mm3,-3) + 46./385.*pow(mm1,7)*pow(mm3,-4) + 59./910.
         *pow(mm1,8)*pow(mm3,-5) + 25./1287.*pow(mm1,9)*pow(mm3,-6) + 5.
         /2002.*pow(mm1,10)*pow(mm3,-7) );

      h210px8 +=  + ln1 * (  - 2./105.*pow(mm1,4)*pow(mm3,-1) - 8./105.
         *pow(mm1,5)*pow(mm3,-2) - 446./3465.*pow(mm1,6)*pow(mm3,-3) - 
         46./385.*pow(mm1,7)*pow(mm3,-4) - 59./910.*pow(mm1,8)*pow(
         mm3,-5) - 25./1287.*pow(mm1,9)*pow(mm3,-6) - 5./2002.*pow(
         mm1,10)*pow(mm3,-7) );

      h210px8 +=  - 4./429.*pow(mm1,-4)*pow(mm3,7) - 32./429.*pow(
         mm1,-3)*pow(mm3,6) - 785./3003.*pow(mm1,-2)*pow(mm3,5) - 75./
         143.*pow(mm1,-1)*pow(mm3,4) - 50581./77220.*pow(mm3,3) - 6511./
         12870.*mm1*pow(mm3,2) - 38639./180180.*pow(mm1,2)*mm3 - 1240./
         27027.*pow(mm1,3) - 85./1092.*pow(mm1,4)*pow(mm3,-1) - 5179./
         30030.*pow(mm1,5)*pow(mm3,-2) - 100447./540540.*pow(mm1,6)*
         pow(mm3,-3) - 5041./45045.*pow(mm1,7)*pow(mm3,-4) - 655./18018.
         *pow(mm1,8)*pow(mm3,-5) - 5./1001.*pow(mm1,9)*pow(mm3,-6);
      return h210px8*pi162/pow(mm2,9);
    default:
      std::cout << "h210psing wrong iprop ="<<iprop<<'\n';
    }
  }
  //case 3
  if (fabs(mm1+mm2-mm3) < 1e-6*fabs(mm1+mm2+mm3)){
    switch(iprop){
    case 1:
      h210px1 =
       + ln3 * ( 1./28.*pow(mm1,-1)*pow(mm2,4) + 1./6.*pow(mm2,3) + 3./
         10.*mm1*pow(mm2,2) + 1./4.*pow(mm1,2)*mm2 + 1./12.*pow(mm1,3)
          + 1./420.*pow(mm1,6)*pow(mm2,-3) );

      h210px1 +=  + ln2 * (  - 1./28.*pow(mm1,-1)*pow(mm2,4) - 1./12.*
         pow(mm2,3) - 1./20.*mm1*pow(mm2,2) );

      h210px1 +=  + ln1 * (  - 1./420.*pow(mm1,6)*pow(mm2,-3) );

      h210px1 +=  + 143./2016.*pow(mm2,3) + 143./672.*mm1*pow(mm2,2) + 
         739./3360.*pow(mm1,2)*mm2 + 859./10080.*pow(mm1,3) + 1./420.*
         pow(mm1,4)*pow(mm2,-1) - 1./210.*pow(mm1,5)*pow(mm2,-2);
      return h210px1*pi162/pow(mm3,3);
    case 2:
   h210px2 =
       + ln3 * (  - 1./252.*pow(mm1,-3)*pow(mm2,6) - 11./504.*pow(
         mm1,-2)*pow(mm2,5) - 13./280.*pow(mm1,-1)*pow(mm2,4) - 19./420.
         *pow(mm2,3) - 1./60.*mm1*pow(mm2,2) + 1./140.*pow(mm1,4)*pow(
         mm2,-1) + 1./70.*pow(mm1,5)*pow(mm2,-2) + 23./2520.*pow(mm1,6)
         *pow(mm2,-3) + 1./504.*pow(mm1,7)*pow(mm2,-4) );

      h210px2 +=  + ln2 * ( 1./252.*pow(mm1,-3)*pow(mm2,6) + 11./504.*
         pow(mm1,-2)*pow(mm2,5) + 13./280.*pow(mm1,-1)*pow(mm2,4) + 19./
         420.*pow(mm2,3) + 1./60.*mm1*pow(mm2,2) );

      h210px2 +=  + ln1 * (  - 1./140.*pow(mm1,4)*pow(mm2,-1) - 1./70.*
         pow(mm1,5)*pow(mm2,-2) - 23./2520.*pow(mm1,6)*pow(mm2,-3) - 1./
         504.*pow(mm1,7)*pow(mm2,-4) );

      h210px2 +=  + 1./126.*pow(mm1,-2)*pow(mm2,5) + 5./126.*pow(
         mm1,-1)*pow(mm2,4) + 331./3780.*pow(mm2,3) + 106./945.*mm1*
         pow(mm2,2) + 29./360.*pow(mm1,2)*mm2 + 17./1080.*pow(mm1,3) - 
         157./7560.*pow(mm1,4)*pow(mm2,-1) - 41./2520.*pow(mm1,5)*pow(
         mm2,-2) - 1./252.*pow(mm1,6)*pow(mm2,-3);
      return h210px2*pi162/pow(mm3,5);
    case 3:
      h210px3 =
       + ln3 * ( 1./72.*pow(mm1,-2)*pow(mm2,5) + 43./504.*pow(mm1,-1)*
         pow(mm2,4) + 3./14.*pow(mm2,3) + 29./105.*mm1*pow(mm2,2) + 11./
         60.*pow(mm1,2)*mm2 + 1./20.*pow(mm1,3) - 1./280.*pow(mm1,6)*
         pow(mm2,-3) - 13./2520.*pow(mm1,7)*pow(mm2,-4) - 1./630.*pow(
         mm1,8)*pow(mm2,-5) );

      h210px3 +=  + ln2 * (  - 1./72.*pow(mm1,-2)*pow(mm2,5) - 43./504.
         *pow(mm1,-1)*pow(mm2,4) - 3./14.*pow(mm2,3) - 29./105.*mm1*
         pow(mm2,2) - 11./60.*pow(mm1,2)*mm2 - 1./20.*pow(mm1,3) );

      h210px3 +=  + ln1 * ( 1./280.*pow(mm1,6)*pow(mm2,-3) + 13./2520.*
         pow(mm1,7)*pow(mm2,-4) + 1./630.*pow(mm1,8)*pow(mm2,-5) );

      h210px3 +=  - 1./36.*pow(mm1,-1)*pow(mm2,4) - 29./252.*pow(mm2,3)
          - 65./378.*mm1*pow(mm2,2) - 781./7560.*pow(mm1,2)*mm2 - 43./
         2520.*pow(mm1,3) - 1./1080.*pow(mm1,4)*pow(mm2,-1) + 23./7560.
         *pow(mm1,5)*pow(mm2,-2) + 11./1260.*pow(mm1,6)*pow(mm2,-3) + 1.
         /315.*pow(mm1,7)*pow(mm2,-4);
      return h210px3*pi162/pow(mm3,5);
    case 4:
      h210px4 =
       + ln3 * (  - 1./72.*pow(mm1,-2)*pow(mm2,5) - 3./56.*pow(mm1,-1)*
         pow(mm2,4) - 1./14.*pow(mm2,3) - 1./30.*mm1*pow(mm2,2) - 1./
         280.*pow(mm1,6)*pow(mm2,-3) - 1./504.*pow(mm1,7)*pow(mm2,-4) )
         ;

      h210px4 +=  + ln2 * ( 1./72.*pow(mm1,-2)*pow(mm2,5) + 3./56.*pow(
         mm1,-1)*pow(mm2,4) + 1./14.*pow(mm2,3) + 1./30.*mm1*pow(mm2,2)
          );

      h210px4 +=  + ln1 * ( 1./280.*pow(mm1,6)*pow(mm2,-3) + 1./504.*
         pow(mm1,7)*pow(mm2,-4) );

      h210px4 +=  + 1./36.*pow(mm1,-1)*pow(mm2,4) + 17./126.*pow(mm2,3)
          + 95./378.*mm1*pow(mm2,2) + 533./2520.*pow(mm1,2)*mm2 + 23./
         360.*pow(mm1,3) - 17./7560.*pow(mm1,4)*pow(mm2,-1) + 13./2520.
         *pow(mm1,5)*pow(mm2,-2) + 1./252.*pow(mm1,6)*pow(mm2,-3);
      return h210px4*pi162/pow(mm3,5);
    case 5:
      h210px5 =
       + ln3 * (  - 1./264.*pow(mm1,-4)*pow(mm2,7) - 163./5544.*pow(
         mm1,-3)*pow(mm2,6) - 61./616.*pow(mm1,-2)*pow(mm2,5) - 53./280.
         *pow(mm1,-1)*pow(mm2,4) - 71./315.*pow(mm2,3) - 6./35.*mm1*
         pow(mm2,2) - 11./140.*pow(mm1,2)*mm2 - 1./60.*pow(mm1,3) - 3./
         280.*pow(mm1,4)*pow(mm2,-1) - 31./840.*pow(mm1,5)*pow(mm2,-2)
          - 127./2520.*pow(mm1,6)*pow(mm2,-3) - 107./3080.*pow(mm1,7)*
         pow(mm2,-4) - 19./1540.*pow(mm1,8)*pow(mm2,-5) - 5./2772.*pow(
         mm1,9)*pow(mm2,-6) );

      h210px5 +=  + ln2 * ( 1./264.*pow(mm1,-4)*pow(mm2,7) + 163./5544.
         *pow(mm1,-3)*pow(mm2,6) + 61./616.*pow(mm1,-2)*pow(mm2,5) + 53.
         /280.*pow(mm1,-1)*pow(mm2,4) + 71./315.*pow(mm2,3) + 6./35.*
         mm1*pow(mm2,2) + 11./140.*pow(mm1,2)*mm2 + 1./60.*pow(mm1,3) )
         ;

      h210px5 +=  + ln1 * ( 3./280.*pow(mm1,4)*pow(mm2,-1) + 31./840.*
         pow(mm1,5)*pow(mm2,-2) + 127./2520.*pow(mm1,6)*pow(mm2,-3) + 
         107./3080.*pow(mm1,7)*pow(mm2,-4) + 19./1540.*pow(mm1,8)*pow(
         mm2,-5) + 5./2772.*pow(mm1,9)*pow(mm2,-6) );

      h210px5 +=  + 1./132.*pow(mm1,-3)*pow(mm2,6) + 305./5544.*pow(
         mm1,-2)*pow(mm2,5) + 949./5544.*pow(mm1,-1)*pow(mm2,4) + 12071.
         /41580.*pow(mm2,3) + 1279./4620.*mm1*pow(mm2,2) + 3709./27720.
         *pow(mm1,2)*mm2 + 505./16632.*pow(mm1,3) + 19./462.*pow(mm1,4)
         *pow(mm2,-1) + 113./1540.*pow(mm1,5)*pow(mm2,-2) + 1213./20790.
         *pow(mm1,6)*pow(mm2,-3) + 317./13860.*pow(mm1,7)*pow(mm2,-4)
          + 5./1386.*pow(mm1,8)*pow(mm2,-5);
      return h210px5*pi162/pow(mm3,7);
    case 6:
      h210px6 =
       + ln3 * ( 1./264.*pow(mm1,-4)*pow(mm2,7) + 131./5544.*pow(
         mm1,-3)*pow(mm2,6) + 31./504.*pow(mm1,-2)*pow(mm2,5) + 71./840.
         *pow(mm1,-1)*pow(mm2,4) + 13./210.*pow(mm2,3) + 2./105.*mm1*
         pow(mm2,2) - 3./280.*pow(mm1,4)*pow(mm2,-1) - 23./840.*pow(
         mm1,5)*pow(mm2,-2) - 67./2520.*pow(mm1,6)*pow(mm2,-3) - 67./
         5544.*pow(mm1,7)*pow(mm2,-4) - 1./462.*pow(mm1,8)*pow(mm2,-5)
          );

      h210px6 +=  + ln2 * (  - 1./264.*pow(mm1,-4)*pow(mm2,7) - 131./
         5544.*pow(mm1,-3)*pow(mm2,6) - 31./504.*pow(mm1,-2)*pow(mm2,5)
          - 71./840.*pow(mm1,-1)*pow(mm2,4) - 13./210.*pow(mm2,3) - 2./
         105.*mm1*pow(mm2,2) );

      h210px6 +=  + ln1 * ( 3./280.*pow(mm1,4)*pow(mm2,-1) + 23./840.*
         pow(mm1,5)*pow(mm2,-2) + 67./2520.*pow(mm1,6)*pow(mm2,-3) + 67.
         /5544.*pow(mm1,7)*pow(mm2,-4) + 1./462.*pow(mm1,8)*pow(mm2,-5)
          );

      h210px6 +=  - 1./132.*pow(mm1,-3)*pow(mm2,6) - 241./5544.*pow(
         mm1,-2)*pow(mm2,5) - 565./5544.*pow(mm1,-1)*pow(mm2,4) - 10673.
         /83160.*pow(mm2,3) - 1211./11880.*mm1*pow(mm2,2) - 221./3960.*
         pow(mm1,2)*mm2 - 491./83160.*pow(mm1,3) + 2923./83160.*pow(
         mm1,4)*pow(mm2,-1) + 131./3080.*pow(mm1,5)*pow(mm2,-2) + 61./
         2772.*pow(mm1,6)*pow(mm2,-3) + 1./231.*pow(mm1,7)*pow(mm2,-4);
      return h210px6*pi162/pow(mm3,7);
    case 7:
      h210px7 =
       + ln3 * (  - 1./99.*pow(mm1,-3)*pow(mm2,6) - 13./198.*pow(
         mm1,-2)*pow(mm2,5) - 23./126.*pow(mm1,-1)*pow(mm2,4) - 71./252.
         *pow(mm2,3) - 109./420.*mm1*pow(mm2,2) - 29./210.*pow(mm1,2)*
         mm2 - 1./30.*pow(mm1,3) + 2./315.*pow(mm1,6)*pow(mm2,-3) + 4./
         315.*pow(mm1,7)*pow(mm2,-4) + 113./13860.*pow(mm1,8)*pow(
         mm2,-5) + 5./2772.*pow(mm1,9)*pow(mm2,-6) );

      h210px7 +=  + ln2 * ( 1./99.*pow(mm1,-3)*pow(mm2,6) + 13./198.*
         pow(mm1,-2)*pow(mm2,5) + 23./126.*pow(mm1,-1)*pow(mm2,4) + 71./
         252.*pow(mm2,3) + 109./420.*mm1*pow(mm2,2) + 29./210.*pow(
         mm1,2)*mm2 + 1./30.*pow(mm1,3) );

      h210px7 +=  + ln1 * (  - 2./315.*pow(mm1,6)*pow(mm2,-3) - 4./315.
         *pow(mm1,7)*pow(mm2,-4) - 113./13860.*pow(mm1,8)*pow(mm2,-5)
          - 5./2772.*pow(mm1,9)*pow(mm2,-6) );

      h210px7 +=  + 2./99.*pow(mm1,-2)*pow(mm2,5) + 4./33.*pow(mm1,-1)*
         pow(mm2,4) + 4861./16632.*pow(mm2,3) + 5825./16632.*mm1*pow(
         mm2,2) + 8417./41580.*pow(mm1,2)*mm2 + 1733./41580.*pow(mm1,3)
          + 103./83160.*pow(mm1,4)*pow(mm2,-1) - 377./83160.*pow(mm1,5)
         *pow(mm2,-2) - 767./41580.*pow(mm1,6)*pow(mm2,-3) - 67./4620.*
         pow(mm1,7)*pow(mm2,-4) - 5./1386.*pow(mm1,8)*pow(mm2,-5);
      return h210px7*pi162/pow(mm3,7);
    case 8:
      h210px8 =
       + ln3 * ( 2./429.*pow(mm1,-5)*pow(mm2,8) + 17./429.*pow(mm1,-4)*
         pow(mm2,7) + 122./819.*pow(mm1,-3)*pow(mm2,6) + 25./77.*pow(
         mm1,-2)*pow(mm2,5) + 174./385.*pow(mm1,-1)*pow(mm2,4) + 263./
         630.*pow(mm2,3) + 9./35.*mm1*pow(mm2,2) + 1./10.*pow(mm1,2)*
         mm2 + 2./105.*pow(mm1,3) + 2./105.*pow(mm1,4)*pow(mm2,-1) + 8./
         105.*pow(mm1,5)*pow(mm2,-2) + 446./3465.*pow(mm1,6)*pow(
         mm2,-3) + 46./385.*pow(mm1,7)*pow(mm2,-4) + 59./910.*pow(
         mm1,8)*pow(mm2,-5) + 25./1287.*pow(mm1,9)*pow(mm2,-6) + 5./
         2002.*pow(mm1,10)*pow(mm2,-7) );

      h210px8 +=  + ln2 * (  - 2./429.*pow(mm1,-5)*pow(mm2,8) - 17./429.
         *pow(mm1,-4)*pow(mm2,7) - 122./819.*pow(mm1,-3)*pow(mm2,6) - 
         25./77.*pow(mm1,-2)*pow(mm2,5) - 174./385.*pow(mm1,-1)*pow(
         mm2,4) - 263./630.*pow(mm2,3) - 9./35.*mm1*pow(mm2,2) - 1./10.
         *pow(mm1,2)*mm2 - 2./105.*pow(mm1,3) );

      h210px8 +=  + ln1 * (  - 2./105.*pow(mm1,4)*pow(mm2,-1) - 8./105.
         *pow(mm1,5)*pow(mm2,-2) - 446./3465.*pow(mm1,6)*pow(mm2,-3) - 
         46./385.*pow(mm1,7)*pow(mm2,-4) - 59./910.*pow(mm1,8)*pow(
         mm2,-5) - 25./1287.*pow(mm1,9)*pow(mm2,-6) - 5./2002.*pow(
         mm1,10)*pow(mm2,-7) );

      h210px8 +=  - 4./429.*pow(mm1,-4)*pow(mm2,7) - 32./429.*pow(
         mm1,-3)*pow(mm2,6) - 785./3003.*pow(mm1,-2)*pow(mm2,5) - 75./
         143.*pow(mm1,-1)*pow(mm2,4) - 50581./77220.*pow(mm2,3) - 6511./
         12870.*mm1*pow(mm2,2) - 38639./180180.*pow(mm1,2)*mm2 - 1240./
         27027.*pow(mm1,3) - 85./1092.*pow(mm1,4)*pow(mm2,-1) - 5179./
         30030.*pow(mm1,5)*pow(mm2,-2) - 100447./540540.*pow(mm1,6)*
         pow(mm2,-3) - 5041./45045.*pow(mm1,7)*pow(mm2,-4) - 655./18018.
         *pow(mm1,8)*pow(mm2,-5) - 5./1001.*pow(mm1,9)*pow(mm2,-6);
      return h210px8*pi162/pow(mm3,9);
    default:
      std::cout << "h210psing, wrong iprop = "<<iprop<<'\n';
      return 0.;
     }
  }
  std::cout << "funny masses in h210psing\n";
  return 0.;
}


//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// clausen's function in terms of li2
double  claus2q(const double x){
  const double dli2one =  pi2/6.;//1.64493406684823;
  dcomplex expt=exp(dcomplex(0.,x));
  return real(dcomplex(0.,-1.)*(jbdli2(expt)-dli2one)
	      +dcomplex(0.,-pi/2.*x+x*x/4.));
}

dcomplex zphiq(const double xm1,const double xm2,const double xm3){
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
    return 2./xlam/xmass2*(claus2q(2.*acos(x1))+claus2q(2.*acos(x2))
			   +claus2q(2.*acos(x3)));}
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
          std::cout << "crazy mass combination in phiq" << std::endl;
          assert(0);}//makes program stop
    }
    dcomplex zx1 = 0.5*(1.+xx-yy-xlam);
    dcomplex zx2 = 0.5*(1.-xx+yy-xlam);     
    return 1./xlam*(2.*log(zx1)*log(zx2)-2.*jbdli2(zx1)-2.*jbdli2(zx2)
		    -log(xx)*log(yy)+pi2/3.)/xmass2;
  }
};

double psiq(const double xm1,const double xm2,const double xm3){
  double xm12 = xm1*xm1;
  double xm22 = xm2*xm2;
  double xm32 = xm3*xm3;
  double xlam2 = pow(xm12-xm22-xm32,2)-4.*xm22*xm32;
  return  -xlam2*real(zphiq(xm1,xm2,xm3));
}

double h0(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq, const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h0sing(iprop,m1sq,m2sq,m3sq,xmu2);

  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double hh0;
  switch(iprop){
  case 1:
    hh0 = -0.5*psiq(m1,m2,m3)
      + m1sq*(pi2/12.+1.5-ln1+0.5*(-ln2*ln3+ln1*ln4))
      + m2sq*(pi2/12.+1.5-ln2+0.5*(-ln1*ln3+ln2*ln4))
      + m3sq*(pi2/12.+1.5-ln3+0.5*(-ln1*ln2+ln3*ln4));
    break;
  case 2:
   hh0 =
     + psiq(m1,m2,m3) * ( 1./2.*m3sq*lm1 + 1./2.*m2sq*lm1 - 1./2.*m1sq
         *lm1 );
   hh0 +=  + 1./2. - 1./2.*ln2*ln3 + ln1 + 1./2.*ln1*ln3 + 1./2.*
         ln1*ln2 + 1./2.*pow(ln1,2) + 1./12.*pi2;
   break;
  case 3:
   hh0 =
       + psiq(m1,m2,m3) * ( 1./2.*m3sq*lm1 - 1./2.*m2sq*lm1 + 1./2.*m1sq
         *lm1 );

      hh0 +=  + 1./2. + ln2 + 1./2.*ln2*ln3 + 1./2.*pow(ln2,2) - 1./2.
         *ln1*ln3 + 1./2.*ln1*ln2 + 1./12.*pi2;
      break;
  case 4:
   hh0 =
       + psiq(m1,m2,m3) * (  - 1./2.*m3sq*lm1 + 1./2.*m2sq*lm1 + 1./2.*
         m1sq*lm1 );

      hh0 +=  + 1./2. + ln3 + 1./2.*pow(ln3,2) + 1./2.*ln2*ln3 + 1./2.
         *ln1*ln3 - 1./2.*ln1*ln2 + 1./12.*pi2;
      break;
  case 5:
   hh0 =
       + psiq(m1,m2,m3) * ( 1./2.*lm1 + 1./2.*pow(m3sq,2)*pow(lm1,2) - 1.
         /2.*pow(m2sq,2)*pow(lm1,2) + m1sq*m2sq*pow(lm1,2) - 1./2.*pow(
         m1sq,2)*pow(lm1,2) );

      hh0 +=  - 1./2.*pow(m2sq,-1)*ln3 + 1./2.*pow(m2sq,-1)*ln1 + 1./2.
         *pow(m2sq,-1)*pow(m3sq,2)*ln3*lm1 - 1./2.*pow(m2sq,-1)*pow(
         m3sq,2)*ln1*lm1 + m3sq*ln3*lm1 - m3sq*ln2*lm1 + 1./2.*m2sq*ln3
         *lm1 - m2sq*ln2*lm1 + 1./2.*m2sq*ln1*lm1 - m1sq*pow(m2sq,-1)*
         m3sq*ln3*lm1 + m1sq*pow(m2sq,-1)*m3sq*ln1*lm1 - m1sq*ln3*lm1
          + m1sq*ln2*lm1 + 1./2.*pow(m1sq,2)*pow(m2sq,-1)*ln3*lm1 - 1./
         2.*pow(m1sq,2)*pow(m2sq,-1)*ln1*lm1;
      break;
  case 6:
   hh0 =
       + psiq(m1,m2,m3) * ( 1./2.*lm1 - 1./2.*pow(m3sq,2)*pow(lm1,2) + 1.
         /2.*pow(m2sq,2)*pow(lm1,2) + m1sq*m3sq*pow(lm1,2) - 1./2.*pow(
         m1sq,2)*pow(lm1,2) );

      hh0 +=  - 1./2.*pow(m3sq,-1)*ln2 + 1./2.*pow(m3sq,-1)*ln1 - m3sq
         *ln3*lm1 + 1./2.*m3sq*ln2*lm1 + 1./2.*m3sq*ln1*lm1 - m2sq*ln3*
         lm1 + m2sq*ln2*lm1 + 1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm1 - 
         1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm1 + m1sq*ln3*lm1 - m1sq*
         ln2*lm1 - m1sq*m2sq*pow(m3sq,-1)*ln2*lm1 + m1sq*m2sq*pow(
         m3sq,-1)*ln1*lm1 + 1./2.*pow(m1sq,2)*pow(m3sq,-1)*ln2*lm1 - 1./
         2.*pow(m1sq,2)*pow(m3sq,-1)*ln1*lm1;
      break;
  case 7:
   hh0 =
       + psiq(m1,m2,m3) * ( 1./2.*lm1 - 1./2.*pow(m3sq,2)*pow(lm1,2) + 
         m2sq*m3sq*pow(lm1,2) - 1./2.*pow(m2sq,2)*pow(lm1,2) + 1./2.*
         pow(m1sq,2)*pow(lm1,2) );

      hh0 +=  + 1./2.*pow(m3sq,-1)*ln2 - 1./2.*pow(m3sq,-1)*ln1 - m3sq
         *ln3*lm1 + 1./2.*m3sq*ln2*lm1 + 1./2.*m3sq*ln1*lm1 + m2sq*ln3*
         lm1 - m2sq*ln1*lm1 - 1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm1 + 
         1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm1 - m1sq*ln3*lm1 + m1sq*
         ln1*lm1 + m1sq*m2sq*pow(m3sq,-1)*ln2*lm1 - m1sq*m2sq*pow(
         m3sq,-1)*ln1*lm1 - 1./2.*pow(m1sq,2)*pow(m3sq,-1)*ln2*lm1 + 1./
         2.*pow(m1sq,2)*pow(m3sq,-1)*ln1*lm1;
      break;
  case 8:
   hh0 =
       + psiq(m1,m2,m3) * ( 1./2.*m3sq*pow(lm1,2) - 3./2.*pow(m3sq,3)*
         pow(lm1,3) + 1./2.*m2sq*pow(lm1,2) + 3./2.*m2sq*pow(m3sq,2)*
         pow(lm1,3) + 3./2.*pow(m2sq,2)*m3sq*pow(lm1,3) - 3./2.*pow(
         m2sq,3)*pow(lm1,3) + 1./2.*m1sq*pow(lm1,2) + 3./2.*m1sq*pow(
         m3sq,2)*pow(lm1,3) - 3*m1sq*m2sq*m3sq*pow(lm1,3) + 3./2.*m1sq*
         pow(m2sq,2)*pow(lm1,3) + 3./2.*pow(m1sq,2)*m3sq*pow(lm1,3) + 3.
         /2.*pow(m1sq,2)*m2sq*pow(lm1,3) - 3./2.*pow(m1sq,3)*pow(lm1,3)
          );

      hh0 +=  - 1./2.*pow(m2sq,-1)*pow(m3sq,-1) + 1./2.*pow(m2sq,-1)*
         m3sq*lm1 + pow(m2sq,-1)*m3sq*ln3*lm1 - pow(m2sq,-1)*m3sq*ln1*
         lm1 - pow(m2sq,-1)*pow(m3sq,3)*ln3*pow(lm1,2) + pow(m2sq,-1)*
         pow(m3sq,3)*ln1*pow(lm1,2) + lm1 - 1./2.*ln2*lm1 + 1./2.*ln1*
         lm1 - 2*pow(m3sq,2)*ln3*pow(lm1,2) + 5./2.*pow(m3sq,2)*ln2*
         pow(lm1,2) - 1./2.*pow(m3sq,2)*ln1*pow(lm1,2) + 1./2.*m2sq*
         pow(m3sq,-1)*lm1 + 1./2.*m2sq*pow(m3sq,-1)*ln2*lm1 - 1./2.*
         m2sq*pow(m3sq,-1)*ln1*lm1 + m2sq*m3sq*ln3*pow(lm1,2) + 1./2.*
         m2sq*m3sq*ln2*pow(lm1,2) - 3./2.*m2sq*m3sq*ln1*pow(lm1,2) + 2*
         pow(m2sq,2)*ln3*pow(lm1,2) - 5./2.*pow(m2sq,2)*ln2*pow(lm1,2)
          + 1./2.*pow(m2sq,2)*ln1*pow(lm1,2) - 1./2.*pow(m2sq,3)*pow(
         m3sq,-1)*ln2*pow(lm1,2) + 1./2.*pow(m2sq,3)*pow(m3sq,-1)*ln1*
         pow(lm1,2) - m1sq*pow(m2sq,-1)*lm1 - m1sq*pow(m2sq,-1)*ln3*lm1
          + m1sq*pow(m2sq,-1)*ln1*lm1 + 3*m1sq*pow(m2sq,-1)*pow(m3sq,2)
         *ln3*pow(lm1,2) - 3*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(
         lm1,2);
      hh0 +=  - m1sq*pow(m3sq,-1)*lm1 - 1./2.*m1sq*pow(m3sq,-1)*ln2*
         lm1 + 1./2.*m1sq*pow(m3sq,-1)*ln1*lm1 + 2*m1sq*m3sq*ln3*pow(
         lm1,2) - 9./2.*m1sq*m3sq*ln2*pow(lm1,2) + 5./2.*m1sq*m3sq*ln1*
         pow(lm1,2) - 3*m1sq*m2sq*ln3*pow(lm1,2) + m1sq*m2sq*ln2*pow(
         lm1,2) + 2*m1sq*m2sq*ln1*pow(lm1,2) + 3./2.*m1sq*pow(m2sq,2)*
         pow(m3sq,-1)*ln2*pow(lm1,2) - 3./2.*m1sq*pow(m2sq,2)*pow(
         m3sq,-1)*ln1*pow(lm1,2) + 1./2.*pow(m1sq,2)*pow(m2sq,-1)*pow(
         m3sq,-1)*lm1 - 3*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,2)
          + 3*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*pow(lm1,2) + 3./2.*pow(
         m1sq,2)*ln2*pow(lm1,2) - 3./2.*pow(m1sq,2)*ln1*pow(lm1,2) - 3./
         2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*pow(lm1,2) + 3./2.*pow(
         m1sq,2)*m2sq*pow(m3sq,-1)*ln1*pow(lm1,2) + pow(m1sq,3)*pow(
         m2sq,-1)*ln3*pow(lm1,2) - pow(m1sq,3)*pow(m2sq,-1)*ln1*pow(
         lm1,2) + 1./2.*pow(m1sq,3)*pow(m3sq,-1)*ln2*pow(lm1,2) - 1./2.
         *pow(m1sq,3)*pow(m3sq,-1)*ln1*pow(lm1,2);
      break;
     default:
    std::cout <<"wrong iprop in h0 (quenched case)\n";
    hh0=0;
    break;
  }
  return pi162*hh0;
}


double h0p(const int iprop, const double m1sq, const double m2sq,
	   const double m3sq, const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = pow(m1sq-m2sq-m3sq,2)-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h0psing(iprop,m1sq,m2sq,m3sq,xmu2);
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double hh0p=0.;
  switch(iprop){
  case 1:
   hh0p = m1sq*m2sq*m3sq/pow(lm,2)*psiq(m1,m2,m3)+1./8.
    + m1sq*dm*ln1/2.*lm1
    + m2sq/2.*lm1*(m2sq-m1sq-m3sq)*ln2
    + m3sq/2.*lm1*(m3sq-m1sq-m2sq)*ln3;
  break;
  case 2:
   hh0p =
       + psiq(m1,m2,m3) * ( m2sq*m3sq*pow(lm1,2) + 3*m1sq*m2sq*pow(
         m3sq,2)*pow(lm1,3) + 3*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) - 3*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) );

      hh0p +=  - 1./2.*m3sq*lm1 - 1./2.*m3sq*ln3*lm1 - 1./2.*m3sq*ln1*
         lm1 + pow(m3sq,3)*ln3*pow(lm1,2) - 1./2.*m2sq*lm1 - 1./2.*m2sq
         *ln2*lm1 - 1./2.*m2sq*ln1*lm1 + m2sq*pow(m3sq,2)*ln3*pow(
         lm1,2) - 2*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) - 2*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,2) + pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + pow(
         m2sq,3)*ln2*pow(lm1,2) + 1./2.*m1sq*lm1 + m1sq*ln1*lm1 - 2*
         m1sq*pow(m3sq,2)*ln3*pow(lm1,2) - m1sq*pow(m3sq,2)*ln1*pow(
         lm1,2) + m1sq*m2sq*m3sq*ln3*pow(lm1,2) + m1sq*m2sq*m3sq*ln2*
         pow(lm1,2) - 4*m1sq*m2sq*m3sq*ln1*pow(lm1,2) - 2*m1sq*pow(
         m2sq,2)*ln2*pow(lm1,2) - m1sq*pow(m2sq,2)*ln1*pow(lm1,2) + 
         pow(m1sq,2)*m3sq*ln3*pow(lm1,2) + 2*pow(m1sq,2)*m3sq*ln1*pow(
         lm1,2) + pow(m1sq,2)*m2sq*ln2*pow(lm1,2) + 2*pow(m1sq,2)*m2sq*
         ln1*pow(lm1,2) - pow(m1sq,3)*ln1*pow(lm1,2);
      break;
  case 3:
   hh0p =
       + psiq(m1,m2,m3) * ( m1sq*m3sq*pow(lm1,2) + 3*m1sq*m2sq*pow(
         m3sq,2)*pow(lm1,3) - 3*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) + 3*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) );

      hh0p +=  - 1./2.*m3sq*lm1 - 1./2.*m3sq*ln3*lm1 - 1./2.*m3sq*ln2*
         lm1 + pow(m3sq,3)*ln3*pow(lm1,2) + 1./2.*m2sq*lm1 + m2sq*ln2*
         lm1 - 2*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) - m2sq*pow(m3sq,2)*ln2
         *pow(lm1,2) + pow(m2sq,2)*m3sq*ln3*pow(lm1,2) + 2*pow(m2sq,2)*
         m3sq*ln2*pow(lm1,2) - pow(m2sq,3)*ln2*pow(lm1,2) - 1./2.*m1sq*
         lm1 - 1./2.*m1sq*ln2*lm1 - 1./2.*m1sq*ln1*lm1 + m1sq*pow(
         m3sq,2)*ln3*pow(lm1,2) - 2*m1sq*pow(m3sq,2)*ln1*pow(lm1,2) + 
         m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 4*m1sq*m2sq*m3sq*ln2*pow(
         lm1,2) + m1sq*m2sq*m3sq*ln1*pow(lm1,2) + 2*m1sq*pow(m2sq,2)*
         ln2*pow(lm1,2) + m1sq*pow(m2sq,2)*ln1*pow(lm1,2) - 2*pow(
         m1sq,2)*m3sq*ln3*pow(lm1,2) + pow(m1sq,2)*m3sq*ln1*pow(lm1,2)
          - pow(m1sq,2)*m2sq*ln2*pow(lm1,2) - 2*pow(m1sq,2)*m2sq*ln1*
         pow(lm1,2) + pow(m1sq,3)*ln1*pow(lm1,2);
      break;
  case 4:
   hh0p =
       + psiq(m1,m2,m3) * ( m1sq*m2sq*pow(lm1,2) - 3*m1sq*m2sq*pow(
         m3sq,2)*pow(lm1,3) + 3*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) + 3*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) );

      hh0p +=  + 1./2.*m3sq*lm1 + m3sq*ln3*lm1 - pow(m3sq,3)*ln3*pow(
         lm1,2) - 1./2.*m2sq*lm1 - 1./2.*m2sq*ln3*lm1 - 1./2.*m2sq*ln2*
         lm1 + 2*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) + m2sq*pow(m3sq,2)*ln2
         *pow(lm1,2) - pow(m2sq,2)*m3sq*ln3*pow(lm1,2) - 2*pow(m2sq,2)*
         m3sq*ln2*pow(lm1,2) + pow(m2sq,3)*ln2*pow(lm1,2) - 1./2.*m1sq*
         lm1 - 1./2.*m1sq*ln3*lm1 - 1./2.*m1sq*ln1*lm1 + 2*m1sq*pow(
         m3sq,2)*ln3*pow(lm1,2) + m1sq*pow(m3sq,2)*ln1*pow(lm1,2) - 4*
         m1sq*m2sq*m3sq*ln3*pow(lm1,2) + m1sq*m2sq*m3sq*ln2*pow(lm1,2)
          + m1sq*m2sq*m3sq*ln1*pow(lm1,2) + m1sq*pow(m2sq,2)*ln2*pow(
         lm1,2) - 2*m1sq*pow(m2sq,2)*ln1*pow(lm1,2) - pow(m1sq,2)*m3sq*
         ln3*pow(lm1,2) - 2*pow(m1sq,2)*m3sq*ln1*pow(lm1,2) - 2*pow(
         m1sq,2)*m2sq*ln2*pow(lm1,2) + pow(m1sq,2)*m2sq*ln1*pow(lm1,2)
          + pow(m1sq,3)*ln1*pow(lm1,2);
      break;
  case 5:
   hh0p =
       + psiq(m1,m2,m3) * ( m3sq*pow(lm1,2) + 3*m2sq*pow(m3sq,2)*pow(
         lm1,3) - 3*pow(m2sq,2)*m3sq*pow(lm1,3) + 3*m1sq*pow(m3sq,2)*
         pow(lm1,3) + 9*m1sq*m2sq*m3sq*pow(lm1,3) + 15*m1sq*m2sq*pow(
         m3sq,3)*pow(lm1,4) - 15*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) - 3*
         pow(m1sq,2)*m3sq*pow(lm1,3) + 30*pow(m1sq,2)*pow(m2sq,2)*m3sq*
         pow(lm1,4) - 15*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) );

      hh0p +=  - lm1 - 1./2.*ln2*lm1 - 1./2.*ln1*lm1 - 3*pow(m3sq,2)*
         pow(lm1,2) + pow(m3sq,2)*ln3*pow(lm1,2) - 2*pow(m3sq,2)*ln2*
         pow(lm1,2) - 2*pow(m3sq,2)*ln1*pow(lm1,2) + 4*pow(m3sq,4)*ln3*
         pow(lm1,3) + m2sq*m3sq*pow(lm1,2) - 2*m2sq*m3sq*ln3*pow(lm1,2)
          - m2sq*m3sq*ln2*pow(lm1,2) + m2sq*m3sq*ln1*pow(lm1,2) - 8*
         m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 2*pow(m2sq,2)*pow(lm1,2) + 4
         *pow(m2sq,2)*ln2*pow(lm1,2) + pow(m2sq,2)*ln1*pow(lm1,2) - 12*
         pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 12*pow(m2sq,2)*pow(
         m3sq,2)*ln2*pow(lm1,3) + 8*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) - 4
         *pow(m2sq,4)*ln2*pow(lm1,3) + m1sq*m3sq*pow(lm1,2) - m1sq*m3sq
         *ln3*pow(lm1,2) + m1sq*m3sq*ln2*pow(lm1,2) - 2*m1sq*m3sq*ln1*
         pow(lm1,2) - m1sq*pow(m3sq,3)*ln3*pow(lm1,3) - 7*m1sq*pow(
         m3sq,3)*ln1*pow(lm1,3) - 4*m1sq*m2sq*pow(lm1,2) - 5*m1sq*m2sq*
         ln2*pow(lm1,2) - 5*m1sq*m2sq*ln1*pow(lm1,2) + 22*m1sq*m2sq*
         pow(m3sq,2)*ln3*pow(lm1,3) - 10*m1sq*m2sq*pow(m3sq,2)*ln2*pow(
         lm1,3);
      hh0p +=  - 12*m1sq*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) - 9*m1sq*pow(
         m2sq,2)*m3sq*ln3*pow(lm1,3) - 14*m1sq*pow(m2sq,2)*m3sq*ln2*
         pow(lm1,3) + 15*m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + 12*m1sq
         *pow(m2sq,3)*ln2*pow(lm1,3) + 4*m1sq*pow(m2sq,3)*ln1*pow(
         lm1,3) + 2*pow(m1sq,2)*pow(lm1,2) + pow(m1sq,2)*ln2*pow(lm1,2)
          + 4*pow(m1sq,2)*ln1*pow(lm1,2) - 10*pow(m1sq,2)*pow(m3sq,2)*
         ln3*pow(lm1,3) + 10*pow(m1sq,2)*pow(m3sq,2)*ln1*pow(lm1,3) - 6
         *pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) + 14*pow(m1sq,2)*m2sq*
         m3sq*ln2*pow(lm1,3) - 16*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3)
          - 12*pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) - 12*pow(m1sq,2)*
         pow(m2sq,2)*ln1*pow(lm1,3) + 7*pow(m1sq,3)*m3sq*ln3*pow(lm1,3)
          + pow(m1sq,3)*m3sq*ln1*pow(lm1,3) + 4*pow(m1sq,3)*m2sq*ln2*
         pow(lm1,3) + 12*pow(m1sq,3)*m2sq*ln1*pow(lm1,3) - 4*pow(
         m1sq,4)*ln1*pow(lm1,3);
      break;
  case 6:
   hh0p =
       + psiq(m1,m2,m3) * ( m2sq*pow(lm1,2) - 3*m2sq*pow(m3sq,2)*pow(
         lm1,3) + 3*pow(m2sq,2)*m3sq*pow(lm1,3) + 9*m1sq*m2sq*m3sq*pow(
         lm1,3) - 15*m1sq*m2sq*pow(m3sq,3)*pow(lm1,4) + 3*m1sq*pow(
         m2sq,2)*pow(lm1,3) + 15*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) - 3*
         pow(m1sq,2)*m2sq*pow(lm1,3) + 30*pow(m1sq,2)*m2sq*pow(m3sq,2)*
         pow(lm1,4) - 15*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) );

      hh0p +=  - lm1 - 1./2.*ln3*lm1 - 1./2.*ln1*lm1 + 2*pow(m3sq,2)*
         pow(lm1,2) + 4*pow(m3sq,2)*ln3*pow(lm1,2) + pow(m3sq,2)*ln1*
         pow(lm1,2) - 4*pow(m3sq,4)*ln3*pow(lm1,3) + m2sq*m3sq*pow(
         lm1,2) - m2sq*m3sq*ln3*pow(lm1,2) - 2*m2sq*m3sq*ln2*pow(lm1,2)
          + m2sq*m3sq*ln1*pow(lm1,2) + 8*m2sq*pow(m3sq,3)*ln2*pow(
         lm1,3) - 3*pow(m2sq,2)*pow(lm1,2) - 2*pow(m2sq,2)*ln3*pow(
         lm1,2) + pow(m2sq,2)*ln2*pow(lm1,2) - 2*pow(m2sq,2)*ln1*pow(
         lm1,2) + 12*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 12*pow(
         m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - 8*pow(m2sq,3)*m3sq*ln3*
         pow(lm1,3) + 4*pow(m2sq,4)*ln2*pow(lm1,3) - 4*m1sq*m3sq*pow(
         lm1,2) - 5*m1sq*m3sq*ln3*pow(lm1,2) - 5*m1sq*m3sq*ln1*pow(
         lm1,2) + 12*m1sq*pow(m3sq,3)*ln3*pow(lm1,3) + 4*m1sq*pow(
         m3sq,3)*ln1*pow(lm1,3) + m1sq*m2sq*pow(lm1,2) + m1sq*m2sq*ln3*
         pow(lm1,2) - m1sq*m2sq*ln2*pow(lm1,2) - 2*m1sq*m2sq*ln1*pow(
         lm1,2) - 14*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) - 9*m1sq*m2sq
         *pow(m3sq,2)*ln2*pow(lm1,3);
      hh0p +=  + 15*m1sq*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) - 10*m1sq*
         pow(m2sq,2)*m3sq*ln3*pow(lm1,3) + 22*m1sq*pow(m2sq,2)*m3sq*ln2
         *pow(lm1,3) - 12*m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) - m1sq*
         pow(m2sq,3)*ln2*pow(lm1,3) - 7*m1sq*pow(m2sq,3)*ln1*pow(lm1,3)
          + 2*pow(m1sq,2)*pow(lm1,2) + pow(m1sq,2)*ln3*pow(lm1,2) + 4*
         pow(m1sq,2)*ln1*pow(lm1,2) - 12*pow(m1sq,2)*pow(m3sq,2)*ln3*
         pow(lm1,3) - 12*pow(m1sq,2)*pow(m3sq,2)*ln1*pow(lm1,3) + 14*
         pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) - 6*pow(m1sq,2)*m2sq*m3sq
         *ln2*pow(lm1,3) - 16*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) - 10
         *pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) + 10*pow(m1sq,2)*pow(
         m2sq,2)*ln1*pow(lm1,3) + 4*pow(m1sq,3)*m3sq*ln3*pow(lm1,3) + 
         12*pow(m1sq,3)*m3sq*ln1*pow(lm1,3) + 7*pow(m1sq,3)*m2sq*ln2*
         pow(lm1,3) + pow(m1sq,3)*m2sq*ln1*pow(lm1,3) - 4*pow(m1sq,4)*
         ln1*pow(lm1,3);
      break;
  case 7:
   hh0p =
       + psiq(m1,m2,m3) * ( m1sq*pow(lm1,2) - 3*m1sq*pow(m3sq,2)*pow(
         lm1,3) + 9*m1sq*m2sq*m3sq*pow(lm1,3) - 15*m1sq*m2sq*pow(
         m3sq,3)*pow(lm1,4) - 3*m1sq*pow(m2sq,2)*pow(lm1,3) + 30*m1sq*
         pow(m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 15*m1sq*pow(m2sq,3)*m3sq*
         pow(lm1,4) + 3*pow(m1sq,2)*m3sq*pow(lm1,3) + 3*pow(m1sq,2)*
         m2sq*pow(lm1,3) + 15*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) );

      hh0p +=  - lm1 - 1./2.*ln3*lm1 - 1./2.*ln2*lm1 + 2*pow(m3sq,2)*
         pow(lm1,2) + 4*pow(m3sq,2)*ln3*pow(lm1,2) + pow(m3sq,2)*ln2*
         pow(lm1,2) - 4*pow(m3sq,4)*ln3*pow(lm1,3) - 4*m2sq*m3sq*pow(
         lm1,2) - 5*m2sq*m3sq*ln3*pow(lm1,2) - 5*m2sq*m3sq*ln2*pow(
         lm1,2) + 12*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) + 4*m2sq*pow(
         m3sq,3)*ln2*pow(lm1,3) + 2*pow(m2sq,2)*pow(lm1,2) + pow(
         m2sq,2)*ln3*pow(lm1,2) + 4*pow(m2sq,2)*ln2*pow(lm1,2) - 12*
         pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 12*pow(m2sq,2)*pow(
         m3sq,2)*ln2*pow(lm1,3) + 4*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) + 
         12*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 4*pow(m2sq,4)*ln2*pow(
         lm1,3) + m1sq*m3sq*pow(lm1,2) - m1sq*m3sq*ln3*pow(lm1,2) + 
         m1sq*m3sq*ln2*pow(lm1,2) - 2*m1sq*m3sq*ln1*pow(lm1,2) + 8*m1sq
         *pow(m3sq,3)*ln1*pow(lm1,3) + m1sq*m2sq*pow(lm1,2) + m1sq*m2sq
         *ln3*pow(lm1,2) - 2*m1sq*m2sq*ln2*pow(lm1,2) - m1sq*m2sq*ln1*
         pow(lm1,2) - 14*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 15*m1sq
         *m2sq*pow(m3sq,2)*ln2*pow(lm1,3);
      hh0p +=  - 9*m1sq*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 14*m1sq*pow(
         m2sq,2)*m3sq*ln3*pow(lm1,3) - 16*m1sq*pow(m2sq,2)*m3sq*ln2*
         pow(lm1,3) - 6*m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + m1sq*
         pow(m2sq,3)*ln2*pow(lm1,3) + 7*m1sq*pow(m2sq,3)*ln1*pow(lm1,3)
          - 3*pow(m1sq,2)*pow(lm1,2) - 2*pow(m1sq,2)*ln3*pow(lm1,2) - 2
         *pow(m1sq,2)*ln2*pow(lm1,2) + pow(m1sq,2)*ln1*pow(lm1,2) + 12*
         pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 12*pow(m1sq,2)*pow(
         m3sq,2)*ln1*pow(lm1,3) - 10*pow(m1sq,2)*m2sq*m3sq*ln3*pow(
         lm1,3) - 12*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3) + 22*pow(
         m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) + 10*pow(m1sq,2)*pow(m2sq,2)*
         ln2*pow(lm1,3) - 10*pow(m1sq,2)*pow(m2sq,2)*ln1*pow(lm1,3) - 8
         *pow(m1sq,3)*m3sq*ln3*pow(lm1,3) - 7*pow(m1sq,3)*m2sq*ln2*pow(
         lm1,3) - pow(m1sq,3)*m2sq*ln1*pow(lm1,3) + 4*pow(m1sq,4)*ln1*
         pow(lm1,3);
      break;
  case 8:
   hh0p =
       + psiq(m1,m2,m3) * ( pow(lm1,2) - 3*pow(m3sq,2)*pow(lm1,3) + 9*
         m2sq*m3sq*pow(lm1,3) - 15*m2sq*pow(m3sq,3)*pow(lm1,4) - 3*pow(
         m2sq,2)*pow(lm1,3) + 30*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 
         15*pow(m2sq,3)*m3sq*pow(lm1,4) + 9*m1sq*m3sq*pow(lm1,3) - 15*
         m1sq*pow(m3sq,3)*pow(lm1,4) + 9*m1sq*m2sq*pow(lm1,3) + 30*m1sq
         *m2sq*pow(m3sq,2)*pow(lm1,4) - 105*m1sq*m2sq*pow(m3sq,4)*pow(
         lm1,5) + 30*m1sq*pow(m2sq,2)*m3sq*pow(lm1,4) + 105*m1sq*pow(
         m2sq,2)*pow(m3sq,3)*pow(lm1,5) - 15*m1sq*pow(m2sq,3)*pow(
         lm1,4) + 105*m1sq*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) - 105*
         m1sq*pow(m2sq,4)*m3sq*pow(lm1,5) - 3*pow(m1sq,2)*pow(lm1,3) + 
         30*pow(m1sq,2)*pow(m3sq,2)*pow(lm1,4) + 30*pow(m1sq,2)*m2sq*
         m3sq*pow(lm1,4) + 105*pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(lm1,5)
          + 30*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,4) - 210*pow(m1sq,2)*
         pow(m2sq,2)*pow(m3sq,2)*pow(lm1,5) + 105*pow(m1sq,2)*pow(
         m2sq,3)*m3sq*pow(lm1,5) - 15*pow(m1sq,3)*m3sq*pow(lm1,4) - 15*
         pow(m1sq,3)*m2sq*pow(lm1,4) );

      hh0p +=  + psiq(m1,m2,m3) * ( 105*pow(m1sq,3)*m2sq*pow(m3sq,2)*
         pow(lm1,5) + 105*pow(m1sq,3)*pow(m2sq,2)*m3sq*pow(lm1,5) - 105
         *pow(m1sq,4)*m2sq*m3sq*pow(lm1,5) );

      hh0p +=  - 3*m3sq*pow(lm1,2) - 2*m3sq*ln2*pow(lm1,2) - 2*m3sq*
         ln1*pow(lm1,2) + 16*pow(m3sq,3)*pow(lm1,3) + 12*pow(m3sq,3)*
         ln3*pow(lm1,3) + 8*pow(m3sq,3)*ln2*pow(lm1,3) + 8*pow(m3sq,3)*
         ln1*pow(lm1,3) - 24*pow(m3sq,5)*ln3*pow(lm1,4) - 3*m2sq*pow(
         lm1,2) - 2*m2sq*ln3*pow(lm1,2) - m2sq*ln2*pow(lm1,2) - m2sq*
         ln1*pow(lm1,2) - 16*m2sq*pow(m3sq,2)*pow(lm1,3) + 6*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,3) - 25*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) - 
         9*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 24*m2sq*pow(m3sq,4)*ln3*
         pow(lm1,4) + 48*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) - 16*pow(
         m2sq,2)*m3sq*pow(lm1,3) - 26*pow(m2sq,2)*m3sq*ln3*pow(lm1,3)
          + 4*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 6*pow(m2sq,2)*m3sq*ln1*
         pow(lm1,3) + 72*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) - 120*
         pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) + 16*pow(m2sq,3)*pow(
         lm1,3) + 8*pow(m2sq,3)*ln3*pow(lm1,3) + 13*pow(m2sq,3)*ln2*
         pow(lm1,3) + 7*pow(m2sq,3)*ln1*pow(lm1,3) - 120*pow(m2sq,3)*
         pow(m3sq,2)*ln3*pow(lm1,4);
      hh0p +=  + 72*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) + 48*pow(
         m2sq,4)*m3sq*ln3*pow(lm1,4) + 24*pow(m2sq,4)*m3sq*ln2*pow(
         lm1,4) - 24*pow(m2sq,5)*ln2*pow(lm1,4) - 2*m1sq*pow(lm1,2) - 
         m1sq*ln3*pow(lm1,2) - m1sq*ln2*pow(lm1,2) - 2*m1sq*ln1*pow(
         lm1,2) - 17*m1sq*pow(m3sq,2)*pow(lm1,3) - m1sq*pow(m3sq,2)*ln3
         *pow(lm1,3) - 9*m1sq*pow(m3sq,2)*ln2*pow(lm1,3) - 18*m1sq*pow(
         m3sq,2)*ln1*pow(lm1,3) + 30*m1sq*pow(m3sq,4)*ln3*pow(lm1,4) + 
         42*m1sq*pow(m3sq,4)*ln1*pow(lm1,4) + 46*m1sq*m2sq*m3sq*pow(
         lm1,3) + 14*m1sq*m2sq*m3sq*ln3*pow(lm1,3) + 9*m1sq*m2sq*m3sq*
         ln2*pow(lm1,3) + m1sq*m2sq*m3sq*ln1*pow(lm1,3) - 168*m1sq*m2sq
         *pow(m3sq,3)*ln3*pow(lm1,4) + 27*m1sq*m2sq*pow(m3sq,3)*ln2*
         pow(lm1,4) + 45*m1sq*m2sq*pow(m3sq,3)*ln1*pow(lm1,4) - 17*m1sq
         *pow(m2sq,2)*pow(lm1,3) - 9*m1sq*pow(m2sq,2)*ln3*pow(lm1,3) - 
         6*m1sq*pow(m2sq,2)*ln2*pow(lm1,3) - 13*m1sq*pow(m2sq,2)*ln1*
         pow(lm1,3) + 114*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,4)
          + 111*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4);
      hh0p +=  - 177*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) + 24*
         m1sq*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) - 171*m1sq*pow(m2sq,3)*
         m3sq*ln2*pow(lm1,4) + 51*m1sq*pow(m2sq,3)*m3sq*ln1*pow(lm1,4)
          + 33*m1sq*pow(m2sq,4)*ln2*pow(lm1,4) + 39*m1sq*pow(m2sq,4)*
         ln1*pow(lm1,4) - 14*pow(m1sq,2)*m3sq*pow(lm1,3) - 18*pow(
         m1sq,2)*m3sq*ln3*pow(lm1,3) - 6*pow(m1sq,2)*m3sq*ln2*pow(
         lm1,3) - 4*pow(m1sq,2)*m3sq*ln1*pow(lm1,3) + 54*pow(m1sq,2)*
         pow(m3sq,3)*ln3*pow(lm1,4) - 102*pow(m1sq,2)*pow(m3sq,3)*ln1*
         pow(lm1,4) - 14*pow(m1sq,2)*m2sq*pow(lm1,3) - 6*pow(m1sq,2)*
         m2sq*ln3*pow(lm1,3) - 14*pow(m1sq,2)*m2sq*ln2*pow(lm1,3) - 8*
         pow(m1sq,2)*m2sq*ln1*pow(lm1,3) + 108*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,4) - 159*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,4) + 99*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) - 
         150*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) + 102*pow(
         m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 96*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*ln1*pow(lm1,4);
      hh0p +=  + 45*pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,4) - 93*pow(
         m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,4) + 15*pow(m1sq,3)*pow(lm1,3)
          + 7*pow(m1sq,3)*ln3*pow(lm1,3) + 7*pow(m1sq,3)*ln2*pow(lm1,3)
          + 14*pow(m1sq,3)*ln1*pow(lm1,3) - 102*pow(m1sq,3)*pow(m3sq,2)
         *ln3*pow(lm1,4) + 54*pow(m1sq,3)*pow(m3sq,2)*ln1*pow(lm1,4) + 
         36*pow(m1sq,3)*m2sq*m3sq*ln3*pow(lm1,4) + 45*pow(m1sq,3)*m2sq*
         m3sq*ln2*pow(lm1,4) - 177*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,4)
          - 93*pow(m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,4) + 45*pow(m1sq,3)*
         pow(m2sq,2)*ln1*pow(lm1,4) + 42*pow(m1sq,4)*m3sq*ln3*pow(
         lm1,4) + 30*pow(m1sq,4)*m3sq*ln1*pow(lm1,4) + 39*pow(m1sq,4)*
         m2sq*ln2*pow(lm1,4) + 33*pow(m1sq,4)*m2sq*ln1*pow(lm1,4) - 24*
         pow(m1sq,5)*ln1*pow(lm1,4);
      break;
     default:
    std::cout <<"wrong iprop in h0 (quenched case)\n";
    hh0p=0;
    break;
  }

  return pi162*hh0p;
}

double h10(const int iprop, const double m1sq, const double m2sq,
	   const double m3sq, const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h10sing(iprop,m1sq,m2sq,m3sq,xmu2);
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double h10;
  switch(iprop){
  case 1:
    h10 = 0.25*(-1.+m1sq*dm*lm1)*psiq(m1,m2,m3)
	   +m1sq*(3./8.-0.5*ln1)
	   +m2sq*(pi2/24.+9./16.-0.25*ln2+0.25*(ln2*ln4-ln1*ln3))
	   +m3sq*(pi2/24.+9./16.-0.25*ln3+0.25*(ln3*ln4-ln1*ln2));
    break;
  case 2:
    h10 =
       + psiq(m1,m2,m3) * ( 1./4.*m1sq*lm1 - 1./4.*m1sq*pow(m3sq,2)*pow(
         lm1,2) - 1./2.*m1sq*m2sq*m3sq*pow(lm1,2) - 1./4.*m1sq*pow(
         m2sq,2)*pow(lm1,2) + 1./2.*pow(m1sq,2)*m3sq*pow(lm1,2) + 1./2.
         *pow(m1sq,2)*m2sq*pow(lm1,2) - 1./4.*pow(m1sq,3)*pow(lm1,2) );

    h10 +=  - 1./8. - 1./4.*ln3 - 1./4.*ln2 - 1./4.*pow(m3sq,2)*ln3
         *lm1 + 1./4.*pow(m3sq,2)*ln2*lm1 + 1./4.*pow(m2sq,2)*ln3*lm1
          - 1./4.*pow(m2sq,2)*ln2*lm1 - 1./2.*m1sq*m3sq*ln2*lm1 + 1./2.
         *m1sq*m3sq*ln1*lm1 - 1./2.*m1sq*m2sq*ln3*lm1 + 1./2.*m1sq*m2sq
         *ln1*lm1 + 1./4.*pow(m1sq,2)*ln3*lm1 + 1./4.*pow(m1sq,2)*ln2*
         lm1 - 1./2.*pow(m1sq,2)*ln1*lm1;
    break;
  case 3:
    h10 =
       + psiq(m1,m2,m3) * ( 1./4.*m3sq*lm1 - 1./4.*m2sq*lm1 - 1./4.*m1sq
         *pow(m3sq,2)*pow(lm1,2) + 1./4.*m1sq*pow(m2sq,2)*pow(lm1,2) - 
         1./2.*pow(m1sq,2)*m2sq*pow(lm1,2) + 1./4.*pow(m1sq,3)*pow(
         lm1,2) );

    h10 +=  + 5./16. + 3./4.*ln2 + 1./4.*ln2*ln3 + 1./4.*pow(ln2,2)
          - 1./4.*ln1*ln3 + 1./4.*ln1*ln2 + 1./4.*m1sq*pow(m2sq,-1)*ln3
          - 1./4.*m1sq*pow(m2sq,-1)*ln1 - 1./4.*m1sq*pow(m2sq,-1)*pow(
         m3sq,2)*ln3*lm1 + 1./4.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm1
          - 1./2.*m1sq*m3sq*ln3*lm1 + 1./2.*m1sq*m3sq*ln2*lm1 - 1./4.*
         m1sq*m2sq*ln3*lm1 + 1./2.*m1sq*m2sq*ln2*lm1 - 1./4.*m1sq*m2sq*
         ln1*lm1 + 1./2.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln3*lm1 - 1./2.*
         pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*lm1 + 1./2.*pow(m1sq,2)*ln3*
         lm1 - 1./2.*pow(m1sq,2)*ln2*lm1 - 1./4.*pow(m1sq,3)*pow(
         m2sq,-1)*ln3*lm1 + 1./4.*pow(m1sq,3)*pow(m2sq,-1)*ln1*lm1 + 1./
         24.*pi2;
    break;
  case 4:
    h10 =
       + psiq(m1,m2,m3) * (  - 1./4.*m3sq*lm1 + 1./4.*m2sq*lm1 + 1./4.*
         m1sq*pow(m3sq,2)*pow(lm1,2) - 1./4.*m1sq*pow(m2sq,2)*pow(
         lm1,2) - 1./2.*pow(m1sq,2)*m3sq*pow(lm1,2) + 1./4.*pow(m1sq,3)
         *pow(lm1,2) );

    h10 +=  + 5./16. + 3./4.*ln3 + 1./4.*pow(ln3,2) + 1./4.*ln2*ln3
          + 1./4.*ln1*ln3 - 1./4.*ln1*ln2 + 1./4.*m1sq*pow(m3sq,-1)*ln2
          - 1./4.*m1sq*pow(m3sq,-1)*ln1 + 1./2.*m1sq*m3sq*ln3*lm1 - 1./
         4.*m1sq*m3sq*ln2*lm1 - 1./4.*m1sq*m3sq*ln1*lm1 + 1./2.*m1sq*
         m2sq*ln3*lm1 - 1./2.*m1sq*m2sq*ln2*lm1 - 1./4.*m1sq*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*lm1 + 1./4.*m1sq*pow(m2sq,2)*pow(
         m3sq,-1)*ln1*lm1 - 1./2.*pow(m1sq,2)*ln3*lm1 + 1./2.*pow(
         m1sq,2)*ln2*lm1 + 1./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*lm1
          - 1./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln1*lm1 - 1./4.*pow(
         m1sq,3)*pow(m3sq,-1)*ln2*lm1 + 1./4.*pow(m1sq,3)*pow(m3sq,-1)*
         ln1*lm1 + 1./24.*pi2;
    break;
  case 5:
    h10 =
       + psiq(m1,m2,m3) * (  - 1./4.*m1sq*m3sq*pow(lm1,2) - 3./4.*m1sq*
         pow(m3sq,3)*pow(lm1,3) - 3./4.*m1sq*m2sq*pow(lm1,2) - 3./4.*
         m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) + 3./4.*m1sq*pow(m2sq,2)*m3sq
         *pow(lm1,3) + 3./4.*m1sq*pow(m2sq,3)*pow(lm1,3) + 3./4.*pow(
         m1sq,2)*pow(lm1,2) + 3./4.*pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3)
          - 3./2.*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 9./4.*pow(m1sq,2)*
         pow(m2sq,2)*pow(lm1,3) + 3./4.*pow(m1sq,3)*m3sq*pow(lm1,3) + 9.
         /4.*pow(m1sq,3)*m2sq*pow(lm1,3) - 3./4.*pow(m1sq,4)*pow(lm1,3)
          );

    h10 +=  - 1./4.*pow(m2sq,-1) + 1./4.*pow(m2sq,-1)*pow(m3sq,2)*
         lm1 - 1./2.*pow(m3sq,3)*ln3*pow(lm1,2) + 1./2.*pow(m3sq,3)*ln2
         *pow(lm1,2) - 1./4.*m2sq*lm1 + 1./2.*m2sq*ln3*lm1 - 1./2.*m2sq
         *ln2*lm1 + 1./2.*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) - 1./2.*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,2) + 1./2.*pow(m2sq,2)*m3sq*ln3*pow(
         lm1,2) - 1./2.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) - 1./2.*pow(
         m2sq,3)*ln3*pow(lm1,2) + 1./2.*pow(m2sq,3)*ln2*pow(lm1,2) - 1./
         2.*m1sq*pow(m2sq,-1)*m3sq*lm1 + 1./4.*m1sq*pow(m2sq,-1)*m3sq*
         ln3*lm1 - 1./4.*m1sq*pow(m2sq,-1)*m3sq*ln1*lm1 - 1./4.*m1sq*
         pow(m2sq,-1)*pow(m3sq,3)*ln3*pow(lm1,2) + 1./4.*m1sq*pow(
         m2sq,-1)*pow(m3sq,3)*ln1*pow(lm1,2) - 1./4.*m1sq*ln3*lm1 - 1./
         2.*m1sq*ln2*lm1 + 3./4.*m1sq*ln1*lm1 - 5./4.*m1sq*pow(m3sq,2)*
         ln3*pow(lm1,2) + 5./4.*m1sq*pow(m3sq,2)*ln1*pow(lm1,2) - 7./4.
         *m1sq*m2sq*m3sq*ln3*pow(lm1,2) + 2*m1sq*m2sq*m3sq*ln2*pow(
         lm1,2) - 1./4.*m1sq*m2sq*m3sq*ln1*pow(lm1,2) + 5./4.*m1sq*pow(
         m2sq,2)*ln3*pow(lm1,2);
    h10 +=  - 5./4.*m1sq*pow(m2sq,2)*ln1*pow(lm1,2) + 1./4.*pow(
         m1sq,2)*pow(m2sq,-1)*lm1 - 1./4.*pow(m1sq,2)*pow(m2sq,-1)*ln3*
         lm1 + 1./4.*pow(m1sq,2)*pow(m2sq,-1)*ln1*lm1 + 3./4.*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(lm1,2) - 3./4.*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm1,2) + 2*pow(
         m1sq,2)*m3sq*ln3*pow(lm1,2) - 3./2.*pow(m1sq,2)*m3sq*ln2*pow(
         lm1,2) - 1./2.*pow(m1sq,2)*m3sq*ln1*pow(lm1,2) - 3./4.*pow(
         m1sq,2)*m2sq*ln3*pow(lm1,2) - 3./2.*pow(m1sq,2)*m2sq*ln2*pow(
         lm1,2) + 9./4.*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) - 3./4.*pow(
         m1sq,3)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,2) + 3./4.*pow(m1sq,3)*
         pow(m2sq,-1)*m3sq*ln1*pow(lm1,2) - 1./4.*pow(m1sq,3)*ln3*pow(
         lm1,2) + pow(m1sq,3)*ln2*pow(lm1,2) - 3./4.*pow(m1sq,3)*ln1*
         pow(lm1,2) + 1./4.*pow(m1sq,4)*pow(m2sq,-1)*ln3*pow(lm1,2) - 1.
         /4.*pow(m1sq,4)*pow(m2sq,-1)*ln1*pow(lm1,2);
    break;
  case 6:
    h10 =
       + psiq(m1,m2,m3) * (  - 3./4.*m1sq*m3sq*pow(lm1,2) + 3./4.*m1sq*
         pow(m3sq,3)*pow(lm1,3) - 1./4.*m1sq*m2sq*pow(lm1,2) + 3./4.*
         m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 3./4.*m1sq*pow(m2sq,2)*m3sq
         *pow(lm1,3) - 3./4.*m1sq*pow(m2sq,3)*pow(lm1,3) + 3./4.*pow(
         m1sq,2)*pow(lm1,2) - 9./4.*pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3)
          - 3./2.*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) + 3./4.*pow(m1sq,2)*
         pow(m2sq,2)*pow(lm1,3) + 9./4.*pow(m1sq,3)*m3sq*pow(lm1,3) + 3.
         /4.*pow(m1sq,3)*m2sq*pow(lm1,3) - 3./4.*pow(m1sq,4)*pow(lm1,3)
          );

    h10 +=  - 1./4.*pow(m3sq,-1) - 1./4.*m3sq*lm1 - 1./2.*m3sq*ln3*
         lm1 + 1./2.*m3sq*ln2*lm1 + 1./2.*pow(m3sq,3)*ln3*pow(lm1,2) - 
         1./2.*pow(m3sq,3)*ln2*pow(lm1,2) - 1./2.*m2sq*pow(m3sq,2)*ln3*
         pow(lm1,2) + 1./2.*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./4.*
         pow(m2sq,2)*pow(m3sq,-1)*lm1 - 1./2.*pow(m2sq,2)*m3sq*ln3*pow(
         lm1,2) + 1./2.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + 1./2.*pow(
         m2sq,3)*ln3*pow(lm1,2) - 1./2.*pow(m2sq,3)*ln2*pow(lm1,2) - 1./
         2.*m1sq*ln3*lm1 - 1./4.*m1sq*ln2*lm1 + 3./4.*m1sq*ln1*lm1 + 5./
         4.*m1sq*pow(m3sq,2)*ln2*pow(lm1,2) - 5./4.*m1sq*pow(m3sq,2)*
         ln1*pow(lm1,2) - 1./2.*m1sq*m2sq*pow(m3sq,-1)*lm1 + 1./4.*m1sq
         *m2sq*pow(m3sq,-1)*ln2*lm1 - 1./4.*m1sq*m2sq*pow(m3sq,-1)*ln1*
         lm1 + 2*m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 7./4.*m1sq*m2sq*m3sq*
         ln2*pow(lm1,2) - 1./4.*m1sq*m2sq*m3sq*ln1*pow(lm1,2) - 5./4.*
         m1sq*pow(m2sq,2)*ln2*pow(lm1,2) + 5./4.*m1sq*pow(m2sq,2)*ln1*
         pow(lm1,2) - 1./4.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln2*pow(
         lm1,2);
    h10 +=  + 1./4.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln1*pow(lm1,2) + 
         1./4.*pow(m1sq,2)*pow(m3sq,-1)*lm1 - 1./4.*pow(m1sq,2)*pow(
         m3sq,-1)*ln2*lm1 + 1./4.*pow(m1sq,2)*pow(m3sq,-1)*ln1*lm1 - 3./
         2.*pow(m1sq,2)*m3sq*ln3*pow(lm1,2) - 3./4.*pow(m1sq,2)*m3sq*
         ln2*pow(lm1,2) + 9./4.*pow(m1sq,2)*m3sq*ln1*pow(lm1,2) - 3./2.
         *pow(m1sq,2)*m2sq*ln3*pow(lm1,2) + 2*pow(m1sq,2)*m2sq*ln2*pow(
         lm1,2) - 1./2.*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) + 3./4.*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(lm1,2) - 3./4.*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(lm1,2) + pow(m1sq,3)*
         ln3*pow(lm1,2) - 1./4.*pow(m1sq,3)*ln2*pow(lm1,2) - 3./4.*pow(
         m1sq,3)*ln1*pow(lm1,2) - 3./4.*pow(m1sq,3)*m2sq*pow(m3sq,-1)*
         ln2*pow(lm1,2) + 3./4.*pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln1*pow(
         lm1,2) + 1./4.*pow(m1sq,4)*pow(m3sq,-1)*ln2*pow(lm1,2) - 1./4.
         *pow(m1sq,4)*pow(m3sq,-1)*ln1*pow(lm1,2);
    break;
  case 7:
    h10 =
       + psiq(m1,m2,m3) * ( 1./4.*lm1 - 1./4.*pow(m3sq,2)*pow(lm1,2) + 1.
         /2.*m2sq*m3sq*pow(lm1,2) - 1./4.*pow(m2sq,2)*pow(lm1,2) - 1./4.
         *m1sq*m3sq*pow(lm1,2) + 3./4.*m1sq*pow(m3sq,3)*pow(lm1,3) - 1./
         4.*m1sq*m2sq*pow(lm1,2) - 3./4.*m1sq*m2sq*pow(m3sq,2)*pow(
         lm1,3) - 3./4.*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) + 3./4.*m1sq*
         pow(m2sq,3)*pow(lm1,3) - 3./4.*pow(m1sq,2)*pow(m3sq,2)*pow(
         lm1,3) + 3./2.*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 3./4.*pow(
         m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 3./4.*pow(m1sq,3)*m3sq*pow(
         lm1,3) - 3./4.*pow(m1sq,3)*m2sq*pow(lm1,3) + 3./4.*pow(m1sq,4)
         *pow(lm1,3) );

    h10 +=  + 1./4.*pow(m3sq,-1)*ln2 - 1./4.*pow(m3sq,-1)*ln1 - 1./
         2.*m3sq*ln3*lm1 + 1./4.*m3sq*ln2*lm1 + 1./4.*m3sq*ln1*lm1 + 1./
         2.*m2sq*ln3*lm1 - 1./2.*m2sq*ln1*lm1 - 1./4.*pow(m2sq,2)*pow(
         m3sq,-1)*ln2*lm1 + 1./4.*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm1 + 1./
         4.*m1sq*pow(m2sq,-1)*pow(m3sq,-1) - 1./4.*m1sq*pow(m2sq,-1)*
         m3sq*lm1 - 1./2.*m1sq*pow(m2sq,-1)*m3sq*ln3*lm1 + 1./2.*m1sq*
         pow(m2sq,-1)*m3sq*ln1*lm1 + 1./2.*m1sq*pow(m2sq,-1)*pow(
         m3sq,3)*ln3*pow(lm1,2) - 1./2.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*
         ln1*pow(lm1,2) - 1./2.*m1sq*lm1 - 1./2.*m1sq*ln3*lm1 + 1./4.*
         m1sq*ln2*lm1 + 1./4.*m1sq*ln1*lm1 + m1sq*pow(m3sq,2)*ln3*pow(
         lm1,2) - 5./4.*m1sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./4.*m1sq*
         pow(m3sq,2)*ln1*pow(lm1,2) - 1./4.*m1sq*m2sq*pow(m3sq,-1)*lm1
          + 1./4.*m1sq*m2sq*pow(m3sq,-1)*ln2*lm1 - 1./4.*m1sq*m2sq*pow(
         m3sq,-1)*ln1*lm1 - 1./2.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 1./4.
         *m1sq*m2sq*m3sq*ln2*pow(lm1,2) + 3./4.*m1sq*m2sq*m3sq*ln1*pow(
         lm1,2);
    h10 +=  - m1sq*pow(m2sq,2)*ln3*pow(lm1,2) + 5./4.*m1sq*pow(
         m2sq,2)*ln2*pow(lm1,2) - 1./4.*m1sq*pow(m2sq,2)*ln1*pow(lm1,2)
          + 1./4.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln2*pow(lm1,2) - 1./4.*
         m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln1*pow(lm1,2) + 1./2.*pow(
         m1sq,2)*pow(m2sq,-1)*lm1 + 1./2.*pow(m1sq,2)*pow(m2sq,-1)*ln3*
         lm1 - 1./2.*pow(m1sq,2)*pow(m2sq,-1)*ln1*lm1 - 3./2.*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(lm1,2) + 3./2.*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm1,2) + 1./2.*pow(
         m1sq,2)*pow(m3sq,-1)*lm1 - pow(m1sq,2)*m3sq*ln3*pow(lm1,2) + 9.
         /4.*pow(m1sq,2)*m3sq*ln2*pow(lm1,2) - 5./4.*pow(m1sq,2)*m3sq*
         ln1*pow(lm1,2) + 3./2.*pow(m1sq,2)*m2sq*ln3*pow(lm1,2) - 1./2.
         *pow(m1sq,2)*m2sq*ln2*pow(lm1,2) - pow(m1sq,2)*m2sq*ln1*pow(
         lm1,2) - 3./4.*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(
         lm1,2) + 3./4.*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(
         lm1,2) - 1./4.*pow(m1sq,3)*pow(m2sq,-1)*pow(m3sq,-1)*lm1 + 3./
         2.*pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,2);
    h10 +=  - 3./2.*pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln1*pow(lm1,2) - 
         3./4.*pow(m1sq,3)*ln2*pow(lm1,2) + 3./4.*pow(m1sq,3)*ln1*pow(
         lm1,2) + 3./4.*pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln2*pow(lm1,2) - 
         3./4.*pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln1*pow(lm1,2) - 1./2.*
         pow(m1sq,4)*pow(m2sq,-1)*ln3*pow(lm1,2) + 1./2.*pow(m1sq,4)*
         pow(m2sq,-1)*ln1*pow(lm1,2) - 1./4.*pow(m1sq,4)*pow(m3sq,-1)*
         ln2*pow(lm1,2) + 1./4.*pow(m1sq,4)*pow(m3sq,-1)*ln1*pow(lm1,2)
         ;
    break;
  case 8:
    h10 =
       + psiq(m1,m2,m3) * (  - 1./4.*m1sq*pow(lm1,2) - 3./2.*m1sq*pow(
         m3sq,2)*pow(lm1,3) + 15./4.*m1sq*pow(m3sq,4)*pow(lm1,4) - 3./2.
         *m1sq*pow(m2sq,2)*pow(lm1,3) - 15./2.*m1sq*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,4) + 15./4.*m1sq*pow(m2sq,4)*pow(lm1,4) - 3./2.
         *pow(m1sq,2)*m3sq*pow(lm1,3) - 15./2.*pow(m1sq,2)*pow(m3sq,3)*
         pow(lm1,4) - 3./2.*pow(m1sq,2)*m2sq*pow(lm1,3) + 15./2.*pow(
         m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) + 15./2.*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,4) - 15./2.*pow(m1sq,2)*pow(m2sq,3)*pow(
         lm1,4) + 3*pow(m1sq,3)*pow(lm1,3) - 15*pow(m1sq,3)*m2sq*m3sq*
         pow(lm1,4) + 15./2.*pow(m1sq,4)*m3sq*pow(lm1,4) + 15./2.*pow(
         m1sq,4)*m2sq*pow(lm1,4) - 15./4.*pow(m1sq,5)*pow(lm1,4) );

    h10 +=  + 1./2.*pow(m2sq,-1)*m3sq*lm1 - 1./2.*pow(m2sq,-1)*pow(
         m3sq,3)*pow(lm1,2) - 3./2.*pow(m3sq,2)*ln3*pow(lm1,2) + 3./2.*
         pow(m3sq,2)*ln2*pow(lm1,2) + 2*pow(m3sq,4)*ln3*pow(lm1,3) - 2*
         pow(m3sq,4)*ln2*pow(lm1,3) + 1./2.*m2sq*pow(m3sq,-1)*lm1 + 
         m2sq*m3sq*pow(lm1,2) - 4*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) + 4*
         m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 3./2.*pow(m2sq,2)*ln3*pow(
         lm1,2) - 3./2.*pow(m2sq,2)*ln2*pow(lm1,2) - 1./2.*pow(m2sq,3)*
         pow(m3sq,-1)*pow(lm1,2) + 4*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) - 
         4*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 2*pow(m2sq,4)*ln3*pow(
         lm1,3) + 2*pow(m2sq,4)*ln2*pow(lm1,3) - 1./4.*m1sq*pow(
         m2sq,-1)*lm1 + 1./4.*m1sq*pow(m2sq,-1)*ln3*lm1 - 1./4.*m1sq*
         pow(m2sq,-1)*ln1*lm1 + 5./4.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*
         pow(lm1,2) - 5./4.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(
         lm1,2) + 5./4.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm1,2) + 
         m1sq*pow(m2sq,-1)*pow(m3sq,4)*ln3*pow(lm1,3) - m1sq*pow(
         m2sq,-1)*pow(m3sq,4)*ln1*pow(lm1,3);
    h10 +=  - 1./4.*m1sq*pow(m3sq,-1)*lm1 - 9./4.*m1sq*m3sq*pow(
         lm1,2) - m1sq*m3sq*ln3*pow(lm1,2) + 3./4.*m1sq*m3sq*ln2*pow(
         lm1,2) + 1./4.*m1sq*m3sq*ln1*pow(lm1,2) + 7./2.*m1sq*pow(
         m3sq,3)*ln3*pow(lm1,3) + 5./4.*m1sq*pow(m3sq,3)*ln2*pow(lm1,3)
          - 19./4.*m1sq*pow(m3sq,3)*ln1*pow(lm1,3) - 9./4.*m1sq*m2sq*
         pow(lm1,2) + 1./4.*m1sq*m2sq*ln3*pow(lm1,2) - m1sq*m2sq*ln2*
         pow(lm1,2) + 3./4.*m1sq*m2sq*ln1*pow(lm1,2) + 11./2.*m1sq*m2sq
         *pow(m3sq,2)*ln3*pow(lm1,3) - 23./2.*m1sq*m2sq*pow(m3sq,2)*ln2
         *pow(lm1,3) + 6*m1sq*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 5./4.*
         m1sq*pow(m2sq,2)*pow(m3sq,-1)*pow(lm1,2) - 3./4.*m1sq*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*pow(lm1,2) + 3./4.*m1sq*pow(m2sq,2)*
         pow(m3sq,-1)*ln1*pow(lm1,2) - 23./2.*m1sq*pow(m2sq,2)*m3sq*ln3
         *pow(lm1,3) + 6*m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) + 11./2.*
         m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + 3./2.*m1sq*pow(m2sq,3)*
         ln3*pow(lm1,3) + 7./2.*m1sq*pow(m2sq,3)*ln2*pow(lm1,3) - 5*
         m1sq*pow(m2sq,3)*ln1*pow(lm1,3);
    h10 +=  + 3./4.*m1sq*pow(m2sq,4)*pow(m3sq,-1)*ln2*pow(lm1,3) - 
         3./4.*m1sq*pow(m2sq,4)*pow(m3sq,-1)*ln1*pow(lm1,3) - 1./4.*
         pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,-1)*lm1 - 3./4.*pow(m1sq,2)*
         pow(m2sq,-1)*m3sq*pow(lm1,2) + 5./2.*pow(m1sq,2)*pow(m2sq,-1)*
         m3sq*ln3*pow(lm1,2) - 5./2.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*
         pow(lm1,2) - 4*pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln3*pow(
         lm1,3) + 4*pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln1*pow(lm1,3)
          + 5./2.*pow(m1sq,2)*pow(lm1,2) - 1./2.*pow(m1sq,2)*ln3*pow(
         lm1,2) - 3./2.*pow(m1sq,2)*ln2*pow(lm1,2) + 2*pow(m1sq,2)*ln1*
         pow(lm1,2) - 23./2.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 
         15./2.*pow(m1sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 4*pow(m1sq,2)*
         pow(m3sq,2)*ln1*pow(lm1,3) - 3./4.*pow(m1sq,2)*m2sq*pow(
         m3sq,-1)*pow(lm1,2) + 3./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*
         pow(lm1,2) - 3./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln1*pow(
         lm1,2) + 7*pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) + 8*pow(
         m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3);
    h10 +=  - 15*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) + 13./2.*pow(
         m1sq,2)*pow(m2sq,2)*ln3*pow(lm1,3) - 21./2.*pow(m1sq,2)*pow(
         m2sq,2)*ln2*pow(lm1,3) + 4*pow(m1sq,2)*pow(m2sq,2)*ln1*pow(
         lm1,3) - 3*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,-1)*ln2*pow(lm1,3)
          + 3*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,-1)*ln1*pow(lm1,3) - 1./
         4.*pow(m1sq,3)*pow(m2sq,-1)*pow(lm1,2) - 5./4.*pow(m1sq,3)*
         pow(m2sq,-1)*ln3*pow(lm1,2) + 5./4.*pow(m1sq,3)*pow(m2sq,-1)*
         ln1*pow(lm1,2) + 6*pow(m1sq,3)*pow(m2sq,-1)*pow(m3sq,2)*ln3*
         pow(lm1,3) - 6*pow(m1sq,3)*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(
         lm1,3) - 1./4.*pow(m1sq,3)*pow(m3sq,-1)*pow(lm1,2) - 3./4.*
         pow(m1sq,3)*pow(m3sq,-1)*ln2*pow(lm1,2) + 3./4.*pow(m1sq,3)*
         pow(m3sq,-1)*ln1*pow(lm1,2) + 9./2.*pow(m1sq,3)*m3sq*ln3*pow(
         lm1,3) - 10*pow(m1sq,3)*m3sq*ln2*pow(lm1,3) + 11./2.*pow(
         m1sq,3)*m3sq*ln1*pow(lm1,3) - 17./2.*pow(m1sq,3)*m2sq*ln3*pow(
         lm1,3) + 5./2.*pow(m1sq,3)*m2sq*ln2*pow(lm1,3) + 6*pow(m1sq,3)
         *m2sq*ln1*pow(lm1,3);
    h10 +=  + 9./2.*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(
      lm1,3) - 9./2.*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(
         lm1,3) + 1./4.*pow(m1sq,4)*pow(m2sq,-1)*pow(m3sq,-1)*pow(
         lm1,2) - 4*pow(m1sq,4)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,3) + 4*
         pow(m1sq,4)*pow(m2sq,-1)*m3sq*ln1*pow(lm1,3) + 3./2.*pow(
         m1sq,4)*ln3*pow(lm1,3) + 5./2.*pow(m1sq,4)*ln2*pow(lm1,3) - 4*
         pow(m1sq,4)*ln1*pow(lm1,3) - 3*pow(m1sq,4)*m2sq*pow(m3sq,-1)*
         ln2*pow(lm1,3) + 3*pow(m1sq,4)*m2sq*pow(m3sq,-1)*ln1*pow(
         lm1,3) + pow(m1sq,5)*pow(m2sq,-1)*ln3*pow(lm1,3) - pow(m1sq,5)
         *pow(m2sq,-1)*ln1*pow(lm1,3) + 3./4.*pow(m1sq,5)*pow(m3sq,-1)*
         ln2*pow(lm1,3) - 3./4.*pow(m1sq,5)*pow(m3sq,-1)*ln1*pow(lm1,3)
         ;
    break;
  default:
    std::cout <<"h10 (quenched) called with wromg iprop ="<<iprop<<'\n';
    h10 = 0.;
  }
  return h10*pi162;
}

double h10p(const int iprop, const double m1sq, const double m2sq,
	    const double m3sq, const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h10psing(iprop,m1sq,m2sq,m3sq,xmu2);
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double h10p;
  switch(iprop){
  case 1:
    h10p = pow(m1sq,2)*m2sq*m3sq*dm*pow(lm1,3)*psiq(m1,m2,m3)
     -m1sq*dm/(6.*lm)+pow(m1sq,2)*ln1/(6.*lm)*(1.+12.*m2sq*m3sq*lm1)
     +m2sq*ln2/(6.*lm)*
           (m2sq-m3sq-2.*m1sq-6.*m1sq*m3sq*lm1*(dm+2.*m2sq))
      +m3sq*ln3/(6.*lm)*
	  (m3sq-m2sq-2.*m1sq-6.*m1sq*m2sq*lm1*(dm+2.*m3sq))
	  + 7./72.;
    break;
  case 2:
    h10p =
       + psiq(m1,m2,m3) * (  - 2*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 2*
         m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) + 3*pow(m1sq,2)*m2sq*m3sq*
         pow(lm1,3) - 5*pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(lm1,4) - 10*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 5*pow(m1sq,2)
         *pow(m2sq,3)*m3sq*pow(lm1,4) + 10*pow(m1sq,3)*m2sq*pow(m3sq,2)
         *pow(lm1,4) + 10*pow(m1sq,3)*pow(m2sq,2)*m3sq*pow(lm1,4) - 5*
         pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) );

    h10p +=  + 1./6.*m3sq*lm1 - 1./3.*m3sq*ln3*lm1 + 1./3.*pow(
         m3sq,3)*ln3*pow(lm1,2) + 1./6.*m2sq*lm1 - 1./3.*m2sq*ln2*lm1
          - m2sq*pow(m3sq,2)*ln3*pow(lm1,2) + 2./3.*m2sq*pow(m3sq,2)*
         ln2*pow(lm1,2) + 2./3.*pow(m2sq,2)*m3sq*ln3*pow(lm1,2) - pow(
         m2sq,2)*m3sq*ln2*pow(lm1,2) + 1./3.*pow(m2sq,3)*ln2*pow(lm1,2)
          - 1./6.*m1sq*lm1 + 1./3.*m1sq*ln1*lm1 + 1./3.*m1sq*pow(
         m3sq,2)*pow(lm1,2) - m1sq*pow(m3sq,2)*ln3*pow(lm1,2) + 8./3.*
         m1sq*m2sq*m3sq*pow(lm1,2) - 7./3.*m1sq*m2sq*m3sq*ln3*pow(
         lm1,2) - 7./3.*m1sq*m2sq*m3sq*ln2*pow(lm1,2) + 4*m1sq*m2sq*
         m3sq*ln1*pow(lm1,2) - 5*m1sq*m2sq*pow(m3sq,3)*ln3*pow(lm1,3)
          + 5*m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 1./3.*m1sq*pow(
         m2sq,2)*pow(lm1,2) - m1sq*pow(m2sq,2)*ln2*pow(lm1,2) + 5*m1sq*
         pow(m2sq,3)*m3sq*ln3*pow(lm1,3) - 5*m1sq*pow(m2sq,3)*m3sq*ln2*
         pow(lm1,3) - 2./3.*pow(m1sq,2)*m3sq*pow(lm1,2) + 2./3.*pow(
         m1sq,2)*m3sq*ln3*pow(lm1,2) + 1./3.*pow(m1sq,2)*m3sq*ln1*pow(
         lm1,2);
    h10p +=  - 2./3.*pow(m1sq,2)*m2sq*pow(lm1,2) + 2./3.*pow(
         m1sq,2)*m2sq*ln2*pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*ln1*pow(
         lm1,2) - 10*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) + 10*
         pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) - 10*pow(m1sq,2)*
         pow(m2sq,2)*m3sq*ln3*pow(lm1,3) + 10*pow(m1sq,2)*pow(m2sq,2)*
         m3sq*ln1*pow(lm1,3) + 1./3.*pow(m1sq,3)*pow(lm1,2) - 1./3.*
         pow(m1sq,3)*ln1*pow(lm1,2) + 5*pow(m1sq,3)*m2sq*m3sq*ln3*pow(
         lm1,3) + 5*pow(m1sq,3)*m2sq*m3sq*ln2*pow(lm1,3) - 10*pow(
         m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3);
    break;
  case 3:
    h10p =
       + psiq(m1,m2,m3) * (  - pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3) - 2*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 5*pow(m1sq,2)*m2sq*pow(
         m3sq,3)*pow(lm1,4) + 5*pow(m1sq,2)*pow(m2sq,3)*m3sq*pow(lm1,4)
          + pow(m1sq,3)*m3sq*pow(lm1,3) - 10*pow(m1sq,3)*pow(m2sq,2)*
         m3sq*pow(lm1,4) + 5*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) );

    h10p +=  - 1./6.*m3sq*lm1 - 1./6.*m3sq*ln3*lm1 - 1./6.*m3sq*ln2
         *lm1 + 1./3.*pow(m3sq,3)*ln3*pow(lm1,2) + 1./6.*m2sq*lm1 + 1./
         3.*m2sq*ln2*lm1 - 2./3.*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) - 1./3.
         *m2sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./3.*pow(m2sq,2)*m3sq*ln3*
         pow(lm1,2) + 2./3.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) - 1./3.*
         pow(m2sq,3)*ln2*pow(lm1,2) - 1./6.*m1sq*lm1 - 1./3.*m1sq*ln2*
         lm1 + 4./3.*m1sq*pow(m3sq,2)*pow(lm1,2) - 4./3.*m1sq*pow(
         m3sq,2)*ln3*pow(lm1,2) + m1sq*pow(m3sq,2)*ln2*pow(lm1,2) - 
         m1sq*m2sq*m3sq*pow(lm1,2) + 7./3.*m1sq*m2sq*m3sq*ln3*pow(
         lm1,2) - 3*m1sq*m2sq*m3sq*ln2*pow(lm1,2) - 4*m1sq*m2sq*pow(
         m3sq,3)*ln3*pow(lm1,3) + 4*m1sq*m2sq*pow(m3sq,3)*ln2*pow(
         lm1,3) - 1./3.*m1sq*pow(m2sq,2)*pow(lm1,2) + m1sq*pow(m2sq,2)*
         ln2*pow(lm1,2) + 8*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3)
          - 8*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - 4*m1sq*pow(
         m2sq,3)*m3sq*ln3*pow(lm1,3) + 4*m1sq*pow(m2sq,3)*m3sq*ln2*pow(
         lm1,3);
    h10p +=  - pow(m1sq,2)*m3sq*pow(lm1,2) - 5./3.*pow(m1sq,2)*m3sq
         *ln3*pow(lm1,2) - pow(m1sq,2)*m3sq*ln2*pow(lm1,2) + 7./3.*pow(
         m1sq,2)*m3sq*ln1*pow(lm1,2) - pow(m1sq,2)*pow(m3sq,3)*ln3*pow(
         lm1,3) + pow(m1sq,2)*pow(m3sq,3)*ln1*pow(lm1,3) + 2./3.*pow(
         m1sq,2)*m2sq*pow(lm1,2) - 2./3.*pow(m1sq,2)*m2sq*ln2*pow(
         lm1,2) - 1./3.*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) - 10*pow(
         m1sq,2)*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 2*pow(m1sq,2)*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,3) + 8*pow(m1sq,2)*m2sq*pow(m3sq,2)*
         ln1*pow(lm1,3) + 7*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,3)
          + 2*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 9*pow(
         m1sq,2)*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) - 1./3.*pow(m1sq,3)*
         pow(lm1,2) + 1./3.*pow(m1sq,3)*ln1*pow(lm1,2) + 2*pow(m1sq,3)*
         pow(m3sq,2)*ln3*pow(lm1,3) - 2*pow(m1sq,3)*pow(m3sq,2)*ln1*
         pow(lm1,3) - 2*pow(m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) - 6*pow(
         m1sq,3)*m2sq*m3sq*ln2*pow(lm1,3) + 8*pow(m1sq,3)*m2sq*m3sq*ln1
         *pow(lm1,3);
    h10p +=  - pow(m1sq,4)*m3sq*ln3*pow(lm1,3) + pow(m1sq,4)*m3sq*
         ln1*pow(lm1,3);
    break;
  case 4:
    h10p =
       + psiq(m1,m2,m3) * (  - 2*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) + 5*
         pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(lm1,4) - pow(m1sq,2)*pow(
         m2sq,2)*pow(lm1,3) - 5*pow(m1sq,2)*pow(m2sq,3)*m3sq*pow(lm1,4)
          + pow(m1sq,3)*m2sq*pow(lm1,3) - 10*pow(m1sq,3)*m2sq*pow(
         m3sq,2)*pow(lm1,4) + 5*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) );

    h10p +=  + 1./6.*m3sq*lm1 + 1./3.*m3sq*ln3*lm1 - 1./3.*pow(
         m3sq,3)*ln3*pow(lm1,2) - 1./6.*m2sq*lm1 - 1./6.*m2sq*ln3*lm1
          - 1./6.*m2sq*ln2*lm1 + 2./3.*m2sq*pow(m3sq,2)*ln3*pow(lm1,2)
          + 1./3.*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) - 1./3.*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,2) - 2./3.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + 
         1./3.*pow(m2sq,3)*ln2*pow(lm1,2) - 1./6.*m1sq*lm1 - 1./3.*m1sq
         *ln3*lm1 - 1./3.*m1sq*pow(m3sq,2)*pow(lm1,2) + m1sq*pow(
         m3sq,2)*ln3*pow(lm1,2) - m1sq*m2sq*m3sq*pow(lm1,2) - 3*m1sq*
         m2sq*m3sq*ln3*pow(lm1,2) + 7./3.*m1sq*m2sq*m3sq*ln2*pow(lm1,2)
          + 4*m1sq*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) - 4*m1sq*m2sq*pow(
         m3sq,3)*ln2*pow(lm1,3) + 4./3.*m1sq*pow(m2sq,2)*pow(lm1,2) + 
         m1sq*pow(m2sq,2)*ln3*pow(lm1,2) - 4./3.*m1sq*pow(m2sq,2)*ln2*
         pow(lm1,2) - 8*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 8
         *m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 4*m1sq*pow(
         m2sq,3)*m3sq*ln3*pow(lm1,3) - 4*m1sq*pow(m2sq,3)*m3sq*ln2*pow(
         lm1,3);
    h10p +=  + 2./3.*pow(m1sq,2)*m3sq*pow(lm1,2) - 2./3.*pow(
         m1sq,2)*m3sq*ln3*pow(lm1,2) - 1./3.*pow(m1sq,2)*m3sq*ln1*pow(
         lm1,2) - pow(m1sq,2)*m2sq*pow(lm1,2) - pow(m1sq,2)*m2sq*ln3*
         pow(lm1,2) - 5./3.*pow(m1sq,2)*m2sq*ln2*pow(lm1,2) + 7./3.*
         pow(m1sq,2)*m2sq*ln1*pow(lm1,2) + 2*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,3) + 7*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,3) - 9*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 2
         *pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 10*pow(m1sq,2)*
         pow(m2sq,2)*m3sq*ln2*pow(lm1,3) + 8*pow(m1sq,2)*pow(m2sq,2)*
         m3sq*ln1*pow(lm1,3) - pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,3)
          + pow(m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,3) - 1./3.*pow(m1sq,3)*
         pow(lm1,2) + 1./3.*pow(m1sq,3)*ln1*pow(lm1,2) - 6*pow(m1sq,3)*
         m2sq*m3sq*ln3*pow(lm1,3) - 2*pow(m1sq,3)*m2sq*m3sq*ln2*pow(
         lm1,3) + 8*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3) + 2*pow(
         m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,3) - 2*pow(m1sq,3)*pow(m2sq,2)
         *ln1*pow(lm1,3);
    h10p +=  - pow(m1sq,4)*m2sq*ln2*pow(lm1,3) + pow(m1sq,4)*m2sq*
         ln1*pow(lm1,3);
    break;
  case 5:
    h10p =
       + psiq(m1,m2,m3) * (  - 2*m1sq*pow(m3sq,2)*pow(lm1,3) - 4*m1sq*
         m2sq*m3sq*pow(lm1,3) - 10*m1sq*m2sq*pow(m3sq,3)*pow(lm1,4) + 
         10*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) + 3*pow(m1sq,2)*m3sq*pow(
         lm1,3) - 5*pow(m1sq,2)*pow(m3sq,3)*pow(lm1,4) - 15*pow(m1sq,2)
         *m2sq*pow(m3sq,2)*pow(lm1,4) - 35*pow(m1sq,2)*m2sq*pow(m3sq,4)
         *pow(lm1,5) - 40*pow(m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,4) - 35*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,5) + 35*pow(
         m1sq,2)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) + 35*pow(m1sq,2)*
         pow(m2sq,4)*m3sq*pow(lm1,5) + 10*pow(m1sq,3)*pow(m3sq,2)*pow(
         lm1,4) + 35*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) + 35*pow(m1sq,3)*
         m2sq*pow(m3sq,3)*pow(lm1,5) - 70*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,5) - 105*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(
         lm1,5) - 5*pow(m1sq,4)*m3sq*pow(lm1,4) + 35*pow(m1sq,4)*m2sq*
         pow(m3sq,2)*pow(lm1,5) + 105*pow(m1sq,4)*pow(m2sq,2)*m3sq*pow(
         lm1,5) - 35*pow(m1sq,5)*m2sq*m3sq*pow(lm1,5) );

    h10p +=  - 1./6.*lm1 - 1./3.*ln2*lm1 + pow(m3sq,2)*pow(lm1,2)
          - 5./3.*pow(m3sq,2)*ln3*pow(lm1,2) + 2./3.*pow(m3sq,2)*ln2*
         pow(lm1,2) + 4./3.*pow(m3sq,4)*ln3*pow(lm1,3) - m2sq*m3sq*pow(
         lm1,2) + 2*m2sq*m3sq*ln3*pow(lm1,2) - 8./3.*m2sq*m3sq*ln2*pow(
         lm1,2) - 16./3.*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) + 8./3.*m2sq*
         pow(m3sq,3)*ln2*pow(lm1,3) + 5./3.*pow(m2sq,2)*ln2*pow(lm1,2)
          + 20./3.*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 20./3.*pow(
         m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - 8./3.*pow(m2sq,3)*m3sq*
         ln3*pow(lm1,3) + 16./3.*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 4./3.
         *pow(m2sq,4)*ln2*pow(lm1,3) + 1./3.*m1sq*m3sq*pow(lm1,2) - 3*
         m1sq*m3sq*ln3*pow(lm1,2) - 7./3.*m1sq*m3sq*ln2*pow(lm1,2) + 14.
         /3.*m1sq*m3sq*ln1*pow(lm1,2) + 19./3.*m1sq*pow(m3sq,3)*pow(
         lm1,3) - 29./3.*m1sq*pow(m3sq,3)*ln3*pow(lm1,3) + 5*m1sq*pow(
         m3sq,3)*ln2*pow(lm1,3) + 2*m1sq*pow(m3sq,3)*ln1*pow(lm1,3) + 1.
         /3.*m1sq*m2sq*pow(lm1,2) - 8./3.*m1sq*m2sq*ln2*pow(lm1,2) - 2./
         3.*m1sq*m2sq*ln1*pow(lm1,2);
    h10p +=  + 28./3.*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 40./3.*
         m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) - 8./3.*m1sq*m2sq*pow(
         m3sq,2)*ln2*pow(lm1,3) + 16*m1sq*m2sq*pow(m3sq,2)*ln1*pow(
         lm1,3) - 30*m1sq*m2sq*pow(m3sq,4)*ln3*pow(lm1,4) + 30*m1sq*
         m2sq*pow(m3sq,4)*ln2*pow(lm1,4) - 43./3.*m1sq*pow(m2sq,2)*m3sq
         *pow(lm1,3) + 25*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 29./3.
         *m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 18*m1sq*pow(m2sq,2)*
         m3sq*ln1*pow(lm1,3) + 30*m1sq*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(
         lm1,4) - 30*m1sq*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) - 4./3.
         *m1sq*pow(m2sq,3)*pow(lm1,3) + 16./3.*m1sq*pow(m2sq,3)*ln2*
         pow(lm1,3) + 30*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 
         30*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) - 30*m1sq*pow(
         m2sq,4)*m3sq*ln3*pow(lm1,4) + 30*m1sq*pow(m2sq,4)*m3sq*ln2*
         pow(lm1,4) - 1./3.*pow(m1sq,2)*pow(lm1,2) + 2./3.*pow(m1sq,2)*
         ln2*pow(lm1,2) + pow(m1sq,2)*ln1*pow(lm1,2) - 34./3.*pow(
         m1sq,2)*pow(m3sq,2)*pow(lm1,3);
    h10p +=  + 11./3.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 10*
         pow(m1sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 19./3.*pow(m1sq,2)*
         pow(m3sq,2)*ln1*pow(lm1,3) - 5*pow(m1sq,2)*pow(m3sq,4)*ln3*
         pow(lm1,4) + 5*pow(m1sq,2)*pow(m3sq,4)*ln1*pow(lm1,4) + 32./3.
         *pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 27*pow(m1sq,2)*m2sq*m3sq*
         ln3*pow(lm1,3) - 38./3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3)
          + 37*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) - 45*pow(m1sq,2)*
         m2sq*pow(m3sq,3)*ln3*pow(lm1,4) - 20*pow(m1sq,2)*m2sq*pow(
         m3sq,3)*ln2*pow(lm1,4) + 65*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln1*
         pow(lm1,4) + 4*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 20./3.*
         pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) - 4./3.*pow(m1sq,2)*
         pow(m2sq,2)*ln1*pow(lm1,3) - 75*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,4) + 80*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*ln2*pow(lm1,4) - 5*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)
         *ln1*pow(lm1,4) + 85*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(
         lm1,4);
    h10p +=  - 20*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) - 65*
         pow(m1sq,2)*pow(m2sq,3)*m3sq*ln1*pow(lm1,4) + 11./3.*pow(
         m1sq,3)*m3sq*pow(lm1,3) + 14./3.*pow(m1sq,3)*m3sq*ln3*pow(
         lm1,3) + 5*pow(m1sq,3)*m3sq*ln2*pow(lm1,3) - 7*pow(m1sq,3)*
         m3sq*ln1*pow(lm1,3) + 15*pow(m1sq,3)*pow(m3sq,3)*ln3*pow(
         lm1,4) - 15*pow(m1sq,3)*pow(m3sq,3)*ln1*pow(lm1,4) - 4*pow(
         m1sq,3)*m2sq*pow(lm1,3) + 8./3.*pow(m1sq,3)*m2sq*ln2*pow(
         lm1,3) + 8./3.*pow(m1sq,3)*m2sq*ln1*pow(lm1,3) + 60*pow(
         m1sq,3)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) - 50*pow(m1sq,3)*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,4) - 10*pow(m1sq,3)*m2sq*pow(m3sq,2)*
         ln1*pow(lm1,4) - 75*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(
         lm1,4) - 50*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 125*
         pow(m1sq,3)*pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 4./3.*pow(
         m1sq,4)*pow(lm1,3) - 4./3.*pow(m1sq,4)*ln1*pow(lm1,3) - 15*
         pow(m1sq,4)*pow(m3sq,2)*ln3*pow(lm1,4) + 15*pow(m1sq,4)*pow(
         m3sq,2)*ln1*pow(lm1,4);
    h10p +=  + 15*pow(m1sq,4)*m2sq*m3sq*ln3*pow(lm1,4) + 40*pow(
         m1sq,4)*m2sq*m3sq*ln2*pow(lm1,4) - 55*pow(m1sq,4)*m2sq*m3sq*
         ln1*pow(lm1,4) + 5*pow(m1sq,5)*m3sq*ln3*pow(lm1,4) - 5*pow(
         m1sq,5)*m3sq*ln1*pow(lm1,4);
    break;
  case 6:
    h10p =
       + psiq(m1,m2,m3) * (  - 4*m1sq*m2sq*m3sq*pow(lm1,3) + 10*m1sq*
         m2sq*pow(m3sq,3)*pow(lm1,4) - 2*m1sq*pow(m2sq,2)*pow(lm1,3) - 
         10*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) + 3*pow(m1sq,2)*m2sq*pow(
         lm1,3) - 40*pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) + 35*pow(
         m1sq,2)*m2sq*pow(m3sq,4)*pow(lm1,5) - 15*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,4) + 35*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,3)*pow(lm1,5) - 5*pow(m1sq,2)*pow(m2sq,3)*pow(lm1,4) - 35
         *pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) - 35*pow(
         m1sq,2)*pow(m2sq,4)*m3sq*pow(lm1,5) + 35*pow(m1sq,3)*m2sq*m3sq
         *pow(lm1,4) - 105*pow(m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5) + 10
         *pow(m1sq,3)*pow(m2sq,2)*pow(lm1,4) - 70*pow(m1sq,3)*pow(
         m2sq,2)*pow(m3sq,2)*pow(lm1,5) + 35*pow(m1sq,3)*pow(m2sq,3)*
         m3sq*pow(lm1,5) - 5*pow(m1sq,4)*m2sq*pow(lm1,4) + 105*pow(
         m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5) + 35*pow(m1sq,4)*pow(
         m2sq,2)*m3sq*pow(lm1,5) - 35*pow(m1sq,5)*m2sq*m3sq*pow(lm1,5)
          );

    h10p +=  - 1./6.*lm1 - 1./3.*ln3*lm1 + 5./3.*pow(m3sq,2)*ln3*
         pow(lm1,2) - 4./3.*pow(m3sq,4)*ln3*pow(lm1,3) - m2sq*m3sq*pow(
         lm1,2) - 8./3.*m2sq*m3sq*ln3*pow(lm1,2) + 2*m2sq*m3sq*ln2*pow(
         lm1,2) + 16./3.*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) - 8./3.*m2sq*
         pow(m3sq,3)*ln2*pow(lm1,3) + pow(m2sq,2)*pow(lm1,2) + 2./3.*
         pow(m2sq,2)*ln3*pow(lm1,2) - 5./3.*pow(m2sq,2)*ln2*pow(lm1,2)
          - 20./3.*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 20./3.*pow(
         m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 8./3.*pow(m2sq,3)*m3sq*
         ln3*pow(lm1,3) - 16./3.*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) + 4./3.
         *pow(m2sq,4)*ln2*pow(lm1,3) + 1./3.*m1sq*m3sq*pow(lm1,2) - 8./
         3.*m1sq*m3sq*ln3*pow(lm1,2) - 2./3.*m1sq*m3sq*ln1*pow(lm1,2)
          - 4./3.*m1sq*pow(m3sq,3)*pow(lm1,3) + 16./3.*m1sq*pow(m3sq,3)
         *ln3*pow(lm1,3) + 1./3.*m1sq*m2sq*pow(lm1,2) - 7./3.*m1sq*m2sq
         *ln3*pow(lm1,2) - 3*m1sq*m2sq*ln2*pow(lm1,2) + 14./3.*m1sq*
         m2sq*ln1*pow(lm1,2) - 43./3.*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3)
          - 29./3.*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3);
    h10p +=  + 25*m1sq*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) - 18*m1sq*
         m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 30*m1sq*m2sq*pow(m3sq,4)*ln3
         *pow(lm1,4) - 30*m1sq*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) + 28./3.
         *m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) - 8./3.*m1sq*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,3) - 40./3.*m1sq*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,3) + 16*m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) - 30*m1sq*
         pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 30*m1sq*pow(m2sq,2)*
         pow(m3sq,3)*ln2*pow(lm1,4) + 19./3.*m1sq*pow(m2sq,3)*pow(
         lm1,3) + 5*m1sq*pow(m2sq,3)*ln3*pow(lm1,3) - 29./3.*m1sq*pow(
         m2sq,3)*ln2*pow(lm1,3) + 2*m1sq*pow(m2sq,3)*ln1*pow(lm1,3) - 
         30*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) + 30*m1sq*pow(
         m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) + 30*m1sq*pow(m2sq,4)*m3sq*
         ln3*pow(lm1,4) - 30*m1sq*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 1./
         3.*pow(m1sq,2)*pow(lm1,2) + 2./3.*pow(m1sq,2)*ln3*pow(lm1,2)
          + pow(m1sq,2)*ln1*pow(lm1,2) + 4*pow(m1sq,2)*pow(m3sq,2)*pow(
         lm1,3);
    h10p +=  - 20./3.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 4./3.
         *pow(m1sq,2)*pow(m3sq,2)*ln1*pow(lm1,3) + 32./3.*pow(m1sq,2)*
         m2sq*m3sq*pow(lm1,3) - 38./3.*pow(m1sq,2)*m2sq*m3sq*ln3*pow(
         lm1,3) - 27*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3) + 37*pow(
         m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) - 20*pow(m1sq,2)*m2sq*pow(
         m3sq,3)*ln3*pow(lm1,4) + 85*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln2*
         pow(lm1,4) - 65*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln1*pow(lm1,4) - 
         34./3.*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 10*pow(m1sq,2)*
         pow(m2sq,2)*ln3*pow(lm1,3) + 11./3.*pow(m1sq,2)*pow(m2sq,2)*
         ln2*pow(lm1,3) + 19./3.*pow(m1sq,2)*pow(m2sq,2)*ln1*pow(lm1,3)
          + 80*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,4) - 75*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4) - 5*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) - 20*pow(
         m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) - 45*pow(m1sq,2)*pow(
         m2sq,3)*m3sq*ln2*pow(lm1,4) + 65*pow(m1sq,2)*pow(m2sq,3)*m3sq*
         ln1*pow(lm1,4);
    h10p +=  - 5*pow(m1sq,2)*pow(m2sq,4)*ln2*pow(lm1,4) + 5*pow(
         m1sq,2)*pow(m2sq,4)*ln1*pow(lm1,4) - 4*pow(m1sq,3)*m3sq*pow(
         lm1,3) + 8./3.*pow(m1sq,3)*m3sq*ln3*pow(lm1,3) + 8./3.*pow(
         m1sq,3)*m3sq*ln1*pow(lm1,3) + 11./3.*pow(m1sq,3)*m2sq*pow(
         lm1,3) + 5*pow(m1sq,3)*m2sq*ln3*pow(lm1,3) + 14./3.*pow(
         m1sq,3)*m2sq*ln2*pow(lm1,3) - 7*pow(m1sq,3)*m2sq*ln1*pow(
         lm1,3) - 50*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) - 75*
         pow(m1sq,3)*m2sq*pow(m3sq,2)*ln2*pow(lm1,4) + 125*pow(m1sq,3)*
         m2sq*pow(m3sq,2)*ln1*pow(lm1,4) - 50*pow(m1sq,3)*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,4) + 60*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,4) - 10*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 15*
         pow(m1sq,3)*pow(m2sq,3)*ln2*pow(lm1,4) - 15*pow(m1sq,3)*pow(
         m2sq,3)*ln1*pow(lm1,4) + 4./3.*pow(m1sq,4)*pow(lm1,3) - 4./3.*
         pow(m1sq,4)*ln1*pow(lm1,3) + 40*pow(m1sq,4)*m2sq*m3sq*ln3*pow(
         lm1,4) + 15*pow(m1sq,4)*m2sq*m3sq*ln2*pow(lm1,4) - 55*pow(
         m1sq,4)*m2sq*m3sq*ln1*pow(lm1,4);
    h10p +=  - 15*pow(m1sq,4)*pow(m2sq,2)*ln2*pow(lm1,4) + 15*pow(
         m1sq,4)*pow(m2sq,2)*ln1*pow(lm1,4) + 5*pow(m1sq,5)*m2sq*ln2*
         pow(lm1,4) - 5*pow(m1sq,5)*m2sq*ln1*pow(lm1,4);
    break;
  case 7:

    h10p =
       + psiq(m1,m2,m3) * (  - 2*pow(m1sq,2)*m3sq*pow(lm1,3) + 5*pow(
         m1sq,2)*pow(m3sq,3)*pow(lm1,4) - 2*pow(m1sq,2)*m2sq*pow(lm1,3)
          - 10*pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) + 35*pow(m1sq,2)
         *m2sq*pow(m3sq,4)*pow(lm1,5) - 10*pow(m1sq,2)*pow(m2sq,2)*m3sq
         *pow(lm1,4) - 35*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*pow(
         lm1,5) + 5*pow(m1sq,2)*pow(m2sq,3)*pow(lm1,4) - 35*pow(m1sq,2)
         *pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) + 35*pow(m1sq,2)*pow(
         m2sq,4)*m3sq*pow(lm1,5) + pow(m1sq,3)*pow(lm1,3) - 10*pow(
         m1sq,3)*pow(m3sq,2)*pow(lm1,4) - 5*pow(m1sq,3)*m2sq*m3sq*pow(
         lm1,4) - 35*pow(m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5) - 10*pow(
         m1sq,3)*pow(m2sq,2)*pow(lm1,4) + 70*pow(m1sq,3)*pow(m2sq,2)*
         pow(m3sq,2)*pow(lm1,5) - 35*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(
         lm1,5) + 5*pow(m1sq,4)*m3sq*pow(lm1,4) + 5*pow(m1sq,4)*m2sq*
         pow(lm1,4) - 35*pow(m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5) - 35*
         pow(m1sq,4)*pow(m2sq,2)*m3sq*pow(lm1,5) + 35*pow(m1sq,5)*m2sq*
         m3sq*pow(lm1,5) );

    h10p +=  - 1./3.*lm1 - 1./6.*ln3*lm1 - 1./6.*ln2*lm1 + 2./3.*
         pow(m3sq,2)*pow(lm1,2) + 4./3.*pow(m3sq,2)*ln3*pow(lm1,2) + 1./
         3.*pow(m3sq,2)*ln2*pow(lm1,2) - 4./3.*pow(m3sq,4)*ln3*pow(
         lm1,3) - 4./3.*m2sq*m3sq*pow(lm1,2) - 5./3.*m2sq*m3sq*ln3*pow(
         lm1,2) - 5./3.*m2sq*m3sq*ln2*pow(lm1,2) + 4*m2sq*pow(m3sq,3)*
         ln3*pow(lm1,3) + 4./3.*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 2./3.
         *pow(m2sq,2)*pow(lm1,2) + 1./3.*pow(m2sq,2)*ln3*pow(lm1,2) + 4.
         /3.*pow(m2sq,2)*ln2*pow(lm1,2) - 4*pow(m2sq,2)*pow(m3sq,2)*ln3
         *pow(lm1,3) - 4*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 4./3.
         *pow(m2sq,3)*m3sq*ln3*pow(lm1,3) + 4*pow(m2sq,3)*m3sq*ln2*pow(
         lm1,3) - 4./3.*pow(m2sq,4)*ln2*pow(lm1,3) + 4./3.*m1sq*m3sq*
         pow(lm1,2) - 3*m1sq*m3sq*ln3*pow(lm1,2) + 7./3.*m1sq*m3sq*ln2*
         pow(lm1,2) - 16./3.*m1sq*pow(m3sq,3)*pow(lm1,3) + 20./3.*m1sq*
         pow(m3sq,3)*ln3*pow(lm1,3) - 4*m1sq*pow(m3sq,3)*ln2*pow(lm1,3)
          + 4./3.*m1sq*m2sq*pow(lm1,2) + 7./3.*m1sq*m2sq*ln3*pow(lm1,2)
          - 3*m1sq*m2sq*ln2*pow(lm1,2);
    h10p +=  + 16./3.*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 88./3.*
         m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 80./3.*m1sq*m2sq*pow(
         m3sq,2)*ln2*pow(lm1,3) + 24*m1sq*m2sq*pow(m3sq,4)*ln3*pow(
         lm1,4) - 24*m1sq*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) + 16./3.*m1sq
         *pow(m2sq,2)*m3sq*pow(lm1,3) + 80./3.*m1sq*pow(m2sq,2)*m3sq*
         ln3*pow(lm1,3) - 88./3.*m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3)
          - 72*m1sq*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 72*m1sq*
         pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) - 16./3.*m1sq*pow(
         m2sq,3)*pow(lm1,3) - 4*m1sq*pow(m2sq,3)*ln3*pow(lm1,3) + 20./3.
         *m1sq*pow(m2sq,3)*ln2*pow(lm1,3) + 72*m1sq*pow(m2sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,4) - 72*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln2*
         pow(lm1,4) - 24*m1sq*pow(m2sq,4)*m3sq*ln3*pow(lm1,4) + 24*m1sq
         *pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 3*pow(m1sq,2)*pow(lm1,2) - 
         5./3.*pow(m1sq,2)*ln3*pow(lm1,2) - 5./3.*pow(m1sq,2)*ln2*pow(
         lm1,2) + 7./3.*pow(m1sq,2)*ln1*pow(lm1,2) + 25./3.*pow(m1sq,2)
         *pow(m3sq,2)*pow(lm1,3);
    h10p +=  + 1./3.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 7*
         pow(m1sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - 22./3.*pow(m1sq,2)*
         pow(m3sq,2)*ln1*pow(lm1,3) + 6*pow(m1sq,2)*pow(m3sq,4)*ln3*
         pow(lm1,4) - 6*pow(m1sq,2)*pow(m3sq,4)*ln1*pow(lm1,4) - 62./3.
         *pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 40./3.*pow(m1sq,2)*m2sq*
         m3sq*ln3*pow(lm1,3) - 37./3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(
         lm1,3) + 77./3.*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) + 40*pow(
         m1sq,2)*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) + 7*pow(m1sq,2)*m2sq*
         pow(m3sq,3)*ln2*pow(lm1,4) - 47*pow(m1sq,2)*m2sq*pow(m3sq,3)*
         ln1*pow(lm1,4) + 25./3.*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) + 7
         *pow(m1sq,2)*pow(m2sq,2)*ln3*pow(lm1,3) + 4./3.*pow(m1sq,2)*
         pow(m2sq,2)*ln2*pow(lm1,3) - 25./3.*pow(m1sq,2)*pow(m2sq,2)*
         ln1*pow(lm1,3) - 54*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln3*
         pow(lm1,4) - 53*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(
         lm1,4) + 107*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(
         lm1,4);
    h10p +=  + 8*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) + 41*
         pow(m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) - 49*pow(m1sq,2)*
         pow(m2sq,3)*m3sq*ln1*pow(lm1,4) + 5*pow(m1sq,2)*pow(m2sq,4)*
         ln2*pow(lm1,4) - 5*pow(m1sq,2)*pow(m2sq,4)*ln1*pow(lm1,4) - 2./
         3.*pow(m1sq,3)*m3sq*pow(lm1,3) - 14./3.*pow(m1sq,3)*m3sq*ln3*
         pow(lm1,3) - 2*pow(m1sq,3)*m3sq*ln2*pow(lm1,3) + 4*pow(m1sq,3)
         *m3sq*ln1*pow(lm1,3) - 18*pow(m1sq,3)*pow(m3sq,3)*ln3*pow(
         lm1,4) + 18*pow(m1sq,3)*pow(m3sq,3)*ln1*pow(lm1,4) - 2./3.*
         pow(m1sq,3)*m2sq*pow(lm1,3) - 2*pow(m1sq,3)*m2sq*ln3*pow(
         lm1,3) - 17./3.*pow(m1sq,3)*m2sq*ln2*pow(lm1,3) + 5*pow(
         m1sq,3)*m2sq*ln1*pow(lm1,3) - 36*pow(m1sq,3)*m2sq*pow(m3sq,2)*
         ln3*pow(lm1,4) + 53*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln2*pow(
         lm1,4) - 17*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) + 50*
         pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) - 34*pow(m1sq,3)*
         pow(m2sq,2)*m3sq*ln2*pow(lm1,4) - 16*pow(m1sq,3)*pow(m2sq,2)*
         m3sq*ln1*pow(lm1,4);
    h10p +=  - 15*pow(m1sq,3)*pow(m2sq,3)*ln2*pow(lm1,4) + 15*pow(
         m1sq,3)*pow(m2sq,3)*ln1*pow(lm1,4) - 7./3.*pow(m1sq,4)*pow(
         lm1,3) - pow(m1sq,4)*ln3*pow(lm1,3) - pow(m1sq,4)*ln2*pow(
         lm1,3) + 10./3.*pow(m1sq,4)*ln1*pow(lm1,3) + 18*pow(m1sq,4)*
         pow(m3sq,2)*ln3*pow(lm1,4) - 18*pow(m1sq,4)*pow(m3sq,2)*ln1*
         pow(lm1,4) - 28*pow(m1sq,4)*m2sq*m3sq*ln3*pow(lm1,4) - 31*pow(
         m1sq,4)*m2sq*m3sq*ln2*pow(lm1,4) + 59*pow(m1sq,4)*m2sq*m3sq*
         ln1*pow(lm1,4) + 15*pow(m1sq,4)*pow(m2sq,2)*ln2*pow(lm1,4) - 
         15*pow(m1sq,4)*pow(m2sq,2)*ln1*pow(lm1,4) - 6*pow(m1sq,5)*m3sq
         *ln3*pow(lm1,4) + 6*pow(m1sq,5)*m3sq*ln1*pow(lm1,4) - 5*pow(
         m1sq,5)*m2sq*ln2*pow(lm1,4) + 5*pow(m1sq,5)*m2sq*ln1*pow(
         lm1,4);
    break;
  case 8:
    h10p =
       + psiq(m1,m2,m3) * (  - 4*m1sq*m3sq*pow(lm1,3) + 10*m1sq*pow(
         m3sq,3)*pow(lm1,4) - 4*m1sq*m2sq*pow(lm1,3) - 20*m1sq*m2sq*
         pow(m3sq,2)*pow(lm1,4) + 70*m1sq*m2sq*pow(m3sq,4)*pow(lm1,5)
          - 20*m1sq*pow(m2sq,2)*m3sq*pow(lm1,4) - 70*m1sq*pow(m2sq,2)*
         pow(m3sq,3)*pow(lm1,5) + 10*m1sq*pow(m2sq,3)*pow(lm1,4) - 70*
         m1sq*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) + 70*m1sq*pow(m2sq,4)*
         m3sq*pow(lm1,5) + 3*pow(m1sq,2)*pow(lm1,3) - 40*pow(m1sq,2)*
         pow(m3sq,2)*pow(lm1,4) + 35*pow(m1sq,2)*pow(m3sq,4)*pow(lm1,5)
          - 35*pow(m1sq,2)*m2sq*m3sq*pow(lm1,4) - 140*pow(m1sq,2)*m2sq*
         pow(m3sq,3)*pow(lm1,5) + 315*pow(m1sq,2)*m2sq*pow(m3sq,5)*pow(
         lm1,6) - 40*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,4) + 70*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,5) - 140*pow(m1sq,2)*
         pow(m2sq,3)*m3sq*pow(lm1,5) - 630*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,3)*pow(lm1,6) + 35*pow(m1sq,2)*pow(m2sq,4)*pow(lm1,5) + 
         315*pow(m1sq,2)*pow(m2sq,5)*m3sq*pow(lm1,6) + 35*pow(m1sq,3)*
         m3sq*pow(lm1,4) );

    h10p +=  + psiq(m1,m2,m3) * (  - 105*pow(m1sq,3)*pow(m3sq,3)*
         pow(lm1,5) + 35*pow(m1sq,3)*m2sq*pow(lm1,4) - 175*pow(m1sq,3)*
         m2sq*pow(m3sq,2)*pow(lm1,5) - 630*pow(m1sq,3)*m2sq*pow(m3sq,4)
         *pow(lm1,6) - 175*pow(m1sq,3)*pow(m2sq,2)*m3sq*pow(lm1,5) + 
         630*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,6) - 105*pow(
         m1sq,3)*pow(m2sq,3)*pow(lm1,5) + 630*pow(m1sq,3)*pow(m2sq,3)*
         pow(m3sq,2)*pow(lm1,6) - 630*pow(m1sq,3)*pow(m2sq,4)*m3sq*pow(
         lm1,6) - 5*pow(m1sq,4)*pow(lm1,4) + 105*pow(m1sq,4)*pow(
         m3sq,2)*pow(lm1,5) + 280*pow(m1sq,4)*m2sq*m3sq*pow(lm1,5) + 
         105*pow(m1sq,4)*pow(m2sq,2)*pow(lm1,5) - 1260*pow(m1sq,4)*pow(
         m2sq,2)*pow(m3sq,2)*pow(lm1,6) - 35*pow(m1sq,5)*m3sq*pow(
         lm1,5) - 35*pow(m1sq,5)*m2sq*pow(lm1,5) + 630*pow(m1sq,5)*m2sq
         *pow(m3sq,2)*pow(lm1,6) + 630*pow(m1sq,5)*pow(m2sq,2)*m3sq*
         pow(lm1,6) - 315*pow(m1sq,6)*m2sq*m3sq*pow(lm1,6) );

    h10p +=  + 2./3.*m3sq*pow(lm1,2) - 10./3.*m3sq*ln3*pow(lm1,2)
          + 2*m3sq*ln2*pow(lm1,2) - 8./3.*pow(m3sq,3)*pow(lm1,3) + 12*
         pow(m3sq,3)*ln3*pow(lm1,3) - 8./3.*pow(m3sq,3)*ln2*pow(lm1,3)
          - 8*pow(m3sq,5)*ln3*pow(lm1,4) + 2./3.*m2sq*pow(lm1,2) + 2*
         m2sq*ln3*pow(lm1,2) - 10./3.*m2sq*ln2*pow(lm1,2) + 8./3.*m2sq*
         pow(m3sq,2)*pow(lm1,3) - 92./3.*m2sq*pow(m3sq,2)*ln3*pow(
         lm1,3) + 64./3.*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) + 40*m2sq*pow(
         m3sq,4)*ln3*pow(lm1,4) - 16*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) + 
         8./3.*pow(m2sq,2)*m3sq*pow(lm1,3) + 64./3.*pow(m2sq,2)*m3sq*
         ln3*pow(lm1,3) - 92./3.*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 72*
         pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 56*pow(m2sq,2)*pow(
         m3sq,3)*ln2*pow(lm1,4) - 8./3.*pow(m2sq,3)*pow(lm1,3) - 8./3.*
         pow(m2sq,3)*ln3*pow(lm1,3) + 12*pow(m2sq,3)*ln2*pow(lm1,3) + 
         56*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 72*pow(m2sq,3)*
         pow(m3sq,2)*ln2*pow(lm1,4) - 16*pow(m2sq,4)*m3sq*ln3*pow(
         lm1,4);
    h10p +=  + 40*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 8*pow(m2sq,5)*
         ln2*pow(lm1,4) - 3*m1sq*pow(lm1,2) - 3*m1sq*ln3*pow(lm1,2) - 3
         *m1sq*ln2*pow(lm1,2) + 14./3.*m1sq*ln1*pow(lm1,2) + 12*m1sq*
         pow(m3sq,2)*pow(lm1,3) - 59./3.*m1sq*pow(m3sq,2)*ln3*pow(
         lm1,3) + 25*m1sq*pow(m3sq,2)*ln2*pow(lm1,3) - 44./3.*m1sq*pow(
         m3sq,2)*ln1*pow(lm1,3) - 38*m1sq*pow(m3sq,4)*pow(lm1,4) + 66*
         m1sq*pow(m3sq,4)*ln3*pow(lm1,4) - 30*m1sq*pow(m3sq,4)*ln2*pow(
         lm1,4) - 12*m1sq*pow(m3sq,4)*ln1*pow(lm1,4) + 4./3.*m1sq*m2sq*
         m3sq*pow(lm1,3) - 68./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,3) - 62./3.
         *m1sq*m2sq*m3sq*ln2*pow(lm1,3) + 154./3.*m1sq*m2sq*m3sq*ln1*
         pow(lm1,3) - 48*m1sq*m2sq*pow(m3sq,3)*pow(lm1,4) - 110*m1sq*
         m2sq*pow(m3sq,3)*ln3*pow(lm1,4) + 172*m1sq*m2sq*pow(m3sq,3)*
         ln2*pow(lm1,4) - 94*m1sq*m2sq*pow(m3sq,3)*ln1*pow(lm1,4) + 240
         *m1sq*m2sq*pow(m3sq,5)*ln3*pow(lm1,5) - 240*m1sq*m2sq*pow(
         m3sq,5)*ln2*pow(lm1,5) + 12*m1sq*pow(m2sq,2)*pow(lm1,3) + 25*
         m1sq*pow(m2sq,2)*ln3*pow(lm1,3);
    h10p +=  - 53./3.*m1sq*pow(m2sq,2)*ln2*pow(lm1,3) - 50./3.*m1sq
         *pow(m2sq,2)*ln1*pow(lm1,3) + 172*m1sq*pow(m2sq,2)*pow(m3sq,2)
         *pow(lm1,4) - 100*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,4)
          - 98*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4) + 214*m1sq*
         pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) - 480*m1sq*pow(m2sq,2)*
         pow(m3sq,4)*ln3*pow(lm1,5) + 480*m1sq*pow(m2sq,2)*pow(m3sq,4)*
         ln2*pow(lm1,5) - 48*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) + 174*
         m1sq*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) - 108*m1sq*pow(m2sq,3)*
         m3sq*ln2*pow(lm1,4) - 98*m1sq*pow(m2sq,3)*m3sq*ln1*pow(lm1,4)
          - 38*m1sq*pow(m2sq,4)*pow(lm1,4) - 30*m1sq*pow(m2sq,4)*ln3*
         pow(lm1,4) + 64*m1sq*pow(m2sq,4)*ln2*pow(lm1,4) - 10*m1sq*pow(
         m2sq,4)*ln1*pow(lm1,4) + 480*m1sq*pow(m2sq,4)*pow(m3sq,2)*ln3*
         pow(lm1,5) - 480*m1sq*pow(m2sq,4)*pow(m3sq,2)*ln2*pow(lm1,5)
          - 240*m1sq*pow(m2sq,5)*m3sq*ln3*pow(lm1,5) + 240*m1sq*pow(
         m2sq,5)*m3sq*ln2*pow(lm1,5) - 49./3.*pow(m1sq,2)*m3sq*pow(
         lm1,3);
    h10p +=  - 32./3.*pow(m1sq,2)*m3sq*ln3*pow(lm1,3) - 27*pow(
         m1sq,2)*m3sq*ln2*pow(lm1,3) + 85./3.*pow(m1sq,2)*m3sq*ln1*pow(
         lm1,3) + 101*pow(m1sq,2)*pow(m3sq,3)*pow(lm1,4) - 90*pow(
         m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 85*pow(m1sq,2)*pow(
         m3sq,3)*ln2*pow(lm1,4) - 11*pow(m1sq,2)*pow(m3sq,3)*ln1*pow(
         lm1,4) + 40*pow(m1sq,2)*pow(m3sq,5)*ln3*pow(lm1,5) - 40*pow(
         m1sq,2)*pow(m3sq,5)*ln1*pow(lm1,5) - 49./3.*pow(m1sq,2)*m2sq*
         pow(lm1,3) - 27*pow(m1sq,2)*m2sq*ln3*pow(lm1,3) - 41./3.*pow(
         m1sq,2)*m2sq*ln2*pow(lm1,3) + 94./3.*pow(m1sq,2)*m2sq*ln1*pow(
         lm1,3) - 121*pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) - pow(
         m1sq,2)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) - 70*pow(m1sq,2)*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,4) + 87*pow(m1sq,2)*m2sq*pow(m3sq,2)*
         ln1*pow(lm1,4) + 150*pow(m1sq,2)*m2sq*pow(m3sq,4)*ln3*pow(
         lm1,5) + 365*pow(m1sq,2)*m2sq*pow(m3sq,4)*ln2*pow(lm1,5) - 515
         *pow(m1sq,2)*m2sq*pow(m3sq,4)*ln1*pow(lm1,5) - 121*pow(m1sq,2)
         *pow(m2sq,2)*m3sq*pow(lm1,4);
    h10p +=  - 82*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) + 11*
         pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 87*pow(m1sq,2)*
         pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 550*pow(m1sq,2)*pow(m2sq,2)*
         pow(m3sq,3)*ln3*pow(lm1,5) - 1110*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,3)*ln2*pow(lm1,5) + 560*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,3)*ln1*pow(lm1,5) + 101*pow(m1sq,2)*pow(m2sq,3)*pow(
         lm1,4) + 85*pow(m1sq,2)*pow(m2sq,3)*ln3*pow(lm1,4) - 78*pow(
         m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,4) - 23*pow(m1sq,2)*pow(
         m2sq,3)*ln1*pow(lm1,4) - 1110*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,5) + 560*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,2)*ln2*pow(lm1,5) + 550*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,2)*ln1*pow(lm1,5) + 370*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln3*
         pow(lm1,5) + 150*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln2*pow(lm1,5)
          - 520*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln1*pow(lm1,5) + 35*pow(
         m1sq,2)*pow(m2sq,5)*ln2*pow(lm1,5) - 35*pow(m1sq,2)*pow(
         m2sq,5)*ln1*pow(lm1,5);
    h10p +=  + 7*pow(m1sq,3)*pow(lm1,3) + 14./3.*pow(m1sq,3)*ln3*
         pow(lm1,3) + 14./3.*pow(m1sq,3)*ln2*pow(lm1,3) - 75*pow(
         m1sq,3)*pow(m3sq,2)*pow(lm1,4) + 19*pow(m1sq,3)*pow(m3sq,2)*
         ln3*pow(lm1,4) - 75*pow(m1sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) + 
         40*pow(m1sq,3)*pow(m3sq,2)*ln1*pow(lm1,4) - 160*pow(m1sq,3)*
         pow(m3sq,4)*ln3*pow(lm1,5) + 160*pow(m1sq,3)*pow(m3sq,4)*ln1*
         pow(lm1,5) + 170*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) - 84*pow(
         m1sq,3)*m2sq*m3sq*ln3*pow(lm1,4) - 102*pow(m1sq,3)*m2sq*m3sq*
         ln2*pow(lm1,4) + 154*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,4) - 
         790*pow(m1sq,3)*m2sq*pow(m3sq,3)*ln3*pow(lm1,5) + 310*pow(
         m1sq,3)*m2sq*pow(m3sq,3)*ln2*pow(lm1,5) + 480*pow(m1sq,3)*m2sq
         *pow(m3sq,3)*ln1*pow(lm1,5) - 75*pow(m1sq,3)*pow(m2sq,2)*pow(
         lm1,4) - 75*pow(m1sq,3)*pow(m2sq,2)*ln3*pow(lm1,4) + pow(
         m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,4) + 58*pow(m1sq,3)*pow(
         m2sq,2)*ln1*pow(lm1,4) + 620*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,5);
    h10p +=  + 640*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,5)
          - 1260*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,5)
          + 290*pow(m1sq,3)*pow(m2sq,3)*m3sq*ln3*pow(lm1,5) - 770*pow(
         m1sq,3)*pow(m2sq,3)*m3sq*ln2*pow(lm1,5) + 480*pow(m1sq,3)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,5) - 140*pow(m1sq,3)*pow(m2sq,4)*ln2*
         pow(lm1,5) + 140*pow(m1sq,3)*pow(m2sq,4)*ln1*pow(lm1,5) - pow(
         m1sq,4)*m3sq*pow(lm1,4) + 8*pow(m1sq,4)*m3sq*ln3*pow(lm1,4) + 
         15*pow(m1sq,4)*m3sq*ln2*pow(lm1,4) + pow(m1sq,4)*m3sq*ln1*pow(
         lm1,4) + 240*pow(m1sq,4)*pow(m3sq,3)*ln3*pow(lm1,5) - 240*pow(
         m1sq,4)*pow(m3sq,3)*ln1*pow(lm1,5) - pow(m1sq,4)*m2sq*pow(
         lm1,4) + 15*pow(m1sq,4)*m2sq*ln3*pow(lm1,4) + 16*pow(m1sq,4)*
         m2sq*ln2*pow(lm1,4) - 7*pow(m1sq,4)*m2sq*ln1*pow(lm1,4) + 170*
         pow(m1sq,4)*m2sq*pow(m3sq,2)*ln3*pow(lm1,5) - 720*pow(m1sq,4)*
         m2sq*pow(m3sq,2)*ln2*pow(lm1,5) + 550*pow(m1sq,4)*m2sq*pow(
         m3sq,2)*ln1*pow(lm1,5) - 690*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln3*
         pow(lm1,5);
    h10p +=  + 130*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln2*pow(lm1,5) + 
         560*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln1*pow(lm1,5) + 210*pow(
         m1sq,4)*pow(m2sq,3)*ln2*pow(lm1,5) - 210*pow(m1sq,4)*pow(
         m2sq,3)*ln1*pow(lm1,5) + 13*pow(m1sq,5)*pow(lm1,4) + 5*pow(
         m1sq,5)*ln3*pow(lm1,4) + 5*pow(m1sq,5)*ln2*pow(lm1,4) - 18*
         pow(m1sq,5)*ln1*pow(lm1,4) - 160*pow(m1sq,5)*pow(m3sq,2)*ln3*
         pow(lm1,5) + 160*pow(m1sq,5)*pow(m3sq,2)*ln1*pow(lm1,5) + 230*
         pow(m1sq,5)*m2sq*m3sq*ln3*pow(lm1,5) + 250*pow(m1sq,5)*m2sq*
         m3sq*ln2*pow(lm1,5) - 480*pow(m1sq,5)*m2sq*m3sq*ln1*pow(lm1,5)
          - 140*pow(m1sq,5)*pow(m2sq,2)*ln2*pow(lm1,5) + 140*pow(
         m1sq,5)*pow(m2sq,2)*ln1*pow(lm1,5) + 40*pow(m1sq,6)*m3sq*ln3*
         pow(lm1,5) - 40*pow(m1sq,6)*m3sq*ln1*pow(lm1,5) + 35*pow(
         m1sq,6)*m2sq*ln2*pow(lm1,5) - 35*pow(m1sq,6)*m2sq*ln1*pow(
         lm1,5);
      break;
  default:
    std::cout <<"h10p (quenched) called with wrong iprop ="<<iprop<<'\n';
    h10p = 0.;
  }

  return h10p*pi162;
}

double h210(const int iprop,const double m1sq, const double m2sq,
	    const double m3sq, const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h210sing(iprop,m1sq,m2sq,m3sq,xmu2);
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double h210;
  switch(iprop){
  case 1:
    h210 = 1./(6.*lm*lm)*
       (-lm*lm+lm*m1sq*dm+2.*pow(m1sq,2)*m2sq*m3sq)*psiq(m1,m2,m3)
     + m1sq*(17./72.-ln1/3.+m1sq*dm*ln1/(6.*lm))
     + m2sq*(pi2/36.+19./54.-ln2/9.
            -m1sq/(6.*lm)*(dm+2.*m3sq)*ln2+1./6.*(ln2*ln4-ln1*ln3))
      + m3sq*(pi2/36.+19./54.-ln3/9.
	      -m1sq/(6.*lm)*(dm+2.*m2sq)*ln3+1./6.*(ln3*ln4-ln1*ln2));
    break;
  case 2:
    h210 =
       + psiq(m1,m2,m3) * ( 1./2.*m3sq*lm*pow(lm1,2) - 1./2.*m3sq*pow(
         lm,2)*pow(lm1,3) + 1./2.*m2sq*lm*pow(lm1,2) - 1./2.*m2sq*pow(
         lm,2)*pow(lm1,3) - 1./3.*m1sq*lm*pow(lm1,2) + 1./2.*m1sq*pow(
         lm,2)*pow(lm1,3) + 1./3.*m1sq*pow(m3sq,2)*pow(lm1,2) - 1./2.*
         m1sq*pow(m3sq,2)*lm*pow(lm1,3) + 4./3.*m1sq*m2sq*m3sq*pow(
         lm1,2) - m1sq*m2sq*m3sq*lm*pow(lm1,3) + 1./3.*m1sq*pow(m2sq,2)
         *pow(lm1,2) - 1./2.*m1sq*pow(m2sq,2)*lm*pow(lm1,3) - 2./3.*
         pow(m1sq,2)*m3sq*pow(lm1,2) + pow(m1sq,2)*m3sq*lm*pow(lm1,3)
          - 2./3.*pow(m1sq,2)*m2sq*pow(lm1,2) + pow(m1sq,2)*m2sq*lm*
         pow(lm1,3) + pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,3) + pow(
         m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,3) + 1./3.*pow(m1sq,3)*pow(
         lm1,2) - 1./2.*pow(m1sq,3)*lm*pow(lm1,3) - pow(m1sq,3)*m2sq*
         m3sq*pow(lm1,3) );

    h210 +=  - 7./72. + 1./6.*pow(m1sq,-1)*m3sq*ln3 - 1./6.*pow(
         m1sq,-1)*m3sq*ln3*pow(lm,2)*pow(lm1,2) - 1./6.*pow(m1sq,-1)*
         m3sq*ln2 + 1./6.*pow(m1sq,-1)*m3sq*ln2*pow(lm,2)*pow(lm1,2) - 
         1./6.*pow(m1sq,-1)*m2sq*ln3 + 1./6.*pow(m1sq,-1)*m2sq*ln3*pow(
         lm,2)*pow(lm1,2) + 1./6.*pow(m1sq,-1)*m2sq*ln2 - 1./6.*pow(
         m1sq,-1)*m2sq*ln2*pow(lm,2)*pow(lm1,2) - 1./6.*ln3*pow(lm,2)*
         pow(lm1,2) - 1./6.*ln2*pow(lm,2)*pow(lm1,2) - 1./3.*ln1 + 1./3.
         *ln1*pow(lm,2)*pow(lm1,2) + 1./6.*pow(m3sq,2)*ln3*lm1 - 1./6.*
         pow(m3sq,2)*ln3*lm*pow(lm1,2) + 1./6.*pow(m3sq,2)*ln2*lm*pow(
         lm1,2) - 1./6.*m2sq*m3sq*ln3*lm1 - 1./6.*m2sq*m3sq*ln2*lm1 + 1.
         /6.*pow(m2sq,2)*ln3*lm*pow(lm1,2) + 1./6.*pow(m2sq,2)*ln2*lm1
          - 1./6.*pow(m2sq,2)*ln2*lm*pow(lm1,2) - 1./6.*m1sq*m3sq*lm1
          - 1./3.*m1sq*m3sq*ln3*lm1 - 1./3.*m1sq*m3sq*ln2*lm*pow(lm1,2)
          - 1./3.*m1sq*m3sq*ln1*lm1 + 1./3.*m1sq*m3sq*ln1*lm*pow(lm1,2)
          + 1./3.*m1sq*pow(m3sq,3)*ln3*pow(lm1,2) - 1./6.*m1sq*m2sq*lm1
          - 1./3.*m1sq*m2sq*ln3*lm*pow(lm1,2);
    h210 +=  - 1./3.*m1sq*m2sq*ln2*lm1 - 1./3.*m1sq*m2sq*ln1*lm1 + 
         1./3.*m1sq*m2sq*ln1*lm*pow(lm1,2) + 1./3.*m1sq*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,2) - 2./3.*m1sq*m2sq*pow(m3sq,2)*ln2*pow(
         lm1,2) - 2./3.*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,2) + 1./3.*
         m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + 1./3.*m1sq*pow(m2sq,3)*
         ln2*pow(lm1,2) + 1./6.*pow(m1sq,2)*lm1 + 1./6.*pow(m1sq,2)*ln3
         *lm*pow(lm1,2) + 1./6.*pow(m1sq,2)*ln2*lm*pow(lm1,2) + 1./2.*
         pow(m1sq,2)*ln1*lm1 - 1./3.*pow(m1sq,2)*ln1*lm*pow(lm1,2) - 2./
         3.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,2) - 1./3.*pow(m1sq,2)*
         pow(m3sq,2)*ln1*pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*m3sq*ln3*
         pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,2) - 4./3.
         *pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,2) - 2./3.*pow(m1sq,2)*pow(
         m2sq,2)*ln2*pow(lm1,2) - 1./3.*pow(m1sq,2)*pow(m2sq,2)*ln1*
         pow(lm1,2) + 1./3.*pow(m1sq,3)*m3sq*ln3*pow(lm1,2) + 2./3.*
         pow(m1sq,3)*m3sq*ln1*pow(lm1,2) + 1./3.*pow(m1sq,3)*m2sq*ln2*
         pow(lm1,2);
    h210 +=  + 2./3.*pow(m1sq,3)*m2sq*ln1*pow(lm1,2) - 1./3.*pow(
         m1sq,4)*ln1*pow(lm1,2);
    break;
  case 3:
    h210 =
       + psiq(m1,m2,m3) * ( 2./3.*m3sq*lm*pow(lm1,2) - 1./2.*m3sq*pow(
         lm,2)*pow(lm1,3) - 2./3.*m2sq*lm*pow(lm1,2) + 1./2.*m2sq*pow(
         lm,2)*pow(lm1,3) + 1./2.*m1sq*lm*pow(lm1,2) - 1./2.*m1sq*pow(
         lm,2)*pow(lm1,3) + 1./3.*m1sq*pow(m3sq,2)*pow(lm1,2) - 1./2.*
         m1sq*pow(m3sq,2)*lm*pow(lm1,3) - 1./3.*m1sq*pow(m2sq,2)*pow(
         lm1,2) + 1./2.*m1sq*pow(m2sq,2)*lm*pow(lm1,3) + 1./3.*pow(
         m1sq,2)*m3sq*pow(lm1,2) + 2./3.*pow(m1sq,2)*m2sq*pow(lm1,2) - 
         pow(m1sq,2)*m2sq*lm*pow(lm1,3) + pow(m1sq,2)*m2sq*pow(m3sq,2)*
         pow(lm1,3) - pow(m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,3) - 1./3.*
         pow(m1sq,3)*pow(lm1,2) + 1./2.*pow(m1sq,3)*lm*pow(lm1,3) + 
         pow(m1sq,3)*m2sq*m3sq*pow(lm1,3) );

    h210 +=  + 13./54. + 1./6.*pow(m2sq,-1)*m3sq*ln3 - 1./6.*pow(
         m2sq,-1)*m3sq*ln3*pow(lm,2)*pow(lm1,2) - 1./6.*pow(m2sq,-1)*
         m3sq*ln1 + 1./6.*pow(m2sq,-1)*m3sq*ln1*pow(lm,2)*pow(lm1,2) + 
         1./6.*ln3 - 1./6.*ln3*pow(lm,2)*pow(lm1,2) + 2./9.*ln2 + 1./3.
         *ln2*pow(lm,2)*pow(lm1,2) + 1./6.*ln2*ln3 + 1./6.*pow(ln2,2)
          + 1./6.*ln1 - 1./6.*ln1*pow(lm,2)*pow(lm1,2) - 1./6.*ln1*ln3
          + 1./6.*ln1*ln2 + 1./6.*m1sq*pow(m2sq,-1)*ln3*pow(lm,2)*pow(
         lm1,2) - 1./6.*m1sq*pow(m2sq,-1)*ln1*pow(lm,2)*pow(lm1,2) - 1./
         6.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln3*lm*pow(lm1,2) + 1./6.*
         m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*pow(lm1,2) - 1./6.*m1sq*
         m3sq*lm1 - 1./6.*m1sq*m3sq*ln3*lm1 - 1./3.*m1sq*m3sq*ln3*lm*
         pow(lm1,2) - 1./6.*m1sq*m3sq*ln2*lm1 + 1./3.*m1sq*m3sq*ln2*lm*
         pow(lm1,2) + 1./3.*m1sq*pow(m3sq,3)*ln3*pow(lm1,2) + 1./6.*
         m1sq*m2sq*lm1 - 1./6.*m1sq*m2sq*ln3*lm*pow(lm1,2) + 1./3.*m1sq
         *m2sq*ln2*lm1 + 1./3.*m1sq*m2sq*ln2*lm*pow(lm1,2) - 1./6.*m1sq
         *m2sq*ln1*lm*pow(lm1,2);
    h210 +=  - 2./3.*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) - 1./3.*
         m1sq*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./3.*m1sq*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,2) + 2./3.*m1sq*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,2) - 1./3.*m1sq*pow(m2sq,3)*ln2*pow(lm1,2) + 1./3.*pow(
         m1sq,2)*pow(m2sq,-1)*m3sq*ln3*lm*pow(lm1,2) - 1./3.*pow(
         m1sq,2)*pow(m2sq,-1)*m3sq*ln1*lm*pow(lm1,2) - 1./6.*pow(
         m1sq,2)*lm1 + 1./3.*pow(m1sq,2)*ln3*lm*pow(lm1,2) - 1./6.*pow(
         m1sq,2)*ln2*lm1 - 1./3.*pow(m1sq,2)*ln2*lm*pow(lm1,2) - 1./6.*
         pow(m1sq,2)*ln1*lm1 + 1./3.*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(
         lm1,2) - 2./3.*pow(m1sq,2)*pow(m3sq,2)*ln1*pow(lm1,2) + 1./3.*
         pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,2) - 4./3.*pow(m1sq,2)*m2sq*
         m3sq*ln2*pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*m3sq*ln1*pow(
         lm1,2) + 2./3.*pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,2) + 1./3.*
         pow(m1sq,2)*pow(m2sq,2)*ln1*pow(lm1,2) - 1./6.*pow(m1sq,3)*
         pow(m2sq,-1)*ln3*lm*pow(lm1,2) + 1./6.*pow(m1sq,3)*pow(
         m2sq,-1)*ln1*lm*pow(lm1,2);
    h210 +=  - 2./3.*pow(m1sq,3)*m3sq*ln3*pow(lm1,2) + 1./3.*pow(
         m1sq,3)*m3sq*ln1*pow(lm1,2) - 1./3.*pow(m1sq,3)*m2sq*ln2*pow(
         lm1,2) - 2./3.*pow(m1sq,3)*m2sq*ln1*pow(lm1,2) + 1./3.*pow(
         m1sq,4)*ln1*pow(lm1,2) + 1./36.*pi2;
    break;
  case 4:
    h210 =
       + psiq(m1,m2,m3) * (  - 2./3.*m3sq*lm*pow(lm1,2) + 1./2.*m3sq*
         pow(lm,2)*pow(lm1,3) + 2./3.*m2sq*lm*pow(lm1,2) - 1./2.*m2sq*
         pow(lm,2)*pow(lm1,3) + 1./2.*m1sq*lm*pow(lm1,2) - 1./2.*m1sq*
         pow(lm,2)*pow(lm1,3) - 1./3.*m1sq*pow(m3sq,2)*pow(lm1,2) + 1./
         2.*m1sq*pow(m3sq,2)*lm*pow(lm1,3) + 1./3.*m1sq*pow(m2sq,2)*
         pow(lm1,2) - 1./2.*m1sq*pow(m2sq,2)*lm*pow(lm1,3) + 2./3.*pow(
         m1sq,2)*m3sq*pow(lm1,2) - pow(m1sq,2)*m3sq*lm*pow(lm1,3) + 1./
         3.*pow(m1sq,2)*m2sq*pow(lm1,2) - pow(m1sq,2)*m2sq*pow(m3sq,2)*
         pow(lm1,3) + pow(m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,3) - 1./3.*
         pow(m1sq,3)*pow(lm1,2) + 1./2.*pow(m1sq,3)*lm*pow(lm1,3) + 
         pow(m1sq,3)*m2sq*m3sq*pow(lm1,3) );

    h210 +=  + 13./54. + 2./9.*ln3 + 1./3.*ln3*pow(lm,2)*pow(lm1,2)
          + 1./6.*pow(ln3,2) + 1./6.*ln2 - 1./6.*ln2*pow(lm,2)*pow(
         lm1,2) + 1./6.*ln2*ln3 + 1./6.*ln1 - 1./6.*ln1*pow(lm,2)*pow(
         lm1,2) + 1./6.*ln1*ln3 - 1./6.*ln1*ln2 + 1./6.*m2sq*pow(
         m3sq,-1)*ln2 - 1./6.*m2sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(
         lm1,2) - 1./6.*m2sq*pow(m3sq,-1)*ln1 + 1./6.*m2sq*pow(m3sq,-1)
         *ln1*pow(lm,2)*pow(lm1,2) + 1./6.*m1sq*pow(m3sq,-1)*ln2*pow(
         lm,2)*pow(lm1,2) - 1./6.*m1sq*pow(m3sq,-1)*ln1*pow(lm,2)*pow(
         lm1,2) + 1./6.*m1sq*m3sq*lm1 + 1./3.*m1sq*m3sq*ln3*lm1 + 1./3.
         *m1sq*m3sq*ln3*lm*pow(lm1,2) - 1./6.*m1sq*m3sq*ln2*lm*pow(
         lm1,2) - 1./6.*m1sq*m3sq*ln1*lm*pow(lm1,2) - 1./3.*m1sq*pow(
         m3sq,3)*ln3*pow(lm1,2) - 1./6.*m1sq*m2sq*lm1 - 1./6.*m1sq*m2sq
         *ln3*lm1 + 1./3.*m1sq*m2sq*ln3*lm*pow(lm1,2) - 1./6.*m1sq*m2sq
         *ln2*lm1 - 1./3.*m1sq*m2sq*ln2*lm*pow(lm1,2) + 2./3.*m1sq*m2sq
         *pow(m3sq,2)*ln3*pow(lm1,2) + 1./3.*m1sq*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,2);
    h210 +=  - 1./6.*m1sq*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm*pow(
      lm1,2) + 1./6.*m1sq*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(lm1,2) - 
         1./3.*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,2) - 2./3.*m1sq*pow(
         m2sq,2)*m3sq*ln2*pow(lm1,2) + 1./3.*m1sq*pow(m2sq,3)*ln2*pow(
         lm1,2) - 1./6.*pow(m1sq,2)*lm1 - 1./6.*pow(m1sq,2)*ln3*lm1 - 1.
         /3.*pow(m1sq,2)*ln3*lm*pow(lm1,2) + 1./3.*pow(m1sq,2)*ln2*lm*
         pow(lm1,2) - 1./6.*pow(m1sq,2)*ln1*lm1 + 2./3.*pow(m1sq,2)*
         pow(m3sq,2)*ln3*pow(lm1,2) + 1./3.*pow(m1sq,2)*pow(m3sq,2)*ln1
         *pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*lm*pow(
         lm1,2) - 1./3.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,2)
          - 4./3.*pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,2) + 1./3.*pow(
         m1sq,2)*m2sq*m3sq*ln2*pow(lm1,2) + 1./3.*pow(m1sq,2)*m2sq*m3sq
         *ln1*pow(lm1,2) + 1./3.*pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,2)
          - 2./3.*pow(m1sq,2)*pow(m2sq,2)*ln1*pow(lm1,2) - 1./6.*pow(
         m1sq,3)*pow(m3sq,-1)*ln2*lm*pow(lm1,2) + 1./6.*pow(m1sq,3)*
         pow(m3sq,-1)*ln1*lm*pow(lm1,2);
    h210 +=  - 1./3.*pow(m1sq,3)*m3sq*ln3*pow(lm1,2) - 2./3.*pow(
         m1sq,3)*m3sq*ln1*pow(lm1,2) - 2./3.*pow(m1sq,3)*m2sq*ln2*pow(
         lm1,2) + 1./3.*pow(m1sq,3)*m2sq*ln1*pow(lm1,2) + 1./3.*pow(
         m1sq,4)*ln1*pow(lm1,2) + 1./36.*pi2;
    break;
  case 5:
    h210 =
       + psiq(m1,m2,m3) * ( 1./2.*lm*pow(lm1,2) - 1./2.*pow(lm,2)*pow(
         lm1,3) - pow(m3sq,2)*pow(lm1,2) + 7./2.*pow(m3sq,2)*lm*pow(
         lm1,3) - 5./2.*pow(m3sq,2)*pow(lm,2)*pow(lm1,4) + pow(m2sq,2)*
         pow(lm1,2) - 7./2.*pow(m2sq,2)*lm*pow(lm1,3) + 5./2.*pow(
         m2sq,2)*pow(lm,2)*pow(lm1,4) + m1sq*m3sq*pow(lm1,2) - 1./2.*
         m1sq*m3sq*lm*pow(lm1,3) + 2*m1sq*pow(m3sq,3)*pow(lm1,3) - 5./2.
         *m1sq*pow(m3sq,3)*lm*pow(lm1,4) - m1sq*m2sq*pow(lm1,2) + 11./2.
         *m1sq*m2sq*lm*pow(lm1,3) - 5*m1sq*m2sq*pow(lm,2)*pow(lm1,4) + 
         4*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 5./2.*m1sq*m2sq*pow(
         m3sq,2)*lm*pow(lm1,4) - 4*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) + 5.
         /2.*m1sq*pow(m2sq,2)*m3sq*lm*pow(lm1,4) - 2*m1sq*pow(m2sq,3)*
         pow(lm1,3) + 5./2.*m1sq*pow(m2sq,3)*lm*pow(lm1,4) - 2*pow(
         m1sq,2)*lm*pow(lm1,3) + 5./2.*pow(m1sq,2)*pow(lm,2)*pow(lm1,4)
          - pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3) + 5./2.*pow(m1sq,2)*pow(
         m3sq,2)*lm*pow(lm1,4) + 8*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 5
         *pow(m1sq,2)*m2sq*m3sq*lm*pow(lm1,4) );

    h210 +=  + psiq(m1,m2,m3) * ( 5*pow(m1sq,2)*m2sq*pow(m3sq,3)*
         pow(lm1,4) + 6*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 15./2.*
         pow(m1sq,2)*pow(m2sq,2)*lm*pow(lm1,4) - 5*pow(m1sq,2)*pow(
         m2sq,3)*m3sq*pow(lm1,4) - 3*pow(m1sq,3)*m3sq*pow(lm1,3) + 5./2.
         *pow(m1sq,3)*m3sq*lm*pow(lm1,4) - 6*pow(m1sq,3)*m2sq*pow(
         lm1,3) + 15./2.*pow(m1sq,3)*m2sq*lm*pow(lm1,4) + 10*pow(
         m1sq,3)*pow(m2sq,2)*m3sq*pow(lm1,4) + 2*pow(m1sq,4)*pow(lm1,3)
          - 5./2.*pow(m1sq,4)*lm*pow(lm1,4) - 5*pow(m1sq,4)*m2sq*m3sq*
         pow(lm1,4) );

    h210 +=  - 1./6.*pow(m1sq,-1)*pow(m2sq,-1)*m3sq + 1./6.*pow(
         m1sq,-1)*pow(m2sq,-1)*m3sq*pow(lm,2)*pow(lm1,2) + 1./6.*pow(
         m1sq,-1) - 1./6.*pow(m1sq,-1)*pow(lm,2)*pow(lm1,2) - 1./6.*
         pow(m1sq,-1)*ln3 + 1./6.*pow(m1sq,-1)*ln3*pow(lm,2)*pow(lm1,2)
          + 1./6.*pow(m1sq,-1)*ln2 - 1./6.*pow(m1sq,-1)*ln2*pow(lm,2)*
         pow(lm1,2) + 2./3.*pow(m1sq,-1)*pow(m3sq,2)*ln3*lm*pow(lm1,2)
          - 2./3.*pow(m1sq,-1)*pow(m3sq,2)*ln3*pow(lm,2)*pow(lm1,3) - 2.
         /3.*pow(m1sq,-1)*pow(m3sq,2)*ln2*lm*pow(lm1,2) + 2./3.*pow(
         m1sq,-1)*pow(m3sq,2)*ln2*pow(lm,2)*pow(lm1,3) - 4./3.*pow(
         m1sq,-1)*m2sq*m3sq*ln3*lm*pow(lm1,2) + 4./3.*pow(m1sq,-1)*m2sq
         *m3sq*ln3*pow(lm,2)*pow(lm1,3) + 4./3.*pow(m1sq,-1)*m2sq*m3sq*
         ln2*lm*pow(lm1,2) - 4./3.*pow(m1sq,-1)*m2sq*m3sq*ln2*pow(lm,2)
         *pow(lm1,3) + 2./3.*pow(m1sq,-1)*pow(m2sq,2)*ln3*lm*pow(lm1,2)
          - 2./3.*pow(m1sq,-1)*pow(m2sq,2)*ln3*pow(lm,2)*pow(lm1,3) - 2.
         /3.*pow(m1sq,-1)*pow(m2sq,2)*ln2*lm*pow(lm1,2) + 2./3.*pow(
         m1sq,-1)*pow(m2sq,2)*ln2*pow(lm,2)*pow(lm1,3);
    h210 +=  - 1./6.*pow(m2sq,-1)*pow(lm,2)*pow(lm1,2) + 1./6.*pow(
         m2sq,-1)*pow(m3sq,2)*lm*pow(lm1,2) + 1./2.*pow(m2sq,-1)*pow(
         m3sq,2)*ln3*lm*pow(lm1,2) - 1./2.*pow(m2sq,-1)*pow(m3sq,2)*ln3
         *pow(lm,2)*pow(lm1,3) - 1./2.*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*
         pow(lm1,2) + 1./2.*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm,2)*pow(
         lm1,3) - 1./6.*m3sq*lm1 - 1./6.*m3sq*ln3*lm1 + 7./3.*m3sq*ln3*
         lm*pow(lm1,2) - 7./3.*m3sq*ln3*pow(lm,2)*pow(lm1,3) - 1./6.*
         m3sq*ln2*lm1 - m3sq*ln2*lm*pow(lm1,2) + m3sq*ln2*pow(lm,2)*
         pow(lm1,3) - 4./3.*m3sq*ln1*lm*pow(lm1,2) + 4./3.*m3sq*ln1*
         pow(lm,2)*pow(lm1,3) + 2./3.*pow(m3sq,3)*ln3*pow(lm1,2) - 2./3.
         *pow(m3sq,3)*ln3*lm*pow(lm1,3) - 1./3.*pow(m3sq,3)*ln2*pow(
         lm1,2) + 2./3.*pow(m3sq,3)*ln2*lm*pow(lm1,3) + 1./6.*m2sq*lm1
          - 1./6.*m2sq*lm*pow(lm1,2) - 1./2.*m2sq*ln3*lm*pow(lm1,2) + 5.
         /6.*m2sq*ln3*pow(lm,2)*pow(lm1,3) + 1./3.*m2sq*ln2*lm1 - 4./3.
         *m2sq*ln2*lm*pow(lm1,2) + m2sq*ln2*pow(lm,2)*pow(lm1,3) + 11./
         6.*m2sq*ln1*lm*pow(lm1,2);
    h210 +=  - 11./6.*m2sq*ln1*pow(lm,2)*pow(lm1,3) - m2sq*pow(
         m3sq,2)*ln3*pow(lm1,2) + 2./3.*m2sq*pow(m3sq,2)*ln3*lm*pow(
         lm1,3) - 2./3.*m2sq*pow(m3sq,2)*ln2*lm*pow(lm1,3) + 2./3.*pow(
         m2sq,2)*m3sq*ln3*lm*pow(lm1,3) + pow(m2sq,2)*m3sq*ln2*pow(
         lm1,2) - 2./3.*pow(m2sq,2)*m3sq*ln2*lm*pow(lm1,3) + 1./3.*pow(
         m2sq,3)*ln3*pow(lm1,2) - 2./3.*pow(m2sq,3)*ln3*lm*pow(lm1,3)
          - 2./3.*pow(m2sq,3)*ln2*pow(lm1,2) + 2./3.*pow(m2sq,3)*ln2*lm
         *pow(lm1,3) - 1./3.*m1sq*pow(m2sq,-1)*m3sq*lm*pow(lm1,2) - 5./
         6.*m1sq*pow(m2sq,-1)*m3sq*ln3*lm*pow(lm1,2) + m1sq*pow(
         m2sq,-1)*m3sq*ln3*pow(lm,2)*pow(lm1,3) + 5./6.*m1sq*pow(
         m2sq,-1)*m3sq*ln1*lm*pow(lm1,2) - m1sq*pow(m2sq,-1)*m3sq*ln1*
         pow(lm,2)*pow(lm1,3) + 1./3.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*ln3
         *pow(lm1,2) - 1./2.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*ln3*lm*pow(
         lm1,3) - 1./3.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*ln1*pow(lm1,2) + 
         1./2.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*ln1*lm*pow(lm1,3) - 1./2.*
         m1sq*lm1;
    h210 +=  - 1./2.*m1sq*ln3*lm*pow(lm1,2) + 1./3.*m1sq*ln3*pow(
         lm,2)*pow(lm1,3) - 1./3.*m1sq*ln2*lm1 + 4./3.*m1sq*ln2*lm*pow(
         lm1,2) - 5./3.*m1sq*ln2*pow(lm,2)*pow(lm1,3) - 1./3.*m1sq*ln1*
         lm1 - 5./6.*m1sq*ln1*lm*pow(lm1,2) + 4./3.*m1sq*ln1*pow(lm,2)*
         pow(lm1,3) - m1sq*pow(m3sq,2)*pow(lm1,2) + 2*m1sq*pow(m3sq,2)*
         ln3*pow(lm1,2) - 13./6.*m1sq*pow(m3sq,2)*ln3*lm*pow(lm1,3) - 
         m1sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./3.*m1sq*pow(m3sq,2)*ln2*
         lm*pow(lm1,3) - 7./3.*m1sq*pow(m3sq,2)*ln1*pow(lm1,2) + 11./6.
         *m1sq*pow(m3sq,2)*ln1*lm*pow(lm1,3) + 4./3.*m1sq*pow(m3sq,4)*
         ln3*pow(lm1,3) + 1./3.*m1sq*m2sq*m3sq*pow(lm1,2) + 4./3.*m1sq*
         m2sq*m3sq*ln3*pow(lm1,2) - 17./6.*m1sq*m2sq*m3sq*ln3*lm*pow(
         lm1,3) - 11./3.*m1sq*m2sq*m3sq*ln2*pow(lm1,2) + 10./3.*m1sq*
         m2sq*m3sq*ln2*lm*pow(lm1,3) + m1sq*m2sq*m3sq*ln1*pow(lm1,2) - 
         1./2.*m1sq*m2sq*m3sq*ln1*lm*pow(lm1,3) - 8./3.*m1sq*m2sq*pow(
         m3sq,3)*ln2*pow(lm1,3) + 2./3.*m1sq*pow(m2sq,2)*pow(lm1,2) - 2.
         /3.*m1sq*pow(m2sq,2)*ln3*pow(lm1,2);
    h210 +=  + 3./2.*m1sq*pow(m2sq,2)*ln3*lm*pow(lm1,3) + 5./3.*
         m1sq*pow(m2sq,2)*ln2*pow(lm1,2) + 1./3.*m1sq*pow(m2sq,2)*ln2*
         lm*pow(lm1,3) + 5./3.*m1sq*pow(m2sq,2)*ln1*pow(lm1,2) - 11./6.
         *m1sq*pow(m2sq,2)*ln1*lm*pow(lm1,3) - 4*m1sq*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,3) + 4*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*
         pow(lm1,3) + 8./3.*m1sq*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) - 4./3.
         *m1sq*pow(m2sq,4)*ln2*pow(lm1,3) + 1./6.*pow(m1sq,2)*pow(
         m2sq,-1)*lm*pow(lm1,2) + 1./3.*pow(m1sq,2)*pow(m2sq,-1)*ln3*lm
         *pow(lm1,2) - 1./2.*pow(m1sq,2)*pow(m2sq,-1)*ln3*pow(lm,2)*
         pow(lm1,3) - 1./3.*pow(m1sq,2)*pow(m2sq,-1)*ln1*lm*pow(lm1,2)
          + 1./2.*pow(m1sq,2)*pow(m2sq,-1)*ln1*pow(lm,2)*pow(lm1,3) - 
         pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(lm1,2) + 3./2.*
         pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*lm*pow(lm1,3) + pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm1,2) - 3./2.*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*pow(lm1,3) + 1./3.*
         pow(m1sq,2)*m3sq*pow(lm1,2);
    h210 +=  - 10./3.*pow(m1sq,2)*m3sq*ln3*pow(lm1,2) + 11./3.*pow(
         m1sq,2)*m3sq*ln3*lm*pow(lm1,3) + 2*pow(m1sq,2)*m3sq*ln2*pow(
         lm1,2) - 8./3.*pow(m1sq,2)*m3sq*ln2*lm*pow(lm1,3) + 1./3.*pow(
         m1sq,2)*m3sq*ln1*pow(lm1,2) - pow(m1sq,2)*m3sq*ln1*lm*pow(
         lm1,3) - 1./3.*pow(m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,3) - 7./3.*
         pow(m1sq,2)*pow(m3sq,3)*ln1*pow(lm1,3) - 4./3.*pow(m1sq,2)*
         m2sq*pow(lm1,2) - 1./2.*pow(m1sq,2)*m2sq*ln3*lm*pow(lm1,3) - 1.
         /3.*pow(m1sq,2)*m2sq*ln2*pow(lm1,2) - 8./3.*pow(m1sq,2)*m2sq*
         ln2*lm*pow(lm1,3) - 4*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) + 19./6.
         *pow(m1sq,2)*m2sq*ln1*lm*pow(lm1,3) + 22./3.*pow(m1sq,2)*m2sq*
         pow(m3sq,2)*ln3*pow(lm1,3) - 10./3.*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*ln2*pow(lm1,3) - 4*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*
         pow(lm1,3) - 3*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 
         14./3.*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) + 5*pow(
         m1sq,2)*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + 4*pow(m1sq,2)*pow(
         m2sq,3)*ln2*pow(lm1,3);
    h210 +=  + 4./3.*pow(m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,3) + pow(
         m1sq,3)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,2) - 3./2.*pow(m1sq,3)*
         pow(m2sq,-1)*m3sq*ln3*lm*pow(lm1,3) - pow(m1sq,3)*pow(m2sq,-1)
         *m3sq*ln1*pow(lm1,2) + 3./2.*pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln1
         *lm*pow(lm1,3) + 2./3.*pow(m1sq,3)*pow(lm1,2) + 2./3.*pow(
         m1sq,3)*ln3*pow(lm1,2) - 5./6.*pow(m1sq,3)*ln3*lm*pow(lm1,3)
          - 2./3.*pow(m1sq,3)*ln2*pow(lm1,2) + 5./3.*pow(m1sq,3)*ln2*lm
         *pow(lm1,3) + 2*pow(m1sq,3)*ln1*pow(lm1,2) - 5./6.*pow(m1sq,3)
         *ln1*lm*pow(lm1,3) - 10./3.*pow(m1sq,3)*pow(m3sq,2)*ln3*pow(
         lm1,3) + 10./3.*pow(m1sq,3)*pow(m3sq,2)*ln1*pow(lm1,3) - 2*
         pow(m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) + 14./3.*pow(m1sq,3)*m2sq
         *m3sq*ln2*pow(lm1,3) - 16./3.*pow(m1sq,3)*m2sq*m3sq*ln1*pow(
         lm1,3) - 4*pow(m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,3) - 4*pow(
         m1sq,3)*pow(m2sq,2)*ln1*pow(lm1,3) - 1./3.*pow(m1sq,4)*pow(
         m2sq,-1)*ln3*pow(lm1,2) + 1./2.*pow(m1sq,4)*pow(m2sq,-1)*ln3*
         lm*pow(lm1,3);
    h210 +=  + 1./3.*pow(m1sq,4)*pow(m2sq,-1)*ln1*pow(lm1,2) - 1./2.
         *pow(m1sq,4)*pow(m2sq,-1)*ln1*lm*pow(lm1,3) + 7./3.*pow(
         m1sq,4)*m3sq*ln3*pow(lm1,3) + 1./3.*pow(m1sq,4)*m3sq*ln1*pow(
         lm1,3) + 4./3.*pow(m1sq,4)*m2sq*ln2*pow(lm1,3) + 4*pow(m1sq,4)
         *m2sq*ln1*pow(lm1,3) - 4./3.*pow(m1sq,5)*ln1*pow(lm1,3);
    break;
  case 6:
    h210 =
       + psiq(m1,m2,m3) * ( 1./2.*lm*pow(lm1,2) - 1./2.*pow(lm,2)*pow(
         lm1,3) + pow(m3sq,2)*pow(lm1,2) - 7./2.*pow(m3sq,2)*lm*pow(
         lm1,3) + 5./2.*pow(m3sq,2)*pow(lm,2)*pow(lm1,4) - pow(m2sq,2)*
         pow(lm1,2) + 7./2.*pow(m2sq,2)*lm*pow(lm1,3) - 5./2.*pow(
         m2sq,2)*pow(lm,2)*pow(lm1,4) - m1sq*m3sq*pow(lm1,2) + 11./2.*
         m1sq*m3sq*lm*pow(lm1,3) - 5*m1sq*m3sq*pow(lm,2)*pow(lm1,4) - 2
         *m1sq*pow(m3sq,3)*pow(lm1,3) + 5./2.*m1sq*pow(m3sq,3)*lm*pow(
         lm1,4) + m1sq*m2sq*pow(lm1,2) - 1./2.*m1sq*m2sq*lm*pow(lm1,3)
          - 4*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) + 5./2.*m1sq*m2sq*pow(
         m3sq,2)*lm*pow(lm1,4) + 4*m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) - 5.
         /2.*m1sq*pow(m2sq,2)*m3sq*lm*pow(lm1,4) + 2*m1sq*pow(m2sq,3)*
         pow(lm1,3) - 5./2.*m1sq*pow(m2sq,3)*lm*pow(lm1,4) - 2*pow(
         m1sq,2)*lm*pow(lm1,3) + 5./2.*pow(m1sq,2)*pow(lm,2)*pow(lm1,4)
          + 6*pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3) - 15./2.*pow(m1sq,2)*
         pow(m3sq,2)*lm*pow(lm1,4) + 8*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3)
          - 5*pow(m1sq,2)*m2sq*m3sq*lm*pow(lm1,4) );

    h210 +=  + psiq(m1,m2,m3) * (  - 5*pow(m1sq,2)*m2sq*pow(m3sq,3)*
         pow(lm1,4) - pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) + 5./2.*pow(
         m1sq,2)*pow(m2sq,2)*lm*pow(lm1,4) + 5*pow(m1sq,2)*pow(m2sq,3)*
         m3sq*pow(lm1,4) - 6*pow(m1sq,3)*m3sq*pow(lm1,3) + 15./2.*pow(
         m1sq,3)*m3sq*lm*pow(lm1,4) - 3*pow(m1sq,3)*m2sq*pow(lm1,3) + 5.
         /2.*pow(m1sq,3)*m2sq*lm*pow(lm1,4) + 10*pow(m1sq,3)*m2sq*pow(
         m3sq,2)*pow(lm1,4) + 2*pow(m1sq,4)*pow(lm1,3) - 5./2.*pow(
         m1sq,4)*lm*pow(lm1,4) - 5*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) );

    h210 +=  + 1./6.*pow(m1sq,-1) - 1./6.*pow(m1sq,-1)*pow(lm,2)*
         pow(lm1,2) + 1./6.*pow(m1sq,-1)*ln3 - 1./6.*pow(m1sq,-1)*ln3*
         pow(lm,2)*pow(lm1,2) - 1./6.*pow(m1sq,-1)*ln2 + 1./6.*pow(
         m1sq,-1)*ln2*pow(lm,2)*pow(lm1,2) - 2./3.*pow(m1sq,-1)*pow(
         m3sq,2)*ln3*lm*pow(lm1,2) + 2./3.*pow(m1sq,-1)*pow(m3sq,2)*ln3
         *pow(lm,2)*pow(lm1,3) + 2./3.*pow(m1sq,-1)*pow(m3sq,2)*ln2*lm*
         pow(lm1,2) - 2./3.*pow(m1sq,-1)*pow(m3sq,2)*ln2*pow(lm,2)*pow(
         lm1,3) - 1./6.*pow(m1sq,-1)*m2sq*pow(m3sq,-1) + 1./6.*pow(
         m1sq,-1)*m2sq*pow(m3sq,-1)*pow(lm,2)*pow(lm1,2) + 4./3.*pow(
         m1sq,-1)*m2sq*m3sq*ln3*lm*pow(lm1,2) - 4./3.*pow(m1sq,-1)*m2sq
         *m3sq*ln3*pow(lm,2)*pow(lm1,3) - 4./3.*pow(m1sq,-1)*m2sq*m3sq*
         ln2*lm*pow(lm1,2) + 4./3.*pow(m1sq,-1)*m2sq*m3sq*ln2*pow(lm,2)
         *pow(lm1,3) - 2./3.*pow(m1sq,-1)*pow(m2sq,2)*ln3*lm*pow(lm1,2)
          + 2./3.*pow(m1sq,-1)*pow(m2sq,2)*ln3*pow(lm,2)*pow(lm1,3) + 2.
         /3.*pow(m1sq,-1)*pow(m2sq,2)*ln2*lm*pow(lm1,2) - 2./3.*pow(
         m1sq,-1)*pow(m2sq,2)*ln2*pow(lm,2)*pow(lm1,3);
    h210 +=  - 1./6.*pow(m3sq,-1)*pow(lm,2)*pow(lm1,2) + 1./6.*m3sq
         *lm1 - 1./6.*m3sq*lm*pow(lm1,2) + 1./3.*m3sq*ln3*lm1 - 4./3.*
         m3sq*ln3*lm*pow(lm1,2) + m3sq*ln3*pow(lm,2)*pow(lm1,3) - 1./2.
         *m3sq*ln2*lm*pow(lm1,2) + 5./6.*m3sq*ln2*pow(lm,2)*pow(lm1,3)
          + 11./6.*m3sq*ln1*lm*pow(lm1,2) - 11./6.*m3sq*ln1*pow(lm,2)*
         pow(lm1,3) - 2./3.*pow(m3sq,3)*ln3*pow(lm1,2) + 2./3.*pow(
         m3sq,3)*ln3*lm*pow(lm1,3) + 1./3.*pow(m3sq,3)*ln2*pow(lm1,2)
          - 2./3.*pow(m3sq,3)*ln2*lm*pow(lm1,3) - 1./6.*m2sq*lm1 - 1./6.
         *m2sq*ln3*lm1 - m2sq*ln3*lm*pow(lm1,2) + m2sq*ln3*pow(lm,2)*
         pow(lm1,3) - 1./6.*m2sq*ln2*lm1 + 7./3.*m2sq*ln2*lm*pow(lm1,2)
          - 7./3.*m2sq*ln2*pow(lm,2)*pow(lm1,3) - 4./3.*m2sq*ln1*lm*
         pow(lm1,2) + 4./3.*m2sq*ln1*pow(lm,2)*pow(lm1,3) + m2sq*pow(
         m3sq,2)*ln3*pow(lm1,2) - 2./3.*m2sq*pow(m3sq,2)*ln3*lm*pow(
         lm1,3) + 2./3.*m2sq*pow(m3sq,2)*ln2*lm*pow(lm1,3) + 1./6.*pow(
         m2sq,2)*pow(m3sq,-1)*lm*pow(lm1,2) + 1./2.*pow(m2sq,2)*pow(
         m3sq,-1)*ln2*lm*pow(lm1,2);
    h210 +=  - 1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(lm,2)*pow(
      lm1,3) - 1./2.*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(lm1,2) + 1./2.
         *pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(lm,2)*pow(lm1,3) - 2./3.*
         pow(m2sq,2)*m3sq*ln3*lm*pow(lm1,3) - pow(m2sq,2)*m3sq*ln2*pow(
         lm1,2) + 2./3.*pow(m2sq,2)*m3sq*ln2*lm*pow(lm1,3) - 1./3.*pow(
         m2sq,3)*ln3*pow(lm1,2) + 2./3.*pow(m2sq,3)*ln3*lm*pow(lm1,3)
          + 2./3.*pow(m2sq,3)*ln2*pow(lm1,2) - 2./3.*pow(m2sq,3)*ln2*lm
         *pow(lm1,3) - 1./2.*m1sq*lm1 - 1./3.*m1sq*ln3*lm1 + 4./3.*m1sq
         *ln3*lm*pow(lm1,2) - 5./3.*m1sq*ln3*pow(lm,2)*pow(lm1,3) - 1./
         2.*m1sq*ln2*lm*pow(lm1,2) + 1./3.*m1sq*ln2*pow(lm,2)*pow(
         lm1,3) - 1./3.*m1sq*ln1*lm1 - 5./6.*m1sq*ln1*lm*pow(lm1,2) + 4.
         /3.*m1sq*ln1*pow(lm,2)*pow(lm1,3) + 2./3.*m1sq*pow(m3sq,2)*
         pow(lm1,2) + 5./3.*m1sq*pow(m3sq,2)*ln3*pow(lm1,2) + 1./3.*
         m1sq*pow(m3sq,2)*ln3*lm*pow(lm1,3) - 2./3.*m1sq*pow(m3sq,2)*
         ln2*pow(lm1,2) + 3./2.*m1sq*pow(m3sq,2)*ln2*lm*pow(lm1,3) + 5./
         3.*m1sq*pow(m3sq,2)*ln1*pow(lm1,2);
    h210 +=  - 11./6.*m1sq*pow(m3sq,2)*ln1*lm*pow(lm1,3) - 4./3.*
         m1sq*pow(m3sq,4)*ln3*pow(lm1,3) - 1./3.*m1sq*m2sq*pow(m3sq,-1)
         *lm*pow(lm1,2) - 5./6.*m1sq*m2sq*pow(m3sq,-1)*ln2*lm*pow(
         lm1,2) + m1sq*m2sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3) + 5./
         6.*m1sq*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,2) - m1sq*m2sq*pow(
         m3sq,-1)*ln1*pow(lm,2)*pow(lm1,3) + 1./3.*m1sq*m2sq*m3sq*pow(
         lm1,2) - 11./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) + 10./3.*m1sq*
         m2sq*m3sq*ln3*lm*pow(lm1,3) + 4./3.*m1sq*m2sq*m3sq*ln2*pow(
         lm1,2) - 17./6.*m1sq*m2sq*m3sq*ln2*lm*pow(lm1,3) + m1sq*m2sq*
         m3sq*ln1*pow(lm1,2) - 1./2.*m1sq*m2sq*m3sq*ln1*lm*pow(lm1,3)
          + 8./3.*m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) - m1sq*pow(
         m2sq,2)*pow(lm1,2) - m1sq*pow(m2sq,2)*ln3*pow(lm1,2) + 1./3.*
         m1sq*pow(m2sq,2)*ln3*lm*pow(lm1,3) + 2*m1sq*pow(m2sq,2)*ln2*
         pow(lm1,2) - 13./6.*m1sq*pow(m2sq,2)*ln2*lm*pow(lm1,3) - 7./3.
         *m1sq*pow(m2sq,2)*ln1*pow(lm1,2) + 11./6.*m1sq*pow(m2sq,2)*ln1
         *lm*pow(lm1,3);
    h210 +=  + 4*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 4*
         m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 1./3.*m1sq*pow(
         m2sq,3)*pow(m3sq,-1)*ln2*pow(lm1,2) - 1./2.*m1sq*pow(m2sq,3)*
         pow(m3sq,-1)*ln2*lm*pow(lm1,3) - 1./3.*m1sq*pow(m2sq,3)*pow(
         m3sq,-1)*ln1*pow(lm1,2) + 1./2.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*
         ln1*lm*pow(lm1,3) - 8./3.*m1sq*pow(m2sq,3)*m3sq*ln3*pow(lm1,3)
          + 4./3.*m1sq*pow(m2sq,4)*ln2*pow(lm1,3) + 1./6.*pow(m1sq,2)*
         pow(m3sq,-1)*lm*pow(lm1,2) + 1./3.*pow(m1sq,2)*pow(m3sq,-1)*
         ln2*lm*pow(lm1,2) - 1./2.*pow(m1sq,2)*pow(m3sq,-1)*ln2*pow(
         lm,2)*pow(lm1,3) - 1./3.*pow(m1sq,2)*pow(m3sq,-1)*ln1*lm*pow(
         lm1,2) + 1./2.*pow(m1sq,2)*pow(m3sq,-1)*ln1*pow(lm,2)*pow(
         lm1,3) - 4./3.*pow(m1sq,2)*m3sq*pow(lm1,2) - 1./3.*pow(m1sq,2)
         *m3sq*ln3*pow(lm1,2) - 8./3.*pow(m1sq,2)*m3sq*ln3*lm*pow(
         lm1,3) - 1./2.*pow(m1sq,2)*m3sq*ln2*lm*pow(lm1,3) - 4*pow(
         m1sq,2)*m3sq*ln1*pow(lm1,2) + 19./6.*pow(m1sq,2)*m3sq*ln1*lm*
         pow(lm1,3);
    h210 +=  + 4*pow(m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,3) + 4./3.*
         pow(m1sq,2)*pow(m3sq,3)*ln1*pow(lm1,3) + 1./3.*pow(m1sq,2)*
         m2sq*pow(lm1,2) + 2*pow(m1sq,2)*m2sq*ln3*pow(lm1,2) - 8./3.*
         pow(m1sq,2)*m2sq*ln3*lm*pow(lm1,3) - 10./3.*pow(m1sq,2)*m2sq*
         ln2*pow(lm1,2) + 11./3.*pow(m1sq,2)*m2sq*ln2*lm*pow(lm1,3) + 1.
         /3.*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) - pow(m1sq,2)*m2sq*ln1*lm*
         pow(lm1,3) - 14./3.*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln3*pow(
         lm1,3) - 3*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) + 5*
         pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) - pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*pow(lm1,2) + 3./2.*pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*lm*pow(lm1,3) + pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,-1)*ln1*pow(lm1,2) - 3./2.*pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(lm1,3) - 10./3.*pow(m1sq,2)*
         pow(m2sq,2)*m3sq*ln3*pow(lm1,3) + 22./3.*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*ln2*pow(lm1,3) - 4*pow(m1sq,2)*pow(m2sq,2)*m3sq*
         ln1*pow(lm1,3);
    h210 +=  - 1./3.*pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,3) - 7./3.
         *pow(m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,3) + 2./3.*pow(m1sq,3)*
         pow(lm1,2) - 2./3.*pow(m1sq,3)*ln3*pow(lm1,2) + 5./3.*pow(
         m1sq,3)*ln3*lm*pow(lm1,3) + 2./3.*pow(m1sq,3)*ln2*pow(lm1,2)
          - 5./6.*pow(m1sq,3)*ln2*lm*pow(lm1,3) + 2*pow(m1sq,3)*ln1*
         pow(lm1,2) - 5./6.*pow(m1sq,3)*ln1*lm*pow(lm1,3) - 4*pow(
         m1sq,3)*pow(m3sq,2)*ln3*pow(lm1,3) - 4*pow(m1sq,3)*pow(m3sq,2)
         *ln1*pow(lm1,3) + pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln2*pow(lm1,2)
          - 3./2.*pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln2*lm*pow(lm1,3) - 
         pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln1*pow(lm1,2) + 3./2.*pow(
         m1sq,3)*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,3) + 14./3.*pow(
         m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) - 2*pow(m1sq,3)*m2sq*m3sq*ln2
         *pow(lm1,3) - 16./3.*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3) - 10.
         /3.*pow(m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,3) + 10./3.*pow(
         m1sq,3)*pow(m2sq,2)*ln1*pow(lm1,3) - 1./3.*pow(m1sq,4)*pow(
         m3sq,-1)*ln2*pow(lm1,2);
    h210 +=  + 1./2.*pow(m1sq,4)*pow(m3sq,-1)*ln2*lm*pow(lm1,3) + 1.
         /3.*pow(m1sq,4)*pow(m3sq,-1)*ln1*pow(lm1,2) - 1./2.*pow(
         m1sq,4)*pow(m3sq,-1)*ln1*lm*pow(lm1,3) + 4./3.*pow(m1sq,4)*
         m3sq*ln3*pow(lm1,3) + 4*pow(m1sq,4)*m3sq*ln1*pow(lm1,3) + 7./3.
         *pow(m1sq,4)*m2sq*ln2*pow(lm1,3) + 1./3.*pow(m1sq,4)*m2sq*ln1*
         pow(lm1,3) - 4./3.*pow(m1sq,5)*ln1*pow(lm1,3);
    break;
  case 7:
    h210 =
       + psiq(m1,m2,m3) * ( 2./3.*lm*pow(lm1,2) - 1./2.*pow(lm,2)*pow(
         lm1,3) + 4./3.*pow(m3sq,2)*pow(lm1,2) - 4*pow(m3sq,2)*lm*pow(
         lm1,3) + 5./2.*pow(m3sq,2)*pow(lm,2)*pow(lm1,4) - 8./3.*m2sq*
         m3sq*pow(lm1,2) + 8*m2sq*m3sq*lm*pow(lm1,3) - 5*m2sq*m3sq*pow(
         lm,2)*pow(lm1,4) + 4./3.*pow(m2sq,2)*pow(lm1,2) - 4*pow(
         m2sq,2)*lm*pow(lm1,3) + 5./2.*pow(m2sq,2)*pow(lm,2)*pow(lm1,4)
          + 1./3.*m1sq*m3sq*pow(lm1,2) - 1./2.*m1sq*m3sq*lm*pow(lm1,3)
          - 2*m1sq*pow(m3sq,3)*pow(lm1,3) + 5./2.*m1sq*pow(m3sq,3)*lm*
         pow(lm1,4) + 1./3.*m1sq*m2sq*pow(lm1,2) - 1./2.*m1sq*m2sq*lm*
         pow(lm1,3) + 2*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 5./2.*m1sq*
         m2sq*pow(m3sq,2)*lm*pow(lm1,4) + 2*m1sq*pow(m2sq,2)*m3sq*pow(
         lm1,3) - 5./2.*m1sq*pow(m2sq,2)*m3sq*lm*pow(lm1,4) - 2*m1sq*
         pow(m2sq,3)*pow(lm1,3) + 5./2.*m1sq*pow(m2sq,3)*lm*pow(lm1,4)
          - 2./3.*pow(m1sq,2)*pow(lm1,2) + 7./2.*pow(m1sq,2)*lm*pow(
         lm1,3) - 5./2.*pow(m1sq,2)*pow(lm,2)*pow(lm1,4) + pow(m1sq,2)*
         pow(m3sq,2)*pow(lm1,3) );

    h210 +=  + psiq(m1,m2,m3) * (  - 5./2.*pow(m1sq,2)*pow(m3sq,2)*
         lm*pow(lm1,4) - pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) + 5*pow(
         m1sq,2)*m2sq*m3sq*lm*pow(lm1,4) - 5*pow(m1sq,2)*m2sq*pow(
         m3sq,3)*pow(lm1,4) + pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 5./2.
         *pow(m1sq,2)*pow(m2sq,2)*lm*pow(lm1,4) + 10*pow(m1sq,2)*pow(
         m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 5*pow(m1sq,2)*pow(m2sq,3)*
         m3sq*pow(lm1,4) + 3*pow(m1sq,3)*m3sq*pow(lm1,3) - 5./2.*pow(
         m1sq,3)*m3sq*lm*pow(lm1,4) + 3*pow(m1sq,3)*m2sq*pow(lm1,3) - 5.
         /2.*pow(m1sq,3)*m2sq*lm*pow(lm1,4) - 2*pow(m1sq,4)*pow(lm1,3)
          + 5./2.*pow(m1sq,4)*lm*pow(lm1,4) + 5*pow(m1sq,4)*m2sq*m3sq*
         pow(lm1,4) );

    h210 +=  + 1./6.*pow(m2sq,-1) - 1./6.*pow(m2sq,-1)*pow(lm,2)*
         pow(lm1,2) + 1./6.*pow(m2sq,-1)*ln3 - 1./6.*pow(m2sq,-1)*ln3*
         pow(lm,2)*pow(lm1,2) - 1./6.*pow(m2sq,-1)*ln1 + 1./6.*pow(
         m2sq,-1)*ln1*pow(lm,2)*pow(lm1,2) - 2./3.*pow(m2sq,-1)*pow(
         m3sq,2)*ln3*lm*pow(lm1,2) + 2./3.*pow(m2sq,-1)*pow(m3sq,2)*ln3
         *pow(lm,2)*pow(lm1,3) + 2./3.*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*
         pow(lm1,2) - 2./3.*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm,2)*pow(
         lm1,3) + 1./6.*pow(m3sq,-1) - 1./6.*pow(m3sq,-1)*pow(lm,2)*
         pow(lm1,2) + 1./6.*pow(m3sq,-1)*ln2 - 1./6.*pow(m3sq,-1)*ln1
          - 4./3.*m3sq*ln3*lm*pow(lm1,2) + m3sq*ln3*pow(lm,2)*pow(
         lm1,3) + 2*m3sq*ln2*lm*pow(lm1,2) - 11./6.*m3sq*ln2*pow(lm,2)*
         pow(lm1,3) - 2./3.*m3sq*ln1*lm*pow(lm1,2) + 5./6.*m3sq*ln1*
         pow(lm,2)*pow(lm1,3) + 2*m2sq*ln3*lm*pow(lm1,2) - 5./3.*m2sq*
         ln3*pow(lm,2)*pow(lm1,3) - 4./3.*m2sq*ln2*lm*pow(lm1,2) + 4./3.
         *m2sq*ln2*pow(lm,2)*pow(lm1,3) - 2./3.*m2sq*ln1*lm*pow(lm1,2)
          + 1./3.*m2sq*ln1*pow(lm,2)*pow(lm1,3);
    h210 +=  - 2./3.*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm*pow(lm1,2) + 1.
         /2.*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3) + 2./3.*
         pow(m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(lm1,2) - 1./2.*pow(m2sq,2)
         *pow(m3sq,-1)*ln1*pow(lm,2)*pow(lm1,3) + 1./6.*m1sq*pow(
         m2sq,-1)*pow(m3sq,-1)*pow(lm,2)*pow(lm1,2) - 1./6.*m1sq*pow(
         m2sq,-1)*m3sq*lm*pow(lm1,2) + m1sq*pow(m2sq,-1)*m3sq*ln3*lm*
         pow(lm1,2) - 4./3.*m1sq*pow(m2sq,-1)*m3sq*ln3*pow(lm,2)*pow(
         lm1,3) - m1sq*pow(m2sq,-1)*m3sq*ln1*lm*pow(lm1,2) + 4./3.*m1sq
         *pow(m2sq,-1)*m3sq*ln1*pow(lm,2)*pow(lm1,3) - 1./3.*m1sq*pow(
         m2sq,-1)*pow(m3sq,3)*ln3*pow(lm1,2) + 2./3.*m1sq*pow(m2sq,-1)*
         pow(m3sq,3)*ln3*lm*pow(lm1,3) + 1./3.*m1sq*pow(m2sq,-1)*pow(
         m3sq,3)*ln1*pow(lm1,2) - 2./3.*m1sq*pow(m2sq,-1)*pow(m3sq,3)*
         ln1*lm*pow(lm1,3) - 1./3.*m1sq*lm1 - 1./3.*m1sq*lm*pow(lm1,2)
          - 1./6.*m1sq*ln3*lm1 - 4./3.*m1sq*ln3*lm*pow(lm1,2) + m1sq*
         ln3*pow(lm,2)*pow(lm1,3) - 1./6.*m1sq*ln2*lm1 - 7./6.*m1sq*ln2
         *lm*pow(lm1,2);
      h210 +=  + 4./3.*m1sq*ln2*pow(lm,2)*pow(lm1,3) + 5./2.*m1sq*ln1
         *lm*pow(lm1,2) - 7./3.*m1sq*ln1*pow(lm,2)*pow(lm1,3) + 2./3.*
         m1sq*pow(m3sq,2)*pow(lm1,2) + 1./3.*m1sq*pow(m3sq,2)*ln3*pow(
         lm1,2) + 5./3.*m1sq*pow(m3sq,2)*ln3*lm*pow(lm1,3) + 4./3.*m1sq
         *pow(m3sq,2)*ln2*pow(lm1,2) - 11./6.*m1sq*pow(m3sq,2)*ln2*lm*
         pow(lm1,3) + 1./6.*m1sq*pow(m3sq,2)*ln1*lm*pow(lm1,3) - 4./3.*
         m1sq*pow(m3sq,4)*ln3*pow(lm1,3) - 1./6.*m1sq*m2sq*pow(m3sq,-1)
         *lm*pow(lm1,2) + 7./6.*m1sq*m2sq*pow(m3sq,-1)*ln2*lm*pow(
         lm1,2) - m1sq*m2sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3) - 7./
         6.*m1sq*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,2) + m1sq*m2sq*pow(
         m3sq,-1)*ln1*pow(lm,2)*pow(lm1,3) - 4./3.*m1sq*m2sq*m3sq*pow(
         lm1,2) - 4./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 2./3.*m1sq*m2sq
         *m3sq*ln3*lm*pow(lm1,3) - 4./3.*m1sq*m2sq*m3sq*ln2*pow(lm1,2)
          - 1./2.*m1sq*m2sq*m3sq*ln2*lm*pow(lm1,3) - 2./3.*m1sq*m2sq*
         m3sq*ln1*pow(lm1,2) + 7./6.*m1sq*m2sq*m3sq*ln1*lm*pow(lm1,3)
          + 4*m1sq*m2sq*pow(m3sq,3)*ln3*pow(lm1,3);
      h210 +=  + 4./3.*m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 2./3.*
         m1sq*pow(m2sq,2)*pow(lm1,2) + 4./3.*m1sq*pow(m2sq,2)*ln3*pow(
         lm1,2) - 5./3.*m1sq*pow(m2sq,2)*ln3*lm*pow(lm1,3) + 1./3.*m1sq
         *pow(m2sq,2)*ln2*pow(lm1,2) + 11./6.*m1sq*pow(m2sq,2)*ln2*lm*
         pow(lm1,3) - 1./6.*m1sq*pow(m2sq,2)*ln1*lm*pow(lm1,3) - 4*m1sq
         *pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 4*m1sq*pow(m2sq,2)*
         pow(m3sq,2)*ln2*pow(lm1,3) - 1./3.*m1sq*pow(m2sq,3)*pow(
         m3sq,-1)*ln2*pow(lm1,2) + 1./2.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*
         ln2*lm*pow(lm1,3) + 1./3.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln1*
         pow(lm1,2) - 1./2.*m1sq*pow(m2sq,3)*pow(m3sq,-1)*ln1*lm*pow(
         lm1,3) + 4./3.*m1sq*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) + 4*m1sq*
         pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 4./3.*m1sq*pow(m2sq,4)*ln2*
         pow(lm1,3) + 1./3.*pow(m1sq,2)*pow(m2sq,-1)*lm*pow(lm1,2) - 1./
         3.*pow(m1sq,2)*pow(m2sq,-1)*ln3*lm*pow(lm1,2) + 2./3.*pow(
         m1sq,2)*pow(m2sq,-1)*ln3*pow(lm,2)*pow(lm1,3) + 1./3.*pow(
         m1sq,2)*pow(m2sq,-1)*ln1*lm*pow(lm1,2);
      h210 +=  - 2./3.*pow(m1sq,2)*pow(m2sq,-1)*ln1*pow(lm,2)*pow(
      lm1,3) + pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(lm1,2) - 2*
         pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln3*lm*pow(lm1,3) - pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm1,2) + 2*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*pow(lm1,3) + 1./3.*
         pow(m1sq,2)*pow(m3sq,-1)*lm*pow(lm1,2) - 1./2.*pow(m1sq,2)*
         pow(m3sq,-1)*ln2*lm*pow(lm1,2) + 1./2.*pow(m1sq,2)*pow(
         m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3) + 1./2.*pow(m1sq,2)*pow(
         m3sq,-1)*ln1*lm*pow(lm1,2) - 1./2.*pow(m1sq,2)*pow(m3sq,-1)*
         ln1*pow(lm,2)*pow(lm1,3) + 1./3.*pow(m1sq,2)*m3sq*pow(lm1,2)
          + 1./3.*pow(m1sq,2)*m3sq*ln3*pow(lm1,2) - 4./3.*pow(m1sq,2)*
         m3sq*ln3*lm*pow(lm1,3) - 4./3.*pow(m1sq,2)*m3sq*ln2*pow(lm1,2)
          + 19./6.*pow(m1sq,2)*m3sq*ln2*lm*pow(lm1,3) + 1./3.*pow(
         m1sq,2)*m3sq*ln1*pow(lm1,2) - 11./6.*pow(m1sq,2)*m3sq*ln1*lm*
         pow(lm1,3) + 8./3.*pow(m1sq,2)*pow(m3sq,3)*ln1*pow(lm1,3) + 1./
         3.*pow(m1sq,2)*m2sq*pow(lm1,2);
      h210 +=  - 4./3.*pow(m1sq,2)*m2sq*ln3*pow(lm1,2) + 8./3.*pow(
         m1sq,2)*m2sq*ln3*lm*pow(lm1,3) - pow(m1sq,2)*m2sq*ln2*lm*pow(
         lm1,3) + 2./3.*pow(m1sq,2)*m2sq*ln1*pow(lm1,2) - 5./3.*pow(
         m1sq,2)*m2sq*ln1*lm*pow(lm1,3) - 14./3.*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,3) + 5*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,3) - 3*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(lm1,2) - 3./2.*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm*pow(lm1,3) - pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(lm1,2) + 3./2.*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(lm1,3) + 14./3.*
         pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 16./3.*pow(
         m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 2*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*ln1*pow(lm1,3) + 1./3.*pow(m1sq,2)*pow(m2sq,3)*
         ln2*pow(lm1,3) + 7./3.*pow(m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,3)
          - 1./6.*pow(m1sq,3)*pow(m2sq,-1)*pow(m3sq,-1)*lm*pow(lm1,2)
          - pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,2);
      h210 +=  + 2*pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln3*lm*pow(lm1,3) + 
         pow(m1sq,3)*pow(m2sq,-1)*m3sq*ln1*pow(lm1,2) - 2*pow(m1sq,3)*
         pow(m2sq,-1)*m3sq*ln1*lm*pow(lm1,3) - pow(m1sq,3)*pow(lm1,2)
          - 1./3.*pow(m1sq,3)*ln3*pow(lm1,2) - 1./3.*pow(m1sq,3)*ln3*lm
         *pow(lm1,3) - 1./3.*pow(m1sq,3)*ln2*pow(lm1,2) - 5./6.*pow(
         m1sq,3)*ln2*lm*pow(lm1,3) - 1./3.*pow(m1sq,3)*ln1*pow(lm1,2)
          + 7./6.*pow(m1sq,3)*ln1*lm*pow(lm1,3) + 4*pow(m1sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,3) - 4*pow(m1sq,3)*pow(m3sq,2)*ln1*pow(
         lm1,3) - pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln2*pow(lm1,2) + 3./2.*
         pow(m1sq,3)*m2sq*pow(m3sq,-1)*ln2*lm*pow(lm1,3) + pow(m1sq,3)*
         m2sq*pow(m3sq,-1)*ln1*pow(lm1,2) - 3./2.*pow(m1sq,3)*m2sq*pow(
         m3sq,-1)*ln1*lm*pow(lm1,3) - 10./3.*pow(m1sq,3)*m2sq*m3sq*ln3*
         pow(lm1,3) - 4*pow(m1sq,3)*m2sq*m3sq*ln2*pow(lm1,3) + 22./3.*
         pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3) + 10./3.*pow(m1sq,3)*pow(
         m2sq,2)*ln2*pow(lm1,3) - 10./3.*pow(m1sq,3)*pow(m2sq,2)*ln1*
         pow(lm1,3);
      h210 +=  + 1./3.*pow(m1sq,4)*pow(m2sq,-1)*ln3*pow(lm1,2) - 2./3.
         *pow(m1sq,4)*pow(m2sq,-1)*ln3*lm*pow(lm1,3) - 1./3.*pow(
         m1sq,4)*pow(m2sq,-1)*ln1*pow(lm1,2) + 2./3.*pow(m1sq,4)*pow(
         m2sq,-1)*ln1*lm*pow(lm1,3) + 1./3.*pow(m1sq,4)*pow(m3sq,-1)*
         ln2*pow(lm1,2) - 1./2.*pow(m1sq,4)*pow(m3sq,-1)*ln2*lm*pow(
         lm1,3) - 1./3.*pow(m1sq,4)*pow(m3sq,-1)*ln1*pow(lm1,2) + 1./2.
         *pow(m1sq,4)*pow(m3sq,-1)*ln1*lm*pow(lm1,3) - 8./3.*pow(
         m1sq,4)*m3sq*ln3*pow(lm1,3) - 7./3.*pow(m1sq,4)*m2sq*ln2*pow(
         lm1,3) - 1./3.*pow(m1sq,4)*m2sq*ln1*pow(lm1,3) + 4./3.*pow(
         m1sq,5)*ln1*pow(lm1,3);
      break;
  case 8:
    h210 =
       + psiq(m1,m2,m3) * (  - m3sq*pow(lm1,2) + 7./2.*m3sq*lm*pow(
         lm1,3) - 5./2.*m3sq*pow(lm,2)*pow(lm1,4) + 10*pow(m3sq,3)*pow(
         lm1,3) - 55./2.*pow(m3sq,3)*lm*pow(lm1,4) + 35./2.*pow(m3sq,3)
         *pow(lm,2)*pow(lm1,5) - m2sq*pow(lm1,2) + 7./2.*m2sq*lm*pow(
         lm1,3) - 5./2.*m2sq*pow(lm,2)*pow(lm1,4) - 10*m2sq*pow(m3sq,2)
         *pow(lm1,3) + 55./2.*m2sq*pow(m3sq,2)*lm*pow(lm1,4) - 35./2.*
         m2sq*pow(m3sq,2)*pow(lm,2)*pow(lm1,5) - 10*pow(m2sq,2)*m3sq*
         pow(lm1,3) + 55./2.*pow(m2sq,2)*m3sq*lm*pow(lm1,4) - 35./2.*
         pow(m2sq,2)*m3sq*pow(lm,2)*pow(lm1,5) + 10*pow(m2sq,3)*pow(
         lm1,3) - 55./2.*pow(m2sq,3)*lm*pow(lm1,4) + 35./2.*pow(m2sq,3)
         *pow(lm,2)*pow(lm1,5) + 3*m1sq*lm*pow(lm1,3) - 5./2.*m1sq*pow(
         lm,2)*pow(lm1,4) - 8*m1sq*pow(m3sq,2)*pow(lm1,3) + 45./2.*m1sq
         *pow(m3sq,2)*lm*pow(lm1,4) - 35./2.*m1sq*pow(m3sq,2)*pow(lm,2)
         *pow(lm1,5) - 15*m1sq*pow(m3sq,4)*pow(lm1,4) + 35./2.*m1sq*
         pow(m3sq,4)*lm*pow(lm1,5) + 26*m1sq*m2sq*m3sq*pow(lm1,3) - 55*
         m1sq*m2sq*m3sq*lm*pow(lm1,4) );

      h210 +=  + psiq(m1,m2,m3) * ( 35*m1sq*m2sq*m3sq*pow(lm,2)*pow(
         lm1,5) - 10*m1sq*m2sq*pow(m3sq,3)*pow(lm1,4) - 8*m1sq*pow(
         m2sq,2)*pow(lm1,3) + 45./2.*m1sq*pow(m2sq,2)*lm*pow(lm1,4) - 
         35./2.*m1sq*pow(m2sq,2)*pow(lm,2)*pow(lm1,5) + 50*m1sq*pow(
         m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 35*m1sq*pow(m2sq,2)*pow(
         m3sq,2)*lm*pow(lm1,5) - 10*m1sq*pow(m2sq,3)*m3sq*pow(lm1,4) - 
         15*m1sq*pow(m2sq,4)*pow(lm1,4) + 35./2.*m1sq*pow(m2sq,4)*lm*
         pow(lm1,5) - 2*pow(m1sq,2)*m3sq*pow(lm1,3) + 45./2.*pow(
         m1sq,2)*m3sq*lm*pow(lm1,4) - 35./2.*pow(m1sq,2)*m3sq*pow(lm,2)
         *pow(lm1,5) + 25*pow(m1sq,2)*pow(m3sq,3)*pow(lm1,4) - 35*pow(
         m1sq,2)*pow(m3sq,3)*lm*pow(lm1,5) - 2*pow(m1sq,2)*m2sq*pow(
         lm1,3) + 45./2.*pow(m1sq,2)*m2sq*lm*pow(lm1,4) - 35./2.*pow(
         m1sq,2)*m2sq*pow(lm,2)*pow(lm1,5) - 20*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*pow(lm1,4) + 35*pow(m1sq,2)*m2sq*pow(m3sq,2)*lm*pow(
         lm1,5) - 35*pow(m1sq,2)*m2sq*pow(m3sq,4)*pow(lm1,5) - 20*pow(
         m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,4) );

      h210 +=  + psiq(m1,m2,m3) * ( 35*pow(m1sq,2)*pow(m2sq,2)*m3sq*lm
         *pow(lm1,5) + 35*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*pow(
         lm1,5) + 25*pow(m1sq,2)*pow(m2sq,3)*pow(lm1,4) - 35*pow(
         m1sq,2)*pow(m2sq,3)*lm*pow(lm1,5) + 35*pow(m1sq,2)*pow(m2sq,3)
         *pow(m3sq,2)*pow(lm1,5) - 35*pow(m1sq,2)*pow(m2sq,4)*m3sq*pow(
         lm1,5) + pow(m1sq,3)*pow(lm1,3) - 35./2.*pow(m1sq,3)*lm*pow(
         lm1,4) + 35./2.*pow(m1sq,3)*pow(lm,2)*pow(lm1,5) + 10*pow(
         m1sq,3)*pow(m3sq,2)*pow(lm1,4) + 75*pow(m1sq,3)*m2sq*m3sq*pow(
         lm1,4) - 70*pow(m1sq,3)*m2sq*m3sq*lm*pow(lm1,5) + 35*pow(
         m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5) + 10*pow(m1sq,3)*pow(
         m2sq,2)*pow(lm1,4) - 70*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*
         pow(lm1,5) + 35*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(lm1,5) - 35*
         pow(m1sq,4)*m3sq*pow(lm1,4) + 35*pow(m1sq,4)*m3sq*lm*pow(
         lm1,5) - 35*pow(m1sq,4)*m2sq*pow(lm1,4) + 35*pow(m1sq,4)*m2sq*
         lm*pow(lm1,5) + 35*pow(m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5) + 
         35*pow(m1sq,4)*pow(m2sq,2)*m3sq*pow(lm1,5) );

      h210 +=  + psiq(m1,m2,m3) * ( 15*pow(m1sq,5)*pow(lm1,4) - 35./2.
         *pow(m1sq,5)*lm*pow(lm1,5) - 35*pow(m1sq,5)*m2sq*m3sq*pow(
         lm1,5) );

      h210 +=  - 1./6.*pow(m1sq,-1)*pow(m2sq,-1) + 1./6.*pow(m1sq,-1)
         *pow(m2sq,-1)*pow(lm,2)*pow(lm1,2) + 2./3.*pow(m1sq,-1)*pow(
         m2sq,-1)*pow(m3sq,2)*lm*pow(lm1,2) - 2./3.*pow(m1sq,-1)*pow(
         m2sq,-1)*pow(m3sq,2)*pow(lm,2)*pow(lm1,3) - 1./6.*pow(m1sq,-1)
         *pow(m3sq,-1) + 1./6.*pow(m1sq,-1)*pow(m3sq,-1)*pow(lm,2)*pow(
         lm1,2) - 2./3.*pow(m1sq,-1)*m3sq*lm*pow(lm1,2) + 2./3.*pow(
         m1sq,-1)*m3sq*pow(lm,2)*pow(lm1,3) + 2*pow(m1sq,-1)*m3sq*ln3*
         lm*pow(lm1,2) - 2*pow(m1sq,-1)*m3sq*ln3*pow(lm,2)*pow(lm1,3)
          - 2*pow(m1sq,-1)*m3sq*ln2*lm*pow(lm1,2) + 2*pow(m1sq,-1)*m3sq
         *ln2*pow(lm,2)*pow(lm1,3) + 4./3.*pow(m1sq,-1)*pow(m3sq,3)*ln3
         *pow(lm1,2) - 16./3.*pow(m1sq,-1)*pow(m3sq,3)*ln3*lm*pow(
         lm1,3) + 4*pow(m1sq,-1)*pow(m3sq,3)*ln3*pow(lm,2)*pow(lm1,4)
          - 4./3.*pow(m1sq,-1)*pow(m3sq,3)*ln2*pow(lm1,2) + 16./3.*pow(
         m1sq,-1)*pow(m3sq,3)*ln2*lm*pow(lm1,3) - 4*pow(m1sq,-1)*pow(
         m3sq,3)*ln2*pow(lm,2)*pow(lm1,4) - 2./3.*pow(m1sq,-1)*m2sq*lm*
         pow(lm1,2);
      h210 +=  + 2./3.*pow(m1sq,-1)*m2sq*pow(lm,2)*pow(lm1,3) - 2*
         pow(m1sq,-1)*m2sq*ln3*lm*pow(lm1,2) + 2*pow(m1sq,-1)*m2sq*ln3*
         pow(lm,2)*pow(lm1,3) + 2*pow(m1sq,-1)*m2sq*ln2*lm*pow(lm1,2)
          - 2*pow(m1sq,-1)*m2sq*ln2*pow(lm,2)*pow(lm1,3) - 4*pow(
         m1sq,-1)*m2sq*pow(m3sq,2)*ln3*pow(lm1,2) + 16*pow(m1sq,-1)*
         m2sq*pow(m3sq,2)*ln3*lm*pow(lm1,3) - 12*pow(m1sq,-1)*m2sq*pow(
         m3sq,2)*ln3*pow(lm,2)*pow(lm1,4) + 4*pow(m1sq,-1)*m2sq*pow(
         m3sq,2)*ln2*pow(lm1,2) - 16*pow(m1sq,-1)*m2sq*pow(m3sq,2)*ln2*
         lm*pow(lm1,3) + 12*pow(m1sq,-1)*m2sq*pow(m3sq,2)*ln2*pow(lm,2)
         *pow(lm1,4) + 2./3.*pow(m1sq,-1)*pow(m2sq,2)*pow(m3sq,-1)*lm*
         pow(lm1,2) - 2./3.*pow(m1sq,-1)*pow(m2sq,2)*pow(m3sq,-1)*pow(
         lm,2)*pow(lm1,3) + 4*pow(m1sq,-1)*pow(m2sq,2)*m3sq*ln3*pow(
         lm1,2) - 16*pow(m1sq,-1)*pow(m2sq,2)*m3sq*ln3*lm*pow(lm1,3) + 
         12*pow(m1sq,-1)*pow(m2sq,2)*m3sq*ln3*pow(lm,2)*pow(lm1,4) - 4*
         pow(m1sq,-1)*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + 16*pow(m1sq,-1)
         *pow(m2sq,2)*m3sq*ln2*lm*pow(lm1,3);
      h210 +=  - 12*pow(m1sq,-1)*pow(m2sq,2)*m3sq*ln2*pow(lm,2)*pow(
      lm1,4) - 4./3.*pow(m1sq,-1)*pow(m2sq,3)*ln3*pow(lm1,2) + 16./3.*
         pow(m1sq,-1)*pow(m2sq,3)*ln3*lm*pow(lm1,3) - 4*pow(m1sq,-1)*
         pow(m2sq,3)*ln3*pow(lm,2)*pow(lm1,4) + 4./3.*pow(m1sq,-1)*pow(
         m2sq,3)*ln2*pow(lm1,2) - 16./3.*pow(m1sq,-1)*pow(m2sq,3)*ln2*
         lm*pow(lm1,3) + 4*pow(m1sq,-1)*pow(m2sq,3)*ln2*pow(lm,2)*pow(
         lm1,4) - 1./2.*pow(m2sq,-1)*m3sq*lm*pow(lm1,2) + 5./6.*pow(
         m2sq,-1)*m3sq*pow(lm,2)*pow(lm1,3) + pow(m2sq,-1)*m3sq*ln3*lm*
         pow(lm1,2) - pow(m2sq,-1)*m3sq*ln3*pow(lm,2)*pow(lm1,3) - pow(
         m2sq,-1)*m3sq*ln1*lm*pow(lm1,2) + pow(m2sq,-1)*m3sq*ln1*pow(
         lm,2)*pow(lm1,3) + 1./3.*pow(m2sq,-1)*pow(m3sq,3)*pow(lm1,2)
          - 2./3.*pow(m2sq,-1)*pow(m3sq,3)*lm*pow(lm1,3) + pow(m2sq,-1)
         *pow(m3sq,3)*ln3*pow(lm1,2) - 4*pow(m2sq,-1)*pow(m3sq,3)*ln3*
         lm*pow(lm1,3) + 3*pow(m2sq,-1)*pow(m3sq,3)*ln3*pow(lm,2)*pow(
         lm1,4) - pow(m2sq,-1)*pow(m3sq,3)*ln1*pow(lm1,2) + 4*pow(
         m2sq,-1)*pow(m3sq,3)*ln1*lm*pow(lm1,3);
      h210 +=  - 3*pow(m2sq,-1)*pow(m3sq,3)*ln1*pow(lm,2)*pow(lm1,4)
          - 1./3.*lm1 + 11./3.*lm*pow(lm1,2) - 11./3.*pow(lm,2)*pow(
         lm1,3) - 1./6.*ln3*lm1 + 2./3.*ln3*lm*pow(lm1,2) - 2./3.*ln3*
         pow(lm,2)*pow(lm1,3) - 1./6.*ln2*lm1 + 1./6.*ln2*lm*pow(lm1,2)
          - 1./6.*ln2*pow(lm,2)*pow(lm1,3) - 5./6.*ln1*lm*pow(lm1,2) + 
         5./6.*ln1*pow(lm,2)*pow(lm1,3) + 2./3.*pow(m3sq,2)*pow(lm1,2)
          + 20./3.*pow(m3sq,2)*ln3*pow(lm1,2) - 55./3.*pow(m3sq,2)*ln3*
         lm*pow(lm1,3) + 12*pow(m3sq,2)*ln3*pow(lm,2)*pow(lm1,4) - 7./3.
         *pow(m3sq,2)*ln2*pow(lm1,2) + 49./6.*pow(m3sq,2)*ln2*lm*pow(
         lm1,3) - 9./2.*pow(m3sq,2)*ln2*pow(lm,2)*pow(lm1,4) - 8./3.*
         pow(m3sq,2)*ln1*pow(lm1,2) + 61./6.*pow(m3sq,2)*ln1*lm*pow(
         lm1,3) - 15./2.*pow(m3sq,2)*ln1*pow(lm,2)*pow(lm1,4) - 4*pow(
         m3sq,4)*ln3*pow(lm1,3) + 4*pow(m3sq,4)*ln3*lm*pow(lm1,4) + 8./
         3.*pow(m3sq,4)*ln2*pow(lm1,3) - 4*pow(m3sq,4)*ln2*lm*pow(
         lm1,4) - 1./2.*m2sq*pow(m3sq,-1)*lm*pow(lm1,2) + 5./6.*m2sq*
         pow(m3sq,-1)*pow(lm,2)*pow(lm1,3);
      h210 +=  + 1./2.*m2sq*pow(m3sq,-1)*ln2*lm*pow(lm1,2) - 1./2.*
         m2sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3) - 1./2.*m2sq*pow(
         m3sq,-1)*ln1*lm*pow(lm1,2) + 1./2.*m2sq*pow(m3sq,-1)*ln1*pow(
         lm,2)*pow(lm1,3) - 2*m2sq*m3sq*pow(lm1,2) + 4./3.*m2sq*m3sq*lm
         *pow(lm1,3) - 16./3.*m2sq*m3sq*ln3*pow(lm1,2) + 44./3.*m2sq*
         m3sq*ln3*lm*pow(lm1,3) - 11*m2sq*m3sq*ln3*pow(lm,2)*pow(lm1,4)
          - 16./3.*m2sq*m3sq*ln2*pow(lm1,2) + 85./6.*m2sq*m3sq*ln2*lm*
         pow(lm1,3) - 21./2.*m2sq*m3sq*ln2*pow(lm,2)*pow(lm1,4) + 22./3.
         *m2sq*m3sq*ln1*pow(lm1,2) - 173./6.*m2sq*m3sq*ln1*lm*pow(
         lm1,3) + 43./2.*m2sq*m3sq*ln1*pow(lm,2)*pow(lm1,4) + 28./3.*
         m2sq*pow(m3sq,3)*ln3*pow(lm1,3) - 8*m2sq*pow(m3sq,3)*ln3*lm*
         pow(lm1,4) - 4*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 8*m2sq*pow(
         m3sq,3)*ln2*lm*pow(lm1,4) + 2./3.*pow(m2sq,2)*pow(lm1,2) - 7./
         3.*pow(m2sq,2)*ln3*pow(lm1,2) + 23./3.*pow(m2sq,2)*ln3*lm*pow(
         lm1,3) - 4*pow(m2sq,2)*ln3*pow(lm,2)*pow(lm1,4) + 20./3.*pow(
         m2sq,2)*ln2*pow(lm1,2);
      h210 +=  - 113./6.*pow(m2sq,2)*ln2*lm*pow(lm1,3) + 25./2.*pow(
         m2sq,2)*ln2*pow(lm,2)*pow(lm1,4) - 8./3.*pow(m2sq,2)*ln1*pow(
         lm1,2) + 67./6.*pow(m2sq,2)*ln1*lm*pow(lm1,3) - 17./2.*pow(
         m2sq,2)*ln1*pow(lm,2)*pow(lm1,4) - 4*pow(m2sq,2)*pow(m3sq,2)*
         ln3*pow(lm1,3) - 4*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 1./
         3.*pow(m2sq,3)*pow(m3sq,-1)*pow(lm1,2) - 2./3.*pow(m2sq,3)*
         pow(m3sq,-1)*lm*pow(lm1,3) + pow(m2sq,3)*pow(m3sq,-1)*ln2*pow(
         lm1,2) - 7./2.*pow(m2sq,3)*pow(m3sq,-1)*ln2*lm*pow(lm1,3) + 5./
         2.*pow(m2sq,3)*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,4) - pow(
         m2sq,3)*pow(m3sq,-1)*ln1*pow(lm1,2) + 7./2.*pow(m2sq,3)*pow(
         m3sq,-1)*ln1*lm*pow(lm1,3) - 5./2.*pow(m2sq,3)*pow(m3sq,-1)*
         ln1*pow(lm,2)*pow(lm1,4) - 4*pow(m2sq,3)*m3sq*ln3*pow(lm1,3)
          + 8*pow(m2sq,3)*m3sq*ln3*lm*pow(lm1,4) + 28./3.*pow(m2sq,3)*
         m3sq*ln2*pow(lm1,3) - 8*pow(m2sq,3)*m3sq*ln2*lm*pow(lm1,4) + 8.
         /3.*pow(m2sq,4)*ln3*pow(lm1,3) - 4*pow(m2sq,4)*ln3*lm*pow(
         lm1,4);
      h210 +=  - 4*pow(m2sq,4)*ln2*pow(lm1,3) + 4*pow(m2sq,4)*ln2*lm*
         pow(lm1,4) - 1./2.*m1sq*pow(m2sq,-1)*lm*pow(lm1,2) + 1./3.*
         m1sq*pow(m2sq,-1)*pow(lm,2)*pow(lm1,3) - 5./6.*m1sq*pow(
         m2sq,-1)*ln3*lm*pow(lm1,2) + m1sq*pow(m2sq,-1)*ln3*pow(lm,2)*
         pow(lm1,3) + 5./6.*m1sq*pow(m2sq,-1)*ln1*lm*pow(lm1,2) - m1sq*
         pow(m2sq,-1)*ln1*pow(lm,2)*pow(lm1,3) - 2./3.*m1sq*pow(
         m2sq,-1)*pow(m3sq,2)*pow(lm1,2) + 3./2.*m1sq*pow(m2sq,-1)*pow(
         m3sq,2)*lm*pow(lm1,3) - 5./3.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*
         ln3*pow(lm1,2) + 59./6.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln3*lm*
         pow(lm1,3) - 9*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln3*pow(lm,2)*
         pow(lm1,4) + 5./3.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(
         lm1,2) - 59./6.*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*lm*pow(
         lm1,3) + 9*m1sq*pow(m2sq,-1)*pow(m3sq,2)*ln1*pow(lm,2)*pow(
         lm1,4) - 7./3.*m1sq*pow(m2sq,-1)*pow(m3sq,4)*ln3*pow(lm1,3) + 
         3*m1sq*pow(m2sq,-1)*pow(m3sq,4)*ln3*lm*pow(lm1,4) + 7./3.*m1sq
         *pow(m2sq,-1)*pow(m3sq,4)*ln1*pow(lm1,3);
      h210 +=  - 3*m1sq*pow(m2sq,-1)*pow(m3sq,4)*ln1*lm*pow(lm1,4) - 
         1./2.*m1sq*pow(m3sq,-1)*lm*pow(lm1,2) + 1./3.*m1sq*pow(
         m3sq,-1)*pow(lm,2)*pow(lm1,3) - 1./2.*m1sq*pow(m3sq,-1)*ln2*lm
         *pow(lm1,2) + 1./2.*m1sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,3)
          + 1./2.*m1sq*pow(m3sq,-1)*ln1*lm*pow(lm1,2) - 1./2.*m1sq*pow(
         m3sq,-1)*ln1*pow(lm,2)*pow(lm1,3) + 4./3.*m1sq*m3sq*pow(lm1,2)
          - 7./2.*m1sq*m3sq*lm*pow(lm1,3) - 7./3.*m1sq*m3sq*ln3*pow(
         lm1,2) + 34./3.*m1sq*m3sq*ln3*lm*pow(lm1,3) - 10*m1sq*m3sq*ln3
         *pow(lm,2)*pow(lm1,4) + 5*m1sq*m3sq*ln2*pow(lm1,2) - 70./3.*
         m1sq*m3sq*ln2*lm*pow(lm1,3) + 37./2.*m1sq*m3sq*ln2*pow(lm,2)*
         pow(lm1,4) - 14./3.*m1sq*m3sq*ln1*pow(lm1,2) + 12*m1sq*m3sq*
         ln1*lm*pow(lm1,3) - 17./2.*m1sq*m3sq*ln1*pow(lm,2)*pow(lm1,4)
          + 16./3.*m1sq*pow(m3sq,3)*pow(lm1,3) - 14./3.*m1sq*pow(
         m3sq,3)*ln3*pow(lm1,3) + 11*m1sq*pow(m3sq,3)*ln3*lm*pow(lm1,4)
          + 4*m1sq*pow(m3sq,3)*ln2*pow(lm1,3) - 1./2.*m1sq*pow(m3sq,3)*
         ln2*lm*pow(lm1,4);
      h210 +=  + 38./3.*m1sq*pow(m3sq,3)*ln1*pow(lm1,3) - 21./2.*m1sq
         *pow(m3sq,3)*ln1*lm*pow(lm1,4) - 8*m1sq*pow(m3sq,5)*ln3*pow(
         lm1,4) + 4./3.*m1sq*m2sq*pow(lm1,2) - 7./2.*m1sq*m2sq*lm*pow(
         lm1,3) + 16./3.*m1sq*m2sq*ln3*pow(lm1,2) - 45./2.*m1sq*m2sq*
         ln3*lm*pow(lm1,3) + 17*m1sq*m2sq*ln3*pow(lm,2)*pow(lm1,4) - 11.
         /3.*m1sq*m2sq*ln2*pow(lm1,2) + 11*m1sq*m2sq*ln2*lm*pow(lm1,3)
          - 9*m1sq*m2sq*ln2*pow(lm,2)*pow(lm1,4) - 11./3.*m1sq*m2sq*ln1
         *pow(lm1,2) + 23./2.*m1sq*m2sq*ln1*lm*pow(lm1,3) - 8*m1sq*m2sq
         *ln1*pow(lm,2)*pow(lm1,4) - 16./3.*m1sq*m2sq*pow(m3sq,2)*pow(
         lm1,3) - 12*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 13*m1sq*
         m2sq*pow(m3sq,2)*ln3*lm*pow(lm1,4) + 16*m1sq*m2sq*pow(m3sq,2)*
         ln2*pow(lm1,3) - 27*m1sq*m2sq*pow(m3sq,2)*ln2*lm*pow(lm1,4) - 
         16*m1sq*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 14*m1sq*m2sq*pow(
         m3sq,2)*ln1*lm*pow(lm1,4) + 8*m1sq*m2sq*pow(m3sq,4)*ln3*pow(
         lm1,4) + 16*m1sq*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) - 2./3.*m1sq*
         pow(m2sq,2)*pow(m3sq,-1)*pow(lm1,2);
      h210 +=  + 3./2.*m1sq*pow(m2sq,2)*pow(m3sq,-1)*lm*pow(lm1,3) - 
         2*m1sq*pow(m2sq,2)*pow(m3sq,-1)*ln2*pow(lm1,2) + 9*m1sq*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*lm*pow(lm1,3) - 15./2.*m1sq*pow(
         m2sq,2)*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,4) + 2*m1sq*pow(
         m2sq,2)*pow(m3sq,-1)*ln1*pow(lm1,2) - 9*m1sq*pow(m2sq,2)*pow(
         m3sq,-1)*ln1*lm*pow(lm1,3) + 15./2.*m1sq*pow(m2sq,2)*pow(
         m3sq,-1)*ln1*pow(lm,2)*pow(lm1,4) - 16./3.*m1sq*pow(m2sq,2)*
         m3sq*pow(lm1,3) + 46./3.*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,3)
          - 27*m1sq*pow(m2sq,2)*m3sq*ln3*lm*pow(lm1,4) - 14*m1sq*pow(
         m2sq,2)*m3sq*ln2*pow(lm1,3) + 14*m1sq*pow(m2sq,2)*m3sq*ln2*lm*
         pow(lm1,4) - 40./3.*m1sq*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + 13*
         m1sq*pow(m2sq,2)*m3sq*ln1*lm*pow(lm1,4) + 24*m1sq*pow(m2sq,2)*
         pow(m3sq,3)*ln3*pow(lm1,4) - 40*m1sq*pow(m2sq,2)*pow(m3sq,3)*
         ln2*pow(lm1,4) + 16./3.*m1sq*pow(m2sq,3)*pow(lm1,3) + 11./3.*
         m1sq*pow(m2sq,3)*ln3*pow(lm1,3) - 4*m1sq*pow(m2sq,3)*ln2*pow(
         lm1,3);
      h210 +=  + 11*m1sq*pow(m2sq,3)*ln2*lm*pow(lm1,4) + 37./3.*m1sq*
         pow(m2sq,3)*ln1*pow(lm1,3) - 11*m1sq*pow(m2sq,3)*ln1*lm*pow(
         lm1,4) - 40*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) + 24*
         m1sq*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) - 2*m1sq*pow(
         m2sq,4)*pow(m3sq,-1)*ln2*pow(lm1,3) + 5./2.*m1sq*pow(m2sq,4)*
         pow(m3sq,-1)*ln2*lm*pow(lm1,4) + 2*m1sq*pow(m2sq,4)*pow(
         m3sq,-1)*ln1*pow(lm1,3) - 5./2.*m1sq*pow(m2sq,4)*pow(m3sq,-1)*
         ln1*lm*pow(lm1,4) + 16*m1sq*pow(m2sq,4)*m3sq*ln3*pow(lm1,4) + 
         8*m1sq*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 8*m1sq*pow(m2sq,5)*
         ln2*pow(lm1,4) + 1./3.*pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,-1)*
         lm*pow(lm1,2) - 1./2.*pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,-1)*
         pow(lm,2)*pow(lm1,3) - 1./2.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*lm*
         pow(lm1,3) + 1./3.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln3*pow(
         lm1,2) - 23./3.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln3*lm*pow(
         lm1,3) + 9*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln3*pow(lm,2)*pow(
         lm1,4);
      h210 +=  - 1./3.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*pow(lm1,2)
          + 23./3.*pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*lm*pow(lm1,3) - 9*
         pow(m1sq,2)*pow(m2sq,-1)*m3sq*ln1*pow(lm,2)*pow(lm1,4) + 28./3.
         *pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln3*pow(lm1,3) - 12*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln3*lm*pow(lm1,4) - 28./3.*
         pow(m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln1*pow(lm1,3) + 12*pow(
         m1sq,2)*pow(m2sq,-1)*pow(m3sq,3)*ln1*lm*pow(lm1,4) - 13./3.*
         pow(m1sq,2)*pow(lm1,2) + 13./3.*pow(m1sq,2)*lm*pow(lm1,3) - 3*
         pow(m1sq,2)*ln3*pow(lm1,2) + 23./3.*pow(m1sq,2)*ln3*lm*pow(
         lm1,3) - 6*pow(m1sq,2)*ln3*pow(lm,2)*pow(lm1,4) - 7./3.*pow(
         m1sq,2)*ln2*pow(lm1,2) + 47./6.*pow(m1sq,2)*ln2*lm*pow(lm1,3)
          - 15./2.*pow(m1sq,2)*ln2*pow(lm,2)*pow(lm1,4) + 3*pow(m1sq,2)
         *ln1*pow(lm1,2) - 31./2.*pow(m1sq,2)*ln1*lm*pow(lm1,3) + 27./2.
         *pow(m1sq,2)*ln1*pow(lm,2)*pow(lm1,4) - 17./3.*pow(m1sq,2)*
         pow(m3sq,2)*pow(lm1,3) + 27*pow(m1sq,2)*pow(m3sq,2)*ln3*pow(
         lm1,3);
      h210 +=  - 31*pow(m1sq,2)*pow(m3sq,2)*ln3*lm*pow(lm1,4) - 21*
         pow(m1sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 23*pow(m1sq,2)*pow(
         m3sq,2)*ln2*lm*pow(lm1,4) - 46./3.*pow(m1sq,2)*pow(m3sq,2)*ln1
         *pow(lm1,3) + 8*pow(m1sq,2)*pow(m3sq,2)*ln1*lm*pow(lm1,4) + 10
         *pow(m1sq,2)*pow(m3sq,4)*ln3*pow(lm1,4) + 14*pow(m1sq,2)*pow(
         m3sq,4)*ln1*pow(lm1,4) - 1./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*
         lm*pow(lm1,3) + pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*pow(lm1,2)
          - 15./2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*lm*pow(lm1,3) + 15.
         /2.*pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln2*pow(lm,2)*pow(lm1,4) - 
         pow(m1sq,2)*m2sq*pow(m3sq,-1)*ln1*pow(lm1,2) + 15./2.*pow(
         m1sq,2)*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,3) - 15./2.*pow(
         m1sq,2)*m2sq*pow(m3sq,-1)*ln1*pow(lm,2)*pow(lm1,4) + 46./3.*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 12*pow(m1sq,2)*m2sq*m3sq*
         ln3*pow(lm1,3) + 18*pow(m1sq,2)*m2sq*m3sq*ln3*lm*pow(lm1,4) - 
         47./3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3) + 20*pow(m1sq,2)*
         m2sq*m3sq*ln2*lm*pow(lm1,4);
      h210 +=  + 107./3.*pow(m1sq,2)*m2sq*m3sq*ln1*pow(lm1,3) - 38*
         pow(m1sq,2)*m2sq*m3sq*ln1*lm*pow(lm1,4) - 56*pow(m1sq,2)*m2sq*
         pow(m3sq,3)*ln3*pow(lm1,4) + 9*pow(m1sq,2)*m2sq*pow(m3sq,3)*
         ln2*pow(lm1,4) + 15*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln1*pow(
         lm1,4) - 17./3.*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) - 59./3.*
         pow(m1sq,2)*pow(m2sq,2)*ln3*pow(lm1,3) + 21*pow(m1sq,2)*pow(
         m2sq,2)*ln3*lm*pow(lm1,4) + 70./3.*pow(m1sq,2)*pow(m2sq,2)*ln2
         *pow(lm1,3) - 29*pow(m1sq,2)*pow(m2sq,2)*ln2*lm*pow(lm1,4) - 
         13*pow(m1sq,2)*pow(m2sq,2)*ln1*pow(lm1,3) + 8*pow(m1sq,2)*pow(
         m2sq,2)*ln1*lm*pow(lm1,4) + 38*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,4) + 37*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*ln2*pow(lm1,4) - 59*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*ln1*pow(lm1,4) + 8*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,-1)*ln2*pow(lm1,3) - 10*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,-1)*ln2*lm*pow(lm1,4) - 8*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,-1)*ln1*pow(lm1,3);
      h210 +=  + 10*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,-1)*ln1*lm*pow(
      lm1,4) + 8*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) - 57*pow(
         m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) + 17*pow(m1sq,2)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,4) + 11*pow(m1sq,2)*pow(m2sq,4)*ln2*
         pow(lm1,4) + 13*pow(m1sq,2)*pow(m2sq,4)*ln1*pow(lm1,4) + 2./3.
         *pow(m1sq,3)*pow(m2sq,-1)*pow(lm1,2) - 5./6.*pow(m1sq,3)*pow(
         m2sq,-1)*lm*pow(lm1,3) + 1./3.*pow(m1sq,3)*pow(m2sq,-1)*ln3*
         pow(lm1,2) + 11./6.*pow(m1sq,3)*pow(m2sq,-1)*ln3*lm*pow(lm1,3)
          - 3*pow(m1sq,3)*pow(m2sq,-1)*ln3*pow(lm,2)*pow(lm1,4) - 1./3.
         *pow(m1sq,3)*pow(m2sq,-1)*ln1*pow(lm1,2) - 11./6.*pow(m1sq,3)*
         pow(m2sq,-1)*ln1*lm*pow(lm1,3) + 3*pow(m1sq,3)*pow(m2sq,-1)*
         ln1*pow(lm,2)*pow(lm1,4) - 14*pow(m1sq,3)*pow(m2sq,-1)*pow(
         m3sq,2)*ln3*pow(lm1,3) + 18*pow(m1sq,3)*pow(m2sq,-1)*pow(
         m3sq,2)*ln3*lm*pow(lm1,4) + 14*pow(m1sq,3)*pow(m2sq,-1)*pow(
         m3sq,2)*ln1*pow(lm1,3) - 18*pow(m1sq,3)*pow(m2sq,-1)*pow(
         m3sq,2)*ln1*lm*pow(lm1,4);
      h210 +=  + 2./3.*pow(m1sq,3)*pow(m3sq,-1)*pow(lm1,2) - 5./6.*
         pow(m1sq,3)*pow(m3sq,-1)*lm*pow(lm1,3) + 2*pow(m1sq,3)*pow(
         m3sq,-1)*ln2*lm*pow(lm1,3) - 5./2.*pow(m1sq,3)*pow(m3sq,-1)*
         ln2*pow(lm,2)*pow(lm1,4) - 2*pow(m1sq,3)*pow(m3sq,-1)*ln1*lm*
         pow(lm1,3) + 5./2.*pow(m1sq,3)*pow(m3sq,-1)*ln1*pow(lm,2)*pow(
         lm1,4) - 14./3.*pow(m1sq,3)*m3sq*pow(lm1,3) - 56./3.*pow(
         m1sq,3)*m3sq*ln3*pow(lm1,3) + 13*pow(m1sq,3)*m3sq*ln3*lm*pow(
         lm1,4) + 52./3.*pow(m1sq,3)*m3sq*ln2*pow(lm1,3) - 26*pow(
         m1sq,3)*m3sq*ln2*lm*pow(lm1,4) - 32./3.*pow(m1sq,3)*m3sq*ln1*
         pow(lm1,3) + 13*pow(m1sq,3)*m3sq*ln1*lm*pow(lm1,4) + 18*pow(
         m1sq,3)*pow(m3sq,3)*ln3*pow(lm1,4) - 34*pow(m1sq,3)*pow(
         m3sq,3)*ln1*pow(lm1,4) - 14./3.*pow(m1sq,3)*m2sq*pow(lm1,3) + 
         46./3.*pow(m1sq,3)*m2sq*ln3*pow(lm1,3) - 23*pow(m1sq,3)*m2sq*
         ln3*lm*pow(lm1,4) - 43./3.*pow(m1sq,3)*m2sq*ln2*pow(lm1,3) + 9
         *pow(m1sq,3)*m2sq*ln2*lm*pow(lm1,4) - 13*pow(m1sq,3)*m2sq*ln1*
         pow(lm1,3);
      h210 +=  + 14*pow(m1sq,3)*m2sq*ln1*lm*pow(lm1,4) + 36*pow(
         m1sq,3)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) - 53*pow(m1sq,3)*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,4) + 33*pow(m1sq,3)*m2sq*pow(m3sq,2)*
         ln1*pow(lm1,4) - 12*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln2*
         pow(lm1,3) + 15*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln2*lm*
         pow(lm1,4) + 12*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln1*pow(
         lm1,3) - 15*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,-1)*ln1*lm*pow(
         lm1,4) - 50*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) + 34*
         pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 32*pow(m1sq,3)*
         pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 15*pow(m1sq,3)*pow(m2sq,3)*
         ln2*pow(lm1,4) - 31*pow(m1sq,3)*pow(m2sq,3)*ln1*pow(lm1,4) - 1.
         /3.*pow(m1sq,4)*pow(m2sq,-1)*pow(m3sq,-1)*pow(lm1,2) + 1./2.*
         pow(m1sq,4)*pow(m2sq,-1)*pow(m3sq,-1)*lm*pow(lm1,3) + 28./3.*
         pow(m1sq,4)*pow(m2sq,-1)*m3sq*ln3*pow(lm1,3) - 12*pow(m1sq,4)*
         pow(m2sq,-1)*m3sq*ln3*lm*pow(lm1,4) - 28./3.*pow(m1sq,4)*pow(
         m2sq,-1)*m3sq*ln1*pow(lm1,3);
      h210 +=  + 12*pow(m1sq,4)*pow(m2sq,-1)*m3sq*ln1*lm*pow(lm1,4)
          + 5*pow(m1sq,4)*pow(lm1,3) + 1./3.*pow(m1sq,4)*ln3*pow(lm1,3)
          + 3*pow(m1sq,4)*ln3*lm*pow(lm1,4) - pow(m1sq,4)*ln2*pow(
         lm1,3) + 5*pow(m1sq,4)*ln2*lm*pow(lm1,4) + 34./3.*pow(m1sq,4)*
         ln1*pow(lm1,3) - 8*pow(m1sq,4)*ln1*lm*pow(lm1,4) - 34*pow(
         m1sq,4)*pow(m3sq,2)*ln3*pow(lm1,4) + 18*pow(m1sq,4)*pow(
         m3sq,2)*ln1*pow(lm1,4) + 8*pow(m1sq,4)*m2sq*pow(m3sq,-1)*ln2*
         pow(lm1,3) - 10*pow(m1sq,4)*m2sq*pow(m3sq,-1)*ln2*lm*pow(
         lm1,4) - 8*pow(m1sq,4)*m2sq*pow(m3sq,-1)*ln1*pow(lm1,3) + 10*
         pow(m1sq,4)*m2sq*pow(m3sq,-1)*ln1*lm*pow(lm1,4) + 12*pow(
         m1sq,4)*m2sq*m3sq*ln3*pow(lm1,4) + 15*pow(m1sq,4)*m2sq*m3sq*
         ln2*pow(lm1,4) - 59*pow(m1sq,4)*m2sq*m3sq*ln1*pow(lm1,4) - 31*
         pow(m1sq,4)*pow(m2sq,2)*ln2*pow(lm1,4) + 15*pow(m1sq,4)*pow(
         m2sq,2)*ln1*pow(lm1,4) - 7./3.*pow(m1sq,5)*pow(m2sq,-1)*ln3*
         pow(lm1,3) + 3*pow(m1sq,5)*pow(m2sq,-1)*ln3*lm*pow(lm1,4) + 7./
         3.*pow(m1sq,5)*pow(m2sq,-1)*ln1*pow(lm1,3);
      h210 +=  - 3*pow(m1sq,5)*pow(m2sq,-1)*ln1*lm*pow(lm1,4) - 2*
         pow(m1sq,5)*pow(m3sq,-1)*ln2*pow(lm1,3) + 5./2.*pow(m1sq,5)*
         pow(m3sq,-1)*ln2*lm*pow(lm1,4) + 2*pow(m1sq,5)*pow(m3sq,-1)*
         ln1*pow(lm1,3) - 5./2.*pow(m1sq,5)*pow(m3sq,-1)*ln1*lm*pow(
         lm1,4) + 14*pow(m1sq,5)*m3sq*ln3*pow(lm1,4) + 10*pow(m1sq,5)*
         m3sq*ln1*pow(lm1,4) + 13*pow(m1sq,5)*m2sq*ln2*pow(lm1,4) + 11*
         pow(m1sq,5)*m2sq*ln1*pow(lm1,4) - 8*pow(m1sq,6)*ln1*pow(lm1,4)
         ;
      break;
  default:
    std::cout <<"h210 (quenched) called with wrong iprop ="<<iprop<<'\n';
    h210 = 0.;
  }
 return h210*pi162;
}

double h210p(const int iprop, const double m1sq, const double m2sq,
	     const double m3sq,const double xmu2){
  double dm = m1sq-m2sq-m3sq;
  double lm = dm*dm-4.*m2sq*m3sq;
  if (fabs(lm/pow(m1sq+m2sq+m3sq,2)) < 1e-8)
    return h210psing(iprop,m1sq,m2sq,m3sq,xmu2);
  double m1 = sqrt(m1sq);
  double m2 = sqrt(m2sq);
  double m3 = sqrt(m3sq);
  double ln1 = log(m1sq/xmu2);
  double ln2 = log(m2sq/xmu2);
  double ln3 = log(m3sq/xmu2);
  //double ln4 = ln1+ln2+ln3;
  double lm1 = 1./lm;
  double h210p;
  switch(iprop){
  case 1:
    h210p = 17./288.-m1sq/(24.*lm)*(dm+2.*m1sq)
	   -5.*pow(m1sq,2)*m2sq*m3sq/(6.*lm*lm)
	   +pow(m1sq,3)*m2sq*m3sq/pow(lm,3)*(1.+5.*m2sq*m3sq/lm)*psiq(m1,m2,m3)
	   +pow(m1sq,3)*dm/(12.*lm*lm)*ln1*(1.+30.*m2sq*m3sq/lm)
	   +pow(m1sq,2)*m2sq*pow(m3sq,2)/pow(lm,3) *ln3*5./2.*(m3sq-m1sq-m2sq)
	   +m3sq/(12.*lm)*ln3*(m3sq-m2sq)
	   +m1sq*m3sq/(12.*pow(lm,2))*ln3*(-3.*pow(m1sq,2)-12.*m1sq*m2sq
		    +5.*m1sq*m3sq+3.*pow(m2sq,2)-m2sq*m3sq-2.*pow(m3,4))
	   +pow(m1sq,2)*m3sq*pow(m2sq,2)/pow(lm,3)*ln2*5./2.*(m2sq-m1sq-m3sq)
	   +m2sq/(12.*lm)*ln2*(m2sq-m3sq)
	   +m1sq*m2sq/(12.*pow(lm,2))*ln2*(-3.*pow(m1sq,2)-12.*m1sq*m3sq
	   +5.*m1sq*m2sq+3.*pow(m3sq,2)-m2sq*m3sq-2.*pow(m2,4));
    break;
  case 2:
    h210p =
       + psiq(m1,m2,m3) * ( 3*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) + 15*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,4) + 5*pow(m1sq,3)*
         m2sq*pow(m3sq,2)*pow(lm1,4) + 5*pow(m1sq,3)*pow(m2sq,2)*m3sq*
         pow(lm1,4) + 35*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,5)
          + 35*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) - 5*pow(
         m1sq,4)*m2sq*m3sq*pow(lm1,4) - 35*pow(m1sq,4)*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,5) );

      h210p +=  + 1./24.*m3sq*lm1 + 1./24.*m2sq*lm1 - 1./12.*m2sq*
         pow(m3sq,2)*ln3*pow(lm1,2) + 1./12.*m2sq*pow(m3sq,2)*ln2*pow(
         lm1,2) + 1./12.*pow(m2sq,2)*m3sq*ln3*pow(lm1,2) - 1./12.*pow(
         m2sq,2)*m3sq*ln2*pow(lm1,2) - 1./4.*m1sq*lm1 + 1./12.*m1sq*
         pow(m3sq,2)*pow(lm1,2) + 2./3.*m1sq*pow(m3sq,2)*ln3*pow(lm1,2)
          - 2./3.*m1sq*pow(m3sq,4)*ln3*pow(lm1,3) - 3./2.*m1sq*m2sq*
         m3sq*pow(lm1,2) - 11./6.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 11./6.
         *m1sq*m2sq*m3sq*ln2*pow(lm1,2) + 4*m1sq*m2sq*pow(m3sq,3)*ln3*
         pow(lm1,3) + m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 1./12.*
         m1sq*pow(m2sq,2)*pow(lm1,2) + 2./3.*m1sq*pow(m2sq,2)*ln2*pow(
         lm1,2) - 13./3.*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 
         13./3.*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + m1sq*pow(
         m2sq,3)*m3sq*ln3*pow(lm1,3) + 4*m1sq*pow(m2sq,3)*m3sq*ln2*pow(
         lm1,3) - 2./3.*m1sq*pow(m2sq,4)*ln2*pow(lm1,3) - 5./12.*pow(
         m1sq,2)*m3sq*pow(lm1,2) - 3./4.*pow(m1sq,2)*m3sq*ln3*pow(
         lm1,2);
      h210p +=  - 1./4.*pow(m1sq,2)*m3sq*ln1*pow(lm1,2) + 7./3.*pow(
         m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,3) - 5./12.*pow(m1sq,2)*m2sq*
         pow(lm1,2) - 3./4.*pow(m1sq,2)*m2sq*ln2*pow(lm1,2) - 1./4.*
         pow(m1sq,2)*m2sq*ln1*pow(lm1,2) - 35./6.*pow(m1sq,2)*m2sq*pow(
         m3sq,2)*pow(lm1,3) - 17./2.*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln3*
         pow(lm1,3) - 6*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) - 
         15./2.*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(lm1,3) + 15*pow(
         m1sq,2)*m2sq*pow(m3sq,4)*ln3*pow(lm1,4) - 35./6.*pow(m1sq,2)*
         pow(m2sq,2)*m3sq*pow(lm1,3) - 6*pow(m1sq,2)*pow(m2sq,2)*m3sq*
         ln3*pow(lm1,3) - 17./2.*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,3) - 15./2.*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln1*pow(lm1,3) + 
         5*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) - 20*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) + 7./3.*pow(
         m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,3) - 20*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) + 5*pow(m1sq,2)*pow(m2sq,3)
         *pow(m3sq,2)*ln2*pow(lm1,4);
      h210p +=  + 15*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) + 1./
         3.*pow(m1sq,3)*pow(lm1,2) + 1./3.*pow(m1sq,3)*ln1*pow(lm1,2)
          - 8./3.*pow(m1sq,3)*pow(m3sq,2)*ln3*pow(lm1,3) - 1./3.*pow(
         m1sq,3)*pow(m3sq,2)*ln1*pow(lm1,3) + 35./6.*pow(m1sq,3)*m2sq*
         m3sq*pow(lm1,3) + 4*pow(m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) + 4*
         pow(m1sq,3)*m2sq*m3sq*ln2*pow(lm1,3) + 22./3.*pow(m1sq,3)*m2sq
         *m3sq*ln1*pow(lm1,3) - 30*pow(m1sq,3)*m2sq*pow(m3sq,3)*ln3*
         pow(lm1,4) - 15*pow(m1sq,3)*m2sq*pow(m3sq,3)*ln1*pow(lm1,4) - 
         8./3.*pow(m1sq,3)*pow(m2sq,2)*ln2*pow(lm1,3) - 1./3.*pow(
         m1sq,3)*pow(m2sq,2)*ln1*pow(lm1,3) + 5*pow(m1sq,3)*pow(m2sq,2)
         *pow(m3sq,2)*ln3*pow(lm1,4) + 5*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,2)*ln2*pow(lm1,4) - 40*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,2)*ln1*pow(lm1,4) - 30*pow(m1sq,3)*pow(m2sq,3)*m3sq*ln2*
         pow(lm1,4) - 15*pow(m1sq,3)*pow(m2sq,3)*m3sq*ln1*pow(lm1,4) + 
         pow(m1sq,4)*m3sq*ln3*pow(lm1,3) + 2./3.*pow(m1sq,4)*m3sq*ln1*
         pow(lm1,3);
      h210p +=  + pow(m1sq,4)*m2sq*ln2*pow(lm1,3) + 2./3.*pow(m1sq,4)
         *m2sq*ln1*pow(lm1,3) + 15*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln3*
         pow(lm1,4) + 30*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) + 
         15*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 30*pow(
         m1sq,4)*pow(m2sq,2)*m3sq*ln1*pow(lm1,4) - 1./3.*pow(m1sq,5)*
         ln1*pow(lm1,3) - 15*pow(m1sq,5)*m2sq*m3sq*ln1*pow(lm1,4);
      break;
  case 3:
   h210p =
       + psiq(m1,m2,m3) * ( pow(m1sq,3)*m3sq*pow(lm1,3) + 15*pow(m1sq,3)
         *m2sq*pow(m3sq,2)*pow(lm1,4) - 5*pow(m1sq,3)*pow(m2sq,2)*m3sq*
         pow(lm1,4) + 35*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,5)
          - 35*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) + 5*pow(
         m1sq,4)*m2sq*m3sq*pow(lm1,4) + 35*pow(m1sq,4)*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,5) );

      h210p +=  - 1./12.*m3sq*lm1 - 1./12.*m3sq*ln3*lm1 - 1./12.*m3sq
         *ln2*lm1 + 1./6.*pow(m3sq,3)*ln3*pow(lm1,2) + 1./12.*m2sq*lm1
          + 1./6.*m2sq*ln2*lm1 - 1./3.*m2sq*pow(m3sq,2)*ln3*pow(lm1,2)
          - 1./6.*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) + 1./6.*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,2) + 1./3.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) - 
         1./6.*pow(m2sq,3)*ln2*pow(lm1,2) + 1./24.*m1sq*lm1 + 1./3.*
         m1sq*pow(m3sq,2)*pow(lm1,2) + 1./12.*m1sq*pow(m3sq,2)*ln3*pow(
         lm1,2) + 1./4.*m1sq*pow(m3sq,2)*ln2*pow(lm1,2) - 2./3.*m1sq*
         pow(m3sq,4)*ln3*pow(lm1,3) - 1./12.*m1sq*m2sq*m3sq*pow(lm1,2)
          + 1./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) - 1./3.*m1sq*m2sq*m3sq*
         ln2*pow(lm1,2) + 1./3.*m1sq*m2sq*pow(m3sq,3)*ln3*pow(lm1,3) + 
         m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) - 1./4.*m1sq*pow(m2sq,2)*
         pow(lm1,2) - 1./3.*m1sq*pow(m2sq,2)*ln2*pow(lm1,2) + 4./3.*
         m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) - 4./3.*m1sq*pow(
         m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - m1sq*pow(m2sq,3)*m3sq*ln3
         *pow(lm1,3);
      h210p +=  - 1./3.*m1sq*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) + 2./3.*
         m1sq*pow(m2sq,4)*ln2*pow(lm1,3) - 2*pow(m1sq,2)*m3sq*pow(
         lm1,2) - pow(m1sq,2)*m3sq*ln3*pow(lm1,2) - pow(m1sq,2)*m3sq*
         ln2*pow(lm1,2) + 7./2.*pow(m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,3)
          + 3./4.*pow(m1sq,2)*m2sq*pow(lm1,2) + 5./6.*pow(m1sq,2)*m2sq*
         ln2*pow(lm1,2) - 35./6.*pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(
         lm1,3) - 11*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) - 8*
         pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) + 15*pow(m1sq,2)*
         m2sq*pow(m3sq,4)*ln3*pow(lm1,4) + 35./6.*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,3) + 5*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*
         pow(lm1,3) + 77./6.*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,3) - 30*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4)
          - 15*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) - 7./
         3.*pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,3) + 15*pow(m1sq,2)*
         pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) + 30*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4);
      h210p +=  - 15*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 1./
         2.*pow(m1sq,3)*pow(lm1,2) - 1./4.*pow(m1sq,3)*ln2*pow(lm1,2)
          - 1./12.*pow(m1sq,3)*ln1*pow(lm1,2) - 5./6.*pow(m1sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,3) - 23./6.*pow(m1sq,3)*pow(m3sq,2)*ln1*
         pow(lm1,3) - 35./6.*pow(m1sq,3)*m2sq*m3sq*pow(lm1,3) - 2*pow(
         m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) - 12*pow(m1sq,3)*m2sq*m3sq*
         ln2*pow(lm1,3) - 4*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3) + 5*
         pow(m1sq,3)*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) - 20*pow(m1sq,3)*
         m2sq*pow(m3sq,3)*ln1*pow(lm1,4) + 8./3.*pow(m1sq,3)*pow(
         m2sq,2)*ln2*pow(lm1,3) + 1./3.*pow(m1sq,3)*pow(m2sq,2)*ln1*
         pow(lm1,3) + 5*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(
         lm1,4) - 40*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4)
          + 5*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) + 30*
         pow(m1sq,3)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) + 15*pow(m1sq,3)*
         pow(m2sq,3)*m3sq*ln1*pow(lm1,4) - 2*pow(m1sq,4)*m3sq*ln3*pow(
         lm1,3);
      h210p +=  + 7./2.*pow(m1sq,4)*m3sq*ln1*pow(lm1,3) - pow(m1sq,4)
         *m2sq*ln2*pow(lm1,3) - 2./3.*pow(m1sq,4)*m2sq*ln1*pow(lm1,3)
          - 20*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) + 5*pow(
         m1sq,4)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) - 15*pow(m1sq,4)*pow(
         m2sq,2)*m3sq*ln2*pow(lm1,4) - 30*pow(m1sq,4)*pow(m2sq,2)*m3sq*
         ln1*pow(lm1,4) + 1./3.*pow(m1sq,5)*ln1*pow(lm1,3) + 15*pow(
         m1sq,5)*m2sq*m3sq*ln1*pow(lm1,4);
      break;
  case 4:
   h210p =
       + psiq(m1,m2,m3) * ( pow(m1sq,3)*m2sq*pow(lm1,3) - 5*pow(m1sq,3)*
         m2sq*pow(m3sq,2)*pow(lm1,4) + 15*pow(m1sq,3)*pow(m2sq,2)*m3sq*
         pow(lm1,4) - 35*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,5)
          + 35*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) + 5*pow(
         m1sq,4)*m2sq*m3sq*pow(lm1,4) + 35*pow(m1sq,4)*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,5) );

      h210p +=  + 1./12.*m3sq*lm1 + 1./6.*m3sq*ln3*lm1 - 1./6.*pow(
         m3sq,3)*ln3*pow(lm1,2) - 1./12.*m2sq*lm1 - 1./12.*m2sq*ln3*lm1
          - 1./12.*m2sq*ln2*lm1 + 1./3.*m2sq*pow(m3sq,2)*ln3*pow(lm1,2)
          + 1./6.*m2sq*pow(m3sq,2)*ln2*pow(lm1,2) - 1./6.*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,2) - 1./3.*pow(m2sq,2)*m3sq*ln2*pow(lm1,2) + 
         1./6.*pow(m2sq,3)*ln2*pow(lm1,2) + 1./24.*m1sq*lm1 - 1./4.*
         m1sq*pow(m3sq,2)*pow(lm1,2) - 1./3.*m1sq*pow(m3sq,2)*ln3*pow(
         lm1,2) + 2./3.*m1sq*pow(m3sq,4)*ln3*pow(lm1,3) - 1./12.*m1sq*
         m2sq*m3sq*pow(lm1,2) - 1./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,2) + 1.
         /3.*m1sq*m2sq*m3sq*ln2*pow(lm1,2) - 1./3.*m1sq*m2sq*pow(
         m3sq,3)*ln3*pow(lm1,3) - m1sq*m2sq*pow(m3sq,3)*ln2*pow(lm1,3)
          + 1./3.*m1sq*pow(m2sq,2)*pow(lm1,2) + 1./4.*m1sq*pow(m2sq,2)*
         ln3*pow(lm1,2) + 1./12.*m1sq*pow(m2sq,2)*ln2*pow(lm1,2) - 4./3.
         *m1sq*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,3) + 4./3.*m1sq*pow(
         m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + m1sq*pow(m2sq,3)*m3sq*ln3
         *pow(lm1,3);
      h210p +=  + 1./3.*m1sq*pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 2./3.*
         m1sq*pow(m2sq,4)*ln2*pow(lm1,3) + 3./4.*pow(m1sq,2)*m3sq*pow(
         lm1,2) + 5./6.*pow(m1sq,2)*m3sq*ln3*pow(lm1,2) - 7./3.*pow(
         m1sq,2)*pow(m3sq,3)*ln3*pow(lm1,3) - 2*pow(m1sq,2)*m2sq*pow(
         lm1,2) - pow(m1sq,2)*m2sq*ln3*pow(lm1,2) - pow(m1sq,2)*m2sq*
         ln2*pow(lm1,2) + 35./6.*pow(m1sq,2)*m2sq*pow(m3sq,2)*pow(
         lm1,3) + 77./6.*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 
         5*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) - 15*pow(m1sq,2)
         *m2sq*pow(m3sq,4)*ln3*pow(lm1,4) - 35./6.*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,3) - 8*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*
         pow(lm1,3) - 11*pow(m1sq,2)*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) + 
         30*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 15*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) + 7./2.*
         pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,3) - 15*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 30*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4);
      h210p +=  + 15*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 1./
         2.*pow(m1sq,3)*pow(lm1,2) - 1./4.*pow(m1sq,3)*ln3*pow(lm1,2)
          - 1./12.*pow(m1sq,3)*ln1*pow(lm1,2) + 8./3.*pow(m1sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,3) + 1./3.*pow(m1sq,3)*pow(m3sq,2)*ln1*
         pow(lm1,3) - 35./6.*pow(m1sq,3)*m2sq*m3sq*pow(lm1,3) - 12*pow(
         m1sq,3)*m2sq*m3sq*ln3*pow(lm1,3) - 2*pow(m1sq,3)*m2sq*m3sq*ln2
         *pow(lm1,3) - 4*pow(m1sq,3)*m2sq*m3sq*ln1*pow(lm1,3) + 30*pow(
         m1sq,3)*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) + 15*pow(m1sq,3)*m2sq*
         pow(m3sq,3)*ln1*pow(lm1,4) - 5./6.*pow(m1sq,3)*pow(m2sq,2)*ln2
         *pow(lm1,3) - 23./6.*pow(m1sq,3)*pow(m2sq,2)*ln1*pow(lm1,3) - 
         40*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,4) + 5*pow(
         m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4) + 5*pow(m1sq,3)
         *pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) + 5*pow(m1sq,3)*pow(
         m2sq,3)*m3sq*ln2*pow(lm1,4) - 20*pow(m1sq,3)*pow(m2sq,3)*m3sq*
         ln1*pow(lm1,4) - pow(m1sq,4)*m3sq*ln3*pow(lm1,3) - 2./3.*pow(
         m1sq,4)*m3sq*ln1*pow(lm1,3);
      h210p +=  - 2*pow(m1sq,4)*m2sq*ln2*pow(lm1,3) + 7./2.*pow(
         m1sq,4)*m2sq*ln1*pow(lm1,3) - 15*pow(m1sq,4)*m2sq*pow(m3sq,2)*
         ln3*pow(lm1,4) - 30*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln1*pow(
         lm1,4) - 20*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) + 5*
         pow(m1sq,4)*pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 1./3.*pow(
         m1sq,5)*ln1*pow(lm1,3) + 15*pow(m1sq,5)*m2sq*m3sq*ln1*pow(
         lm1,4);
      break;
  case 5:
   h210p =
       + psiq(m1,m2,m3) * ( 3*pow(m1sq,2)*m3sq*pow(lm1,3) + 45*pow(
         m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) - 15*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,4) + 105*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,3)*pow(lm1,5) - 105*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,2)*
         pow(lm1,5) + 5*pow(m1sq,3)*pow(m3sq,2)*pow(lm1,4) + 25*pow(
         m1sq,3)*m2sq*m3sq*pow(lm1,4) + 105*pow(m1sq,3)*m2sq*pow(
         m3sq,3)*pow(lm1,5) + 210*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*
         pow(lm1,5) + 315*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*pow(
         lm1,6) - 35*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(lm1,5) - 315*pow(
         m1sq,3)*pow(m2sq,4)*pow(m3sq,2)*pow(lm1,6) - 5*pow(m1sq,4)*
         m3sq*pow(lm1,4) - 70*pow(m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5)
          + 70*pow(m1sq,4)*pow(m2sq,2)*m3sq*pow(lm1,5) + 630*pow(
         m1sq,4)*pow(m2sq,3)*pow(m3sq,2)*pow(lm1,6) - 35*pow(m1sq,5)*
         m2sq*m3sq*pow(lm1,5) - 315*pow(m1sq,5)*pow(m2sq,2)*pow(m3sq,2)
         *pow(lm1,6) );

      h210p +=  + 1./24.*lm1 + 1./6.*pow(m3sq,2)*pow(lm1,2) - 1./12.*
         pow(m3sq,2)*ln3*pow(lm1,2) + 1./12.*pow(m3sq,2)*ln2*pow(lm1,2)
          - 1./12.*m2sq*m3sq*pow(lm1,2) + 1./6.*m2sq*m3sq*ln3*pow(
         lm1,2) - 1./6.*m2sq*m3sq*ln2*pow(lm1,2) - 1./3.*m2sq*pow(
         m3sq,3)*ln3*pow(lm1,3) + 1./3.*m2sq*pow(m3sq,3)*ln2*pow(lm1,3)
          - 1./12.*pow(m2sq,2)*pow(lm1,2) + 2./3.*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,3) - 2./3.*pow(m2sq,2)*pow(m3sq,2)*ln2*
         pow(lm1,3) - 1./3.*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) + 1./3.*
         pow(m2sq,3)*m3sq*ln2*pow(lm1,3) - 15./4.*m1sq*m3sq*pow(lm1,2)
          - 11./6.*m1sq*m3sq*ln3*pow(lm1,2) - 11./6.*m1sq*m3sq*ln2*pow(
         lm1,2) + 4./3.*m1sq*pow(m3sq,3)*pow(lm1,3) + 20./3.*m1sq*pow(
         m3sq,3)*ln3*pow(lm1,3) + m1sq*pow(m3sq,3)*ln2*pow(lm1,3) - 4*
         m1sq*pow(m3sq,5)*ln3*pow(lm1,4) + 17./12.*m1sq*m2sq*pow(lm1,2)
          + 4./3.*m1sq*m2sq*ln2*pow(lm1,2) - 32./3.*m1sq*m2sq*pow(
         m3sq,2)*pow(lm1,3) - 19*m1sq*m2sq*pow(m3sq,2)*ln3*pow(lm1,3)
          - 47./3.*m1sq*m2sq*pow(m3sq,2)*ln2*pow(lm1,3);
      h210p +=  + 28*m1sq*m2sq*pow(m3sq,4)*ln3*pow(lm1,4) + 6*m1sq*
         m2sq*pow(m3sq,4)*ln2*pow(lm1,4) + 31./3.*m1sq*pow(m2sq,2)*m3sq
         *pow(lm1,3) + 32./3.*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) + 65.
         /3.*m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 50*m1sq*pow(m2sq,2)
         *pow(m3sq,3)*ln3*pow(lm1,4) - 32*m1sq*pow(m2sq,2)*pow(m3sq,3)*
         ln2*pow(lm1,4) - m1sq*pow(m2sq,3)*pow(lm1,3) - 16./3.*m1sq*
         pow(m2sq,3)*ln2*pow(lm1,3) + 32*m1sq*pow(m2sq,3)*pow(m3sq,2)*
         ln3*pow(lm1,4) + 50*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(
         lm1,4) - 6*m1sq*pow(m2sq,4)*m3sq*ln3*pow(lm1,4) - 28*m1sq*pow(
         m2sq,4)*m3sq*ln2*pow(lm1,4) + 4*m1sq*pow(m2sq,5)*ln2*pow(
         lm1,4) - 5./3.*pow(m1sq,2)*pow(lm1,2) - 3./4.*pow(m1sq,2)*ln2*
         pow(lm1,2) - 1./4.*pow(m1sq,2)*ln1*pow(lm1,2) - 79./6.*pow(
         m1sq,2)*pow(m3sq,2)*pow(lm1,3) - 35./6.*pow(m1sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,3) - 6*pow(m1sq,2)*pow(m3sq,2)*ln2*pow(
         lm1,3) - 23./2.*pow(m1sq,2)*pow(m3sq,2)*ln1*pow(lm1,3) + 25*
         pow(m1sq,2)*pow(m3sq,4)*ln3*pow(lm1,4);
      h210p +=  - 157./6.*pow(m1sq,2)*m2sq*m3sq*pow(lm1,3) - 40./3.*
         pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) - 100./3.*pow(m1sq,2)*
         m2sq*m3sq*ln2*pow(lm1,3) - 12*pow(m1sq,2)*m2sq*m3sq*ln1*pow(
         lm1,3) - 55*pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(lm1,4) - 16*pow(
         m1sq,2)*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) - 70*pow(m1sq,2)*m2sq*
         pow(m3sq,3)*ln2*pow(lm1,4) - 60*pow(m1sq,2)*m2sq*pow(m3sq,3)*
         ln1*pow(lm1,4) + 120*pow(m1sq,2)*m2sq*pow(m3sq,5)*ln3*pow(
         lm1,5) + 13./3.*pow(m1sq,2)*pow(m2sq,2)*pow(lm1,3) + 38./3.*
         pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) + pow(m1sq,2)*pow(
         m2sq,2)*ln1*pow(lm1,3) + 5*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)
         *pow(lm1,4) - 56*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(
         lm1,4) - 56*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4)
          + 15*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,4) - 80*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln3*pow(lm1,5) - 160*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln2*pow(lm1,5) + 50*pow(
         m1sq,2)*pow(m2sq,3)*m3sq*pow(lm1,4);
      h210p +=  + 42*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) + 
         149*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) + 45*pow(
         m1sq,2)*pow(m2sq,3)*m3sq*ln1*pow(lm1,4) - 200*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,3)*ln3*pow(lm1,5) + 200*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,3)*ln2*pow(lm1,5) - 18*pow(m1sq,2)*pow(
         m2sq,4)*ln2*pow(lm1,4) + 160*pow(m1sq,2)*pow(m2sq,4)*pow(
         m3sq,2)*ln3*pow(lm1,5) + 80*pow(m1sq,2)*pow(m2sq,4)*pow(
         m3sq,2)*ln2*pow(lm1,5) - 120*pow(m1sq,2)*pow(m2sq,5)*m3sq*ln2*
         pow(lm1,5) + 19./2.*pow(m1sq,3)*m3sq*pow(lm1,3) - 2*pow(
         m1sq,3)*m3sq*ln3*pow(lm1,3) + 4*pow(m1sq,3)*m3sq*ln2*pow(
         lm1,3) + 32./3.*pow(m1sq,3)*m3sq*ln1*pow(lm1,3) - 27*pow(
         m1sq,3)*pow(m3sq,3)*ln3*pow(lm1,4) - 22*pow(m1sq,3)*pow(
         m3sq,3)*ln1*pow(lm1,4) - 17./3.*pow(m1sq,3)*m2sq*pow(lm1,3) - 
         25./3.*pow(m1sq,3)*m2sq*ln2*pow(lm1,3) - 3*pow(m1sq,3)*m2sq*
         ln1*pow(lm1,3) + 5*pow(m1sq,3)*m2sq*pow(m3sq,2)*pow(lm1,4) - 6
         *pow(m1sq,3)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4);
      h210p +=  - 12*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln2*pow(lm1,4) - 64
         *pow(m1sq,3)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) - 85*pow(m1sq,3)*
         m2sq*pow(m3sq,4)*ln3*pow(lm1,5) - 155*pow(m1sq,3)*m2sq*pow(
         m3sq,4)*ln1*pow(lm1,5) - 100*pow(m1sq,3)*pow(m2sq,2)*m3sq*pow(
         lm1,4) - 55*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) - 191*
         pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) - 131*pow(m1sq,3)*
         pow(m2sq,2)*m3sq*ln1*pow(lm1,4) + 390*pow(m1sq,3)*pow(m2sq,2)*
         pow(m3sq,3)*ln3*pow(lm1,5) - 190*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,3)*ln2*pow(lm1,5) - 200*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,3)*ln1*pow(lm1,5) + 30*pow(m1sq,3)*pow(m2sq,3)*ln2*pow(
         lm1,4) + 2*pow(m1sq,3)*pow(m2sq,3)*ln1*pow(lm1,4) - 165*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,5) - 310*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,5) + 235*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln1*pow(lm1,5) + 360*pow(
         m1sq,3)*pow(m2sq,4)*m3sq*ln2*pow(lm1,5) + 120*pow(m1sq,3)*pow(
         m2sq,4)*m3sq*ln1*pow(lm1,5);
      h210p +=  + 7./3.*pow(m1sq,4)*pow(lm1,3) + pow(m1sq,4)*ln2*pow(
         lm1,3) + 2*pow(m1sq,4)*ln1*pow(lm1,3) - 5*pow(m1sq,4)*pow(
         m3sq,2)*ln3*pow(lm1,4) + 42*pow(m1sq,4)*pow(m3sq,2)*ln1*pow(
         lm1,4) + 50*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) + 8*pow(m1sq,4)*
         m2sq*m3sq*ln3*pow(lm1,4) + 70*pow(m1sq,4)*m2sq*m3sq*ln2*pow(
         lm1,4) + 104*pow(m1sq,4)*m2sq*m3sq*ln1*pow(lm1,4) - 190*pow(
         m1sq,4)*m2sq*pow(m3sq,3)*ln3*pow(lm1,5) + 190*pow(m1sq,4)*m2sq
         *pow(m3sq,3)*ln1*pow(lm1,5) - 22*pow(m1sq,4)*pow(m2sq,2)*ln2*
         pow(lm1,4) - 6*pow(m1sq,4)*pow(m2sq,2)*ln1*pow(lm1,4) - 150*
         pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,5) + 230*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,5) - 320*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,5) - 360*pow(
         m1sq,4)*pow(m2sq,3)*m3sq*ln2*pow(lm1,5) - 360*pow(m1sq,4)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,5) + 11*pow(m1sq,5)*m3sq*ln3*pow(
         lm1,4) - 18*pow(m1sq,5)*m3sq*ln1*pow(lm1,4) + 6*pow(m1sq,5)*
         m2sq*ln2*pow(lm1,4);
      h210p +=  + 6*pow(m1sq,5)*m2sq*ln1*pow(lm1,4) + 155*pow(m1sq,5)
         *m2sq*pow(m3sq,2)*ln3*pow(lm1,5) + 85*pow(m1sq,5)*m2sq*pow(
         m3sq,2)*ln1*pow(lm1,5) + 120*pow(m1sq,5)*pow(m2sq,2)*m3sq*ln2*
         pow(lm1,5) + 360*pow(m1sq,5)*pow(m2sq,2)*m3sq*ln1*pow(lm1,5)
          - 2*pow(m1sq,6)*ln1*pow(lm1,4) - 120*pow(m1sq,6)*m2sq*m3sq*
         ln1*pow(lm1,5);
      break;
  case 6:
   h210p =
       + psiq(m1,m2,m3) * ( 3*pow(m1sq,2)*m2sq*pow(lm1,3) - 15*pow(
         m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) + 45*pow(m1sq,2)*pow(
         m2sq,2)*m3sq*pow(lm1,4) - 105*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,3)*pow(lm1,5) + 105*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,2)*
         pow(lm1,5) + 25*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) - 35*pow(
         m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5) + 5*pow(m1sq,3)*pow(
         m2sq,2)*pow(lm1,4) + 210*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*
         pow(lm1,5) - 315*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*pow(
         lm1,6) + 105*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(lm1,5) + 315*
         pow(m1sq,3)*pow(m2sq,4)*pow(m3sq,2)*pow(lm1,6) - 5*pow(m1sq,4)
         *m2sq*pow(lm1,4) + 70*pow(m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5)
          - 70*pow(m1sq,4)*pow(m2sq,2)*m3sq*pow(lm1,5) + 630*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,6) - 35*pow(m1sq,5)*
         m2sq*m3sq*pow(lm1,5) - 315*pow(m1sq,5)*pow(m2sq,2)*pow(m3sq,2)
         *pow(lm1,6) );

      h210p +=  + 1./24.*lm1 - 1./12.*pow(m3sq,2)*pow(lm1,2) - 1./12.
         *m2sq*m3sq*pow(lm1,2) - 1./6.*m2sq*m3sq*ln3*pow(lm1,2) + 1./6.
         *m2sq*m3sq*ln2*pow(lm1,2) + 1./3.*m2sq*pow(m3sq,3)*ln3*pow(
         lm1,3) - 1./3.*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 1./6.*pow(
         m2sq,2)*pow(lm1,2) + 1./12.*pow(m2sq,2)*ln3*pow(lm1,2) - 1./12.
         *pow(m2sq,2)*ln2*pow(lm1,2) - 2./3.*pow(m2sq,2)*pow(m3sq,2)*
         ln3*pow(lm1,3) + 2./3.*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3)
          + 1./3.*pow(m2sq,3)*m3sq*ln3*pow(lm1,3) - 1./3.*pow(m2sq,3)*
         m3sq*ln2*pow(lm1,3) + 17./12.*m1sq*m3sq*pow(lm1,2) + 4./3.*
         m1sq*m3sq*ln3*pow(lm1,2) - m1sq*pow(m3sq,3)*pow(lm1,3) - 16./3.
         *m1sq*pow(m3sq,3)*ln3*pow(lm1,3) + 4*m1sq*pow(m3sq,5)*ln3*pow(
         lm1,4) - 15./4.*m1sq*m2sq*pow(lm1,2) - 11./6.*m1sq*m2sq*ln3*
         pow(lm1,2) - 11./6.*m1sq*m2sq*ln2*pow(lm1,2) + 31./3.*m1sq*
         m2sq*pow(m3sq,2)*pow(lm1,3) + 65./3.*m1sq*m2sq*pow(m3sq,2)*ln3
         *pow(lm1,3) + 32./3.*m1sq*m2sq*pow(m3sq,2)*ln2*pow(lm1,3) - 28
         *m1sq*m2sq*pow(m3sq,4)*ln3*pow(lm1,4);
      h210p +=  - 6*m1sq*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) - 32./3.*
         m1sq*pow(m2sq,2)*m3sq*pow(lm1,3) - 47./3.*m1sq*pow(m2sq,2)*
         m3sq*ln3*pow(lm1,3) - 19*m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3)
          + 50*m1sq*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 32*m1sq*
         pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) + 4./3.*m1sq*pow(
         m2sq,3)*pow(lm1,3) + m1sq*pow(m2sq,3)*ln3*pow(lm1,3) + 20./3.*
         m1sq*pow(m2sq,3)*ln2*pow(lm1,3) - 32*m1sq*pow(m2sq,3)*pow(
         m3sq,2)*ln3*pow(lm1,4) - 50*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln2*
         pow(lm1,4) + 6*m1sq*pow(m2sq,4)*m3sq*ln3*pow(lm1,4) + 28*m1sq*
         pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 4*m1sq*pow(m2sq,5)*ln2*pow(
         lm1,4) - 5./3.*pow(m1sq,2)*pow(lm1,2) - 3./4.*pow(m1sq,2)*ln3*
         pow(lm1,2) - 1./4.*pow(m1sq,2)*ln1*pow(lm1,2) + 13./3.*pow(
         m1sq,2)*pow(m3sq,2)*pow(lm1,3) + 38./3.*pow(m1sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,3) + pow(m1sq,2)*pow(m3sq,2)*ln1*pow(
         lm1,3) - 18*pow(m1sq,2)*pow(m3sq,4)*ln3*pow(lm1,4) - 157./6.*
         pow(m1sq,2)*m2sq*m3sq*pow(lm1,3);
      h210p +=  - 100./3.*pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) - 40./
         3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3) - 12*pow(m1sq,2)*m2sq*
         m3sq*ln1*pow(lm1,3) + 50*pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(
         lm1,4) + 149*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) + 42*
         pow(m1sq,2)*m2sq*pow(m3sq,3)*ln2*pow(lm1,4) + 45*pow(m1sq,2)*
         m2sq*pow(m3sq,3)*ln1*pow(lm1,4) - 120*pow(m1sq,2)*m2sq*pow(
         m3sq,5)*ln3*pow(lm1,5) - 79./6.*pow(m1sq,2)*pow(m2sq,2)*pow(
         lm1,3) - 6*pow(m1sq,2)*pow(m2sq,2)*ln3*pow(lm1,3) - 35./6.*
         pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) - 23./2.*pow(m1sq,2)*
         pow(m2sq,2)*ln1*pow(lm1,3) + 5*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,2)*pow(lm1,4) - 56*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*
         ln3*pow(lm1,4) - 56*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln2*
         pow(lm1,4) + 15*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(
         lm1,4) + 80*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln3*pow(lm1,5)
          + 160*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln2*pow(lm1,5) - 55
         *pow(m1sq,2)*pow(m2sq,3)*m3sq*pow(lm1,4);
      h210p +=  - 70*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln3*pow(lm1,4) - 16
         *pow(m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) - 60*pow(m1sq,2)*
         pow(m2sq,3)*m3sq*ln1*pow(lm1,4) + 200*pow(m1sq,2)*pow(m2sq,3)*
         pow(m3sq,3)*ln3*pow(lm1,5) - 200*pow(m1sq,2)*pow(m2sq,3)*pow(
         m3sq,3)*ln2*pow(lm1,5) + 25*pow(m1sq,2)*pow(m2sq,4)*ln2*pow(
         lm1,4) - 160*pow(m1sq,2)*pow(m2sq,4)*pow(m3sq,2)*ln3*pow(
         lm1,5) - 80*pow(m1sq,2)*pow(m2sq,4)*pow(m3sq,2)*ln2*pow(lm1,5)
          + 120*pow(m1sq,2)*pow(m2sq,5)*m3sq*ln2*pow(lm1,5) - 17./3.*
         pow(m1sq,3)*m3sq*pow(lm1,3) - 25./3.*pow(m1sq,3)*m3sq*ln3*pow(
         lm1,3) - 3*pow(m1sq,3)*m3sq*ln1*pow(lm1,3) + 30*pow(m1sq,3)*
         pow(m3sq,3)*ln3*pow(lm1,4) + 2*pow(m1sq,3)*pow(m3sq,3)*ln1*
         pow(lm1,4) + 19./2.*pow(m1sq,3)*m2sq*pow(lm1,3) + 4*pow(
         m1sq,3)*m2sq*ln3*pow(lm1,3) - 2*pow(m1sq,3)*m2sq*ln2*pow(
         lm1,3) + 32./3.*pow(m1sq,3)*m2sq*ln1*pow(lm1,3) - 100*pow(
         m1sq,3)*m2sq*pow(m3sq,2)*pow(lm1,4) - 191*pow(m1sq,3)*m2sq*
         pow(m3sq,2)*ln3*pow(lm1,4);
      h210p +=  - 55*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln2*pow(lm1,4) - 
         131*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) + 360*pow(
         m1sq,3)*m2sq*pow(m3sq,4)*ln3*pow(lm1,5) + 120*pow(m1sq,3)*m2sq
         *pow(m3sq,4)*ln1*pow(lm1,5) + 5*pow(m1sq,3)*pow(m2sq,2)*m3sq*
         pow(lm1,4) - 12*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) - 
         6*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) - 64*pow(m1sq,3)
         *pow(m2sq,2)*m3sq*ln1*pow(lm1,4) - 310*pow(m1sq,3)*pow(m2sq,2)
         *pow(m3sq,3)*ln3*pow(lm1,5) - 165*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,3)*ln2*pow(lm1,5) + 235*pow(m1sq,3)*pow(m2sq,2)*pow(
         m3sq,3)*ln1*pow(lm1,5) - 27*pow(m1sq,3)*pow(m2sq,3)*ln2*pow(
         lm1,4) - 22*pow(m1sq,3)*pow(m2sq,3)*ln1*pow(lm1,4) - 190*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,5) + 390*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,5) - 200*pow(
         m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln1*pow(lm1,5) - 85*pow(
         m1sq,3)*pow(m2sq,4)*m3sq*ln2*pow(lm1,5) - 155*pow(m1sq,3)*pow(
         m2sq,4)*m3sq*ln1*pow(lm1,5);
      h210p +=  + 7./3.*pow(m1sq,4)*pow(lm1,3) + pow(m1sq,4)*ln3*pow(
         lm1,3) + 2*pow(m1sq,4)*ln1*pow(lm1,3) - 22*pow(m1sq,4)*pow(
         m3sq,2)*ln3*pow(lm1,4) - 6*pow(m1sq,4)*pow(m3sq,2)*ln1*pow(
         lm1,4) + 50*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) + 70*pow(m1sq,4)*
         m2sq*m3sq*ln3*pow(lm1,4) + 8*pow(m1sq,4)*m2sq*m3sq*ln2*pow(
         lm1,4) + 104*pow(m1sq,4)*m2sq*m3sq*ln1*pow(lm1,4) - 360*pow(
         m1sq,4)*m2sq*pow(m3sq,3)*ln3*pow(lm1,5) - 360*pow(m1sq,4)*m2sq
         *pow(m3sq,3)*ln1*pow(lm1,5) - 5*pow(m1sq,4)*pow(m2sq,2)*ln2*
         pow(lm1,4) + 42*pow(m1sq,4)*pow(m2sq,2)*ln1*pow(lm1,4) + 230*
         pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,5) - 150*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,5) - 320*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,5) - 190*pow(
         m1sq,4)*pow(m2sq,3)*m3sq*ln2*pow(lm1,5) + 190*pow(m1sq,4)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,5) + 6*pow(m1sq,5)*m3sq*ln3*pow(
         lm1,4) + 6*pow(m1sq,5)*m3sq*ln1*pow(lm1,4) + 11*pow(m1sq,5)*
         m2sq*ln2*pow(lm1,4);
      h210p +=  - 18*pow(m1sq,5)*m2sq*ln1*pow(lm1,4) + 120*pow(
         m1sq,5)*m2sq*pow(m3sq,2)*ln3*pow(lm1,5) + 360*pow(m1sq,5)*m2sq
         *pow(m3sq,2)*ln1*pow(lm1,5) + 155*pow(m1sq,5)*pow(m2sq,2)*m3sq
         *ln2*pow(lm1,5) + 85*pow(m1sq,5)*pow(m2sq,2)*m3sq*ln1*pow(
         lm1,5) - 2*pow(m1sq,6)*ln1*pow(lm1,4) - 120*pow(m1sq,6)*m2sq*
         m3sq*ln1*pow(lm1,5);
      break;
  case 7:
   h210p =
       + psiq(m1,m2,m3) * ( pow(m1sq,3)*pow(lm1,3) - 5*pow(m1sq,3)*pow(
         m3sq,2)*pow(lm1,4) + 35*pow(m1sq,3)*m2sq*m3sq*pow(lm1,4) - 105
         *pow(m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5) - 5*pow(m1sq,3)*pow(
         m2sq,2)*pow(lm1,4) + 245*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*
         pow(lm1,5) - 315*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*pow(
         lm1,6) - 105*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(lm1,5) + 630*
         pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,3)*pow(lm1,6) - 315*pow(
         m1sq,3)*pow(m2sq,4)*pow(m3sq,2)*pow(lm1,6) + 5*pow(m1sq,4)*
         m3sq*pow(lm1,4) + 5*pow(m1sq,4)*m2sq*pow(lm1,4) + 70*pow(
         m1sq,4)*m2sq*pow(m3sq,2)*pow(lm1,5) + 70*pow(m1sq,4)*pow(
         m2sq,2)*m3sq*pow(lm1,5) + 35*pow(m1sq,5)*m2sq*m3sq*pow(lm1,5)
          + 315*pow(m1sq,5)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,6) );

      h210p +=  - 1./6.*lm1 - 1./12.*ln3*lm1 - 1./12.*ln2*lm1 + 1./3.
         *pow(m3sq,2)*pow(lm1,2) + 2./3.*pow(m3sq,2)*ln3*pow(lm1,2) + 1.
         /6.*pow(m3sq,2)*ln2*pow(lm1,2) - 2./3.*pow(m3sq,4)*ln3*pow(
         lm1,3) - 2./3.*m2sq*m3sq*pow(lm1,2) - 5./6.*m2sq*m3sq*ln3*pow(
         lm1,2) - 5./6.*m2sq*m3sq*ln2*pow(lm1,2) + 2*m2sq*pow(m3sq,3)*
         ln3*pow(lm1,3) + 2./3.*m2sq*pow(m3sq,3)*ln2*pow(lm1,3) + 1./3.
         *pow(m2sq,2)*pow(lm1,2) + 1./6.*pow(m2sq,2)*ln3*pow(lm1,2) + 2.
         /3.*pow(m2sq,2)*ln2*pow(lm1,2) - 2*pow(m2sq,2)*pow(m3sq,2)*ln3
         *pow(lm1,3) - 2*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) + 2./3.
         *pow(m2sq,3)*m3sq*ln3*pow(lm1,3) + 2*pow(m2sq,3)*m3sq*ln2*pow(
         lm1,3) - 2./3.*pow(m2sq,4)*ln2*pow(lm1,3) + 1./2.*m1sq*m3sq*
         pow(lm1,2) + 1./3.*m1sq*m3sq*ln2*pow(lm1,2) - 2*m1sq*pow(
         m3sq,3)*pow(lm1,3) - 7./3.*m1sq*pow(m3sq,3)*ln3*pow(lm1,3) - 
         m1sq*pow(m3sq,3)*ln2*pow(lm1,3) + 4*m1sq*pow(m3sq,5)*ln3*pow(
         lm1,4) + 1./2.*m1sq*m2sq*pow(lm1,2) + 1./3.*m1sq*m2sq*ln3*pow(
         lm1,2);
      h210p +=  + 2*m1sq*m2sq*pow(m3sq,2)*pow(lm1,3) - 4./3.*m1sq*
         m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 14./3.*m1sq*m2sq*pow(m3sq,2)
         *ln2*pow(lm1,3) - 6*m1sq*m2sq*pow(m3sq,4)*ln3*pow(lm1,4) - 6*
         m1sq*m2sq*pow(m3sq,4)*ln2*pow(lm1,4) + 2*m1sq*pow(m2sq,2)*m3sq
         *pow(lm1,3) + 14./3.*m1sq*pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 4./
         3.*m1sq*pow(m2sq,2)*m3sq*ln2*pow(lm1,3) - 6*m1sq*pow(m2sq,2)*
         pow(m3sq,3)*ln3*pow(lm1,4) + 14*m1sq*pow(m2sq,2)*pow(m3sq,3)*
         ln2*pow(lm1,4) - 2*m1sq*pow(m2sq,3)*pow(lm1,3) - m1sq*pow(
         m2sq,3)*ln3*pow(lm1,3) - 7./3.*m1sq*pow(m2sq,3)*ln2*pow(lm1,3)
          + 14*m1sq*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 6*m1sq*
         pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) - 6*m1sq*pow(m2sq,4)*
         m3sq*ln3*pow(lm1,4) - 6*m1sq*pow(m2sq,4)*m3sq*ln2*pow(lm1,4)
          + 4*m1sq*pow(m2sq,5)*ln2*pow(lm1,4) - 35./12.*pow(m1sq,2)*
         pow(lm1,2) - pow(m1sq,2)*ln3*pow(lm1,2) - pow(m1sq,2)*ln2*pow(
         lm1,2) + 77./6.*pow(m1sq,2)*pow(m3sq,2)*pow(lm1,3) + 89./6.*
         pow(m1sq,2)*pow(m3sq,2)*ln3*pow(lm1,3);
      h210p +=  + 5*pow(m1sq,2)*pow(m3sq,2)*ln2*pow(lm1,3) - 25*pow(
         m1sq,2)*pow(m3sq,4)*ln3*pow(lm1,4) - 34*pow(m1sq,2)*m2sq*m3sq*
         pow(lm1,3) - 74./3.*pow(m1sq,2)*m2sq*m3sq*ln3*pow(lm1,3) - 74./
         3.*pow(m1sq,2)*m2sq*m3sq*ln2*pow(lm1,3) + 50*pow(m1sq,2)*m2sq*
         pow(m3sq,3)*pow(lm1,4) + 149*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln3*
         pow(lm1,4) + 54*pow(m1sq,2)*m2sq*pow(m3sq,3)*ln2*pow(lm1,4) - 
         120*pow(m1sq,2)*m2sq*pow(m3sq,5)*ln3*pow(lm1,5) + 77./6.*pow(
         m1sq,2)*pow(m2sq,2)*pow(lm1,3) + 5*pow(m1sq,2)*pow(m2sq,2)*ln3
         *pow(lm1,3) + 89./6.*pow(m1sq,2)*pow(m2sq,2)*ln2*pow(lm1,3) - 
         100*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 178*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,4) - 178*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,4) + 360*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln3*pow(lm1,5) + 120*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,4)*ln2*pow(lm1,5) + 50*pow(
         m1sq,2)*pow(m2sq,3)*m3sq*pow(lm1,4) + 54*pow(m1sq,2)*pow(
         m2sq,3)*m3sq*ln3*pow(lm1,4);
      h210p +=  + 149*pow(m1sq,2)*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) - 
         360*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,3)*ln3*pow(lm1,5) - 360*
         pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,3)*ln2*pow(lm1,5) - 25*pow(
         m1sq,2)*pow(m2sq,4)*ln2*pow(lm1,4) + 120*pow(m1sq,2)*pow(
         m2sq,4)*pow(m3sq,2)*ln3*pow(lm1,5) + 360*pow(m1sq,2)*pow(
         m2sq,4)*pow(m3sq,2)*ln2*pow(lm1,5) - 120*pow(m1sq,2)*pow(
         m2sq,5)*m3sq*ln2*pow(lm1,5) - 41./6.*pow(m1sq,3)*m3sq*pow(
         lm1,3) - 23./3.*pow(m1sq,3)*m3sq*ln3*pow(lm1,3) - 2*pow(
         m1sq,3)*m3sq*ln2*pow(lm1,3) - 19./3.*pow(m1sq,3)*m3sq*ln1*pow(
         lm1,3) + 26*pow(m1sq,3)*pow(m3sq,3)*ln3*pow(lm1,4) + 23*pow(
         m1sq,3)*pow(m3sq,3)*ln1*pow(lm1,4) - 41./6.*pow(m1sq,3)*m2sq*
         pow(lm1,3) - 2*pow(m1sq,3)*m2sq*ln3*pow(lm1,3) - 26./3.*pow(
         m1sq,3)*m2sq*ln2*pow(lm1,3) - 16./3.*pow(m1sq,3)*m2sq*ln1*pow(
         lm1,3) + 5*pow(m1sq,3)*m2sq*pow(m3sq,2)*pow(lm1,4) - 74*pow(
         m1sq,3)*m2sq*pow(m3sq,2)*ln3*pow(lm1,4) + 39*pow(m1sq,3)*m2sq*
         pow(m3sq,2)*ln2*pow(lm1,4);
      h210p +=  - 44*pow(m1sq,3)*m2sq*pow(m3sq,2)*ln1*pow(lm1,4) + 80
         *pow(m1sq,3)*m2sq*pow(m3sq,4)*ln3*pow(lm1,5) + 160*pow(m1sq,3)
         *m2sq*pow(m3sq,4)*ln1*pow(lm1,5) + 5*pow(m1sq,3)*pow(m2sq,2)*
         m3sq*pow(lm1,4) + 38*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln3*pow(
         lm1,4) - 81*pow(m1sq,3)*pow(m2sq,2)*m3sq*ln2*pow(lm1,4) - 36*
         pow(m1sq,3)*pow(m2sq,2)*m3sq*ln1*pow(lm1,4) - 310*pow(m1sq,3)*
         pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,5) + 235*pow(m1sq,3)*pow(
         m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,5) - 165*pow(m1sq,3)*pow(
         m2sq,2)*pow(m3sq,3)*ln1*pow(lm1,5) + 27*pow(m1sq,3)*pow(
         m2sq,3)*ln2*pow(lm1,4) + 22*pow(m1sq,3)*pow(m2sq,3)*ln1*pow(
         lm1,4) + 230*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(
         lm1,5) - 320*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln2*pow(
         lm1,5) - 150*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,2)*ln1*pow(
         lm1,5) + 85*pow(m1sq,3)*pow(m2sq,4)*m3sq*ln2*pow(lm1,5) + 155*
         pow(m1sq,3)*pow(m2sq,4)*m3sq*ln1*pow(lm1,5) - 4*pow(m1sq,4)*
         pow(lm1,3);
      h210p +=  - 2*pow(m1sq,4)*ln3*pow(lm1,3) - 2*pow(m1sq,4)*ln2*
         pow(lm1,3) + 25./6.*pow(m1sq,4)*ln1*pow(lm1,3) + 7*pow(m1sq,4)
         *pow(m3sq,2)*ln3*pow(lm1,4) - 44*pow(m1sq,4)*pow(m3sq,2)*ln1*
         pow(lm1,4) - 55*pow(m1sq,4)*m2sq*m3sq*pow(lm1,4) - 74*pow(
         m1sq,4)*m2sq*m3sq*ln3*pow(lm1,4) - 76*pow(m1sq,4)*m2sq*m3sq*
         ln2*pow(lm1,4) + 31*pow(m1sq,4)*m2sq*m3sq*ln1*pow(lm1,4) + 200
         *pow(m1sq,4)*m2sq*pow(m3sq,3)*ln3*pow(lm1,5) - 200*pow(m1sq,4)
         *m2sq*pow(m3sq,3)*ln1*pow(lm1,5) + 5*pow(m1sq,4)*pow(m2sq,2)*
         ln2*pow(lm1,4) - 42*pow(m1sq,4)*pow(m2sq,2)*ln1*pow(lm1,4) - 
         190*pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,5) - 200*
         pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,5) + 390*pow(
         m1sq,4)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,5) + 190*pow(
         m1sq,4)*pow(m2sq,3)*m3sq*ln2*pow(lm1,5) - 190*pow(m1sq,4)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,5) - 12*pow(m1sq,5)*m3sq*ln3*pow(
         lm1,4) + 19*pow(m1sq,5)*m3sq*ln1*pow(lm1,4) - 11*pow(m1sq,5)*
         m2sq*ln2*pow(lm1,4);
      h210p +=  + 18*pow(m1sq,5)*m2sq*ln1*pow(lm1,4) - 160*pow(
         m1sq,5)*m2sq*pow(m3sq,2)*ln3*pow(lm1,5) - 80*pow(m1sq,5)*m2sq*
         pow(m3sq,2)*ln1*pow(lm1,5) - 155*pow(m1sq,5)*pow(m2sq,2)*m3sq*
         ln2*pow(lm1,5) - 85*pow(m1sq,5)*pow(m2sq,2)*m3sq*ln1*pow(
         lm1,5) + 2*pow(m1sq,6)*ln1*pow(lm1,4) + 120*pow(m1sq,6)*m2sq*
         m3sq*ln1*pow(lm1,5);
      break;
  case 8:
   h210p =
       + psiq(m1,m2,m3) * ( 3*pow(m1sq,2)*pow(lm1,3) - 15*pow(m1sq,2)*
         pow(m3sq,2)*pow(lm1,4) + 105*pow(m1sq,2)*m2sq*m3sq*pow(lm1,4)
          - 315*pow(m1sq,2)*m2sq*pow(m3sq,3)*pow(lm1,5) - 15*pow(
         m1sq,2)*pow(m2sq,2)*pow(lm1,4) + 735*pow(m1sq,2)*pow(m2sq,2)*
         pow(m3sq,2)*pow(lm1,5) - 945*pow(m1sq,2)*pow(m2sq,2)*pow(
         m3sq,4)*pow(lm1,6) - 315*pow(m1sq,2)*pow(m2sq,3)*m3sq*pow(
         lm1,5) + 1890*pow(m1sq,2)*pow(m2sq,3)*pow(m3sq,3)*pow(lm1,6)
          - 945*pow(m1sq,2)*pow(m2sq,4)*pow(m3sq,2)*pow(lm1,6) + 25*
         pow(m1sq,3)*m3sq*pow(lm1,4) - 35*pow(m1sq,3)*pow(m3sq,3)*pow(
         lm1,5) + 25*pow(m1sq,3)*m2sq*pow(lm1,4) + 490*pow(m1sq,3)*m2sq
         *pow(m3sq,2)*pow(lm1,5) - 945*pow(m1sq,3)*m2sq*pow(m3sq,4)*
         pow(lm1,6) + 490*pow(m1sq,3)*pow(m2sq,2)*m3sq*pow(lm1,5) + 
         1260*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,6) - 3465*
         pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,5)*pow(lm1,7) - 35*pow(
         m1sq,3)*pow(m2sq,3)*pow(lm1,5) + 1260*pow(m1sq,3)*pow(m2sq,3)*
         pow(m3sq,2)*pow(lm1,6) );

      h210p +=  + psiq(m1,m2,m3) * ( 3465*pow(m1sq,3)*pow(m2sq,3)*pow(
         m3sq,4)*pow(lm1,7) - 945*pow(m1sq,3)*pow(m2sq,4)*m3sq*pow(
         lm1,6) + 3465*pow(m1sq,3)*pow(m2sq,4)*pow(m3sq,3)*pow(lm1,7)
          - 3465*pow(m1sq,3)*pow(m2sq,5)*pow(m3sq,2)*pow(lm1,7) - 5*
         pow(m1sq,4)*pow(lm1,4) + 70*pow(m1sq,4)*pow(m3sq,2)*pow(lm1,5)
          + 1575*pow(m1sq,4)*m2sq*pow(m3sq,3)*pow(lm1,6) + 70*pow(
         m1sq,4)*pow(m2sq,2)*pow(lm1,5) + 630*pow(m1sq,4)*pow(m2sq,2)*
         pow(m3sq,2)*pow(lm1,6) + 3465*pow(m1sq,4)*pow(m2sq,2)*pow(
         m3sq,4)*pow(lm1,7) + 1575*pow(m1sq,4)*pow(m2sq,3)*m3sq*pow(
         lm1,6) - 6930*pow(m1sq,4)*pow(m2sq,3)*pow(m3sq,3)*pow(lm1,7)
          + 3465*pow(m1sq,4)*pow(m2sq,4)*pow(m3sq,2)*pow(lm1,7) - 35*
         pow(m1sq,5)*m3sq*pow(lm1,5) - 35*pow(m1sq,5)*m2sq*pow(lm1,5)
          - 315*pow(m1sq,5)*m2sq*pow(m3sq,2)*pow(lm1,6) - 315*pow(
         m1sq,5)*pow(m2sq,2)*m3sq*pow(lm1,6) + 3465*pow(m1sq,5)*pow(
         m2sq,2)*pow(m3sq,3)*pow(lm1,7) + 3465*pow(m1sq,5)*pow(m2sq,3)*
         pow(m3sq,2)*pow(lm1,7) );

      h210p +=  + psiq(m1,m2,m3) * (  - 315*pow(m1sq,6)*m2sq*m3sq*pow(
         lm1,6) - 3465*pow(m1sq,6)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,7) )
         ;

      h210p +=  + 1./6.*m3sq*pow(lm1,2) - 1./6.*m3sq*ln3*pow(lm1,2)
          + 1./6.*m3sq*ln2*pow(lm1,2) - 2./3.*pow(m3sq,3)*pow(lm1,3) + 
         1./3.*pow(m3sq,3)*ln3*pow(lm1,3) - 1./3.*pow(m3sq,3)*ln2*pow(
         lm1,3) + 1./6.*m2sq*pow(lm1,2) + 1./6.*m2sq*ln3*pow(lm1,2) - 1.
         /6.*m2sq*ln2*pow(lm1,2) + 2./3.*m2sq*pow(m3sq,2)*pow(lm1,3) - 
         2*m2sq*pow(m3sq,2)*ln3*pow(lm1,3) + 2*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,3) + 2*m2sq*pow(m3sq,4)*ln3*pow(lm1,4) - 2*m2sq*pow(
         m3sq,4)*ln2*pow(lm1,4) + 2./3.*pow(m2sq,2)*m3sq*pow(lm1,3) + 2
         *pow(m2sq,2)*m3sq*ln3*pow(lm1,3) - 2*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,3) - 6*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,4) + 6*pow(
         m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) - 2./3.*pow(m2sq,3)*pow(
         lm1,3) - 1./3.*pow(m2sq,3)*ln3*pow(lm1,3) + 1./3.*pow(m2sq,3)*
         ln2*pow(lm1,3) + 6*pow(m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 6*
         pow(m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) - 2*pow(m2sq,4)*m3sq*
         ln3*pow(lm1,4) + 2*pow(m2sq,4)*m3sq*ln2*pow(lm1,4) - 11./2.*
         m1sq*pow(lm1,2);
      h210p +=  - 11./6.*m1sq*ln3*pow(lm1,2) - 11./6.*m1sq*ln2*pow(
         lm1,2) + 79./3.*m1sq*pow(m3sq,2)*pow(lm1,3) + 27*m1sq*pow(
         m3sq,2)*ln3*pow(lm1,3) + 32./3.*m1sq*pow(m3sq,2)*ln2*pow(
         lm1,3) - 12*m1sq*pow(m3sq,4)*pow(lm1,4) - 60*m1sq*pow(m3sq,4)*
         ln3*pow(lm1,4) - 6*m1sq*pow(m3sq,4)*ln2*pow(lm1,4) + 32*m1sq*
         pow(m3sq,6)*ln3*pow(lm1,5) - 184./3.*m1sq*m2sq*m3sq*pow(lm1,3)
          - 134./3.*m1sq*m2sq*m3sq*ln3*pow(lm1,3) - 134./3.*m1sq*m2sq*
         m3sq*ln2*pow(lm1,3) + 100*m1sq*m2sq*pow(m3sq,3)*pow(lm1,4) + 
         264*m1sq*m2sq*pow(m3sq,3)*ln3*pow(lm1,4) + 126*m1sq*m2sq*pow(
         m3sq,3)*ln2*pow(lm1,4) - 256*m1sq*m2sq*pow(m3sq,5)*ln3*pow(
         lm1,5) - 48*m1sq*m2sq*pow(m3sq,5)*ln2*pow(lm1,5) + 79./3.*m1sq
         *pow(m2sq,2)*pow(lm1,3) + 32./3.*m1sq*pow(m2sq,2)*ln3*pow(
         lm1,3) + 27*m1sq*pow(m2sq,2)*ln2*pow(lm1,3) - 176*m1sq*pow(
         m2sq,2)*pow(m3sq,2)*pow(lm1,4) - 324*m1sq*pow(m2sq,2)*pow(
         m3sq,2)*ln3*pow(lm1,4) - 324*m1sq*pow(m2sq,2)*pow(m3sq,2)*ln2*
         pow(lm1,4);
      h210p +=  + 624*m1sq*pow(m2sq,2)*pow(m3sq,4)*ln3*pow(lm1,5) + 
         304*m1sq*pow(m2sq,2)*pow(m3sq,4)*ln2*pow(lm1,5) + 100*m1sq*
         pow(m2sq,3)*m3sq*pow(lm1,4) + 126*m1sq*pow(m2sq,3)*m3sq*ln3*
         pow(lm1,4) + 264*m1sq*pow(m2sq,3)*m3sq*ln2*pow(lm1,4) - 656*
         m1sq*pow(m2sq,3)*pow(m3sq,3)*ln3*pow(lm1,5) - 656*m1sq*pow(
         m2sq,3)*pow(m3sq,3)*ln2*pow(lm1,5) - 12*m1sq*pow(m2sq,4)*pow(
         lm1,4) - 6*m1sq*pow(m2sq,4)*ln3*pow(lm1,4) - 60*m1sq*pow(
         m2sq,4)*ln2*pow(lm1,4) + 304*m1sq*pow(m2sq,4)*pow(m3sq,2)*ln3*
         pow(lm1,5) + 624*m1sq*pow(m2sq,4)*pow(m3sq,2)*ln2*pow(lm1,5)
          - 48*m1sq*pow(m2sq,5)*m3sq*ln3*pow(lm1,5) - 256*m1sq*pow(
         m2sq,5)*m3sq*ln2*pow(lm1,5) + 32*m1sq*pow(m2sq,6)*ln2*pow(
         lm1,5) - 81./2.*pow(m1sq,2)*m3sq*pow(lm1,3) - 25*pow(m1sq,2)*
         m3sq*ln3*pow(lm1,3) - 40./3.*pow(m1sq,2)*m3sq*ln2*pow(lm1,3)
          - 19*pow(m1sq,2)*m3sq*ln1*pow(lm1,3) + 112*pow(m1sq,2)*pow(
         m3sq,3)*pow(lm1,4) + 175*pow(m1sq,2)*pow(m3sq,3)*ln3*pow(
         lm1,4);
      h210p +=  + 42*pow(m1sq,2)*pow(m3sq,3)*ln2*pow(lm1,4) + 69*pow(
         m1sq,2)*pow(m3sq,3)*ln1*pow(lm1,4) - 232*pow(m1sq,2)*pow(
         m3sq,5)*ln3*pow(lm1,5) - 81./2.*pow(m1sq,2)*m2sq*pow(lm1,3) - 
         40./3.*pow(m1sq,2)*m2sq*ln3*pow(lm1,3) - 28*pow(m1sq,2)*m2sq*
         ln2*pow(lm1,3) - 16*pow(m1sq,2)*m2sq*ln1*pow(lm1,3) - 167*pow(
         m1sq,2)*m2sq*pow(m3sq,2)*pow(lm1,4) - 207*pow(m1sq,2)*m2sq*
         pow(m3sq,2)*ln3*pow(lm1,4) - 95*pow(m1sq,2)*m2sq*pow(m3sq,2)*
         ln2*pow(lm1,4) - 132*pow(m1sq,2)*m2sq*pow(m3sq,2)*ln1*pow(
         lm1,4) + 560*pow(m1sq,2)*m2sq*pow(m3sq,4)*pow(lm1,5) + 1152*
         pow(m1sq,2)*m2sq*pow(m3sq,4)*ln3*pow(lm1,5) + 608*pow(m1sq,2)*
         m2sq*pow(m3sq,4)*ln2*pow(lm1,5) + 480*pow(m1sq,2)*m2sq*pow(
         m3sq,4)*ln1*pow(lm1,5) - 1200*pow(m1sq,2)*m2sq*pow(m3sq,6)*ln3
         *pow(lm1,6) - 167*pow(m1sq,2)*pow(m2sq,2)*m3sq*pow(lm1,4) - 98
         *pow(m1sq,2)*pow(m2sq,2)*m3sq*ln3*pow(lm1,4) - 228*pow(m1sq,2)
         *pow(m2sq,2)*m3sq*ln2*pow(lm1,4) - 108*pow(m1sq,2)*pow(m2sq,2)
         *m3sq*ln1*pow(lm1,4);
      h210p +=  - 560*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*pow(lm1,5)
          - 610*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,5) - 
         903*pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,5) - 495*
         pow(m1sq,2)*pow(m2sq,2)*pow(m3sq,3)*ln1*pow(lm1,5) + 2000*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,5)*ln3*pow(lm1,6) + 1600*pow(
         m1sq,2)*pow(m2sq,2)*pow(m3sq,5)*ln2*pow(lm1,6) + 112*pow(
         m1sq,2)*pow(m2sq,3)*pow(lm1,4) + 42*pow(m1sq,2)*pow(m2sq,3)*
         ln3*pow(lm1,4) + 178*pow(m1sq,2)*pow(m2sq,3)*ln2*pow(lm1,4) + 
         66*pow(m1sq,2)*pow(m2sq,3)*ln1*pow(lm1,4) - 560*pow(m1sq,2)*
         pow(m2sq,3)*pow(m3sq,2)*pow(lm1,5) - 918*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,5) - 640*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,5) - 450*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,2)*ln1*pow(lm1,5) + 1200*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,4)*ln3*pow(lm1,6) - 3600*pow(m1sq,2)*pow(
         m2sq,3)*pow(m3sq,4)*ln2*pow(lm1,6) + 560*pow(m1sq,2)*pow(
         m2sq,4)*m3sq*pow(lm1,5);
      h210p +=  + 608*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln3*pow(lm1,5) + 
         1167*pow(m1sq,2)*pow(m2sq,4)*m3sq*ln2*pow(lm1,5) + 465*pow(
         m1sq,2)*pow(m2sq,4)*m3sq*ln1*pow(lm1,5) - 3600*pow(m1sq,2)*
         pow(m2sq,4)*pow(m3sq,3)*ln3*pow(lm1,6) + 1200*pow(m1sq,2)*pow(
         m2sq,4)*pow(m3sq,3)*ln2*pow(lm1,6) - 232*pow(m1sq,2)*pow(
         m2sq,5)*ln2*pow(lm1,5) + 1600*pow(m1sq,2)*pow(m2sq,5)*pow(
         m3sq,2)*ln3*pow(lm1,6) + 2000*pow(m1sq,2)*pow(m2sq,5)*pow(
         m3sq,2)*ln2*pow(lm1,6) - 1200*pow(m1sq,2)*pow(m2sq,6)*m3sq*ln2
         *pow(lm1,6) + 5./6.*pow(m1sq,3)*pow(lm1,3) - 2*pow(m1sq,3)*ln3
         *pow(lm1,3) - 2*pow(m1sq,3)*ln2*pow(lm1,3) + 38./3.*pow(
         m1sq,3)*ln1*pow(lm1,3) - 163*pow(m1sq,3)*pow(m3sq,2)*pow(
         lm1,4) - 114*pow(m1sq,3)*pow(m3sq,2)*ln3*pow(lm1,4) - 55*pow(
         m1sq,3)*pow(m3sq,2)*ln2*pow(lm1,4) - 194*pow(m1sq,3)*pow(
         m3sq,2)*ln1*pow(lm1,4) + 416*pow(m1sq,3)*pow(m3sq,4)*ln3*pow(
         lm1,5) + 176*pow(m1sq,3)*pow(m3sq,4)*ln1*pow(lm1,5) - 62*pow(
         m1sq,3)*m2sq*m3sq*pow(lm1,4);
      h210p +=  - 154*pow(m1sq,3)*m2sq*m3sq*ln3*pow(lm1,4) - 165*pow(
         m1sq,3)*m2sq*m3sq*ln2*pow(lm1,4) - 53*pow(m1sq,3)*m2sq*m3sq*
         ln1*pow(lm1,4) - 565*pow(m1sq,3)*m2sq*pow(m3sq,3)*pow(lm1,5)
          - 846*pow(m1sq,3)*m2sq*pow(m3sq,3)*ln3*pow(lm1,5) - 359*pow(
         m1sq,3)*m2sq*pow(m3sq,3)*ln2*pow(lm1,5) - 659*pow(m1sq,3)*m2sq
         *pow(m3sq,3)*ln1*pow(lm1,5) + 2050*pow(m1sq,3)*m2sq*pow(
         m3sq,5)*ln3*pow(lm1,6) + 1550*pow(m1sq,3)*m2sq*pow(m3sq,5)*ln1
         *pow(lm1,6) - 163*pow(m1sq,3)*pow(m2sq,2)*pow(lm1,4) - 55*pow(
         m1sq,3)*pow(m2sq,2)*ln3*pow(lm1,4) - 125*pow(m1sq,3)*pow(
         m2sq,2)*ln2*pow(lm1,4) - 183*pow(m1sq,3)*pow(m2sq,2)*ln1*pow(
         lm1,4) + 1270*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*pow(lm1,5)
          + 694*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,5) + 
         624*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,5) + 266*
         pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,5) - 6180*pow(
         m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*ln3*pow(lm1,6) + 615*pow(
         m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*ln2*pow(lm1,6);
      h210p +=  + 765*pow(m1sq,3)*pow(m2sq,2)*pow(m3sq,4)*ln1*pow(
      lm1,6) - 565*pow(m1sq,3)*pow(m2sq,3)*m3sq*pow(lm1,5) - 364*pow(
         m1sq,3)*pow(m2sq,3)*m3sq*ln3*pow(lm1,5) - 916*pow(m1sq,3)*pow(
         m2sq,3)*m3sq*ln2*pow(lm1,5) - 584*pow(m1sq,3)*pow(m2sq,3)*m3sq
         *ln1*pow(lm1,5) + 3550*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,3)*ln3
         *pow(lm1,6) + 3515*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,3)*ln2*
         pow(lm1,6) - 4665*pow(m1sq,3)*pow(m2sq,3)*pow(m3sq,3)*ln1*pow(
         lm1,6) + 421*pow(m1sq,3)*pow(m2sq,4)*ln2*pow(lm1,5) + 171*pow(
         m1sq,3)*pow(m2sq,4)*ln1*pow(lm1,5) + 580*pow(m1sq,3)*pow(
         m2sq,4)*pow(m3sq,2)*ln3*pow(lm1,6) - 6215*pow(m1sq,3)*pow(
         m2sq,4)*pow(m3sq,2)*ln2*pow(lm1,6) + 835*pow(m1sq,3)*pow(
         m2sq,4)*pow(m3sq,2)*ln1*pow(lm1,6) + 2085*pow(m1sq,3)*pow(
         m2sq,5)*m3sq*ln2*pow(lm1,6) + 1515*pow(m1sq,3)*pow(m2sq,5)*
         m3sq*ln1*pow(lm1,6) + 38*pow(m1sq,4)*m3sq*pow(lm1,4) - 12*pow(
         m1sq,4)*m3sq*ln3*pow(lm1,4) + 8*pow(m1sq,4)*m3sq*ln2*pow(
         lm1,4);
      h210p +=  + 136*pow(m1sq,4)*m3sq*ln1*pow(lm1,4) - 176*pow(
         m1sq,4)*pow(m3sq,3)*ln3*pow(lm1,5) - 512*pow(m1sq,4)*pow(
         m3sq,3)*ln1*pow(lm1,5) + 38*pow(m1sq,4)*m2sq*pow(lm1,4) + 8*
         pow(m1sq,4)*m2sq*ln3*pow(lm1,4) - 4*pow(m1sq,4)*m2sq*ln2*pow(
         lm1,4) + 128*pow(m1sq,4)*m2sq*ln1*pow(lm1,4) - 550*pow(m1sq,4)
         *m2sq*pow(m3sq,2)*pow(lm1,5) - 582*pow(m1sq,4)*m2sq*pow(
         m3sq,2)*ln3*pow(lm1,5) - 831*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln2*
         pow(lm1,5) - 403*pow(m1sq,4)*m2sq*pow(m3sq,2)*ln1*pow(lm1,5)
          + 1050*pow(m1sq,4)*m2sq*pow(m3sq,4)*ln3*pow(lm1,6) - 3450*
         pow(m1sq,4)*m2sq*pow(m3sq,4)*ln1*pow(lm1,6) - 550*pow(m1sq,4)*
         pow(m2sq,2)*m3sq*pow(lm1,5) - 816*pow(m1sq,4)*pow(m2sq,2)*m3sq
         *ln3*pow(lm1,5) - 542*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln2*pow(
         lm1,5) - 458*pow(m1sq,4)*pow(m2sq,2)*m3sq*ln1*pow(lm1,5) + 
         3500*pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,3)*ln3*pow(lm1,6) - 4515
         *pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,3)*ln2*pow(lm1,6) + 3415*
         pow(m1sq,4)*pow(m2sq,2)*pow(m3sq,3)*ln1*pow(lm1,6);
      h210p +=  - 191*pow(m1sq,4)*pow(m2sq,3)*ln2*pow(lm1,5) - 497*
         pow(m1sq,4)*pow(m2sq,3)*ln1*pow(lm1,5) - 4410*pow(m1sq,4)*pow(
         m2sq,3)*pow(m3sq,2)*ln3*pow(lm1,6) + 3430*pow(m1sq,4)*pow(
         m2sq,3)*pow(m3sq,2)*ln2*pow(lm1,6) + 3380*pow(m1sq,4)*pow(
         m2sq,3)*pow(m3sq,2)*ln1*pow(lm1,6) + 945*pow(m1sq,4)*pow(
         m2sq,4)*m3sq*ln2*pow(lm1,6) - 3345*pow(m1sq,4)*pow(m2sq,4)*
         m3sq*ln1*pow(lm1,6) + 25*pow(m1sq,5)*pow(lm1,4) + 11*pow(
         m1sq,5)*ln3*pow(lm1,4) + 11*pow(m1sq,5)*ln2*pow(lm1,4) - 11*
         pow(m1sq,5)*ln1*pow(lm1,4) - 128*pow(m1sq,5)*pow(m3sq,2)*ln3*
         pow(lm1,5) + 480*pow(m1sq,5)*pow(m3sq,2)*ln1*pow(lm1,5) + 555*
         pow(m1sq,5)*m2sq*m3sq*pow(lm1,5) + 532*pow(m1sq,5)*m2sq*m3sq*
         ln3*pow(lm1,5) + 547*pow(m1sq,5)*m2sq*m3sq*ln2*pow(lm1,5) + 
         705*pow(m1sq,5)*m2sq*m3sq*ln1*pow(lm1,5) - 3450*pow(m1sq,5)*
         m2sq*pow(m3sq,3)*ln3*pow(lm1,6) + 1050*pow(m1sq,5)*m2sq*pow(
         m3sq,3)*ln1*pow(lm1,6) - 113*pow(m1sq,5)*pow(m2sq,2)*ln2*pow(
         lm1,5);
      h210p +=  + 465*pow(m1sq,5)*pow(m2sq,2)*ln1*pow(lm1,5) + 680*
         pow(m1sq,5)*pow(m2sq,2)*pow(m3sq,2)*ln3*pow(lm1,6) + 785*pow(
         m1sq,5)*pow(m2sq,2)*pow(m3sq,2)*ln2*pow(lm1,6) - 6265*pow(
         m1sq,5)*pow(m2sq,2)*pow(m3sq,2)*ln1*pow(lm1,6) - 3345*pow(
         m1sq,5)*pow(m2sq,3)*m3sq*ln2*pow(lm1,6) + 945*pow(m1sq,5)*pow(
         m2sq,3)*m3sq*ln1*pow(lm1,6) + 88*pow(m1sq,6)*m3sq*ln3*pow(
         lm1,5) - 128*pow(m1sq,6)*m3sq*ln1*pow(lm1,5) + 83*pow(m1sq,6)*
         m2sq*ln2*pow(lm1,5) - 123*pow(m1sq,6)*m2sq*ln1*pow(lm1,5) + 
         1550*pow(m1sq,6)*m2sq*pow(m3sq,2)*ln3*pow(lm1,6) + 2050*pow(
         m1sq,6)*m2sq*pow(m3sq,2)*ln1*pow(lm1,6) + 1515*pow(m1sq,6)*
         pow(m2sq,2)*m3sq*ln2*pow(lm1,6) + 2085*pow(m1sq,6)*pow(m2sq,2)
         *m3sq*ln1*pow(lm1,6) - 16*pow(m1sq,7)*ln1*pow(lm1,5) - 1200*
         pow(m1sq,7)*m2sq*m3sq*ln1*pow(lm1,6);
     break;
  default:
    std::cout <<"h210p (quenched) called with wrong iprop ="<<iprop<<'\n';
    h210p = 0.;
  }
 return h210p*pi162;
}

// outputs functions

double hh(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh (quenched case)"<<std::endl;
    return 0.;}
  return h0(iprop,m1sq,m2sq,m3sq,xmu2)+qsq*h0p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbar(iprop,m1sq,m2sq,m3sq,qsq,0);
}

double hh1(const int iprop, const double m1sq, const double m2sq,
	   const double m3sq, const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh1 (quenched case)"<<std::endl;
    return 0.;}
  return h10(iprop,m1sq,m2sq,m3sq,xmu2)+qsq*h10p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbar(iprop,m1sq,m2sq,m3sq,qsq,1);
}

double hh21(const int iprop, const double m1sq, const double m2sq,
	    const double m3sq,const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh21 (quenched case)"<<std::endl;
    return 0.;}
  return h210(iprop,m1sq,m2sq,m3sq,xmu2)+qsq*h210p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbar(iprop,m1sq,m2sq,m3sq,qsq,2);
}


double hhd(const int iprop, const double m1sq, const double m2sq,
	   const double m3sq, const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hhd (quenched case)"<<std::endl;
    return 0.;}
  return h0p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbard(iprop,m1sq,m2sq,m3sq,qsq,0);
}

double hh1d(const int iprop, const double m1sq, const double m2sq,
	    const double m3sq, const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh1d (quenched case)"<<std::endl;
    return 0.;}
  return h10p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbard(iprop,m1sq,m2sq,m3sq,qsq,1);
}

double hh21d(const int iprop, const double m1sq, const double m2sq,
	     const double m3sq,const double qsq, const double xmu2){
  if (qsq > pow(sqrt(m1sq)+sqrt(m2sq)+sqrt(m3sq),2)){
    std::cout <<"wrong qsq in hh21d (quenched case)"<<std::endl;
    return 0.;}
  return h210p(iprop,m1sq,m2sq,m3sq,xmu2)
    +hbard(iprop,m1sq,m2sq,m3sq,qsq,2);
}

