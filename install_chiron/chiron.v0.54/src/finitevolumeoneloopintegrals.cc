// finitevolumeoneloopintegrals.cc is part of the CHIRON ChPT at two loops
// program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// uses the methods as derived in
//  J.~Bijnens, E.~Boström and T.~A.~Lähde,
//  %``Two-loop Sunset Integrals at Finite Volume,''
//  JHEP {\bf 1401} (2014) 019
//  [arXiv:1311.3531 [hep-lat]].
// and references therein

// finitevolume integrals needed for the masses and decay constants
// the main programs are in the euclidean
// the calling routines AbVt,BbVt,A22bVt,A23bVt
// use Minkowski conventions

#include <cmath>
#include <iostream>
#include <iomanip>
#include "jbnumlib.h"
#include "finitevolumeoneloopintegrals.h"

#include <complex>
typedef std::complex<double> dcomplex;

#ifndef DINTEGRAL
#define DINTEGRAL jbdgauss
#endif

const double pi = M_PI;
const double pi2 = M_PI*M_PI;
const double pi16 = 1./(16.*pi*pi);
const double pi162 = pi16*pi16;

// this contains the long lists of numbers needed for the Bessel method
namespace besseloneloop{
#include "besselonesum.cc"
}
using namespace besseloneloop;
//////////// finite volume integrals Euclidean conventions //////////////
///////// tadpoles theta functions ///////////////////////////////////////
namespace Abvtspace{
  double releps=1e-10;
  double m2l2;
  int nprop,ntype;
}

double Abvtinternal(double x){
  using namespace Abvtspace;
  const double alpha = 0.5;
  double lt = pow(x/(1.-x),alpha);
  double ovfac = alpha*pow(lt,-1.-1./alpha)/pow(1.-x,2);
  switch(nprop){
  case 1:
    break;
  case 2:
    ovfac *= lt;
    break;
  case 3:
    ovfac *= lt*lt;
    break;
  case 4:
    ovfac *= pow(lt,3);
    break;
  default:
    std::cout << "Bad nprop in Abvt, nprop = "<<nprop<<std::endl;
    break;
  }
  switch(ntype){
  case 0:
    ovfac *= (pow(jbdtheta30(exp(-1./lt)),3)-1.);
    break;
 case 22:
    ovfac *= 0.5/lt*(pow(jbdtheta30(exp(-1./lt)),3)-1.);
    break;
 case 23:
   ovfac *= -1./(lt*lt)*jbdtheta32(exp(-1./lt))*pow(jbdtheta30(exp(-1./lt)),2);
    break;
  default:
    std::cout << "Bad ntype in Abvt, ntype = "<<ntype<<std::endl;
    break;
  }
  return ovfac*exp(-lt*m2l2/4.);
}

double Abvt(const int nprop,const int ntype,const double msq,const double xl){
  Abvtspace::nprop = nprop;
  Abvtspace::ntype = ntype;
  Abvtspace::m2l2 = msq*xl*xl;
  double ovfac=nan("");
  switch(nprop){
  case 1:
    ovfac = pi16*4./(xl*xl);
    break;
  case 2:
    ovfac = pi16;
    break;
  case 3:
    ovfac = pi16*(xl*xl)/4./2.;
    break;
  case 4:
    ovfac = pi16*pow(xl/2.,4)/6.;
    break;
  default:
    std::cout << "Bad nprop in Abvt, nprop = "<<nprop<<std::endl;
    break;
  }
  switch(ntype){
  case 0:
    break;
  case 22:
  case 23:
    ovfac *= 4./(xl*xl);
    break;
  default:
    std::cout << "Bad ntype in Abvt, ntype = "<<ntype<<std::endl;
    break;
  }
  double result = DINTEGRAL(Abvtinternal,0.,1.,Abvtspace::releps);
  return ovfac*result;
}
// // defined for later use
namespace Bbvtspace{
  double m1sq,m2sq,qsq,xl,releps=1e-9;
  int nprop,ntype;
}

/* the bessel combinations used are
KK_n(Y,Z) = 2 (Y/Z)^(n/2) K_n(2sqrt(YZ))
one propagator
Y = l_r^2/4 = iL^2/4
Z = m^2
2sqrt(YZ) = sqrt(i)mL = xml below
(Y/Z) = (i L^2/4 m^2) =  (xml^2/4 m^4) (note often negative powers
*/


//////////// finite volume integrals Euclidean conventions //////////////
////////////// with Bessel functions functions//////////////////////////////
///////// tadpoles //////////////////////////////////////////////////////

double Abvb(const int nprop, const int ntype, double msq, double xl){
  double xm = sqrt(msq);
  double xl2 = xl*xl;
  double result = 0.;
  for(int i=maxsumbessel; i>=1;i--){
    double xml = sqrt(double(i))*xm*xl;
    double kk2,kk1,kk0,kkm1,kkm2,kkm3;
    switch(ntype){
    case 0:
      switch(nprop){
      case 1:
	kkm1 = 4.*msq/xml*jbdbesk1(xml);
	result += double(mm[i])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(xml);
	result += double(mm[i])*kk0;
	break;
      case 3:
	kk1 = xml/msq*jbdbesk1(xml);
	result += double(mm[i])*kk1*0.5;
	break;
      case 4:
	kk2 = 0.5*pow(xml/msq,2)*jbdbesk2(xml);
	result += double(mm[i])*kk2/6.;
	break;
     default:
	std::cout << "Bad nprop in Abvb, nprop ntype = "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 22:
      switch(nprop){
      case 1:
	kkm2 = 8.*pow(msq/xml,2)*jbdbesk2(xml);
	result += 0.5*double(mm[i])*kkm2;
	break;
      case 2:
	kkm1 = 4.*msq/xml*jbdbesk1(xml);
	result += 0.5*double(mm[i])*kkm1;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(xml);
	result += 0.25*double(mm[i])*kk0;
	break;
      case 4:
	kk1 = xml/msq*jbdbesk1(xml);
	result += 0.5/6.*double(mm[i])*kk1;
	break;
      default:
	std::cout << "Bad nprop in Abvb, nprop ntype = "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 23:
      switch(nprop){
      case 1:
	kkm3 = 16.*pow(msq/xml,3)*jbdbesk3(xml);
	result += -1./12.*xl2*double(i*mm[i])*kkm3;
	break;
      case 2:
	kkm2 = 8.*pow(msq/xml,2)*jbdbesk2(xml);
	result += -1./12.*xl2*double(i*mm[i])*kkm2;
	break;
      case 3:
	kkm1 = 4.*msq/xml*jbdbesk1(xml);
	result += -0.5/12.*xl2*double(i*mm[i])*kkm1;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(xml);
	result += -1./12./6.*xl2*double(i*mm[i])*kk0;
	break;
      default:
	std::cout << "Bad nprop in Abvb, nprop ntype = "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    default:
      std::cout << "Bad ntype in Abvb, ntype = "<<ntype<<std::endl;
      break;
    }
  }
  return pi16*result;
}
// defined for later use
namespace Bbvbspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps = 1e-8;
  int nprop,ntype;
}


////////////// finite volume abbreviations in Minkowski convention //////////
// theta functions
double AbVt(const double msq, const double xl){
  return -Abvt(1,0,msq,xl);
}
double A22bVt(const double msq, const double xl){
  return Abvt(1,22,msq,xl);
}
double A23bVt(const double msq, const double xl){
  return Abvt(1,23,msq,xl);
}
double BbVt(const double msq, const double xl){
  return  Abvt(2,0,msq,xl);
}
double B22bVt(const double msq, const double xl){
  return -Abvt(2,22,msq,xl);
}
double B23bVt(const double msq, const double xl){
  return -Abvt(2,23,msq,xl);
}
double CbVt(const double msq, const double xl){
  return  -Abvt(3,0,msq,xl);
}
double C22bVt(const double msq, const double xl){
  return Abvt(3,22,msq,xl);
}
double C23bVt(const double msq, const double xl){
  return Abvt(3,23,msq,xl);
}
double DbVt(const double msq, const double xl){
  return  Abvt(4,0,msq,xl);
}
double D22bVt(const double msq, const double xl){
  return -Abvt(4,22,msq,xl);
}
double D23bVt(const double msq, const double xl){
  return -Abvt(4,23,msq,xl);
}

double BbVt(const double m1sq, const double m2sq, const double xl){
  if (m1sq == m2sq) return BbVt(m1sq,xl);
  return (AbVt(m1sq,xl)-AbVt(m2sq,xl))/(m1sq-m2sq);
}
// Bessel functions ////////////////////////////////
double AbVb(const double msq, const double xl){
  return -Abvb(1,0,msq,xl);
}
double A22bVb(const double msq, const double xl){
  return Abvb(1,22,msq,xl);
}
double A23bVb(const double msq, const double xl){
  return Abvb(1,23,msq,xl);
}
double BbVb(const double msq, const double xl){
  return  Abvb(2,0,msq,xl);
}
double B22bVb(const double msq, const double xl){
  return -Abvb(2,22,msq,xl);
}
double B23bVb(const double msq, const double xl){
  return -Abvb(2,23,msq,xl);
}
double CbVb(const double msq, const double xl){
  return  -Abvb(3,0,msq,xl);
}
double C22bVb(const double msq, const double xl){
  return  Abvb(3,22,msq,xl);
}
double C23bVb(const double msq, const double xl){
  return  Abvb(3,23,msq,xl);
}
double DbVb(const double msq, const double xl){
  return  Abvb(4,0,msq,xl);
}
double D22bVb(const double msq, const double xl){
  return -Abvb(4,22,msq,xl);
}
double D23bVb(const double msq, const double xl){
  return -Abvb(4,23,msq,xl);
}

double BbVb(const double m1sq, const double m2sq, const double xl){
  if (m1sq == m2sq) return BbVb(m1sq,xl);
  return (AbVb(m1sq,xl)-AbVb(m2sq,xl))/(m1sq-m2sq);
}
////////////// accuracy setting //////////////////////////////////////
// theta functions
void setprecisionfinitevolumeoneloopt(const double Abacc,const double Bbacc,
bool out){
  Abvtspace::releps = Abacc;
  Bbvtspace::releps = Bbacc;
  if(out){
  std::cout << "#accuracies AbV BbV : "
	    << Abacc <<' '<< Bbacc<<' '<<std::endl;}
}

// Bessel functions
void setprecisionfinitevolumeoneloopb(const int maxonesum, const double Bbacc,bool out){
  maxsumbessel = maxonesum;// see besselonesum.cc for maximum
  Bbvbspace::releps = Bbacc;// accuracy of integral for Bb functions
  if(out){
    std::cout<<"#accuracies maxonesum Bbacc : "
	     <<maxonesum<<' '<<Bbacc<<std::endl;}
}

