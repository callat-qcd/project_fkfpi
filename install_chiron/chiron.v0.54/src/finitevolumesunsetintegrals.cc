// finitevolumesunsetintegrals.cc is part of the CHIRON ChPT at two loops
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

// finitevolume sunset integrals needed for the masses and decay constants
// the main programs are in the euclidean
// the calling routines HHbV,HH1bV,HH21bV,HH27bV
//                      HHdV,HH1dV,HH21dV,HH27dV
// use Minkowski conventions
// or via the functions setprecisionfinitevolumesunset

#include <cmath>
#include <iostream>
#include <iomanip>
#include "jbnumlib.h"
#include "finitevolumesunsetintegrals.h"
// this are defined as an internal functions in finitevolumeoneloopintegrals.cc
// and are called here
double Abvt(const int, const int,const double,const double);
double Abvb(const int, const int,const double,const double);


#include <complex>
typedef std::complex<double> dcomplex; 

const double pi = M_PI;
const double pi2 = M_PI*M_PI;
const double pi16 = 1./(16.*pi*pi);
const double pi162 = pi16*pi16;
/////////////////sunset ////////////////////////////////////////////
//////// rs /////////////////////////////////////////////////////

namespace hhvrstspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype;
}

double hhvrstinternal(double xx[3]){
  using namespace hhvrstspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  const double alpha = 1.;
  // a possible stretching lambda small alpha >1 or large alpha < 1
  double lam = pow(xx[2]/(1.-xx[2]),alpha);
  double ovfac = (1.-x)*alpha/pow(1.-xx[2],2)*pow(lam,1.-1./alpha);
  double sig = x*y+y*z+z*x;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;
  double lam1 = xl*xl/(4.*sig*lam);
  double xb = x*lam1;
  double yb = y*lam1;
  double zb = z*lam1;
  double result = 0.;

  // make sure we get NaN if needed
  double theta0=nan(""),theta02x=nan(""),theta02y=nan(""),theta02z=nan("");
  if(ntype <= 26){
//    theta0 = pow(jbdtheta2d0(yb,xb,zb),3)
//      -pow(jbdtheta30(exp(-yb-zb)),3)-pow(jbdtheta30(exp(-xb-zb)),3)
//      -pow(jbdtheta30(exp(-xb-yb)),3)
//	   +2.;
    double t0m1 = jbdtheta2d0m1(yb,xb,zb);
    double t1m1 = jbdtheta30m1(exp(-yb-zb));
    double t2m1 = jbdtheta30m1(exp(-xb-zb));
    double t3m1 = jbdtheta30m1(exp(-xb-yb));
    theta0 =  3.*t0m1+3.*t0m1*t0m1+t0m1*t0m1*t0m1
             -3.*t1m1-3.*t1m1*t1m1-t1m1*t1m1*t1m1
             -3.*t2m1-3.*t2m1*t2m1-t2m1*t2m1*t2m1
             -3.*t3m1-3.*t3m1*t3m1-t3m1*t3m1*t3m1
          ;
  }
  else{
    double t02y = jbdtheta2d02(yb,xb,zb);
    double t02x = jbdtheta2d02(xb,yb,zb);
    double t02z = jbdtheta2d02(zb,xb,yb);
    double t0 = jbdtheta2d0(xb,yb,zb);
    double t32xy = jbdtheta32(exp(-xb-yb));
    double t32xz = jbdtheta32(exp(-xb-zb));
    double t32yz = jbdtheta32(exp(-yb-zb));
    double t30xy = jbdtheta30(exp(-xb-yb));
    double t30xz = jbdtheta30(exp(-xb-zb));
    double t30yz = jbdtheta30(exp(-yb-zb));
    theta02y = t02y*t0*t0-t32yz*t30yz*t30yz-t32xy*t30xy*t30xy;
    theta02x = t02x*t0*t0-t32xy*t30xy*t30xy-t32xz*t30xz*t30xz;
    theta02z = t02z*t0*t0-t32yz*t30yz*t30yz-t32xz*t30xz*t30xz;
  }
  switch(ntype){
  case 0:
    result = theta0;
    break;
  case 1:
    result = y*z/sig*theta0;
    break;
  case 2:
    result = x*z/sig*theta0;
    break;
  case 21:
    result = pow(y*z/sig,2)*theta0;
    break;
  case 22:
    result = (y+z)/(2.*lam*sig)*theta0;
    break;
  case 23:
    result = x*y*pow(z/sig,2)*theta0;
    break;
  case 24:
    result = (-z)/(2.*lam*sig)*theta0;
    break;
  case 25:
    result = pow(x*z/sig,2)*theta0;
    break;
  case 26:
    result = (x+z)/(2.*lam*sig)*theta0;
    break;
  case 27:
    result = pow(xl/(2.*sig*lam),2)*
      (-y*(y+z)*theta02y+y*z*theta02x-z*(y+z)*theta02z);
    break;
  case 28:
    result = 0.5*pow(xl/(2.*sig*lam),2)*
      ((2.*y*z-sig)*theta02y+(2.*x*z-sig)*theta02x+(2.*z*z+sig)*theta02z);
    break;
  case 29:
    result = pow(xl/(2.*sig*lam),2)*
      (x*z*theta02y-x*(x+z)*theta02x-z*(x+z)*theta02z);
    break;
  default:
    std::cout <<"hhvrst called with wrong ntype: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }

  switch(nprop){
  case 1:
    result /= (lam*lam);
    break;
  case 2:
    result *= x/lam;
    break;
  case 3:
    result *= y/lam;
    break;
  case 4:
    result *= z/lam;
    break;
  case 5:
    result *= x*y;
    break;
  case 6:
    result *= x*z;
    break;
  case 7:
    result *= y*z;
    break;
  case 8:
    result *= x*y*z*lam;
    break;
  default:
    std::cout <<"hhvrst called with wrong nprop: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }
  result = ovfac*pi162*result*exp(-lam*zrs)/pow(sig,2);
  return result;
}

double hhvrst(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrstspace::m1sq = m1sq;
  hhvrstspace::m2sq = m2sq;
  hhvrstspace::m3sq = m3sq;
  hhvrstspace::qsq = qsq;
  hhvrstspace::xl = xl;
  hhvrstspace::nprop = nprop;
  hhvrstspace::ntype = ntype;
  double a[3] = {0.,0.,0.};
  double b[3] = {1.,1.,1.};
  double releps = hhvrstspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad3(hhvrstinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"#precision not reached in hhvrst, nprop, ntype, ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

//////// rs derivative /////////////////////////////////////////////////////

namespace hhvrstdspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype;
}

double hhvrstdinternal(double xx[3]){
  using namespace hhvrstdspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  const double alpha = 1.;
  // a possible stretching lambda small alpha >1 or large alpha < 1
  double lam = pow(xx[2]/(1.-xx[2]),alpha);
  double ovfac = (1.-x)*alpha/pow(1.-xx[2],2)*pow(lam,1.-1./alpha);
  double sig = x*y+y*z+z*x;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;
  double lam1 = xl*xl/(4.*sig*lam);
  double xb = x*lam1;
  double yb = y*lam1;
  double zb = z*lam1;
  double result = 0.;

  // make sure we get NaN if needed
  double theta0=nan(""),theta02x=nan(""),theta02y=nan(""),theta02z=nan("");
  if(ntype <= 26){
//    theta0 = pow(jbdtheta2d0(yb,xb,zb),3)
//      -pow(jbdtheta30(exp(-yb-zb)),3)-pow(jbdtheta30(exp(-xb-zb)),3)
//      -pow(jbdtheta30(exp(-xb-yb)),3)
//	   +2.;
    double t0m1 = jbdtheta2d0m1(yb,xb,zb);
    double t1m1 = jbdtheta30m1(exp(-yb-zb));
    double t2m1 = jbdtheta30m1(exp(-xb-zb));
    double t3m1 = jbdtheta30m1(exp(-xb-yb));
    theta0 =  3.*t0m1+3.*t0m1*t0m1+t0m1*t0m1*t0m1
             -3.*t1m1-3.*t1m1*t1m1-t1m1*t1m1*t1m1
             -3.*t2m1-3.*t2m1*t2m1-t2m1*t2m1*t2m1
             -3.*t3m1-3.*t3m1*t3m1-t3m1*t3m1*t3m1
          ;
  }
  else{
    double t02y = jbdtheta2d02(yb,xb,zb);
    double t02x = jbdtheta2d02(xb,yb,zb);
    double t02z = jbdtheta2d02(zb,xb,yb);
    double t0 = jbdtheta2d0(xb,yb,zb);
    double t32xy = jbdtheta32(exp(-xb-yb));
    double t32xz = jbdtheta32(exp(-xb-zb));
    double t32yz = jbdtheta32(exp(-yb-zb));
    double t30xy = jbdtheta30(exp(-xb-yb));
    double t30xz = jbdtheta30(exp(-xb-zb));
    double t30yz = jbdtheta30(exp(-yb-zb));
    theta02y = t02y*t0*t0-t32yz*t30yz*t30yz-t32xy*t30xy*t30xy;
    theta02x = t02x*t0*t0-t32xy*t30xy*t30xy-t32xz*t30xz*t30xz;
    theta02z = t02z*t0*t0-t32yz*t30yz*t30yz-t32xz*t30xz*t30xz;
  }
  switch(ntype){
  case 0:
    result = theta0;
    break;
  case 1:
    result = y*z/sig*theta0;
    break;
  case 2:
    result = x*z/sig*theta0;
    break;
  case 21:
    result = pow(y*z/sig,2)*theta0;
    break;
  case 22:
    result = (y+z)/(2.*lam*sig)*theta0;
    break;
  case 23:
    result = x*y*pow(z/sig,2)*theta0;
    break;
  case 24:
    result = (-z)/(2.*lam*sig)*theta0;
    break;
  case 25:
    result = pow(x*z/sig,2)*theta0;
    break;
  case 26:
    result = (x+z)/(2.*lam*sig)*theta0;
    break;
  case 27:
    result = pow(xl/(2.*sig*lam),2)*
      (-y*(y+z)*theta02y+y*z*theta02x-z*(y+z)*theta02z);
    break;
  case 28:
    result = 0.5*pow(xl/(2.*sig*lam),2)*
      ((2.*y*z-sig)*theta02y+(2.*x*z-sig)*theta02x+(2.*z*z+sig)*theta02z);
    break;
  case 29:
    result = pow(xl/(2.*sig*lam),2)*
      (x*z*theta02y-x*(x+z)*theta02x-z*(x+z)*theta02z);
    break;
  default:
    std::cout <<"hhvrstd called with wrong ntype: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }

  switch(nprop){
  case 1:
    result /= (lam*lam);
    break;
  case 2:
    result *= x/lam;
    break;
  case 3:
    result *= y/lam;
    break;
  case 4:
    result *= z/lam;
    break;
  case 5:
    result *= x*y;
    break;
  case 6:
    result *= x*z;
    break;
  case 7:
    result *= y*z;
    break;
  case 8:
    result *= x*y*z*lam;
    break;
  default:
    std::cout <<"hhvrstd called with wrong nprop: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }
  result = ovfac*pi162*result*exp(-lam*zrs)/pow(sig,2)*(-lam*x*y*z/sig);
  return result;
}

double hhvrstd(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrstdspace::m1sq = m1sq;
  hhvrstdspace::m2sq = m2sq;
  hhvrstdspace::m3sq = m3sq;
  hhvrstdspace::qsq = qsq;
  hhvrstdspace::xl = xl;
  hhvrstdspace::nprop = nprop;
  hhvrstdspace::ntype = ntype;
  double a[3] = {0.,0.,0.};
  double b[3] = {1.,1.,1.};
  double releps = hhvrstdspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad3(hhvrstdinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"#precision not reached in hhvrstd, nprop, ntype, ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

////////////////////////// rH /////////////////////////////////////////
namespace hhvrHtspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype;
}

double hhvrHtinternal(double xx[3]){
  using namespace hhvrHtspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  const double alpha = 1.;
  double lam = pow(xx[2]/(1.-xx[2]),alpha);
  double ovfac = (1.-x)*alpha/pow(1.-xx[2],2)*pow(lam,1.-1./alpha);

  double sig = x*y+y*z+z*x;
  double rho = (y+z)/sig;
  double del = (y-z)/sig;
  double tau = y*z/sig;
  double lamh = lam*rho*xl*xl/4.;

  double theta30m1 = jbdtheta30m1(exp(-1./lam));
  double theta30 = 1.+theta30m1;
  double theta3 = theta30m1*(3.+theta30m1*(3.+theta30m1));//
  //double theta30 = jbdtheta30(exp(-1./lam));
  //double theta30m1 = theta30m1-1.;
  //double theta3 = pow(theta30,3)-1.;
  double theta32=nan(""),theta34=nan("");
  if((nprop <= 2) || (ntype >= 26) || 
     ((ntype == 26) && (nprop <= 6)))
    theta32 = jbdtheta32(exp(-1./lam));
  if((nprop <= 2) && (ntype >= 26)) theta34 = jbdtheta34(exp(-1./lam));

  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;

  switch(nprop){
  case 1:
    break;
  case 2:
    ovfac *= x*lamh;
    break;
  case 3:
  case 4:
    break;
  case 5:
  case 6:
    ovfac *= x*lamh;
    break;
  case 7:
    break;
  case 8:
    ovfac *= x*lamh;
    break;
  default:
    std::cout <<"hhvrHtinternal called with wrong nprop: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
    break;
  }

  double A=0.,result=0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = z*((A+2.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*theta3;
      break;
    case 4:
    case 6:
      result = z*theta3;
      break;
    case 7:
    case 8:
      result = y*z*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 1:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = z*((tau*A+2.*tau*del/lamh-x*del*rho/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*tau*theta3;
      break;
    case 4:
    case 6:
      result = z*tau*theta3;
      break;
    case 7:
    case 8:
      result = y*z*tau*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 2:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = x*z*z/2./sig*((A+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*y*z/sig*theta3;
      break;
    case 4:
    case 6:
      result = x*z*z/sig*theta3;
      break;
    case 7:
    case 8:
      result = x*y*z*z*lamh/sig*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 21:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = z*tau*((tau*A+(2.*tau*del-2.*x*del*rho)/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*tau*tau*theta3;
      break;
    case 4:
    case 6:
      result = z*tau*tau*theta3;
      break;
    case 7:
    case 8:
      result = y*z*tau*tau*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 22:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = z*rho/2./lamh*((A+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*rho/2./lamh*theta3;
      break;
    case 4:
    case 6:
      result = z*rho/2./lamh*theta3;
      break;
    case 7:
    case 8:
      result = y*z*rho/2.*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 23:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = x*z*z/2./sig*(
		 (tau*A+(3.*tau*del-x*del*rho)/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*y*z*tau/sig*theta3;
      break;
    case 4:
    case 6:
      result = x*z*z*tau/sig*theta3;
      break;
    case 7:
    case 8:
      result = x*y*z*z*tau*lamh/sig*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 24:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = -z*z/(4.*sig*lamh)*((A+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = -y*z/(2.*sig*lamh)*theta3;
      break;
    case 4:
    case 6:
      result = -z*z/(2.*sig*lamh)*theta3;
      break;
    case 7:
    case 8:
      result = -y*z*z/(2.*sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = x*x*z*z*z/(3.*sig*sig)*(
		 (A+4.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*x*y*z*z/(sig*sig)*theta3;
      break;
    case 4:
    case 6:
      result = x*x*z*z*z/(sig*sig)*theta3;
      break;
    case 7:
    case 8:
      result = x*x*y*z*z*z*lamh/(sig*sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result =
        (-z*m2sq/2.*(A+2.*del/lamh)
         +z*z/(6.*rho*sig*sig*lamh)*(5.*z+3.*y)*(m3sq-m2sq)
	 -3.*y*z*z/(2.*rho*sig*sig*lamh)*(A+del/lamh)
         -z*z*A/(12.*rho*sig)*(3.*A+4.*z/(lamh*sig)+4.*x*x*z*rho*qsq/sig))
          *theta3
         +(del/lamh/8.*pow(z/sig,2)*(2.*z+9.*y)
	   +pow(z/sig,2)/8.*((z+3.*y)*A-2.*z*(m3sq-m2sq))+3.*z*del*rho/8.*m2sq)
	*pow(xl/lamh,2)*theta32*theta30*theta30
         -z*z*del*rho/(64.*sig*sig)*(z+3.*y)*pow(xl/lamh,4)*
            (theta34*theta30*theta30+2.*theta32*theta32*theta30);
      break;
    case 3:
    case 5:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = 
	(1./4.*(z+tau/rho)*A+z*del/lamh+pow(z,3)/(2.*rho*sig*sig*lamh))*
	    theta3
        -3.*del*rho/16.*(z+tau/rho)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = 
	(1./4.*(z-tau/rho)*A+z*z/(2.*sig*lamh)-pow(z,3)/(2.*rho*sig*sig*lamh))*
	    theta3
        -3.*del*rho/16.*(z-tau/rho)*pow(xl/lamh,2)*theta32*theta30*theta30;
       break;
    case 7:
    case 8:
      result = tau/2./rho*(1+z*z/sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = -z*rho*rho/4.*(
           (A+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = -y*rho*rho/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = -z*rho*rho/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = -y*z*rho*rho*lamh/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = z*z*rho/8./sig*(
           (A+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = y*z*rho/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = z*z*rho/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = y*z*z*rho*lamh/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 29:
    switch(nprop){
    case 1:
    case 2:
      A = m3sq-m2sq+del*rho*x*x*qsq;
      result = -z*z*z/(12.*sig*sig)*(
           (A+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = -y*z*z/(4.*sig*sig)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = -z*z*z/(4.*sig*sig)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = -y*z*z*z*lamh/(4.*sig*sig)*pow(xl/lamh,2)*
                    theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  default:
    std::cout <<"hhvrHtinternal called with wrong ntype: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }
  result = ovfac*pi162*result*exp(-lamh*zrs)/pow(sig,2)/lam;

  return result;
}

double hhvrHt(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrHtspace::m1sq = m1sq;
  hhvrHtspace::m2sq = m2sq;
  hhvrHtspace::m3sq = m3sq;
  hhvrHtspace::qsq = qsq;
  hhvrHtspace::xl = xl;
  hhvrHtspace::nprop = nprop;
  hhvrHtspace::ntype = ntype;
  double a[3] = {0.,0.,0.};
  double b[3] = {1.,1.,1.};
  double releps = hhvrHtspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad3(hhvrHtinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"#precision not reached in hhvrHt, nprop ntype ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}


double hhvrGt(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  double result = 0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(1,0,m1sq,xl); 
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(2,0,m1sq,xl);
      break;
    default:
    break;
    }
    break;
  case 1:
    break;
  case 2:
    switch(nprop){
    case 1:
      result = -0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(1,0,m1sq,xl); 
      break;
    case 2:
      result = -0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 21:
    break;
  case 22:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(1,22,m1sq,xl);
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(2,22,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(1,23,m1sq,xl); 
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvt(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 23:
    break;
  case 24:
    switch(nprop){
    case 1:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(1,22,m1sq,xl); 
      break;
    case 2:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(2,22,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(1,23,m1sq,xl); 
      break;
    case 2:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvt(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvt(1,0,m1sq,xl);
      break;
    case 2:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvt(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
      result = log(m3sq/xmu2)*(m2sq/4.+m3sq/4.+qsq/12.)*Abvt(1,0,m1sq,xl)
	-1./3.*Abvt(1,22,m1sq,xl)
	+1./4.*log(m3sq/xmu2)*Abvt(1,23,m1sq,xl);
      result *= pi16;
      break;
    case 2:
      result = log(m3sq/xmu2)*(m2sq/4.+m3sq/4.+qsq/12.)*Abvt(2,0,m1sq,xl)
	-1./3.*Abvt(2,22,m1sq,xl)
	+1./4.*log(m3sq/xmu2)*Abvt(2,23,m1sq,xl);
      result *= pi16;
      break;
    case 3:
    case 4:
      result = -1./4.*pi16*(1.+log(m3sq/xmu2))*Abvt(1,0,m1sq,xl);
      break;
    case 5:
    case 6:
      result = -1./4.*pi16*(1.+log(m3sq/xmu2))*Abvt(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
 case 29:
    switch(nprop){
    case 1:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvt(1,23,m1sq,xl);
      break;
    case 2:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvt(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return result;
}

double hhvrFt(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  return hhvrGt(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2)
    +hhvrHt(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
}
////////////////////////// rH derivative ////////////////////////////
namespace hhvrHtdspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype;
}

double hhvrHtdinternal(double xx[3]){
  using namespace hhvrHtdspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  const double alpha = 1.;
  double lam = pow(xx[2]/(1.-xx[2]),alpha);
  double ovfac = (1.-x)*alpha/pow(1.-xx[2],2)*pow(lam,1.-1./alpha);

  double sig = x*y+y*z+z*x;
  double rho = (y+z)/sig;
  double del = (y-z)/sig;
  double tau = y*z/sig;
  double lamh = lam*rho*xl*xl/4.;

  double theta30m1 = jbdtheta30m1(exp(-1./lam));
  double theta30 = 1.+theta30m1;
  double theta3 = theta30m1*(3.+theta30m1*(3.+theta30m1));//
  //double theta30 = jbdtheta30(exp(-1./lam));
  //double theta30m1 = theta30m1-1.;
  //double theta3 = pow(theta30,3)-1.;
  double theta32=nan(""),theta34=nan("");
  if((nprop <= 2) || (ntype >= 26) || 
     ((ntype == 26) && (nprop <= 6)))
    theta32 = jbdtheta32(exp(-1./lam));
  if((nprop <= 2) && (ntype >= 26)) theta34 = jbdtheta34(exp(-1./lam));

  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;

  switch(nprop){
  case 1:
    break;
  case 2:
    ovfac *= x*lamh;
    break;
  case 3:
  case 4:
    break;
  case 5:
  case 6:
    ovfac *= x*lamh;
    break;
  case 7:
    break;
  case 8:
    ovfac *= x*lamh;
    break;
  default:
    std::cout <<"hhvrHtdinternal called with wrong nprop: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
    break;
  }

  double At=0.,Am=0.,As=0.,result=0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = z*((At+2.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*theta3;
      break;
    case 4:
    case 6:
      result = z*theta3;
      break;
    case 7:
    case 8:
      result = y*z*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 1:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = z*((tau*At+2.*tau*del/lamh-x*del*rho/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*tau*theta3;
      break;
    case 4:
    case 6:
      result = z*tau*theta3;
      break;
    case 7:
    case 8:
      result = y*z*tau*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 2:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = x*z*z/2./sig*((At+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*y*z/sig*theta3;
      break;
    case 4:
    case 6:
      result = x*z*z/sig*theta3;
      break;
    case 7:
    case 8:
      result = x*y*z*z*lamh/sig*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 21:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = z*tau*((tau*At+(2.*tau*del-2.*x*del*rho)/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*tau*tau*theta3;
      break;
    case 4:
    case 6:
      result = z*tau*tau*theta3;
      break;
    case 7:
    case 8:
      result = y*z*tau*tau*lamh*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 22:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = z*rho/2./lamh*((At+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = y*rho/2./lamh*theta3;
      break;
    case 4:
    case 6:
      result = z*rho/2./lamh*theta3;
      break;
    case 7:
    case 8:
      result = y*z*rho/2.*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 23:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = x*z*z/2./sig*(
		 (tau*At+(3.*tau*del-x*del*rho)/lamh)*theta3
                  -0.75*del*rho*tau*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*y*z*tau/sig*theta3;
      break;
    case 4:
    case 6:
      result = x*z*z*tau/sig*theta3;
      break;
    case 7:
    case 8:
      result = x*y*z*z*tau*lamh/sig*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 24:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = -z*z/(4.*sig*lamh)*((At+3.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = -y*z/(2.*sig*lamh)*theta3;
      break;
    case 4:
    case 6:
      result = -z*z/(2.*sig*lamh)*theta3;
      break;
    case 7:
    case 8:
      result = -y*z*z/(2.*sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = x*x*z*z*z/(3.*sig*sig)*(
		 (At+4.*del/lamh)*theta3
                  -0.75*del*rho*pow(xl/lamh,2)*theta32*theta30*theta30);
      break;
    case 3:
    case 5:
      result = x*x*y*z*z/(sig*sig)*theta3;
      break;
    case 4:
    case 6:
      result = x*x*z*z*z/(sig*sig)*theta3;
      break;
    case 7:
    case 8:
      result = x*x*y*z*z*z*lamh/(sig*sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
    case 2:
      Am = m3sq-m2sq+del*rho*x*x*qsq;
      At = Am-del*rho*x*sig/(y*z*lamh);
      As = Am-2.*del*rho*x*sig/(y*z*lamh);
      result =
        (-z*m2sq/2.*(At+2.*del/lamh)
         +z*z/(6.*rho*sig*sig*lamh)*(5.*z+3.*y)*(m3sq-m2sq)
	 -3.*y*z*z/(2.*rho*sig*sig*lamh)*(At+del/lamh)
         -z*z/(12.*rho*sig)*(3.*Am*As+At*(4.*z/(lamh*sig))
			     +4.*x*x*z*rho*(At*qsq-sig/(lamh*x*y*z)*Am)/sig)
	   )*theta3
         +(del/lamh/8.*pow(z/sig,2)*(2.*z+9.*y)
	   +pow(z/sig,2)/8.*((z+3.*y)*At-2.*z*(m3sq-m2sq))+3.*z*del*rho/8.*m2sq)
	*pow(xl/lamh,2)*theta32*theta30*theta30
         -z*z*del*rho/(64.*sig*sig)*(z+3.*y)*pow(xl/lamh,4)*
            (theta34*theta30*theta30+2.*theta32*theta32*theta30);
      break;
    case 3:
    case 5:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = 
	(1./4.*(z+tau/rho)*At+z*del/lamh+pow(z,3)/(2.*rho*sig*sig*lamh))*
	    theta3
        -3.*del*rho/16.*(z+tau/rho)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = 
	(1./4.*(z-tau/rho)*At+z*z/(2.*sig*lamh)-pow(z,3)/(2.*rho*sig*sig*lamh))*
	    theta3
        -3.*del*rho/16.*(z-tau/rho)*pow(xl/lamh,2)*theta32*theta30*theta30;
       break;
    case 7:
    case 8:
      result = tau/2./rho*(1+z*z/sig)*theta3;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = -z*rho*rho/4.*(
           (At+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = -y*rho*rho/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = -z*rho*rho/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = -y*z*rho*rho*lamh/4.*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = z*z*rho/8./sig*(
           (At+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = y*z*rho/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = z*z*rho/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = y*z*z*rho*lamh/4./sig*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  case 29:
    switch(nprop){
    case 1:
    case 2:
      At = m3sq-m2sq+del*rho*x*x*qsq-del*rho*x*sig/(y*z*lamh);
      result = -z*z*z/(12.*sig*sig)*(
           (At+4.*del/lamh)*pow(xl/lamh,2)*theta32*theta30*theta30
           -del*rho/4.*pow(xl/lamh,4)*
	         (theta34*theta30+2.*theta32*theta32)*theta30);
      break;
    case 3:
    case 5:
      result = -y*z*z/(4.*sig*sig)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 4:
    case 6:
      result = -z*z*z/(4.*sig*sig)*pow(xl/lamh,2)*theta32*theta30*theta30;
      break;
    case 7:
    case 8:
      result = -y*z*z*z*lamh/(4.*sig*sig)*pow(xl/lamh,2)*
                    theta32*theta30*theta30;
      break;
    default:
      std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
		<<nprop<<' '<<ntype<<std::endl;
      break;
    }
    break;
  default:
    std::cout <<"hhvrHtdinternal called with wrong ntype: nprop, ntype: "
	      <<nprop<<' '<<ntype<<std::endl;
  }
  result = ovfac*pi162*result*exp(-lamh*zrs)/pow(sig,2)/lam*(-lamh*x*y*z/sig);

  return result;
}

double hhvrHtd(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrHtdspace::m1sq = m1sq;
  hhvrHtdspace::m2sq = m2sq;
  hhvrHtdspace::m3sq = m3sq;
  hhvrHtdspace::qsq = qsq;
  hhvrHtdspace::xl = xl;
  hhvrHtdspace::nprop = nprop;
  hhvrHtdspace::ntype = ntype;
  double a[3] = {0.,0.,0.};
  double b[3] = {1.,1.,1.};
  double releps = hhvrHtdspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad3(hhvrHtdinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"#precision not reached in hhvrHtd, nprop ntype ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}


double hhvrGtd(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  double result = 0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
    break;
    }
    break;
  case 1:
    break;
  case 2:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 21:
    break;
  case 22:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 23:
    break;
  case 24:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
      result = log(m3sq/xmu2)*(1./12.)*Abvt(1,0,m1sq,xl);
      result *= pi16;
      break;
    case 2:
      result = log(m3sq/xmu2)*(1./12.)*Abvt(2,0,m1sq,xl);
      result *= pi16;
      break;
    case 3:
    case 4:
      result = 0.;
      break;
    case 5:
    case 6:
      result = 0.;
      break;
    default:
      break;
    }
    break;
 case 29:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return result;
}

double hhvrFtd(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  return hhvrGtd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2)
    +hhvrHtd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
}

//////////// overall sunset finite volume theta ///////////////////////////////
double hhvt(const int nprop, const int ntype,
	    const double m1sq, const double m2sq,const double m3sq,
	   const double qsq, const double xl, const double xmu2){
  int nprops[9]={0,1,3,2,4,5,7,6,8};
  int npropt[9]={0,1,4,3,2,7,6,5,8};
  double result;
  result = hhvrst(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
  double vr,vs,vt; 
  vr = hhvrFt(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2);
  // no simplifications made for equal masses here
  switch(ntype){
  case 0:
    vs = hhvrFt(nprops[nprop],0,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = hhvrFt(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 1:
    vs = hhvrFt(nprops[nprop],2,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFt(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFt(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 21:
    vs = hhvrFt(nprops[nprop],25,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFt(npropt[nprop],21,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFt(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFt(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFt(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 22:
    vs = hhvrFt(nprops[nprop],26,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFt(npropt[nprop],22,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFt(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 27:
    vs = hhvrFt(nprops[nprop],29,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFt(npropt[nprop],27,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFt(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 23:
    vs = hhvrFt(nprops[nprop],23,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFt(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFt(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFt(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 24:
    vs = hhvrFt(nprops[nprop],24,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFt(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFt(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 28:
    vs = hhvrFt(nprops[nprop],28,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFt(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFt(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  default:
    std::cout << "full sunset integral with type "
	      <<ntype<<" not implemented\n";
    return 0.;
    break;
  }
  result += vr+vs+vt;
  return result; 
}

double hhvtd(const int nprop, const int ntype,
	    const double m1sq, const double m2sq,const double m3sq,
	   const double qsq, const double xl, const double xmu2){
  int nprops[9]={0,1,3,2,4,5,7,6,8};
  int npropt[9]={0,1,4,3,2,7,6,5,8};
  double result;
  result = hhvrstd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
  double vr,vs,vt; 
  vr = hhvrFtd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2);
  // no simplifications made for equal masses here
  switch(ntype){
  case 0:
    vs = hhvrFtd(nprops[nprop],0,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = hhvrFtd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 1:
    vs = hhvrFtd(nprops[nprop],2,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFtd(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFtd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 21:
    vs = hhvrFtd(nprops[nprop],25,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFtd(npropt[nprop],21,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFtd(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFtd(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFtd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 22:
    vs = hhvrFtd(nprops[nprop],26,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFtd(npropt[nprop],22,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFtd(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 27:
    vs = hhvrFtd(nprops[nprop],29,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFtd(npropt[nprop],27,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFtd(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 23:
    vs = hhvrFtd(nprops[nprop],23,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFtd(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFtd(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFtd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 24:
    vs = hhvrFtd(nprops[nprop],24,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFtd(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFtd(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 28:
    vs = hhvrFtd(nprops[nprop],28,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFtd(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFtd(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  default:
    std::cout << "full sunset integral with type "
	      <<ntype<<" not implemented\n";
    return 0.;
    break;
  }
  result += vr+vs+vt;
  return result; 
}

/* the bessel combinations used are
KK_n(Y,Z) = 2 (Y/Z)^(n/2) K_n(2sqrt(YZ))
one propagator
Y = l_r^2/4 = iL^2/4
Z = m^2
2sqrt(YZ) = sqrt(i)mL = xml below
(Y/Z) = (i L^2/4 m^2) =  (xml^2/4 m^4) (note often negative powers
*/

// this contains the long lists of numbers
namespace besselonesumsunset{
#include "besselonesum.cc"
}
using namespace besselonesumsunset;
#include "besseltwosum.cc"
////////////////sunset/////////////////////////////////////////////
//////// rs /////////////////////////////////////////////////////

namespace hhvrsbspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps;
  int nprop,ntype,maxl2=40;
}

double hhvrsbinternal(double xx[2]){
  using namespace besselsummationdouble;
  using namespace hhvrsbspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  double sig = x*y+y*z+z*x;
  double sig2 = sig*sig;
  double xl2 = xl*xl;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;
  double result = 0.;
  double kk2,kk1,kk0,kkm1,kkm2,kkm3,fac;
  for (int i=maxll[maxl2]-1; i >=0;i--){
    double yrs = (y*mmall[i][1]+x*mmall[i][2]+z*mmall[i][3])*xl2/sig/4.;
    switch(ntype){
    case 0:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 1:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*z/sig*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*y*z/sig*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 2:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z/sig*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*x*z/sig*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 21:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += pow(y*z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*pow(y*z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 22:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += (1.-x)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*y*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  y*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 27:
      switch(nprop){
      case 1:
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 23:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*pow(z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 24:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += (-z)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*y*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*z*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  y*z*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z*(-z)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 28:
      switch(nprop){
      case 1:
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 25:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += pow(x*z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
 	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*pow(x*z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 26:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += (1.-y)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*y*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  x*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result +=  y*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 29:
      switch(nprop){
      case 1:
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsb nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    default:
	std::cout<<"invalid ntype in hhvrsb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
    }
  }
    return result*pi162*(1.-x)/sig2;
}

double hhvrsb(const int nprop, const int ntype,
		const double m1sq, const double m2sq,const double m3sq,
		const double qsq, const double xl){
  hhvrsbspace::m1sq = m1sq;
  hhvrsbspace::m2sq = m2sq;
  hhvrsbspace::m3sq = m3sq;
  hhvrsbspace::qsq = qsq;
  hhvrsbspace::xl = xl;
  hhvrsbspace::nprop = nprop;
  hhvrsbspace::ntype = ntype;
  double a[2] = {0.,0.};
  double b[2] = {1.,1.};
  double releps = hhvrsbspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad2(hhvrsbinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"precision not reached in hhvrsb, nprop, ntype, ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

//////// rs derivative ////////////////////////////////////////////

namespace hhvrsbdspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps;
  int nprop,ntype,maxl2=40;
}

double hhvrsbdinternal(double xx[2]){
  using namespace besselsummationdouble;
  using namespace hhvrsbdspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  double sig = x*y+y*z+z*x;
  double sig2 = sig*sig;
  double xl2 = xl*xl;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;
  double dzrs = x*y*z/sig;
  double result = 0.;
  double kk2,kk1,kk0,kkm1,kkm2,kkm3,fac;
  for (int i=maxll[maxl2]-1; i >=0;i--){
    double yrs = (y*mmall[i][1]+x*mmall[i][2]+z*mmall[i][3])*xl2/sig/4.;
    switch(ntype){
    case 0:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 1:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*z/sig*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*z/sig*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*y*z/sig*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*y*z/sig*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 2:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z/sig*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*x*z/sig*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*x*z/sig*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*x*z/sig*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 21:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += pow(y*z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*pow(y*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*pow(y*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*pow(y*z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 22:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += (1.-x)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*(1.-x)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*(1.-x)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 27:
      switch(nprop){
      case 1:
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (-y*(y+z)*mmall[i][1]+y*z*mmall[i][2]-z*(y+z)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 23:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*pow(z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*x*y*pow(z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 24:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += (-z)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*(-z)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*(-z)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*(-z)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 28:
      switch(nprop){
      case 1:
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(24.*sig2)*
	  ((2.*y*z-sig)*mmall[i][1]+(2.*x*z-sig)*mmall[i][2]
	   +(2.*z*z+sig)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 25:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += pow(x*z/sig,2)*double(mmall[i][0])*kkm1;
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*pow(x*z/sig,2)*double(mmall[i][0])*kk0;
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*y*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  x*z*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result +=  y*z*pow(x*z/sig,2)*double(mmall[i][0])*kk1;
	break;
      case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*y*z*pow(x*z/sig,2)*double(mmall[i][0])*kk2;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 26:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += (1.-y)/(2.*sig)*double(mmall[i][0])*kkm2;
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 3:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 4:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*(1.-y)/(2.*sig)*double(mmall[i][0])*kkm1;
	break;
      case 5:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*y*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 6:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  x*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result +=  y*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk0;
	break;
      case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*(1.-y)/(2.*sig)*double(mmall[i][0])*kk1;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 29:
      switch(nprop){
      case 1:
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += fac*double(mmall[i][0])*kkm3;
	break;
      case 2:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += x*fac*double(mmall[i][0])*kkm2;
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += y*fac*double(mmall[i][0])*kkm2;
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += z*fac*double(mmall[i][0])*kkm2;
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  x*y*fac*double(mmall[i][0])*kkm1;
	break;
      case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  x*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result +=  y*z*fac*double(mmall[i][0])*kkm1;
	break;
      case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	fac = xl2/(12.*sig2)*
	  (x*z*mmall[i][1]-x*(x+z)*mmall[i][2]-z*(x+z)*mmall[i][3]);
	result += x*y*z*fac*double(mmall[i][0])*kk0;
	break;
      default:
	std::cout<<"invalid nprop in hhvrsbd nprop ntype "
		  <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    default:
	std::cout<<"invalid ntype in hhvrsbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
    }
  }
    return result*pi162*(1.-x)/sig2;
}

double hhvrsbd(const int nprop, const int ntype,
		const double m1sq, const double m2sq,const double m3sq,
		const double qsq, const double xl){
  hhvrsbdspace::m1sq = m1sq;
  hhvrsbdspace::m2sq = m2sq;
  hhvrsbdspace::m3sq = m3sq;
  hhvrsbdspace::qsq = qsq;
  hhvrsbdspace::xl = xl;
  hhvrsbdspace::nprop = nprop;
  hhvrsbdspace::ntype = ntype;
  double a[2] = {0.,0.};
  double b[2] = {1.,1.};
  double releps = hhvrsbdspace::releps;
  double relerr;
  int ifail;
  double result =  jbdad2(hhvrsbdinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"precision not reached in hhvrsbd, nprop, ntype, ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

////////////////////////// rH /////////////////////////////////////////
namespace hhvrHbspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype,maxsumbessels=400;
}

double hhvrHbinternal(double xx[2]){
  using namespace hhvrHbspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  double sig = x*y+y*z+z*x;
  double sig2 = sig*sig;
  double xl2 = xl*xl;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;

  double rho = (y+z)/sig;
  double del = (y-z)/sig;
  double tau = y*z/sig;
  double A = m3sq-m2sq+del*rho*x*x*qsq;
  double C = del*rho*xl2/4.; 

  double result = 0.;
  double kk2,kk1,kk0,kkm1,kkm2,kkm3,kkm4;
  for (int i=maxsumbessels; i >=1;i--){
    double yrs = rho/4.*xl2*double(i);
    switch(ntype){
    case 0:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(A*kk0+2.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(A*kk1+2.*del*kk0-double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*y*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 1:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(tau*A*kk0+(2.*tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(tau*A*kk1+(2.*tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*tau*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*tau*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*tau*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*tau*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*tau*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*y*tau*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 2:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*z/(2.*sig)*
	  (A*kk0+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z/(2.*sig)*
	  (A*kk1+3.*del*kk0-double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*x*z/sig*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*x*z/sig*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*x*z/sig*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*x*z/sig*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*x*z/sig*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*y*x*z/sig*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 21:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*tau*(tau*A*kk0+2.*(tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*tau*(tau*A*kk1+2.*(tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*tau*tau*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*tau*tau*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*tau*tau*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*tau*tau*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*tau*tau*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*y*tau*tau*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 22:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += z*rho/2.*
	  (A*kkm1+3.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*rho/2.*
	  (A*kk0+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 3:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*rho/2.*kkm1*double(mm[i]);
	break;
      case 4:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*rho/2.*kkm1*double(mm[i]);
	break;
      case 5:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*rho/2.*kk0*double(mm[i]);
	break;
       case 6:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*rho/2.*kk0*double(mm[i]);
	break;
       case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*y*rho/2.*kk0*double(mm[i]);
	break;
       case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*y*rho/2.*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 27:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	kkm4 = 2.*pow(zrs/yrs,2)*jbdbesk4(2.*sqrt(yrs*zrs));
	result += z*(-rho*rho*xl2*double(i)/12.)*
	  (A*kkm2+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*(-rho*rho*xl2*double(i)/12.)*
	  (A*kkm1+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += y*(-rho*rho*xl2*double(i)/12.)*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(-rho*rho*xl2*double(i)/12.)*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*y*(-rho*rho*xl2*double(i)/12.)*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 23:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*z/(2.*sig)*(tau*A*kk0+(3.*tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z/(2.*sig)*(tau*A*kk1+(3.*tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*z*tau/sig*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*z*tau/sig*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*y*z*tau/sig*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z*tau/sig*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z*z*tau/sig*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z*tau/sig*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 24:
      switch(nprop){
      case 1:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += -z*z/(4.*sig)*
	  (A*kkm1+3.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 2:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += -x*z*z/(4.*sig)*
	  (A*kk0+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 3:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -y*z/(2.*sig)*kkm1*double(mm[i]);
	break;
      case 4:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -z*z/(2.*sig)*kkm1*double(mm[i]);
	break;
      case 5:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -x*y*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 6:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -x*z*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -y*z*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -x*y*z*z/(2.*sig)*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 28:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	kkm4 = 2.*pow(zrs/yrs,2)*jbdbesk4(2.*sqrt(yrs*zrs));
	result += z*(z*rho*xl2*double(i)/(24.*sig))*
	  (A*kkm2+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*(z*rho*xl2*double(i)/(24.*sig))*
	  (A*kkm1+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += y*(z*rho*xl2*double(i)/(12.*sig))*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(z*rho*xl2*double(i)/(12.*sig))*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*y*(z*rho*xl2*double(i)/(12.*sig))*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 25:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*z*z*z/(3.*sig2)*(A*kk0+4.*del*kkm1
		     -double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*x*z*z*z/(3.*sig2)*(A*kk1+4.*del*kk0
		       -double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z/sig2*kk0*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*x*z*z*z/sig2*kk0*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*x*y*z*z/sig2*kk1*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*x*z*z*z/sig2*kk1*double(mm[i]);
	break;
       case 7:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z*z/sig2*kk1*double(mm[i]);
	break;
       case 8:
	kk2 = 2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*x*y*z*z*z/sig2*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 26:
      switch(nprop){
      case 1:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	kkm4 = 2.*pow(zrs/yrs,2)*jbdbesk4(2.*sqrt(yrs*zrs));
	result += (-z*A/12.*(6*m2sq+3.*z/(rho*sig)*A+pow(2.*x*z/sig,2)*qsq)*kk0
		   +(-z*del*m2sq+z*z/(6.*rho*sig2)*(m3sq-m2sq)*(5.*z+3.*y)
		     -z*z*A/(6.*rho*sig2)*(2.*z+9.*y))*kkm1
		   +(z*z*xl2*double(i)/(24.*sig2)*((z+3.*y)*A-2.*z*(m3sq-m2sq))
		     -3.*y*z*z*del/(2.*rho*sig2)+z/2.*m2sq*C*double(i))
                      *kkm2
                  +z*z*xl2*double(i)*del/(24.*sig2)*(2.*z+9.*y)*kkm3
		   -z*z*C*xl2*double(i*i)/(48.*sig2)*(z+3.*y)*kkm4)
                    *double(mm[i]);
	break;
      case 2:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*
	      (-z*A/12.*(6*m2sq+3.*z/(rho*sig)*A+pow(2.*x*z/sig,2)*qsq)*kk1
		   +(-z*del*m2sq+z*z/(6.*rho*sig2)*(m3sq-m2sq)*(5.*z+3.*y)
		     -z*z*A/(6.*rho*sig2)*(2.*z+9.*y))*kk0
		   +(z*z*xl2*double(i)/(24.*sig2)*((z+3.*y)*A-2.*z*(m3sq-m2sq))
		     -3.*y*z*z*del/(2.*rho*sig2)+z/2.*m2sq*C*double(i))
                      *kkm1
                  +z*z*xl2*double(i)*del/(24.*sig2)*(2.*z+9.*y)*kkm2
		   -z*z*C*xl2*double(i*i)/(48.*sig2)*(z+3.*y)*kkm3)
                    *double(mm[i]);
	break;
      case 3:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += ((tau/rho+z)/4.*(A*kk0-2.*z/sig*kkm1-double(i)*C*kkm2)
                   +tau*kkm1)*double(mm[i]);
	break;
      case 4:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*z/(4.*rho*sig)*
	  (A*kk0+2.*y/sig*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 5:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*((tau/rho+z)/4.*(A*kk1-2.*z/sig*kk0-double(i)*C*kkm1)
                   +tau*kk0)*double(mm[i]);
	break;
       case 6:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*z/(4.*rho*sig)*
	  (A*kk1+2.*y/sig*kk0-double(i)*C*kkm1)*double(mm[i]);
	break;
       case 7:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += tau/(2.*rho)*(1.+z*z/sig)*kk0*double(mm[i]);
	break;
       case 8:
	kk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*tau/(2.*rho)*(1.+z*z/sig)*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 29:
      switch(nprop){
      case 1:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	kkm4 = 2.*pow(zrs/yrs,2)*jbdbesk4(2.*sqrt(yrs*zrs));
	result += z*(-z*z*xl2*double(i)/(36.*sig2))*
	  (A*kkm2+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm3 = 2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*(-z*z*xl2*double(i)/(36.*sig2))*
	  (A*kkm1+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += y*(-z*z*xl2*double(i)/(12.*sig2))*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(-z*z*xl2*double(i)/(12.*sig2))*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 8:
	kk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*y*(-z*z*xl2*double(i)/(12.*sig2))*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    default:
	std::cout<<"invalid ntype in hhvrHb nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
    }
  }
  return result*pi162*(1.-x)/sig2;
}

double hhvrHb(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrHbspace::m1sq = m1sq;
  hhvrHbspace::m2sq = m2sq;
  hhvrHbspace::m3sq = m3sq;
  hhvrHbspace::qsq = qsq;
  hhvrHbspace::xl = xl;
  hhvrHbspace::nprop = nprop;
  hhvrHbspace::ntype = ntype;
  double a[2] = {0.,0.};
  double b[2] = {1.,1.};
  double releps = hhvrHbspace::releps;
  double relerr;
  int ifail;
 double result =  jbdad2(hhvrHbinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"precision not reached in hhvrs, nprop ntype ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

double hhvrGb(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  double result = 0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(1,0,m1sq,xl); 
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(2,0,m1sq,xl);
      break;
    default:
    break;
    }
    break;
  case 1:
    break;
  case 2:
    switch(nprop){
    case 1:
      result = -0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(1,0,m1sq,xl); 
      break;
    case 2:
      result = -0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 21:
    break;
  case 22:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(1,22,m1sq,xl);
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(2,22,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(1,23,m1sq,xl); 
      break;
    case 2:
      result = -pi16*(1.+log(m3sq/xmu2))*Abvb(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 23:
    break;
  case 24:
    switch(nprop){
    case 1:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(1,22,m1sq,xl); 
      break;
    case 2:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(2,22,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(1,23,m1sq,xl); 
      break;
    case 2:
      result = 0.5*pi16*(1.+log(m3sq/xmu2))*Abvb(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvb(1,0,m1sq,xl);
      break;
    case 2:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvb(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
      result = log(m3sq/xmu2)*(m2sq/4.+m3sq/4.+qsq/12.)*Abvb(1,0,m1sq,xl)
	-1./3.*Abvb(1,22,m1sq,xl)
	+1./4.*log(m3sq/xmu2)*Abvb(1,23,m1sq,xl);
      result *= pi16;
      break;
    case 2:
      result = log(m3sq/xmu2)*(m2sq/4.+m3sq/4.+qsq/12.)*Abvb(2,0,m1sq,xl)
	-1./3.*Abvb(2,22,m1sq,xl)
	+1./4.*log(m3sq/xmu2)*Abvb(2,23,m1sq,xl);
      result *= pi16;
      break;
    case 3:
    case 4:
      result = -1./4.*pi16*(1.+log(m3sq/xmu2))*Abvb(1,0,m1sq,xl);
      break;
    case 5:
    case 6:
      result = -1./4.*pi16*(1.+log(m3sq/xmu2))*Abvb(2,0,m1sq,xl);
      break;
    default:
      break;
    }
    break;
 case 29:
    switch(nprop){
    case 1:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvb(1,23,m1sq,xl);
      break;
    case 2:
      result = -1./3.*pi16*(1.+log(m3sq/xmu2))*Abvb(2,23,m1sq,xl);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return result;
}

double hhvrFb(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  return hhvrGb(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2)
    +hhvrHb(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
}

////////////////////////// rH derivative //////////////////////////////////
namespace hhvrHbdspace{
  double m1sq,m2sq,m3sq,qsq,xl,releps=1e-4;
  int nprop,ntype,maxsumbessels=400;
}

double hhvrHbdinternal(double xx[2]){
  using namespace hhvrHbdspace;
  double x = xx[0];
  double y = xx[1]*(1.-x);
  double z = 1.-x-y;
  double sig = x*y+y*z+z*x;
  double sig2 = sig*sig;
  double xl2 = xl*xl;
  double zrs = x*m1sq+y*m2sq+z*m3sq+x*y*z/sig*qsq;
  double dzrs = x*y*z/sig;

  double rho = (y+z)/sig;
  double del = (y-z)/sig;
  double tau = y*z/sig;
  double At = m3sq-m2sq+del*rho*x*x*qsq;
  double dA = del*rho*x*x;
  double C = del*rho*xl2/4.; 

  double result = 0.;
  double kk2,kk1,kk0,kkm1,kkm2,kkm3,kkm4;
  double akk1,akk0,akkm1,akkm2;
  for (int i=maxsumbessels; i >=1;i--){
    double yrs = rho/4.*xl2*double(i);
    switch(ntype){
    case 0:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*(dA*akk0+At*kk0+2.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(dA*akk1+At*kk1+2.*del*kk0-double(i)*C*kkm1)
	  *double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*y*kk1*double(mm[i]);
	break;
       case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*y*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 1:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*(tau*(dA*akk0+At*kk0)+(2.*tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(tau*(dA*akk1+At*kk1)+(2.*tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*tau*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*tau*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*tau*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*tau*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*y*tau*kk1*double(mm[i]);
	break;
       case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*y*tau*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 2:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*z/(2.*sig)*
	  (dA*akk0+At*kk0+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z/(2.*sig)*
	  (dA*akk1+At*kk1+3.*del*kk0-double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*x*z/sig*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*x*z/sig*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*x*z/sig*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*x*z/sig*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*y*x*z/sig*kk1*double(mm[i]);
	break;
       case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*y*x*z/sig*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 21:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*tau*(tau*(dA*akk0+At*kk0)+2.*(tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*tau*(tau*(dA*akk1+At*kk1)+2.*(tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*tau*tau*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*tau*tau*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*tau*tau*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*tau*tau*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*y*tau*tau*kk1*double(mm[i]);
	break;
       case 8:
	kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*z*y*tau*tau*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 22:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*rho/2.*
	  ((dA*akkm1+At*kkm1)+3.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*rho/2.*
	  ((dA*akk0+At*kk0)+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 3:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += y*rho/2.*kkm1*double(mm[i]);
	break;
      case 4:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*rho/2.*kkm1*double(mm[i]);
	break;
      case 5:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*rho/2.*kk0*double(mm[i]);
	break;
       case 6:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*rho/2.*kk0*double(mm[i]);
	break;
       case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*y*rho/2.*kk0*double(mm[i]);
	break;
       case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*z*y*rho/2.*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 27:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm4 = -dzrs*2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	akkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(-rho*rho*xl2*double(i)/12.)*
	  ((dA*akkm2+At*kkm2)+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(-rho*rho*xl2*double(i)/12.)*
	  ((dA*akkm1+At*kkm1)+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(-rho*rho*xl2*double(i)/12.)*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(-rho*rho*xl2*double(i)/12.)*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*y*(-rho*rho*xl2*double(i)/12.)*kkm1*double(mm[i]);
	break;
       case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*y*(-rho*rho*xl2*double(i)/12.)*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 23:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*z/(2.*sig)*(tau*(dA*akk0+At*kk0)+(3.*tau-rho*x)*del*kkm1
		     -tau*double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z/(2.*sig)*(tau*(dA*akk1+At*kk1)+(3.*tau-rho*x)*del*kk0
		       -tau*double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*y*z*tau/sig*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*z*tau/sig*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*y*z*tau/sig*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*z*z*tau/sig*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*y*z*z*tau/sig*kk1*double(mm[i]);
	break;
       case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z*tau/sig*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 24:
      switch(nprop){
      case 1:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -z*z/(4.*sig)*
	  ((dA*akkm1+At*kkm1)+3.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 2:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -x*z*z/(4.*sig)*
	  ((dA*akk0+At*kk0)+3.*del*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 3:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -y*z/(2.*sig)*kkm1*double(mm[i]);
	break;
      case 4:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += -z*z/(2.*sig)*kkm1*double(mm[i]);
	break;
      case 5:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -x*y*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 6:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -x*z*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += -y*z*z/(2.*sig)*kk0*double(mm[i]);
	break;
       case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += -x*y*z*z/(2.*sig)*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 28:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm4 = -dzrs*2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	akkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(z*rho*xl2*double(i)/(24.*sig))*
	  ((dA*akkm2+At*kkm2)+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(z*rho*xl2*double(i)/(24.*sig))*
	  ((dA*akkm1+At*kkm1)+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(z*rho*xl2*double(i)/(12.*sig))*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(z*rho*xl2*double(i)/(12.*sig))*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*y*(z*rho*xl2*double(i)/(12.*sig))*kkm1*double(mm[i]);
	break;
       case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*y*(z*rho*xl2*double(i)/(12.*sig))*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 25:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*x*z*z*z/(3.*sig2)*((dA*akk0+At*kk0)+4.*del*kkm1
		     -double(i)*C*kkm2)*double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*x*z*z*z/(3.*sig2)*((dA*akk1+At*kk1)+4.*del*kk0
		       -double(i)*C*kkm1)*double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z/sig2*kk0*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*x*z*z*z/sig2*kk0*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*x*y*z*z/sig2*kk1*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*x*z*z*z/sig2*kk1*double(mm[i]);
	break;
       case 7:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*x*y*z*z*z/sig2*kk1*double(mm[i]);
	break;
       case 8:
        kk2  = -dzrs*2.*pow(yrs/zrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	result += x*x*x*y*z*z*z/sig2*kk2*double(mm[i]);
	break;
      default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 26:
      switch(nprop){
      case 1:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm4 = -dzrs*2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += (
		   -z/12.*(6*m2sq)*(dA*akk0+At*kk0)
		   -z/12.*(pow(2.*x*z/sig,2))*(dA*qsq*akk0+At*akk0+At*qsq*kk0)
		   -z/12.*(3.*z/(rho*sig))*(2.*At*dA*akk0+At*At*kk0)
		   +(-z*del*m2sq+z*z/(6.*rho*sig2)*(m3sq-m2sq)*(5.*z+3.*y)
		     )*kkm1
		   +(-z*z/(6.*rho*sig2)*(2.*z+9.*y))*(dA*akkm1+At*kkm1)
		   +(z*z*xl2*double(i)/(24.*sig2)*(-2.*z*(m3sq-m2sq))
		     -3.*y*z*z*del/(2.*rho*sig2)+z/2.*m2sq*C*double(i))
                      *kkm2
		   +(z*z*xl2*double(i)/(24.*sig2)*((z+3.*y)))
		     *(dA*akkm2+At*kkm2)
                  +z*z*xl2*double(i)*del/(24.*sig2)*(2.*z+9.*y)*kkm3
		   -z*z*C*xl2*double(i*i)/(48.*sig2)*(z+3.*y)*kkm4)
                    *double(mm[i]);
	break;
      case 2:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*(
		     -z/12.*(6*m2sq)*(dA*akk1+At*kk1)
		     -z/12.*(pow(2.*x*z/sig,2))*(dA*qsq*akk1+At*akk1+At*qsq*kk1)
		     -z/12.*(3.*z/(rho*sig))*(2.*At*dA*akk1+At*At*kk1)
		   +(-z*del*m2sq+z*z/(6.*rho*sig2)*(m3sq-m2sq)*(5.*z+3.*y)
		     )*kk0
		     +(-z*z/(6.*rho*sig2)*(2.*z+9.*y))*(dA*akk0+At*kk0)
		   +(z*z*xl2*double(i)/(24.*sig2)*(-2.*z*(m3sq-m2sq))
		     -3.*y*z*z*del/(2.*rho*sig2)+z/2.*m2sq*C*double(i))
                      *kkm1
		     +(z*z*xl2*double(i)/(24.*sig2)*((z+3.*y)))
		     *(dA*akkm1+At*kkm1)
                  +z*z*xl2*double(i)*del/(24.*sig2)*(2.*z+9.*y)*kkm2
		   -z*z*C*xl2*double(i*i)/(48.*sig2)*(z+3.*y)*kkm3)
                    *double(mm[i]);
	break;
      case 3:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += ((tau/rho+z)/4.*((dA*akk0+At*kk0)
				   -2.*z/sig*kkm1-double(i)*C*kkm2)
                   +tau*kkm1)*double(mm[i]);
	break;
      case 4:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	akk0 = 2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*z/(4.*rho*sig)*
	  ((dA*akk0+At*kk0)+2.*y/sig*kkm1-double(i)*C*kkm2)*double(mm[i]);
	break;
      case 5:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*((tau/rho+z)/4.*((dA*akk1+At*kk1)
				     -2.*z/sig*kk0-double(i)*C*kkm1)
		     +tau*kk0)*double(mm[i]);
	break;
       case 6:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	akk1 = 2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*z/(4.*rho*sig)*
	  ((dA*akk1+At*kk1)+2.*y/sig*kk0-double(i)*C*kkm1)*double(mm[i]);
	break;
       case 7:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += tau/(2.*rho)*(1.+z*z/sig)*kk0*double(mm[i]);
	break;
       case 8:
	kk1  = -dzrs*2.*yrs/zrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += x*tau/(2.*rho)*(1.+z*z/sig)*kk1*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    case 29:
      switch(nprop){
      case 1:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	kkm4 = -dzrs*2.*pow(zrs/yrs,1.5)*jbdbesk3(2.*sqrt(yrs*zrs));
	akkm2 = 2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	result += z*(-z*z*xl2*double(i)/(36.*sig2))*
	  ((dA*akkm2+At*kkm2)+4.*del*kkm3-double(i)*C*kkm4)*double(mm[i]);
	break;
      case 2:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	kkm3 = -dzrs*2.*zrs/yrs*jbdbesk2(2.*sqrt(yrs*zrs));
	akkm1 = 2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*(-z*z*xl2*double(i)/(36.*sig2))*
	  ((dA*akkm1+At*kkm1)+4.*del*kkm2-double(i)*C*kkm3)*double(mm[i]);
	break;
      case 3:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += y*(-z*z*xl2*double(i)/(12.*sig2))*kkm2*double(mm[i]);
	break;
      case 4:
	kkm2 = -dzrs*2.*pow(zrs/yrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += z*(-z*z*xl2*double(i)/(12.*sig2))*kkm2*double(mm[i]);
	break;
      case 5:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*y*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 6:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += x*z*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 7:
	kkm1 = -dzrs*2.*jbdbesk0(2.*sqrt(yrs*zrs));
	result += z*y*(-z*z*xl2*double(i)/(12.*sig2))*kkm1*double(mm[i]);
	break;
       case 8:
	kk0  = -dzrs*2.*pow(yrs/zrs,0.5)*jbdbesk1(2.*sqrt(yrs*zrs));
	result += x*z*y*(-z*z*xl2*double(i)/(12.*sig2))*kk0*double(mm[i]);
	break;
     default:
	std::cout<<"invalid nprop in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
	break;
      }
      break;
    default:
	std::cout<<"invalid ntype in hhvrHbd nprop ntype "
		 <<nprop<<' '<<ntype<<std::endl;
    }
  }
  return result*pi162*(1.-x)/sig2;
}

double hhvrHbd(const int nprop, const int ntype,
	      const double m1sq, const double m2sq,const double m3sq,
	      const double qsq, const double xl){
  hhvrHbdspace::m1sq = m1sq;
  hhvrHbdspace::m2sq = m2sq;
  hhvrHbdspace::m3sq = m3sq;
  hhvrHbdspace::qsq = qsq;
  hhvrHbdspace::xl = xl;
  hhvrHbdspace::nprop = nprop;
  hhvrHbdspace::ntype = ntype;
  double a[2] = {0.,0.};
  double b[2] = {1.,1.};
  double releps = hhvrHbdspace::releps;
  double relerr;
  int ifail;
 double result =  jbdad2(hhvrHbdinternal,a,b,releps,relerr,ifail);
  if(ifail != 0){
    std::cout <<"precision not reached in hhvrHd, nprop ntype ifail relerr: "
	      <<nprop<<' '<<ntype<<' '<<ifail<<' '<<relerr<<std::endl;
  }
  return result;
}

double hhvrGbd(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  double result = 0.;
  switch(ntype){
  case 0:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
    break;
    }
    break;
  case 1:
    break;
  case 2:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 21:
    break;
  case 22:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 27:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 23:
    break;
  case 24:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 28:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 25:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  case 26:
    switch(nprop){
    case 1:
      result = log(m3sq/xmu2)*(1./12.)*Abvb(1,0,m1sq,xl);
      result *= pi16;
      break;
    case 2:
      result = log(m3sq/xmu2)*(1./12.)*Abvb(2,0,m1sq,xl);
      result *= pi16;
      break;
    case 3:
    case 4:
      result = 0.;
      break;
    case 5:
    case 6:
      result = 0.;
      break;
    default:
      break;
    }
    break;
 case 29:
    switch(nprop){
    case 1:
      result = 0.;
      break;
    case 2:
      result = 0.;
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
  return result;
}

double hhvrFbd(const int nprop, const int ntype,
             const double m1sq, const double m2sq,const double m3sq,
	     const double qsq, const double xl, const double xmu2){
  return hhvrGbd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2)
    +hhvrHbd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
}


//////////// overall sunset finite volume Bessel ///////////////////////////////
double hhvb(const int nprop, const int ntype,
	    const double m1sq, const double m2sq,const double m3sq,
	   const double qsq, const double xl, const double xmu2){
  int nprops[9]={0,1,3,2,4,5,7,6,8};
  int npropt[9]={0,1,4,3,2,7,6,5,8};
  double result;
  result = hhvrsb(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
  double vr,vs,vt; 
  vr = hhvrFb(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2);
  // no simplifications made for equal masses here
  switch(ntype){
  case 0:
    vs = hhvrFb(nprops[nprop],0,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = hhvrFb(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 1:
    vs = hhvrFb(nprops[nprop],2,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFb(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFb(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 21:
    vs = hhvrFb(nprops[nprop],25,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFb(npropt[nprop],21,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFb(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFb(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFb(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 22:
    vs = hhvrFb(nprops[nprop],26,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFb(npropt[nprop],22,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFb(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 27:
    vs = hhvrFb(nprops[nprop],29,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFb(npropt[nprop],27,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFb(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 23:
    vs = hhvrFb(nprops[nprop],23,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFb(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFb(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFb(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 24:
    vs = hhvrFb(nprops[nprop],24,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFb(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFb(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 28:
    vs = hhvrFb(nprops[nprop],28,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFb(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFb(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  default:
    std::cout << "full sunset integral with type "
	      <<ntype<<" not implemented\n";
    return 0.;
    break;
  }
  result += vr+vs+vt;
  return result; 
}

//////////// overall sunset finite volume Bessel derivative//////////////
double hhvbd(const int nprop, const int ntype,
	    const double m1sq, const double m2sq,const double m3sq,
	   const double qsq, const double xl, const double xmu2){
  int nprops[9]={0,1,3,2,4,5,7,6,8};
  int npropt[9]={0,1,4,3,2,7,6,5,8};
  double result;
  result = hhvrsbd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl);
  double vr,vs,vt; 
  vr = hhvrFbd(nprop,ntype,m1sq,m2sq,m3sq,qsq,xl,xmu2);
  // no simplifications made for equal masses here
  switch(ntype){
  case 0:
    vs = hhvrFbd(nprops[nprop],0,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = hhvrFbd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 1:
    vs = hhvrFbd(nprops[nprop],2,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFbd(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFbd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 21:
    vs = hhvrFbd(nprops[nprop],25,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFbd(npropt[nprop],21,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFbd(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFbd(npropt[nprop],1,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      -2.*hhvrFbd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],0,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 22:
    vs = hhvrFbd(nprops[nprop],26,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFbd(npropt[nprop],22,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFbd(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 27:
    vs = hhvrFbd(nprops[nprop],29,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = +hhvrFbd(npropt[nprop],27,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2)
      +2.*hhvrFbd(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 23:
    vs = hhvrFbd(nprops[nprop],23,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFbd(npropt[nprop],23,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFbd(npropt[nprop],25,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         +hhvrFbd(npropt[nprop],2,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 24:
    vs = hhvrFbd(nprops[nprop],24,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFbd(npropt[nprop],24,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFbd(npropt[nprop],26,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  case 28:
    vs = hhvrFbd(nprops[nprop],28,m2sq,m1sq,m3sq,qsq,xl,xmu2);
    vt = -hhvrFbd(npropt[nprop],28,m3sq,m2sq,m1sq,qsq,xl,xmu2)
         -hhvrFbd(npropt[nprop],29,m3sq,m2sq,m1sq,qsq,xl,xmu2);
    break;
  default:
    std::cout << "full sunset integral with type "
	      <<ntype<<" not implemented\n";
    return 0.;
    break;
  }
  result += vr+vs+vt;
  return result; 
}


////////////// finite volume abbreviations in Minkowski convention //////////
// normal case
// theta function
double hhVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvt(1,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1Vt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvt(1,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21Vt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvt(1,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22Vt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvt(1,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27Vt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvt(1,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

double hhdVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvtd(1,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1dVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvtd(1,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21dVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvtd(1,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22dVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvtd(1,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27dVt(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvtd(1,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

//bessel function
double hhVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvb(1,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1Vb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvb(1,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21Vb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvb(1,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22Vb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvb(1,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27Vb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvb(1,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

double hhdVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvbd(1,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1dVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvbd(1,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21dVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return hhvbd(1,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22dVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvbd(1,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27dVb(const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  return -hhvbd(1,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
// with other powers of propagators //////////////////////////////////////
// theta function
double hhVt(const int nprop,
	    const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvt(nprop,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1Vt(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvt(nprop,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvt(nprop,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvt(nprop,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvt(nprop,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

double hhdVt(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvtd(nprop,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1dVt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvtd(nprop,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvtd(nprop,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvtd(nprop,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvtd(nprop,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

//bessel function
double hhVb(const int nprop,
	    const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvb(nprop,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1Vb(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvb(nprop,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvb(nprop,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvb(nprop,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvb(nprop,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

double hhdVb(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvbd(nprop,0,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh1dVb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvbd(nprop,1,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh21dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,1.,-1.,-1.,-1.,1.,1.,1.,-1};
  return extrasign[nprop]*hhvbd(nprop,21,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh22dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvbd(nprop,22,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}
double hh27dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2){
  const double extrasign[9]={0.,-1.,1.,1.,1.,-1.,-1.,-1.,1};
  return extrasign[nprop]*hhvbd(nprop,27,m1sq,m2sq,m3sq,-qsq,xl,mu2);
}

////////////// accuracy setting //////////////////////////////////////
void setprecisionfinitevolumesunsett(const double sunsetracc,
				     const double sunsetrsacc,
				     const bool out){
  hhvrHtspace::releps = sunsetracc;
  hhvrstspace::releps = sunsetrsacc;
  hhvrHtdspace::releps = sunsetracc;
  hhvrstdspace::releps = sunsetrsacc;
  if(out){
    std::cout << "#accuracies sunsetr sunsetrs: "
	      << sunsetracc<<' '<<sunsetrsacc<<std::endl;}
}

void setprecisionfinitevolumesunsetb(const int maxonesum,
				     const int maxtwosum,
				     const double sunsetracc,
				     const double sunsetrsacc,
				     const bool out){
  hhvrHbspace::maxsumbessels = maxonesum;// see besselonesum.cc for maximum
  hhvrHbdspace::maxsumbessels = maxonesum;// see besselonesum.cc for maximum
  hhvrsbspace::maxl2 = maxtwosum;// see besseltwosum.cc for maximum
  hhvrsbdspace::maxl2 = maxtwosum;// see besseltwosum.cc for maximum
  hhvrHbspace::releps = sunsetracc;
  hhvrsbspace::releps = sunsetrsacc;
  hhvrHbdspace::releps = sunsetracc;
  hhvrsbdspace::releps = sunsetrsacc;
  if (out){
    std::cout<<"#accuracies maxonesum maxtwosum sunsetracc sunsetrsacc : "
	     <<maxonesum<<' '<<maxtwosum<<' '
	     <<sunsetracc<<' '<<sunsetrsacc<<std::endl;}
}
