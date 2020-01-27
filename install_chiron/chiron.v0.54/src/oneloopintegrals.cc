// oneloopintegrals.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// tadpoles Ab,Bb,Cb
// one-loop integrals Bb: analytical, equal mass, numerical

//Note: commented out functions are probbaly correct but not fully tested


// allow for different numerical integration
// present possible (CHIRON >= v0.51: jbwgauss,jbwgauss2,jbwquad15,jbwquad21)
// the ifndef allows you to specify it while compiling
#ifndef WINTEGRAL
#define WINTEGRAL jbwgauss
#endif

#include <iostream>
#include <cmath>
#include <complex>

#include "oneloopintegrals.h"
#include "jbnumlib.h"

typedef std::complex<double> dcomplex;
double precisiononeloopintegrals = 1e-10;
const double pi = M_PI; //4.*atan(1.0);
const double pi2 = M_PI*M_PI; //4.*atan(1.0);
const double pi16 = 1./(M_PI*M_PI*16.); //1./(16.*pow(pi,2));


// Abar one propagator integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
double Ab(const double msq, const double mu2){
    return  -pi16*msq*log(msq/mu2);
}

double Bb(const double msq, const double mu2){
    return  -pi16*(1.+log(msq/mu2));
}

double Cb(const double msq, const double mu2=pow(0.77,2)){
    return  -pi16*0.5/msq;
}

double Ab(const int n, const double msq, const double mu2){
  switch(n){
  case 1:
    return Ab(msq,mu2);
  case 2:
    return Bb(msq,mu2);
  case 3:
    return Cb(msq,mu2);
  default:
    std::cout << "Ab(n,msq,mu2) called with unimplemented n=" << n <<" \n";
    return 0.;
  }
}
// epsilon versions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

double Abeps(const double msq, const double mu2){
  return pi16*msq*(pi2/12.+0.5+0.5*pow(log(msq/mu2),2));
}

double Bbeps(const double msq, const double mu2){
  double lm = log(msq/mu2);
  return pi16*(pi2/12.+0.5+lm+0.5*lm*lm);
}

double Cbeps(const double msq, const double mu2){
  return pi16*0.5/msq*(1.+log(msq/mu2));
}


// two-propagator integrals xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// used for the numerical cases below
struct bintdata{
  double m1sq,m2sq,qsq,mu2;
};
struct cintdata{
  double m1sq,m2sq,m3sq,qsq,mu2;
};

// set precision for the numerical integrations
void setprecisiononeloopintegrals(const double eps){
  precisiononeloopintegrals = eps;
}
double getprecisiononeloopintegrals(void){
  return precisiononeloopintegrals;
}


// Bbar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
dcomplex Bb(const double m1sq,const double m2sq,const double qsq,
		  const double mu2){
  double d,f1,f0;
  dcomplex b0iana;
  if ((qsq == 0.)){
    if (m1sq == m2sq) return -pi16*(log(m1sq/mu2)+1.);
    return -pi16/(m1sq-m2sq)*(m1sq*log(m1sq/mu2)-m2sq*log(m2sq/mu2));
  }
  double a=qsq/mu2;
  double b=(m2sq-m1sq-qsq)/mu2;
  double c=m1sq/mu2;
  double dsq = b*b/(4.*a*a)-c/a;
  double y1=1.0+b/(2.*a);
  double y0=b/(2.*a);
  if (dsq >= 0.){
    d=sqrt(dsq);
    f1=y1*log(fabs(y1*y1-dsq))-2.0*y1+d*log(fabs((y1+d)/(y1-d)));
    f0=y0*log(fabs(y0*y0-dsq))-2.0*y0+d*log(fabs((y0+d)/(y0-d)));
  }
  else{
    d=sqrt(-dsq);
    f1=y1*log(fabs(y1*y1-dsq))-2.0*y1+2.0*d*atan(y1/d);
    f0=y0*log(fabs(y0*y0-dsq))-2.0*y0+2.0*d*atan(y0/d);
  }
  b0iana = pi16*(-1.0-log(fabs(a))-f1+f0);
  if (qsq >= pow(sqrt(m1sq)+sqrt(m2sq),2)){
    b0iana = b0iana
      +dcomplex(0.0,pi16*pi/qsq*sqrt(pow((qsq-m1sq-m2sq),2)-4.0*m1sq*m2sq));
  }
  return b0iana;
}

// simplified version equal masses
dcomplex Bb(const double msq, const double qsq, const double mu2){
  if ((qsq == 0.)){
    return  -pi16*(1.+log(msq/mu2));
  }
  double b0iana = pi16*(1.-log(msq/mu2));

  double qsqt = qsq/msq;
  double rhop = 1.-4./qsqt;
  dcomplex ff;
  if (rhop>=0.){
    double rho = sqrt(rhop);
    ff = -rho*log(fabs((rho+1.)/(rho-1.)));
    if (qsqt > 4.) ff += dcomplex(0.,pi*rho);
  }
  else{
    double rho = sqrt(-rhop);
    ff = -2.*rho*atan(1./rho);
  }
  return b0iana+pi16*ff;
}

// version with numerical integration over x left in
static bintdata b0intdat;

dcomplex Bbnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b0intdat.m2sq+(1.-x)*b0intdat.m1sq-x*(1.-x)*b0intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b0intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b0intdat.mu2)-dcomplex(0.,pi);}
  return -lmtsq;
}

dcomplex Bbnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b0intdat.m1sq = m1sq;
  b0intdat.m2sq = m2sq;
  b0intdat.qsq  = qsq;
  b0intdat.mu2  = mu2;
  return pi16*(-1.+WINTEGRAL(Bbnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			 precisiononeloopintegrals/pi16));
}

// B1bar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

dcomplex B1b(const double m1sq,const double m2sq,const double qsq,
		  const double mu2){
  if(m1sq == m2sq){
    return 0.5*Bb(m1sq,qsq,mu2);
  }
  if (qsq == 0.){
    return pi16*(-0.5-0.5*log(m2sq/mu2)
		 -0.5*pow(m1sq/(m2sq-m1sq),2)*log(m1sq/m2sq)
		 -(0.75*m1sq-0.25*m2sq)/(m2sq-m1sq));
  }
  return ( -0.5/qsq*(Ab(m1sq,mu2)-Ab(m2sq,mu2)+
		     +(m2sq-m1sq-qsq)*Bb(m1sq,m2sq,qsq,mu2)));
}

// version with numerical integration over x left in
static bintdata b1intdat;
dcomplex B1bnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b1intdat.m2sq+(1.-x)*b1intdat.m1sq-x*(1.-x)*b1intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b1intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b1intdat.mu2)-dcomplex(0.,pi);}
  return -x*lmtsq;
}

dcomplex B1bnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b1intdat.m1sq = m1sq;
  b1intdat.m2sq = m2sq;
  b1intdat.qsq  = qsq;
  b1intdat.mu2  = mu2;
  return pi16*(-0.5+WINTEGRAL(B1bnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			 precisiononeloopintegrals/pi16));
}


// B21bar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

dcomplex B21b(const double m1sq,const double m2sq,const double qsq,
	      const double mu2){
  return (1./qsq*(Ab(m2sq,mu2)+m1sq*Bb(m1sq,m2sq,qsq,mu2)
		  -4.*B22b(m1sq,m2sq,qsq,mu2)
		  +2.*pi16*(m1sq/4.+m2sq/4.-qsq/12.)));
}


// version with numerical integration over x left in
static bintdata b21intdat;
dcomplex B21bnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b21intdat.m2sq+(1.-x)*b21intdat.m1sq-x*(1.-x)*b21intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b21intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b21intdat.mu2)-dcomplex(0.,pi);}
  return -pow(x,2)*lmtsq;
}

dcomplex B21bnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b21intdat.m1sq = m1sq;
  b21intdat.m2sq = m2sq;
  b21intdat.qsq  = qsq;
  b21intdat.mu2  = mu2;
  return pi16*(-1./3.+WINTEGRAL(B21bnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			     precisiononeloopintegrals/pi16));
}

// B22bar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

dcomplex B22b(const double m1sq,const double m2sq,const double qsq,
		   const double mu2){
  if (m1sq == m2sq){
    return B22b(m1sq,qsq,mu2);
  }
  // small qsq
  if((fabs(qsq/m1sq) < 1.e-5) || (fabs(qsq/m2sq) < 1.e-5)){
    double sq3 = m1sq -m2sq;
    return dcomplex(((m1sq + m2sq)/8.
		     -(m1sq*m1sq*log(m1sq/mu2)-m2sq*m2sq*log(m2sq/mu2))/(4.*sq3)
		     +qsq/(4.*pow(sq3,3))*((pow(m1sq,3)-pow(m2sq,3))/18.
					   +m1sq*m2sq*sq3/2. 
					   +m1sq*m1sq*(m1sq/3.-m2sq)*log(m1sq/mu2)
					   -m2sq*m2sq*(m2sq/3.-m1sq)*log(m2sq/mu2)
					   ))*pi16);
  }
  return dcomplex((m1sq-m2sq+qsq)/(12.*qsq)*Ab(m1sq,mu2)
		  -(m1sq-m2sq-qsq)/(12.*qsq)*Ab(m2sq,mu2)
		  +pi16*(3.*(m1sq+m2sq)-qsq)/18.)
    -dcomplex((pow(m1sq-m2sq-qsq,2)- 4.*qsq*m2sq)/(12.*qsq))*Bb(m1sq,m2sq,qsq,mu2);
}

//
dcomplex B22b(const double m1sq,const double qsq,const double mu2){
	if(fabs(qsq/mu2) < 1.e-7){
	  return dcomplex(pi16*(-m1sq/2.*log(m1sq/mu2) 
				+ qsq/12.*(1.+log(m1sq/mu2))));
	}
	return dcomplex(Ab(m1sq,mu2)/6.+pi16*(6.*m1sq - qsq)/18.)
	  + dcomplex((4.*m1sq - qsq)/12.)*Bb(m1sq,m1sq,qsq,mu2);
}

// version with numerical integration over x left in
static bintdata b22intdat;
dcomplex B22bnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b22intdat.m2sq+(1.-x)*b22intdat.m1sq-x*(1.-x)*b22intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b22intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b22intdat.mu2)-dcomplex(0.,pi);}
  return -0.5*mtsq*lmtsq;
}

dcomplex B22bnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b22intdat.m1sq = m1sq;
  b22intdat.m2sq = m2sq;
  b22intdat.qsq  = qsq;
  b22intdat.mu2  = mu2;
  return pi16*WINTEGRAL(B22bnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			 precisiononeloopintegrals/pi16);
}

// B31bar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
static bintdata b31intdat;

dcomplex B31bnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b31intdat.m2sq+(1.-x)*b31intdat.m1sq-x*(1.-x)*b31intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b31intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b31intdat.mu2)-dcomplex(0.,pi);}
  return -pow(x,3)*lmtsq;
}

dcomplex B31bnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b31intdat.m1sq = m1sq;
  b31intdat.m2sq = m2sq;
  b31intdat.qsq  = qsq;
  b31intdat.mu2  = mu2;
  return pi16*(-1./4.+WINTEGRAL(B31bnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			     precisiononeloopintegrals/pi16));
}

// B32bar integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
static bintdata b32intdat;

dcomplex B32bnumint(const dcomplex xx){
  double x = real(xx);
  dcomplex lmtsq;
  double mtsq = x*b32intdat.m2sq+(1.-x)*b32intdat.m1sq-x*(1.-x)*b32intdat.qsq;
  if ( mtsq >= 0.){
    lmtsq = log(mtsq/b32intdat.mu2);}
  else{
    lmtsq = log(-mtsq/b32intdat.mu2)-dcomplex(0.,pi);}
  return -0.5*x*mtsq*lmtsq;
}

dcomplex B32bnum(const double m1sq,const double m2sq,const double qsq,
	       const double mu2){
  b32intdat.m1sq = m1sq;
  b32intdat.m2sq = m2sq;
  b32intdat.qsq  = qsq;
  b32intdat.mu2  = mu2;
  return pi16*(WINTEGRAL(B32bnumint,dcomplex(0.,0.),dcomplex(1.,0.),
			     precisiononeloopintegrals/pi16));
}
//
//
//// The C integrals xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
//
//static bintdata c0intdat;
//
//dcomplex c0wint(const dcomplex xx){
//  dcomplex x = real(xx);
//  dcomplex ovfac = 1.;
//  // going into the complex plane avoiding the poles
//  if(c0intdat.qsq > pow(sqrt(c0intdat.m1sq)+sqrt(c0intdat.m2sq),2))
//  {
//    x = xx-dcomplex(0.,1.)*xx*(1.-xx)*
//      (-c0intdat.m1sq+c0intdat.m2sq-c0intdat.qsq*(1.-2.*x))/c0intdat.qsq;
//    ovfac = ovfac*(1.-dcomplex(0.,1./c0intdat.qsq)*
//		   ((1.-2.*xx)*(-c0intdat.m1sq+c0intdat.m2sq-c0intdat.qsq)
//		    +2.*xx*(2.-3.*xx)*c0intdat.qsq));
//  }
//
//  dcomplex mtsq = 
//    x*c0intdat.m2sq+(1.-x)*c0intdat.m1sq-x*(1.-x)*c0intdat.qsq;
//  return -(1.-x)/mtsq*ovfac;
//}
//
//dcomplex c0int(double m1sq,double m2sq,double m3sq,double qsq, double mu2){
//  if (m1sq != m2sq){
//    return 1./(m1sq-m2sq)*(b0int(m1sq,m3sq,qsq,mu2)
//			   -b0int(m2sq,m3sq,qsq,mu2));
//  }
//  c0intdat.m1sq = m1sq;
//  c0intdat.m2sq = m3sq;
//  c0intdat.qsq  = qsq;
//  c0intdat.mu2  = mu2;
//  return pi16*(WINTEGRAL(c0wint,dcomplex(0.,0.),dcomplex(1.,0.),
//			     precisiononeloopintegrals/pi16));
//}
//
//
//dcomplex c0intana(double m1sq,double m2sq,double m3sq,double qsq, double mu2){
//  if (m1sq == m2sq){
//    if(m1sq == m3sq){
//      if (qsq == 0.){
//	return  -pi16/(2.*m1sq);
//      }
//      return -1./(qsq-4.*m1sq)*(b0intanae(m1sq,qsq,mu2)
//			      -pi16+pi16*log(m1sq/mu2));
//    }
//    else{
//      if (qsq == 0.){
//	return pi16*(m3sq-m1sq-m3sq*log(m3sq/m1sq))/pow(m3sq-m1sq,2);
//      }
//      double b = m3sq-m1sq-qsq;
//      double d = b*b-4.*m1sq*qsq;
//      double bt = b/(2.*qsq);
//      double y1 = 1.+bt;
//      double y0 = bt;
//
//      if (d >0){
//	double dt = sqrt(d/(4.*qsq*qsq));
//	dcomplex c0 = 1./qsq*(
//	  +(1.+bt)/(2.*dt)*log(fabs((dt+y1)/(dt-y1)))
//	  +1./2.*log(fabs(y1*y1-dt*dt))
//	  -(1.+bt)/(2.*dt)*log(fabs((dt+y0)/(dt-y0)))
//	  -1./2.*log(fabs(y0*y0-dt*dt)));
//	if(qsq > pow(sqrt(m1sq)+sqrt(m2sq),2)){
//	  double xp = -bt+dt;
//	  double xm = -bt-dt;
//	  c0 += dcomplex(0.,pi*(-(1.-xp)/fabs(-m1sq+m3sq-(1.-2.*xp)*qsq)
//				-(1.-xm)/fabs(-m1sq+m3sq-(1.-2.*xm)*qsq)));
//	}
//	return pi16*c0;
//      }
//      else{
//	double dt = sqrt(-d/(4.*qsq*qsq));
//	dcomplex c0 = 1./qsq*(
//	  -(1.+bt)/dt*atan(y1/dt)+1./2.*log(fabs(y1*y1+dt*dt))
//	  +(1.+bt)/dt*atan(y0/dt)-1./2.*log(fabs(y0*y0+dt*dt)));
//	return pi16*c0;
//      }
//    }
//
//  }
//  if (m1sq != m2sq){
//    return 1./(m1sq-m2sq)*(b0intana(m1sq,m3sq,qsq,mu2)
//			   -b0intana(m2sq,m3sq,qsq,mu2));
//  }
//  //  std::cout << "mass case not implemented in c0intana yet"<<'\n';
//  return 0;
//}
//
//
//
//
//
//// The epsilon versions
//// Beps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b0epsdat;
//
//dcomplex b0weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b0epsdat.m2sq+(1.-x)*b0epsdat.m1sq-x*(1.-x)*b0epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b0epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b0epsdat.mu2)-dcomplex(0.,pi);}
//  return  0.5*pow(lmtsq,2)+lmtsq;
//}
//
//dcomplex b0eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b0epsdat.m1sq = m1sq;
//  b0epsdat.m2sq = m2sq;
//  b0epsdat.qsq = qsq;
//  b0epsdat.mu2 = mu2;
//  return pi16*(pi*pi/12.+0.5+
//	       WINTEGRAL(b0weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		      precisiononeloopintegrals/pi16));
//}
//
//dcomplex b0epsana(double m1sq,double m2sq,double qsq, double mu2){
//
//  if (qsq <= 0.){
//    double b = -m1sq+m2sq-qsq;
//    double dt = sqrt(b*b-4.*m1sq*qsq);
//    double xm = (-b+dt)/(2.*qsq);  // beware: xp must be > 1
//    double xp = (-b-dt)/(2.*qsq);  // beware: xm must be < 0
//    double xpm = xp-xm;
//    double lqsq = log(-qsq/mu2);
//    double lxp  = log(xp);
//    double lxp1 = log(xp-1.);
//    double lmxm = log(-xm);
//    double l1xm = log(1.-xm);
//    return pi16*(
//		 pi2/12.+0.5+lqsq+(lqsq-1.)*
//		    (xp*lxp-(xp-1.)*lxp1+(1.-xm)*l1xm+xm*lmxm-2.)
//		 +0.5*lqsq*lqsq+0.5*xp*lxp*lxp-0.5*(xp-1.)*lxp1*lxp1
//                 +0.5*(1.-xm)*l1xm*l1xm+0.5*xm*lmxm*lmxm
//                 +xp*lxp*lmxm-(xp-1.)*lxp1*l1xm
//                 +xpm*log(xpm)*(l1xm-lmxm)-xpm*jbdli2((1.-xm)/xpm)
//		 +xpm*jbdli2(-xm/xpm));
//  }
//  if (qsq <= pow(sqrt(m1sq)-sqrt(m2sq),2)){
//    double b = -m1sq+m2sq-qsq;
//    double dt = sqrt(b*b-4.*m1sq*qsq);
//    double xp,xm;
//    if (b >= 0){              // m2sq > m1sq
//      xp = (-b+dt)/(2.*qsq);  // beware:  xm < xp < 0 here
//      xm = (-b-dt)/(2.*qsq);
//    }
//    else{
//      xm = 1.-(-b+dt)/(2.*qsq);  // corresponds to m1sq<->m2sq
//      xp = 1.-(-b-dt)/(2.*qsq);  //  and gets xm < xp < 0
//    }
//    double xpm = xp-xm;
//    double lqsq = log(qsq/mu2);
//    double lmxp = log(-xp);
//    double l1xp = log(1.-xp);
//    double lmxm = log(-xm);
//    double l1xm = log(1.-xm);
//    return pi16*(
//		 pi2/12.+0.5+lqsq+(lqsq-1.)*
//		    (xp*lmxp+(1.-xp)*l1xp+(1.-xm)*l1xm+xm*lmxm-2.)
//		 +0.5*lqsq*lqsq+0.5*xp*lmxp*lmxp+0.5*(1.-xp)*l1xp*l1xp
//                 +0.5*(1.-xm)*l1xm*l1xm+0.5*xm*lmxm*lmxm
//                 +xp*lmxp*lmxm+(1.-xp)*l1xp*l1xm
//		 +0.5*xpm*(l1xm*l1xm-lmxm*lmxm)-xpm*jbdli2(-xpm/xm)
//		 +xpm*jbdli2(xpm/(1.-xm)));
//
//  }
//  if (qsq <= pow(sqrt(m1sq)+sqrt(m2sq),2)){
//    double b = -m1sq+m2sq-qsq;
//    double bt = b/(2.*qsq);
//    double dt = sqrt(-b*b+4.*m1sq*qsq)/(2.*qsq);
//    dcomplex xp = dcomplex(-bt,dt);
//    dcomplex xm = dcomplex(-bt,-dt);
//    dcomplex xpm = dcomplex(0.,2.*dt);
//    dcomplex lqsq = log(qsq/mu2);
//    dcomplex lmxp = log(-xp);
//    dcomplex l1xp = log(1.-xp);
//    dcomplex lmxm = log(-xm);
//    dcomplex l1xm = log(1.-xm);
//    dcomplex c0 =
//		 pi2/12.+0.5+lqsq+(lqsq-1.)*
//		    (xp*lmxp+(1.-xp)*l1xp+(1.-xm)*l1xm+xm*lmxm-2.)
//		 +0.5*lqsq*lqsq+0.5*xp*lmxp*lmxp+0.5*(1.-xp)*l1xp*l1xp
//                 +0.5*(1.-xm)*l1xm*l1xm+0.5*xm*lmxm*lmxm
//		 +xp*lmxp*lmxm+(1.-xp)*l1xp*l1xm
//		 +0.5*xpm*(l1xm*l1xm-lmxm*lmxm)-xpm*jbdli2(-xpm/xm)
//		 +xpm*jbdli2(xpm/(1.-xm));
//    if (bt < 0.){
//      if (bt < -1){
//	c0 += -2.*pi*dcomplex(0.,1.)*xpm*(l1xm-lmxm);
//      }
//      else{
//	c0 += -2.*pi*dcomplex(0.,1.)*xpm*(log(dcomplex(0.,dt))+log(2.)-lmxm);
//      }
//    }
//    return pi16*c0;
//  }
//
//  //  if (qsq >= pow(sqrt(m1sq)+sqrt(m2sq),2))
//
//  double b = m2sq-m1sq-qsq;
//  double d = b*b-4.*m1sq*qsq;
//  double dt = sqrt(d);
//  double xp = (-b+dt)/(2.*qsq);
//  double xm = (-b-dt)/(2.*qsq);
//  double xpm = xp-xm;
//  double lqsq = log(qsq/mu2);
//  double lxp =  log(xp);
//  double l1xp = log(1.-xp);
//  double lxm =  log(xm);
//  double l1xm = log(1.-xm);
//  double lxpm = log(xpm);
//  // part from the discontinuity
//  dcomplex c0 = dcomplex(0.,pi)*xpm*(1.-2.*lxpm-lqsq)-0.5*pi2*xpm;
//  // rest
//  dcomplex b0 = 
//    pi2/12.+0.5+lqsq+(lqsq-1.)*
//        (xp*lxp+(1.-xp)*l1xp+(1.-xm)*l1xm+xm*lxm-2.)
//    +0.5*lqsq*lqsq+0.5*xp*lxp*lxp+0.5*(1.-xp)*l1xp*l1xp
//    +0.5*(1.-xm)*l1xm*l1xm+0.5*xm*lxm*lxm
//    +xp*lxp*lxm+(1.-xp)*l1xp*l1xm
//    +xpm*lxpm*lxpm-xpm*lxpm*lxm+xpm*jbdli2(-xm/xpm)
//    +0.5*xpm*(l1xm*l1xm-lxpm*lxpm)+xpm*jbdli2(xpm/(1.-xm))-xpm*pi2/3.;
//  return pi16*(b0+c0);
//}
//
//// B1eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b1epsdat;
//
//dcomplex b1weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b1epsdat.m2sq+(1.-x)*b1epsdat.m1sq-x*(1.-x)*b1epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b1epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b1epsdat.mu2)-dcomplex(0.,pi);}
//  return  x*(0.5*pow(lmtsq,2)+lmtsq);
//}
//
//dcomplex b1eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b1epsdat.m1sq = m1sq;
//  b1epsdat.m2sq = m2sq;
//  b1epsdat.qsq = qsq;
//  b1epsdat.mu2 = mu2;
//  return  pi16*(pi*pi/24.+0.25+
//      WINTEGRAL(b1weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		      precisiononeloopintegrals/pi16));
//}
//
//// B21eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b21epsdat;
//
//dcomplex b21weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b21epsdat.m2sq+(1.-x)*b21epsdat.m1sq-x*(1.-x)*b21epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b21epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b21epsdat.mu2)-dcomplex(0.,pi);}
//  return  x*x*(0.5*pow(lmtsq,2)+lmtsq);
//}
//
//dcomplex b21eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b21epsdat.m1sq = m1sq;
//  b21epsdat.m2sq = m2sq;
//  b21epsdat.qsq = qsq;
//  b21epsdat.mu2 = mu2;
//  return  pi16*(pi*pi/36.+1./6.+
//      WINTEGRAL(b21weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		      precisiononeloopintegrals/pi16));
//}
//
//// B22eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b22epsdat;
//
//dcomplex b22weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b22epsdat.m2sq+(1.-x)*b22epsdat.m1sq-x*(1.-x)*b22epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b22epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b22epsdat.mu2)-dcomplex(0.,pi);}
//  return  mtsq*(pi*pi/24.+0.25+0.25*pow(lmtsq,2));
//}
//
//dcomplex b22eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b22epsdat.m1sq = m1sq;
//  b22epsdat.m2sq = m2sq;
//  b22epsdat.qsq = qsq;
//  b22epsdat.mu2 = mu2;
//  return  pi16*(WINTEGRAL(b22weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		       precisiononeloopintegrals/pi16));
//}
//
//
//// B31eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b31epsdat;
//
//dcomplex b31weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b31epsdat.m2sq+(1.-x)*b31epsdat.m1sq-x*(1.-x)*b31epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b31epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b31epsdat.mu2)-dcomplex(0.,pi);}
//  return  pow(x,3)*(0.5*pow(lmtsq,2)+lmtsq);
//}
//
//dcomplex b31eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b31epsdat.m1sq = m1sq;
//  b31epsdat.m2sq = m2sq;
//  b31epsdat.qsq = qsq;
//  b31epsdat.mu2 = mu2;
//  return  pi16*(pi*pi/48.+1./8.+
//      WINTEGRAL(b31weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		      precisiononeloopintegrals/pi16));
//}
//
//
//// B32eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//static bintdata b32epsdat;
//
//dcomplex b32weps(const dcomplex xx){
//  double x = real(xx);
//  dcomplex lmtsq;
//  double mtsq = x*b32epsdat.m2sq+(1.-x)*b32epsdat.m1sq-x*(1.-x)*b32epsdat.qsq;
//  if (mtsq >= 0.){
//    lmtsq = log(mtsq/b32epsdat.mu2);}
//  else{
//    lmtsq = log(-mtsq/b32epsdat.mu2)-dcomplex(0.,pi);}
//  return  x*mtsq*(pi*pi/24.+0.25+0.25*pow(lmtsq,2));
//}
//
//dcomplex b32eps(const double m1sq,const double m2sq,const double qsq,
//    const double mu2){
//  b32epsdat.m1sq = m1sq;
//  b32epsdat.m2sq = m2sq;
//  b32epsdat.qsq = qsq;
//  b32epsdat.mu2 = mu2;
//  return  pi16*(WINTEGRAL(b32weps,dcomplex(0.,0.),dcomplex(1.,0.),
//		      precisiononeloopintegrals/pi16));
//}
//
//// C0eps integral xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
//static bintdata c0epsdat;
//
//dcomplex c0weps(const dcomplex xx){
//  dcomplex x = xx;
//  dcomplex ovfac = 1.;
//  // going into the complex plane avoiding the poles
//  if(c0epsdat.qsq > pow(sqrt(c0epsdat.m1sq)+sqrt(c0epsdat.m2sq),2))
//  {
//    x = xx-dcomplex(0.,1.)*xx*(1.-xx)*
//      (-c0epsdat.m1sq+c0epsdat.m2sq-c0epsdat.qsq*(1.-2.*x))/c0epsdat.qsq;
//    ovfac = ovfac*(1.-dcomplex(0.,1./c0epsdat.qsq)*
//		   ((1.-2.*xx)*(-c0epsdat.m1sq+c0epsdat.m2sq-c0epsdat.qsq)
//		    +2.*xx*(2.-3.*xx)*c0epsdat.qsq));
//  }
//
//  dcomplex mtsq = 
//    x*c0epsdat.m2sq+(1.-x)*c0epsdat.m1sq-x*(1.-x)*c0epsdat.qsq;
//  return (1.-x)/mtsq*(1.+log(mtsq/c0epsdat.mu2))*ovfac;
//}
//
//dcomplex c0eps(double m1sq,double m2sq,double m3sq,double qsq, double mu2){
//  if (m1sq != m2sq){
//    return 1./(m1sq-m2sq)*(b0eps(m1sq,m3sq,qsq,mu2)
//			   -b0eps(m2sq,m3sq,qsq,mu2));
//  }
//  c0epsdat.m1sq = m1sq;
//  c0epsdat.m2sq = m3sq;
//  c0epsdat.qsq  = qsq;
//  c0epsdat.mu2  = mu2;
//  return pi16*(WINTEGRAL(c0weps,dcomplex(0.,0.),dcomplex(1.,0.),
//			     precisiononeloopintegrals/pi16));
//}
//
//// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
//
//dcomplex dqsqBb(double m1sq,double m2sq,double qsq,double mu2){
//  double d;
//  dcomplex b0intana;
//  double a=qsq;
//  double b=m1sq-m2sq-qsq;
//  double c=m2sq;
//  double dsq = (b*b-4.*c*a)/(4.*a*a);
//  double ddsq = (qsq-m1sq-m2sq-4.*qsq*dsq)/(2.*qsq*qsq);
//  double y1=1.0+b/(2.*a); // 1+e in the notes
//  double y0=b/(2.*a); //      e in the notes
//  b0intana = 1./qsq-(1.+2.*y0)/(2.*qsq)*(log(fabs((y1*y1-dsq)/(y0*y0-dsq))));
//  if (dsq >= 0){
//    d=sqrt(dsq);
//    b0intana += ddsq/(2.*d)*(log(fabs((y1+d)/(y1-d)))
//			    -log(fabs((y0+d)/(y0-d))));
//  }
//  else{
//    d=sqrt(-dsq);
//    b0intana += (-ddsq)/d*(atan(y1/d)-atan(y0/d));
//  }
//  b0intana *= (-pi16);
//  if (qsq >= pow(sqrt(m1sq)+sqrt(m2sq),2)){
//    double dsqrt = sqrt(pow((qsq-m1sq-m2sq),2)-4.0*m1sq*m2sq);
//    b0intana = b0intana
//	+dcomplex(0.0,pi16*pi/qsq*(-dsqrt/qsq+(qsq-m1sq-m2sq)/dsqrt));
//  }
//  return b0intana;
//}
//
//// incomplete ones xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//
//dcomplex c0epsana(double m1sq,double m2sq,double m3sq,double qsq, double mu2){
//  // contains only imaginary part for m1sq = m2sq
//  if (m1sq != m2sq) return dcomplex(0.,0.);
//  if (qsq <= pow(sqrt(m1sq)+sqrt(m3sq),2)) return dcomplex(0.,0.);
//  double b = m3sq-m1sq-qsq;
//  double d = b*b-4.*m1sq*qsq;
//  //double bt = b/(2.*qsq);
//  //double y1 = 1.+bt;
//  //double y0 = bt;
//  double dt = sqrt(d/(4.*qsq*qsq));
//  //double xp = -bt+dt;
//  //double xm = -bt-dt;
//  dcomplex c0 = dcomplex(0.,pi)*(m1sq-m3sq-qsq)/(qsq*qsq)/(2.*dt)*
//			 (-1.-2.*log(2.*dt)-log(qsq/mu2));
//  return pi16*c0;
//}
