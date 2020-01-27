// massdecayvevnfV.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// finite volume results results for QCD like theories for masses,
// decay constants and // the qbarq vacuume expectation value
// derived in
// J.~Bijnens and T.~Rössler

// in terms of the lowest order mass and the chiral limit
// decay constant

#include "oneloopintegrals.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"
#include "massdecayvevnfV.h"
#include <cmath>

const double pi16 = 1./(16.*M_PI*M_PI);

// the code that produces results with Bessel method or theta method
#ifdef CHIRONTHETA
#define AbV AbVt
#define BbV BbVt
#define A23bV A23bVt
#define hhV hhVt
#define hh1V hh1Vt
#define hh21V hh21Vt
#define hh27V hh27Vt
#define hhdV hhdVt
#define hh1dV hh1dVt
#define hh21dV hh21dVt
#define hh27dV hh27dVt
#define mnfSUNp4V mnfSUNp4Vt
#define mnfSUNp6V mnfSUNp6Vt
#define mnfSUNp6LV mnfSUNp6LVt
#define mnfSUNp6RV mnfSUNp6RVt
#define fnfSUNp4V fnfSUNp4Vt
#define fnfSUNp6V fnfSUNp6Vt
#define fnfSUNp6LV fnfSUNp6LVt
#define fnfSUNp6RV fnfSUNp6RVt
#define qnfSUNp4V qnfSUNp4Vt
#define qnfSUNp6V qnfSUNp6Vt
#define qnfSUNp6LV qnfSUNp6LVt
#define qnfSUNp6RV qnfSUNp6RVt
#define mnfSONp4V mnfSONp4Vt
#define mnfSONp6V mnfSONp6Vt
#define mnfSONp6LV mnfSONp6LVt
#define mnfSONp6RV mnfSONp6RVt
#define fnfSONp4V fnfSONp4Vt
#define fnfSONp6V fnfSONp6Vt
#define fnfSONp6LV fnfSONp6LVt
#define fnfSONp6RV fnfSONp6RVt
#define qnfSONp4V qnfSONp4Vt
#define qnfSONp6V qnfSONp6Vt
#define qnfSONp6LV qnfSONp6LVt
#define qnfSONp6RV qnfSONp6RVt
#define mnfSPNp4V mnfSPNp4Vt
#define mnfSPNp6V mnfSPNp6Vt
#define mnfSPNp6LV mnfSPNp6LVt
#define mnfSPNp6RV mnfSPNp6RVt
#define fnfSPNp4V fnfSPNp4Vt
#define fnfSPNp6V fnfSPNp6Vt
#define fnfSPNp6LV fnfSPNp6LVt
#define fnfSPNp6RV fnfSPNp6RVt
#define qnfSPNp4V qnfSPNp4Vt
#define qnfSPNp6V qnfSPNp6Vt
#define qnfSPNp6LV qnfSPNp6LVt
#define qnfSPNp6RV qnfSPNp6RVt
#else
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef CHIRONBESSEL
#define AbV AbVb
#define BbV BbVb
#define A23bV A23bVb
#define hhV hhVb
#define hh1V hh1Vb
#define hh21V hh21Vb
#define hh27V hh27Vb
#define hhdV hhdVb
#define hh1dV hh1dVb
#define hh21dV hh21dVb
#define hh27dV hh27dVb
#define mnfSUNp4V mnfSUNp4Vb
#define mnfSUNp6V mnfSUNp6Vb
#define mnfSUNp6LV mnfSUNp6LVb
#define mnfSUNp6RV mnfSUNp6RVb
#define fnfSUNp4V fnfSUNp4Vb
#define fnfSUNp6V fnfSUNp6Vb
#define fnfSUNp6LV fnfSUNp6LVb
#define fnfSUNp6RV fnfSUNp6RVb
#define qnfSUNp4V qnfSUNp4Vb
#define qnfSUNp6V qnfSUNp6Vb
#define qnfSUNp6LV qnfSUNp6LVb
#define qnfSUNp6RV qnfSUNp6RVb
#define mnfSONp4V mnfSONp4Vb
#define mnfSONp6V mnfSONp6Vb
#define mnfSONp6LV mnfSONp6LVb
#define mnfSONp6RV mnfSONp6RVb
#define fnfSONp4V fnfSONp4Vb
#define fnfSONp6V fnfSONp6Vb
#define fnfSONp6LV fnfSONp6LVb
#define fnfSONp6RV fnfSONp6RVb
#define qnfSONp4V qnfSONp4Vb
#define qnfSONp6V qnfSONp6Vb
#define qnfSONp6LV qnfSONp6LVb
#define qnfSONp6RV qnfSONp6RVb
#define mnfSPNp4V mnfSPNp4Vb
#define mnfSPNp6V mnfSPNp6Vb
#define mnfSPNp6LV mnfSPNp6LVb
#define mnfSPNp6RV mnfSPNp6RVb
#define fnfSPNp4V fnfSPNp4Vb
#define fnfSPNp6V fnfSPNp6Vb
#define fnfSPNp6LV fnfSPNp6LVb
#define fnfSPNp6RV fnfSPNp6RVb
#define qnfSPNp4V qnfSPNp4Vb
#define qnfSPNp6V qnfSPNp6Vb
#define qnfSPNp6LV qnfSPNp6LVb
#define qnfSPNp6RV qnfSPNp6RVb
#else
// just some garbage to produce an error
x = 1./0.0.;
#endif
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//////////////////////////////////////////////////////////////////////////
//////////////// SUN /////////////////////////////////////////////////////
// mass SUN
double mnfSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double  MM12p4x11VR =
       + abvm11 * (  - m11*nfm );
  return  MM12p4x11VR/pow(mass.getf0(),2);
}

double mnfSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfSUNp6LV(nf,mass,Liin,L)+mnfSUNp6RV(nf,mass,L);
}

double mnfSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MM12p6x11VL =
       + bbvm11 * (  - 16*pow(m11,3)*L8r*nfm - 16*pow(m11,3)*L6r + 8*
         pow(m11,3)*L5r*nfm + 8*pow(m11,3)*L4r );

      MM12p6x11VL +=  + a23bvm11 * ( 24*m11*L3r*nfm - 12*m11*L3r*nf - 
         12*m11*L2r*pow(nf,2) - 24*m11*L1r + 24*m11*L0r*nfm - 24*m11*
         L0r*nf );

      MM12p6x11VL +=  + abvm11 * (  - 64*pow(m11,2)*L8r*nfm + 32*pow(
         m11,2)*L8r*nf + 16*pow(m11,2)*L6r*pow(nf,2) + 32*pow(m11,2)*
         L5r*nfm - 16*pow(m11,2)*L5r*nf - 24*pow(m11,2)*L3r*nfm + 20*
         pow(m11,2)*L3r*nf - 24*pow(m11,2)*L0r*nfm + 8*pow(m11,2)*L0r*
         nf - 16*( - 2 + pow(nf,2))*pow(m11,2)*L4r + 8*( - 1 + 2*pow(
         nf,2))*pow(m11,2)*L1r + 4*(4 + pow(nf,2))*pow(m11,2)*L2r );
 
  return MM12p6x11VL/pow(mass.getf0(),4);
}

double mnfSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MM12p6x11VR =
       + abvm11*bbvm11 * ( pow(m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + pow(abvm11,2) * ( 1./2.*m11*pow(nfm,2) - 1./4.*
         ( - 2 + pow(nf,2))*m11 );

      MM12p6x11VR +=  + pi16*abvm11 * (  - pow(m11,2)*pow(nfm,2) + 1./3.
         *pow(m11,2)*pow(nf,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( pow(m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 2*m11*pow(nfm,2) - 1./2.
         *( - 2 + pow(nf,2))*m11 );

      MM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * ( 2*pow(m11,2)*
         pow(nfm,2) + 1./12.*( - 8 + 3*pow(nf,2))*pow(m11,2) );

      MM12p6x11VR +=  + hh21V(1,m11,m11,m11,m11,L,mu2) * ( 3./4.*pow(
         m11,2)*pow(nf,2) );

      MM12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * (  - 3./4.*m11
         *pow(nf,2) );

  return  MM12p6x11VR/pow(mass.getf0(),4);
}
///////////////////////////////////////////////////////////////////////////
// decay SUN
double fnfSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double  MF12p4x11VR =
       + abvm11 * ( 1./2.*nf );

  return  MF12p4x11VR/pow(mass.getf0(),2);
}

double fnfSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfSUNp6LV(nf,mass,Liin,L)+fnfSUNp6RV(nf,mass,L);
}

double fnfSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MF12p6x11VL =
       + bbvm11 * ( 8*pow(m11,2)*L8r*nf + 8*pow(m11,2)*L6r*pow(nf,2) - 
         4*pow(m11,2)*L5r*nf - 4*pow(m11,2)*L4r*pow(nf,2) );

      MF12p6x11VL +=  + a23bvm11 * (  - 12*L3r*nfm + 6*L3r*nf + 6*L2r*
         pow(nf,2) + 12*L1r - 12*L0r*nfm + 12*L0r*nf );

      MF12p6x11VL +=  + abvm11 * (  - 4*m11*L5r*nfm + 2*m11*L5r*nf + 12
         *m11*L3r*nfm - 10*m11*L3r*nf + 12*m11*L0r*nfm - 4*m11*L0r*nf
          + 2*( - 2 + pow(nf,2))*m11*L4r - 4*( - 1 + 2*pow(nf,2))*m11*
         L1r - 2*(4 + pow(nf,2))*m11*L2r );
 
  return MF12p6x11VL/pow(mass.getf0(),4);
}

double fnfSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MF12p6x11VR =
       + abvm11*bbvm11 * (  - 1./2.*m11 );

      MF12p6x11VR +=  + pi16*abvm11 * (  - 1./24.*( - 12 + 5*pow(nf,2))
         *m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * (  - 1./2.*m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*abvm11 * (  - 1./2. );

      MF12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * (  - 1./8.*m11*
         pow(nf,2) );

      MF12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * ( 3./8.*pow(
         nf,2) );

      MF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * ( pow(m11,2)*
         pow(nfm,2) + 1./24.*( - 8 + 3*pow(nf,2))*pow(m11,2) );

      MF12p6x11VR +=  + hh21dV(1,m11,m11,m11,m11,L,mu2) * ( 3./8.*pow(
         m11,2)*pow(nf,2) );

      MF12p6x11VR +=  + hh27dV(1,m11,m11,m11,m11,L,mu2) * (  - 3./8.*
         m11*pow(nf,2) );

  return  MF12p6x11VR/pow(mass.getf0(),4);
}

///////////////////////////////////////////////////////////////////////////
// vev SUN
double qnfSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double  MQ11p4x11VR =
       + abvm11 * (  - nfm + nf );

  return  MQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfSUNp6LV(nf,mass,Liin,L)+qnfSUNp6RV(nf,mass,L);
}

double qnfSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  //double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MQ11p6x11VL =
       + bbvm11 * (  - 16*pow(m11,2)*L8r*nfm + 16*pow(m11,2)*L8r*nf + 8
         *pow(m11,2)*L5r*nfm - 8*pow(m11,2)*L5r*nf + 16*( - 1 + nf)*(1
          + nf)*pow(m11,2)*L6r - 8*( - 1 + nf)*(1 + nf)*pow(m11,2)*L4r
          );

      MQ11p6x11VL +=  + abvm11 * (  - 32*m11*L8r*nfm + 32*m11*L8r*nf + 
         16*m11*L5r*nfm - 16*m11*L5r*nf + 32*( - 1 + nf)*(1 + nf)*m11*
         L6r - 16*( - 1 + nf)*(1 + nf)*m11*L4r );
 
  return MQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MQ11p6x11VR =
       + abvm11*bbvm11 * (  - m11 + m11*pow(nfm,2) );

      MQ11p6x11VR +=  + pow(abvm11,2) * (  - 1./2. + 1./2.*pow(nfm,2) )
         ;

      MQ11p6x11VR +=  + pi16*abvm11 * ( m11 - m11*pow(nfm,2) );

      MQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * (  - m11 + m11*pow(nfm,2)
          );

      MQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * (  - 2 + 2*pow(nfm,2) );

  return  MQ11p6x11VR/pow(mass.getf0(),4);
}


//////////////// SON /////////////////////////////////////////////////////
// mass SON
double mnfSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double   MM12p4x11VR =
       + abvm11 * ( 1./2.*m11 - m11*nfm );

  return  MM12p4x11VR/pow(mass.getf0(),2);
}

double mnfSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfSONp6LV(nf,mass,Liin,L)+mnfSONp6RV(nf,mass,L);
}

double mnfSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MM12p6x11VL =
       + bbvm11 * ( 8*pow(m11,3)*L8r - 16*pow(m11,3)*L8r*nfm - 4*pow(
         m11,3)*L5r + 8*pow(m11,3)*L5r*nfm + 8*( - 2 + nf)*pow(m11,3)*
         L6r - 4*( - 2 + nf)*pow(m11,3)*L4r );

      MM12p6x11VL +=  + a23bvm11 * ( 24*m11*L3r*nfm - 24*m11*L1r + 24*
         m11*L0r*nfm - 6*(1 + nf)*m11*L2r*nf - 12*(1 + nf)*m11*L0r - 6*
         (2 + nf)*m11*L3r );

      MM12p6x11VL +=  + abvm11 * (  - 64*pow(m11,2)*L8r*nfm + 32*pow(
         m11,2)*L5r*nfm - 24*pow(m11,2)*L3r*nfm - 24*pow(m11,2)*L0r*nfm
          - 8*( - 4 + 2*nf + pow(nf,2))*pow(m11,2)*L4r + 8*( - 1 + nf
          + pow(nf,2))*pow(m11,2)*L1r + 16*(2 + nf)*pow(m11,2)*L8r + 8*
         (2 + nf)*pow(m11,2)*L6r*nf - 8*(2 + nf)*pow(m11,2)*L5r + 4*(3
          + nf)*pow(m11,2)*L0r + 2*(6 + 5*nf)*pow(m11,2)*L3r + 2*(8 + 
         nf + pow(nf,2))*pow(m11,2)*L2r );
 
  return MM12p6x11VL/pow(mass.getf0(),4);
}

double mnfSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MM12p6x11VR =
       + abvm11*bbvm11 * ( 1./4.*pow(m11,2) - pow(m11,2)*nfm + pow(
         m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + pow(abvm11,2) * (  - 1./2.*m11*nfm + 1./2.*m11*
         pow(nfm,2) - 1./16.*( - 6 + 4*nf + pow(nf,2))*m11 );

      MM12p6x11VR +=  + pi16*abvm11 * ( pow(m11,2)*nfm - pow(m11,2)*
         pow(nfm,2) + 1./12.*( - 1 + nf)*(3 + nf)*pow(m11,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 1./4.*pow(m11,2) - pow(
         m11,2)*nfm + pow(m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * (  - 2*m11*nfm + 2*m11*
         pow(nfm,2) - 1./8.*( - 8 + 4*nf + pow(nf,2))*m11 );

      MM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * (  - pow(m11,2)*
         nfm + 2*pow(m11,2)*pow(nfm,2) + 1./48.*( - 8 + 18*nf + 3*pow(
         nf,2))*pow(m11,2) );

      MM12p6x11VR +=  + hh1V(1,m11,m11,m11,m11,L,mu2) * (  - 1./2.*pow(
         m11,2)*nf );

      MM12p6x11VR +=  + hh21V(1,m11,m11,m11,m11,L,mu2) * ( 3./16.*(2 + 
         nf)*pow(m11,2)*nf );

      MM12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * (  - 3./16.*(2
          + nf)*m11*nf );

  return  MM12p6x11VR/pow(mass.getf0(),4);
}
///////////////////////////////////////////////////////////////////////////
// decay SON

double fnfSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double   MF12p4x11VR =
       + abvm11 * ( 1./4.*nf );

  return  MF12p4x11VR/pow(mass.getf0(),2);
}

double fnfSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfSONp6LV(nf,mass,Liin,L)+fnfSONp6RV(nf,mass,L);
}

double fnfSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MF12p6x11VL =
       + bbvm11 * ( 4*pow(m11,2)*L8r*nf + 4*pow(m11,2)*L6r*pow(nf,2) - 
         2*pow(m11,2)*L5r*nf - 2*pow(m11,2)*L4r*pow(nf,2) );

      MF12p6x11VL +=  + a23bvm11 * (  - 12*L3r*nfm + 12*L1r - 12*L0r*
         nfm + 3*(1 + nf)*L2r*nf + 6*(1 + nf)*L0r + 3*(2 + nf)*L3r );

      MF12p6x11VL +=  + abvm11 * (  - 4*m11*L5r*nfm + 12*m11*L3r*nfm + 
         12*m11*L0r*nfm + ( - 4 + 2*nf + pow(nf,2))*m11*L4r - 4*( - 1
          + nf + pow(nf,2))*m11*L1r + (2 + nf)*m11*L5r - 2*(3 + nf)*m11
         *L0r - (6 + 5*nf)*m11*L3r - (8 + nf + pow(nf,2))*m11*L2r );

  return MF12p6x11VL/pow(mass.getf0(),4);
}

double fnfSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MF12p6x11VR =
       + abvm11*bbvm11 * ( 1./8.*( - 2 + nf)*m11 );

      MF12p6x11VR +=  + pow(abvm11,2) * ( 1./16.*nf );

      MF12p6x11VR +=  + pi16*abvm11 * (  - 1./96.*( - 24 + 22*nf + 5*
         pow(nf,2))*m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 1./8.*( - 2 + nf)*m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 1./4.*( - 1 + nf) );

      MF12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * (  - 1./32.*(2
          + nf)*m11*nf );

      MF12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * ( 3./32.*(2 + 
         nf)*nf );

      MF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * (  - 1./2.*pow(
         m11,2)*nfm + pow(m11,2)*pow(nfm,2) + 1./96.*( - 8 + 18*nf + 3*
         pow(nf,2))*pow(m11,2) );

      MF12p6x11VR +=  + hh1dV(1,m11,m11,m11,m11,L,mu2) * (  - 1./4.*
         pow(m11,2)*nf );

      MF12p6x11VR +=  + hh21dV(1,m11,m11,m11,m11,L,mu2) * ( 3./32.*(2
          + nf)*pow(m11,2)*nf );

      MF12p6x11VR +=  + hh27dV(1,m11,m11,m11,m11,L,mu2) * (  - 3./32.*(
         2 + nf)*m11*nf );


  return  MF12p6x11VR/pow(mass.getf0(),4);
}

///////////////////////////////////////////////////////////////////////////
// vev SON

double qnfSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double   MQ11p4x11VR =
       + abvm11 * (  - nfm + 1./2.*(1 + nf) );

  return  MQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfSONp6LV(nf,mass,Liin,L)+qnfSONp6RV(nf,mass,L);
}

double qnfSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  //double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MQ11p6x11VL =
       + bbvm11 * (  - 16*pow(m11,2)*L8r*nfm + 8*pow(m11,2)*L5r*nfm + 8
         *( - 1 + nf)*(2 + nf)*pow(m11,2)*L6r - 4*( - 1 + nf)*(2 + nf)*
         pow(m11,2)*L4r + 8*(1 + nf)*pow(m11,2)*L8r - 4*(1 + nf)*pow(
         m11,2)*L5r );

      MQ11p6x11VL +=  + abvm11 * (  - 32*m11*L8r*nfm + 16*m11*L5r*nfm
          + 16*( - 1 + nf)*(2 + nf)*m11*L6r - 8*( - 1 + nf)*(2 + nf)*
         m11*L4r + 16*(1 + nf)*m11*L8r - 8*(1 + nf)*m11*L5r );
 
  return MQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double  MQ11p6x11VR =
       + abvm11*bbvm11 * (  - m11*nfm + m11*pow(nfm,2) + 1./4.*( - 1 + 
         nf)*m11 );

      MQ11p6x11VR +=  + pow(abvm11,2) * (  - 1./2.*nfm + 1./2.*pow(
         nfm,2) + 1./8.*( - 1 + nf) );

      MQ11p6x11VR +=  + pi16*abvm11 * ( m11*nfm - m11*pow(nfm,2) - 1./4.
         *( - 1 + nf)*m11 );

      MQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * (  - m11*nfm + m11*pow(
         nfm,2) + 1./4.*( - 1 + nf)*m11 );

      MQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * (  - 2*nfm + 2*pow(nfm,2)
          + 1./2.*( - 1 + nf) );

  return  MQ11p6x11VR/pow(mass.getf0(),4);
}


//////////////// SPN /////////////////////////////////////////////////////
// mass SPN
double mnfSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double  MM12p4x11VR =
       + abvm11 * (  - 1./2.*m11 - 1./2.*m11*nfm );

  return  MM12p4x11VR/pow(mass.getf0(),2);
}

double mnfSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfSPNp6LV(nf,mass,Liin,L)+mnfSPNp6RV(nf,mass,L);
}

double mnfSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MM12p6x11VL =
       + bbvm11 * (  - 8*pow(m11,3)*L8r - 8*pow(m11,3)*L8r*nfm + 4*pow(
         m11,3)*L5r + 4*pow(m11,3)*L5r*nfm - 16*(1 + nf)*pow(m11,3)*L6r
          + 8*(1 + nf)*pow(m11,3)*L4r );

      MM12p6x11VL +=  + a23bvm11 * ( 12*m11*L3r*nfm - 24*m11*L1r + 12*
         m11*L0r*nfm - 12*( - 1 + nf)*m11*L3r - 12*( - 1 + 2*nf)*m11*
         L2r*nf - 12*( - 1 + 2*nf)*m11*L0r );

      MM12p6x11VL +=  + abvm11 * (  - 32*pow(m11,2)*L8r*nfm + 16*pow(
         m11,2)*L5r*nfm - 12*pow(m11,2)*L3r*nfm - 12*pow(m11,2)*L0r*nfm
          + 4*( - 3 + 2*nf)*pow(m11,2)*L0r + 4*( - 3 + 5*nf)*pow(m11,2)
         *L3r + 8*( - 1 - 2*nf + 4*pow(nf,2))*pow(m11,2)*L1r - 32*( - 1
          - nf + pow(nf,2))*pow(m11,2)*L4r + 32*( - 1 + nf)*pow(m11,2)*
         L8r + 32*( - 1 + nf)*pow(m11,2)*L6r*nf - 16*( - 1 + nf)*pow(
         m11,2)*L5r + 4*(4 - nf + 2*pow(nf,2))*pow(m11,2)*L2r );
 
  return MM12p6x11VL/pow(mass.getf0(),4);
}

double mnfSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MM12p6x11VR =
       + abvm11*bbvm11 * ( 1./4.*pow(m11,2) + 1./2.*pow(m11,2)*nfm + 1./
         4.*pow(m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + pow(abvm11,2) * ( 1./4.*m11*nfm + 1./8.*m11*
         pow(nfm,2) - 1./8.*( - 3 - 4*nf + 2*pow(nf,2))*m11 );

      MM12p6x11VR +=  + pi16*abvm11 * (  - 1./2.*pow(m11,2)*nfm - 1./4.
         *pow(m11,2)*pow(nfm,2) + 1./12.*(1 + 2*nf)*( - 3 + 2*nf)*pow(
         m11,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 1./4.*pow(m11,2) + 1./2.
         *pow(m11,2)*nfm + 1./4.*pow(m11,2)*pow(nfm,2) );

      MM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( m11*nfm + 1./2.*m11*pow(
         nfm,2) - 1./2.*( - 2 - 2*nf + pow(nf,2))*m11 );

      MM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * ( 1./2.*pow(
         m11,2)*nfm + 1./2.*pow(m11,2)*pow(nfm,2) + 1./12.*( - 2 - 9*nf
          + 3*pow(nf,2))*pow(m11,2) );

      MM12p6x11VR +=  + hh1V(1,m11,m11,m11,m11,L,mu2) * ( pow(m11,2)*nf
          );

      MM12p6x11VR +=  + hh21V(1,m11,m11,m11,m11,L,mu2) * ( 3./4.*( - 1
          + nf)*pow(m11,2)*nf );

      MM12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * (  - 3./4.*(
          - 1 + nf)*m11*nf );

  return  MM12p6x11VR/pow(mass.getf0(),4);
}
///////////////////////////////////////////////////////////////////////////
// decay SPN

double fnfSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double   MF12p4x11VR =
       + abvm11 * ( 1./2.*nf );

  return  MF12p4x11VR/pow(mass.getf0(),2);
}

double fnfSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfSPNp6LV(nf,mass,Liin,L)+fnfSPNp6RV(nf,mass,L);
}

double fnfSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MF12p6x11VL =
       + bbvm11 * ( 8*pow(m11,2)*L8r*nf + 16*pow(m11,2)*L6r*pow(nf,2)
          - 4*pow(m11,2)*L5r*nf - 8*pow(m11,2)*L4r*pow(nf,2) );

      MF12p6x11VL +=  + a23bvm11 * (  - 6*L3r*nfm + 12*L1r - 6*L0r*nfm
          + 6*( - 1 + nf)*L3r + 6*( - 1 + 2*nf)*L2r*nf + 6*( - 1 + 2*nf
         )*L0r );

      MF12p6x11VL +=  + abvm11 * (  - 2*m11*L5r*nfm + 6*m11*L3r*nfm + 6
         *m11*L0r*nfm - 2*( - 3 + 2*nf)*m11*L0r - 2*( - 3 + 5*nf)*m11*
         L3r - 4*( - 1 - 2*nf + 4*pow(nf,2))*m11*L1r + 4*( - 1 - nf + 
         pow(nf,2))*m11*L4r + 2*( - 1 + nf)*m11*L5r - 2*(4 - nf + 2*
         pow(nf,2))*m11*L2r );

  return MF12p6x11VL/pow(mass.getf0(),4);
}

double fnfSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MF12p6x11VR =
       + abvm11*bbvm11 * (  - 1./4.*(1 + nf)*m11 );

      MF12p6x11VR +=  + pow(abvm11,2) * (  - 1./8.*nf );

      MF12p6x11VR +=  + pi16*abvm11 * (  - 1./24.*( - 6 - 11*nf + 5*
         pow(nf,2))*m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * (  - 1./4.*(1 + nf)*m11 );

      MF12p6x11VR +=  + Ab(m11,mu2)*abvm11 * (  - 1./4.*(1 + 2*nf) );

      MF12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * (  - 1./8.*( - 1
          + nf)*m11*nf );

      MF12p6x11VR +=  + hh27V(1,m11,m11,m11,m11,L,mu2) * ( 3./8.*( - 1
          + nf)*nf );

      MF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * ( 1./4.*pow(
         m11,2)*nfm + 1./4.*pow(m11,2)*pow(nfm,2) + 1./24.*( - 2 - 9*nf
          + 3*pow(nf,2))*pow(m11,2) );

      MF12p6x11VR +=  + hh1dV(1,m11,m11,m11,m11,L,mu2) * ( 1./2.*pow(
         m11,2)*nf );

      MF12p6x11VR +=  + hh21dV(1,m11,m11,m11,m11,L,mu2) * ( 3./8.*( - 1
          + nf)*pow(m11,2)*nf );

      MF12p6x11VR +=  + hh27dV(1,m11,m11,m11,m11,L,mu2) * (  - 3./8.*(
          - 1 + nf)*m11*nf );

  return  MF12p6x11VR/pow(mass.getf0(),4);
}

///////////////////////////////////////////////////////////////////////////
// vev SPN

double qnfSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double  MQ11p4x11VR =
       + abvm11 * (  - 1./2.*nfm + 1./2.*( - 1 + 2*nf) );

  return  MQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfSPNp6LV(nf,mass,Liin,L)+qnfSPNp6RV(nf,mass,L);
}

double qnfSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  //double a23bvm11 = A23bV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MQ11p6x11VL =
       + bbvm11 * (  - 8*pow(m11,2)*L8r*nfm + 4*pow(m11,2)*L5r*nfm + 16
         *( - 1 + nf)*(1 + 2*nf)*pow(m11,2)*L6r - 8*( - 1 + nf)*(1 + 2*
         nf)*pow(m11,2)*L4r + 8*( - 1 + 2*nf)*pow(m11,2)*L8r - 4*( - 1
          + 2*nf)*pow(m11,2)*L5r );

      MQ11p6x11VL +=  + abvm11 * (  - 16*m11*L8r*nfm + 8*m11*L5r*nfm + 
         32*( - 1 + nf)*(1 + 2*nf)*m11*L6r - 16*( - 1 + nf)*(1 + 2*nf)*
         m11*L4r + 16*( - 1 + 2*nf)*m11*L8r - 8*( - 1 + 2*nf)*m11*L5r )
         ;
 
  return MQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   MQ11p6x11VR =
       + abvm11*bbvm11 * ( 1./2.*m11*nfm + 1./4.*m11*pow(nfm,2) - 1./4.
         *(1 + 2*nf)*m11 );

      MQ11p6x11VR +=  + pow(abvm11,2) * ( 1./4.*nfm + 1./8.*pow(nfm,2)
          - 1./8.*(1 + 2*nf) );

      MQ11p6x11VR +=  + pi16*abvm11 * (  - 1./2.*m11*nfm - 1./4.*m11*
         pow(nfm,2) + 1./4.*(1 + 2*nf)*m11 );

      MQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 1./2.*m11*nfm + 1./4.*
         m11*pow(nfm,2) - 1./4.*(1 + 2*nf)*m11 );

      MQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( nfm + 1./2.*pow(nfm,2)
          - 1./2.*(1 + 2*nf) );

  return  MQ11p6x11VR/pow(mass.getf0(),4);
}

