// massdecayvevnPQfV.cc is part of the CHIRON ChPT at two loops program collection
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
#include "massdecayvevnfPQV.h"
#include <cmath>

const double pi16 = 1./(16.*M_PI*M_PI);

// the code that produces results with Bessel method or theta method
#ifdef CHIRONTHETA
#define AbV AbVt
#define BbV BbVt
#define CbV CbVt
#define A23bV A23bVt
#define B23bV B23bVt
#define hhV hhVt
#define hh1V hh1Vt
#define hh21V hh21Vt
#define hh27V hh27Vt
#define hhdV hhdVt
#define hh1dV hh1dVt
#define hh21dV hh21dVt
#define hh27dV hh27dVt
#define mnfPQSUNp4V mnfPQSUNp4Vt
#define mnfPQSUNp6V mnfPQSUNp6Vt
#define mnfPQSUNp6LV mnfPQSUNp6LVt
#define mnfPQSUNp6RV mnfPQSUNp6RVt
#define fnfPQSUNp4V fnfPQSUNp4Vt
#define fnfPQSUNp6V fnfPQSUNp6Vt
#define fnfPQSUNp6LV fnfPQSUNp6LVt
#define fnfPQSUNp6RV fnfPQSUNp6RVt
#define qnfPQSUNp4V qnfPQSUNp4Vt
#define qnfPQSUNp6V qnfPQSUNp6Vt
#define qnfPQSUNp6LV qnfPQSUNp6LVt
#define qnfPQSUNp6RV qnfPQSUNp6RVt
#define mnfPQSONp4V mnfPQSONp4Vt
#define mnfPQSONp6V mnfPQSONp6Vt
#define mnfPQSONp6LV mnfPQSONp6LVt
#define mnfPQSONp6RV mnfPQSONp6RVt
#define fnfPQSONp4V fnfPQSONp4Vt
#define fnfPQSONp6V fnfPQSONp6Vt
#define fnfPQSONp6LV fnfPQSONp6LVt
#define fnfPQSONp6RV fnfPQSONp6RVt
#define qnfPQSONp4V qnfPQSONp4Vt
#define qnfPQSONp6V qnfPQSONp6Vt
#define qnfPQSONp6LV qnfPQSONp6LVt
#define qnfPQSONp6RV qnfPQSONp6RVt
#define mnfPQSPNp4V mnfPQSPNp4Vt
#define mnfPQSPNp6V mnfPQSPNp6Vt
#define mnfPQSPNp6LV mnfPQSPNp6LVt
#define mnfPQSPNp6RV mnfPQSPNp6RVt
#define fnfPQSPNp4V fnfPQSPNp4Vt
#define fnfPQSPNp6V fnfPQSPNp6Vt
#define fnfPQSPNp6LV fnfPQSPNp6LVt
#define fnfPQSPNp6RV fnfPQSPNp6RVt
#define qnfPQSPNp4V qnfPQSPNp4Vt
#define qnfPQSPNp6V qnfPQSPNp6Vt
#define qnfPQSPNp6LV qnfPQSPNp6LVt
#define qnfPQSPNp6RV qnfPQSPNp6RVt
#else
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef CHIRONBESSEL
#define AbV AbVb
#define BbV BbVb
#define CbV CbVb
#define A23bV A23bVb
#define B23bV B23bVb
#define hhV hhVb
#define hh1V hh1Vb
#define hh21V hh21Vb
#define hh27V hh27Vb
#define hhdV hhdVb
#define hh1dV hh1dVb
#define hh21dV hh21dVb
#define hh27dV hh27dVb
#define mnfPQSUNp4V mnfPQSUNp4Vb
#define mnfPQSUNp6V mnfPQSUNp6Vb
#define mnfPQSUNp6LV mnfPQSUNp6LVb
#define mnfPQSUNp6RV mnfPQSUNp6RVb
#define fnfPQSUNp4V fnfPQSUNp4Vb
#define fnfPQSUNp6V fnfPQSUNp6Vb
#define fnfPQSUNp6LV fnfPQSUNp6LVb
#define fnfPQSUNp6RV fnfPQSUNp6RVb
#define qnfPQSUNp4V qnfPQSUNp4Vb
#define qnfPQSUNp6V qnfPQSUNp6Vb
#define qnfPQSUNp6LV qnfPQSUNp6LVb
#define qnfPQSUNp6RV qnfPQSUNp6RVb
#define mnfPQSONp4V mnfPQSONp4Vb
#define mnfPQSONp6V mnfPQSONp6Vb
#define mnfPQSONp6LV mnfPQSONp6LVb
#define mnfPQSONp6RV mnfPQSONp6RVb
#define fnfPQSONp4V fnfPQSONp4Vb
#define fnfPQSONp6V fnfPQSONp6Vb
#define fnfPQSONp6LV fnfPQSONp6LVb
#define fnfPQSONp6RV fnfPQSONp6RVb
#define qnfPQSONp4V qnfPQSONp4Vb
#define qnfPQSONp6V qnfPQSONp6Vb
#define qnfPQSONp6LV qnfPQSONp6LVb
#define qnfPQSONp6RV qnfPQSONp6RVb
#define mnfPQSPNp4V mnfPQSPNp4Vb
#define mnfPQSPNp6V mnfPQSPNp6Vb
#define mnfPQSPNp6LV mnfPQSPNp6LVb
#define mnfPQSPNp6RV mnfPQSPNp6RVb
#define fnfPQSPNp4V fnfPQSPNp4Vb
#define fnfPQSPNp6V fnfPQSPNp6Vb
#define fnfPQSPNp6LV fnfPQSPNp6LVb
#define fnfPQSPNp6RV fnfPQSPNp6RVb
#define qnfPQSPNp4V qnfPQSPNp4Vb
#define qnfPQSPNp6V qnfPQSPNp6Vb
#define qnfPQSPNp6LV qnfPQSPNp6LVb
#define qnfPQSPNp6RV qnfPQSPNp6RVb
#else
// just some garbage to produce an error
x = 1./0.0.;
#endif
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

////////////////////////////////////////////////////////////////////////////
///////////////////// SUN //////////////////////////////////////////////////
// mass SUN
double mnfPQSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double PQM12p4x11VR =
       + bbvm11 * ( (m44 - m11)*m11*nfm );

      PQM12p4x11VR +=  + abvm11 * (  - m11*nfm );

  return  PQM12p4x11VR/pow(mass.getf0(),2);
}

double mnfPQSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfPQSUNp6LV(nf,mass,Liin,L)+mnfPQSUNp6RV(nf,mass,L);
}

double mnfPQSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double  PQM12p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*pow(m11,2)*m44*L6r - 16*(m44 - m11)*
         pow(m11,2)*m44*L4r + 32*(m44 - m11)*pow(m11,3)*L8r*nfm - 16*(
         m44 - m11)*pow(m11,3)*L5r*nfm );

      PQM12p6x11VL +=  + b23bvm11 * (  - 24*(m44 - m11)*m11*L3r*nfm - 
         24*(m44 - m11)*m11*L0r*nfm );

      PQM12p6x11VL +=  + bbvm11 * ( 24*(m44 - m11)*pow(m11,2)*L3r*nfm
          + 24*(m44 - m11)*pow(m11,2)*L0r*nfm + 16*(m44 - m11)*(m44 - 
         m11)*m11*L7r + 16*(2*m44 - 3*m11)*m11*m44*L6r - 8*(3*m44 - 4*
         m11)*m11*m44*L4r - 8*(6*m44 - 7*m11)*pow(m11,2)*L5r*nfm + 16*(
         pow(m44,2) + 4*m11*m44 - 6*pow(m11,2))*m11*L8r*nfm );

      PQM12p6x11VL +=  + a23bvm44 * (  - 12*( - 1 + nf)*(1 + nf)*m11*
         L2r );

      PQM12p6x11VL +=  + a23bvm14 * (  - 12*m11*L3r*nf - 24*m11*L0r*nf
          );

      PQM12p6x11VL +=  + a23bvm11 * ( 24*m11*L3r*nfm - 12*m11*L2r - 24*
         m11*L1r + 24*m11*L0r*nfm );

      PQM12p6x11VL +=  + abvm44 * ( 16*( - 1 + nf)*(1 + nf)*m11*m44*L6r
          - 16*( - 1 + nf)*(1 + nf)*m11*m44*L4r + 4*( - 1 + nf)*(1 + nf
         )*m11*m44*L2r + 16*( - 1 + nf)*(1 + nf)*m11*m44*L1r );

      PQM12p6x11VL +=  + abvm14 * ( 16*(m44 + m11)*m11*L8r*nf - 8*(m44
          + m11)*m11*L5r*nf + 10*(m44 + m11)*m11*L3r*nf + 4*(m44 + m11)
         *m11*L0r*nf );

      PQM12p6x11VL +=  + abvm11 * ( 16*m11*m44*L4r - 64*pow(m11,2)*L8r*
         nfm + 20*pow(m11,2)*L2r + 8*pow(m11,2)*L1r - 16*(m44 - 3*m11)*
         m11*L5r*nfm - 16*(m44 - 2*m11)*m11*L6r + 24*(m44 - 2*m11)*m11*
         L3r*nfm + 24*(m44 - 2*m11)*m11*L0r*nfm - 32*(m44 - m11)*m11*
         L7r );
 
  return PQM12p6x11VL/pow(mass.getf0(),4);
}

double mnfPQSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double cbvm11 = CbV(m11,L);
  double  PQM12p6x11VR =
       + bbvm11*cbvm11 * ( 2*(m44 - m11)*(m44 - m11)*pow(m11,2)*pow(
         nfm,2) );

      PQM12p6x11VR +=  + pow(bbvm11,2) * ( 1./2.*(m44 - m11)*(m44 - 5*
         m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm44*bbvm11 * ( m11*m44 - m11*m44*pow(nfm,2)
          );

      PQM12p6x11VR +=  + abvm14*bbvm11 * (  - 2*m11*m44 );

      PQM12p6x11VR +=  + pow(abvm14,2) * (  - 1./4.*m11*pow(nf,2) );

      PQM12p6x11VR +=  + abvm11*cbvm11 * (  - 2*(m44 - m11)*pow(m11,2)*
         pow(nfm,2) );

      PQM12p6x11VR +=  + abvm11*bbvm11 * ( pow(m11,2) - (m44 - 3*m11)*
         m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pow(abvm11,2) * ( 1./2.*m11 + 1./2.*m11*pow(
         nfm,2) );

      PQM12p6x11VR +=  + pi16*cbvm11 * (  - 2*(m44 - m11)*(m44 - m11)*
         pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*bbvm11 * (  - 2*(m44 - m11)*(m44 - 3*m11)
         *m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm44 * (  - m11*m44 + m11*m44*pow(
         nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm14 * ( 1./3.*(6 + pow(nf,2))*m11*m44
          );

      PQM12p6x11VR +=  + pi16*abvm11 * (  - pow(m11,2) + 2*(m44 - 2*m11
         )*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*cbvm11 * ( 2*(m44 - m11)*(m44 - 2*
         m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( pow(m11,2) + (pow(
         m44,2) - 7*m11*m44 + 8*pow(m11,2))*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( m44 - m44*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm14 * (  - 2*m44 );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 2*m11 - (m44 - 4*m11)*
         pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - 2*m11*m44 );

      PQM12p6x11VR +=  + Ab(m14,mu2)*abvm14 * (  - 1./2.*m11*pow(nf,2)
          );

      PQM12p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( m11*m44 - m11*m44*pow(
         nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * ( 1./3.*pow(
         m11,2) + 2*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * ( 1./4.*(m44 - 
         4*m11)*m11 );

      PQM12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * ( 1./4.*( - 1
          + nf)*(1 + nf)*m11*m44 );

      PQM12p6x11VR +=  + hhV(2,m11,m11,m11,m11,L,mu2) * (  - 4*(m44 - 
         m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * ( 5./4.*(m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hhV(5,m11,m11,m11,m11,L,mu2) * ( 2*(m44 - m11)
         *(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * (  - 2*(m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m11,m14,m14,m11,L,mu2) * ( 3./4.*pow(
         m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m44,m14,m14,m11,L,mu2) * ( 3./4.*( - 1
          + nf)*(1 + nf)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(2,m11,m14,m14,m11,L,mu2) * ( 3./4.*(m44
          - m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * (  - 3./4.*
         m11 );

      PQM12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * (  - 3./4.*(
          - 1 + nf)*(1 + nf)*m11 );

      PQM12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * (  - 3./4.*(
         m44 - m11)*m11 );

  return  PQM12p6x11VR/pow(mass.getf0(),4);
}
////////////////////////////////////////////////////////////////////////
// decay SUN
double fnfPQSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double abvm14 = AbV(m14,L);
  //double bbvm11 = BbV(m11,L);
  double PQF12p4x11VR =
       + abvm14 * ( 1./2.*nf );

  return  PQF12p4x11VR/pow(mass.getf0(),2);
}

double fnfPQSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfPQSUNp6LV(nf,mass,Liin,L)+fnfPQSUNp6RV(nf,mass,L);
}

double fnfPQSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double b23bvm11 = B23bV(m11,L);
  //double cbvm11 = CbV(m11,L);
  double  PQF12p6x11VL =
       + b23bvm11 * ( 12*(m44 - m11)*L3r*nfm + 12*(m44 - m11)*L0r*nfm )
         ;

      PQF12p6x11VL +=  + bbvm14 * ( 4*(m44 + m11)*m44*L6r*pow(nf,2) - 2
         *(m44 + m11)*m44*L4r*pow(nf,2) + 2*(m44 + m11)*(m44 + m11)*L8r
         *nf - (m44 + m11)*(m44 + m11)*L5r*nf );

      PQF12p6x11VL +=  + bbvm11 * ( 4*(m44 - m11)*m11*L5r*nfm - 12*(m44
          - m11)*m11*L3r*nfm - 12*(m44 - m11)*m11*L0r*nfm );

      PQF12p6x11VL +=  + a23bvm44 * ( 6*( - 1 + nf)*(1 + nf)*L2r );

      PQF12p6x11VL +=  + a23bvm14 * ( 6*L3r*nf + 12*L0r*nf );

      PQF12p6x11VL +=  + a23bvm11 * (  - 12*L3r*nfm + 6*L2r + 12*L1r - 
         12*L0r*nfm );

      PQF12p6x11VL +=  + abvm44 * ( 4*( - 1 + nf)*(1 + nf)*m44*L4r - 2*
         ( - 1 + nf)*(1 + nf)*m44*L2r - 8*( - 1 + nf)*(1 + nf)*m44*L1r
          );

      PQF12p6x11VL +=  + abvm14 * (  - 2*m44*L4r*pow(nf,2) + 2*m11*L5r*
         nf - 5*(m44 + m11)*L3r*nf - 2*(m44 + m11)*L0r*nf );

      PQF12p6x11VL +=  + abvm11 * (  - 4*m11*L5r*nfm - 10*m11*L2r - 4*
         m11*L1r - 12*(m44 - 2*m11)*L3r*nfm - 12*(m44 - 2*m11)*L0r*nfm
          );
 
  return PQF12p6x11VL/pow(mass.getf0(),4);
}

double fnfPQSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  //double abvm44 = AbV(m44,L);
  //double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double cbvm11 = CbV(m11,L);
  double  PQF12p6x11VR =
       + abvm11*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQF12p6x11VR +=  + pi16*abvm14 * (  - 1./24.*(4*m44 + m11)*pow(
         nf,2) );

      PQF12p6x11VR +=  + pi16*abvm11 * ( 1./4.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./4.*(m44 + m11)*
         pow(m14,-1) );

      PQF12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * (  - 1./8.*m44
          );

      PQF12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * (  - 1./8.*( - 
         1 + nf)*(1 + nf)*m44 );

      PQF12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * (  - 5./36.*(
         m44 - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * ( 1./72.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(3,m14,m11,m14,m11,L,mu2) * ( 1./36.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * ( 3./8. );

      PQF12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * ( 3./8.*( - 1
          + nf)*(1 + nf) );

      PQF12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * ( 3./8.*(m44
          - m11) );

      PQF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * ( 1./6.*pow(
         m11,2) + pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hhdV(1,m11,m14,m14,m11,L,mu2) * ( 1./8.*(m44
          - 4*m11)*m11 );

      PQF12p6x11VR +=  + hhdV(1,m44,m14,m14,m11,L,mu2) * ( 1./8.*( - 1
          + nf)*(1 + nf)*m11*m44 );

      PQF12p6x11VR +=  + hhdV(2,m11,m11,m11,m11,L,mu2) * (  - 2*(m44 - 
         m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hhdV(2,m11,m14,m14,m11,L,mu2) * ( 5./8.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hhdV(5,m11,m11,m11,m11,L,mu2) * ( (m44 - m11)*
         (m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hh1dV(2,m11,m14,m14,m11,L,mu2) * (  - (m44 - 
         m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m11,m14,m14,m11,L,mu2) * ( 3./8.*pow(
         m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m44,m14,m14,m11,L,mu2) * ( 3./8.*( - 
         1 + nf)*(1 + nf)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(2,m11,m14,m14,m11,L,mu2) * ( 3./8.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh27dV(1,m11,m14,m14,m11,L,mu2) * (  - 3./8.*
         m11 );

      PQF12p6x11VR +=  + hh27dV(1,m44,m14,m14,m11,L,mu2) * (  - 3./8.*(
          - 1 + nf)*(1 + nf)*m11 );

      PQF12p6x11VR +=  + hh27dV(2,m11,m14,m14,m11,L,mu2) * (  - 3./8.*(
         m44 - m11)*m11 );

  return  PQF12p6x11VR/pow(mass.getf0(),4);
}

////////////////////////////////////////////////////////////////////////
// vev SUN
double qnfPQSUNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double bbvm11 = BbV(m11,L);
  double   PQQ11p4x11VR =
       + bbvm11 * ( (m44 - m11)*nfm );

      PQQ11p4x11VR +=  + abvm14 * ( nf );

      PQQ11p4x11VR +=  + abvm11 * (  - nfm );

  return  PQQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfPQSUNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfPQSUNp6LV(nf,mass,Liin,L)+qnfPQSUNp6RV(nf,mass,L);
}

double qnfPQSUNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  //double a23bvm11 = A23bV(m11,L);
  //double a23bvm14 = A23bV(m14,L);
  //double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQQ11p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*m11*m44*L6r - 16*(m44 - m11)*m11*m44
         *L4r + 32*(m44 - m11)*pow(m11,2)*L8r*nfm - 16*(m44 - m11)*pow(
         m11,2)*L5r*nfm );

      PQQ11p6x11VL +=  + bbvm14 * ( 8*(m44 + m11)*m44*L6r*pow(nf,2) - 4
         *(m44 + m11)*m44*L4r*pow(nf,2) + 4*(m44 + m11)*(m44 + m11)*L8r
         *nf - 2*(m44 + m11)*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + bbvm11 * ( 16*(m44 - m11)*(m44 - m11)*L7r + 16
         *(2*m44 - 3*m11)*m44*L6r - 8*(2*m44 - 3*m11)*m44*L4r - 8*(4*
         m44 - 5*m11)*m11*L5r*nfm + 16*(pow(m44,2) + 2*m11*m44 - 4*pow(
         m11,2))*L8r*nfm );

      PQQ11p6x11VL +=  + abvm44 * ( 16*( - 1 + nf)*(1 + nf)*m44*L6r - 8
         *( - 1 + nf)*(1 + nf)*m44*L4r );

      PQQ11p6x11VL +=  + abvm14 * ( 16*m44*L6r*pow(nf,2) - 8*m44*L4r*
         pow(nf,2) + 16*(m44 + m11)*L8r*nf - 8*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + abvm11 * (  - 16*m44*L6r + 8*m44*L4r - 32*m11*
         L8r*nfm - 16*(m44 - 2*m11)*L5r*nfm - 32*(m44 - m11)*L7r );
 
  return PQQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfPQSUNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double cbvm11 = CbV(m11,L);
  double  PQQ11p6x11VR =
       + bbvm11*cbvm11 * ( 2*(m44 - m11)*(m44 - m11)*m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + pow(bbvm11,2) * ( 1./2.*(m44 - m11)*(m44 - 5*
         m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm44*bbvm11 * ( m44 - m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm14*bbvm11 * (  - (m44 + m11) );

      PQQ11p6x11VR +=  + abvm11*cbvm11 * (  - 2*(m44 - m11)*m11*pow(
         nfm,2) );

      PQQ11p6x11VR +=  + abvm11*bbvm14 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + abvm11*bbvm11 * ( m11 - (m44 - 3*m11)*pow(
         nfm,2) );

      PQQ11p6x11VR +=  + abvm11*abvm14 * (  - 1 );

      PQQ11p6x11VR +=  + pow(abvm11,2) * ( 1./2. + 1./2.*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*cbvm11 * (  - 2*(m44 - m11)*(m44 - m11)*
         m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*bbvm11 * (  - 2*(m44 - m11)*(m44 - 3*m11)
         *pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm44 * (  - m44 + m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm14 * ( (m44 + m11) );

      PQQ11p6x11VR +=  + pi16*abvm11 * ( 2*(m44 - 2*m11)*pow(nfm,2) + 1.
         /2.*(m44 - m11) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*cbvm11 * ( 2*(m44 - m11)*(m44 - 2*
         m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( pow(m11,-1)*pow(m44,2)*
         pow(nfm,2) + m11 - (7*m44 - 8*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( pow(m11,-1)*m44 - pow(
         m11,-1)*m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm14 * (  - 2 - pow(m11,-1)*m44
          );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 2 - pow(m11,-1)*m44*
         pow(nfm,2) + 4*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - (m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1 - 1./2.*(m44 + m11
         )*pow(m14,-1) );

      PQQ11p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( m44 - m44*pow(nfm,2) );

  return  PQQ11p6x11VR/pow(mass.getf0(),4);
}

///////////////////// SON //////////////////////////////////////////////////
// mass SON
double mnfPQSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   PQM12p4x11VR =
       + bbvm11 * ( (m44 - m11)*m11*nfm );

      PQM12p4x11VR +=  + abvm11 * ( 1./2.*m11 - m11*nfm );

  return  PQM12p4x11VR/pow(mass.getf0(),2);
}

double mnfPQSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfPQSONp6LV(nf,mass,Liin,L)+mnfPQSONp6RV(nf,mass,L);
}

double mnfPQSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQM12p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*pow(m11,2)*m44*L6r - 16*(m44 - m11)*
         pow(m11,2)*m44*L4r + 32*(m44 - m11)*pow(m11,3)*L8r*nfm - 16*(
         m44 - m11)*pow(m11,3)*L5r*nfm );

      PQM12p6x11VL +=  + b23bvm11 * (  - 24*(m44 - m11)*m11*L3r*nfm - 
         24*(m44 - m11)*m11*L0r*nfm );

      PQM12p6x11VL +=  + bbvm11 * ( 8*pow(m11,3)*L8r - 4*pow(m11,3)*L5r
          + 24*(m44 - m11)*pow(m11,2)*L3r*nfm + 24*(m44 - m11)*pow(
         m11,2)*L0r*nfm + 16*(m44 - m11)*(m44 - m11)*m11*L7r + 8*(4*m44
          - 6*m11 + m11*nf)*m11*m44*L6r - 4*(6*m44 - 8*m11 + m11*nf)*
         m11*m44*L4r - 8*(6*m44 - 7*m11)*pow(m11,2)*L5r*nfm + 16*(pow(
         m44,2) + 4*m11*m44 - 6*pow(m11,2))*m11*L8r*nfm );

      PQM12p6x11VL +=  + a23bvm44 * (  - 6*( - 1 + nf)*(2 + nf)*m11*L2r
          );

      PQM12p6x11VL +=  + a23bvm14 * (  - 6*m11*L3r*nf - 12*m11*L0r*nf )
         ;

      PQM12p6x11VL +=  + a23bvm11 * (  - 12*m11*L3r + 24*m11*L3r*nfm - 
         12*m11*L2r - 24*m11*L1r - 12*m11*L0r + 24*m11*L0r*nfm );

      PQM12p6x11VL +=  + abvm44 * ( 8*( - 1 + nf)*(2 + nf)*m11*m44*L6r
          - 8*( - 1 + nf)*(2 + nf)*m11*m44*L4r + 2*( - 1 + nf)*(2 + nf)
         *m11*m44*L2r + 8*( - 1 + nf)*(2 + nf)*m11*m44*L1r );

      PQM12p6x11VL +=  + abvm14 * ( 8*(m44 + m11)*m11*L8r*nf - 4*(m44
          + m11)*m11*L5r*nf + 5*(m44 + m11)*m11*L3r*nf + 2*(m44 + m11)*
         m11*L0r*nf );

      PQM12p6x11VL +=  + abvm11 * ( 32*pow(m11,2)*L8r - 64*pow(m11,2)*
         L8r*nfm - 16*pow(m11,2)*L5r + 12*pow(m11,2)*L3r + 20*pow(
         m11,2)*L2r + 8*pow(m11,2)*L1r + 12*pow(m11,2)*L0r - 8*( - 2 + 
         nf)*m11*m44*L4r + 8*( - 2*m44 + m44*nf + 4*m11)*m11*L6r - 16*(
         m44 - 3*m11)*m11*L5r*nfm + 24*(m44 - 2*m11)*m11*L3r*nfm + 24*(
         m44 - 2*m11)*m11*L0r*nfm - 32*(m44 - m11)*m11*L7r );
 
  return PQM12p6x11VL/pow(mass.getf0(),4);
}

double mnfPQSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQM12p6x11VR =
       + bbvm11*cbvm11 * ( 2*(m44 - m11)*(m44 - m11)*pow(m11,2)*pow(
         nfm,2) );

      PQM12p6x11VR +=  + pow(bbvm11,2) * ( 1./2.*(m44 - m11)*pow(m11,2)
         *nfm + 1./2.*(m44 - m11)*(m44 - 5*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm44*bbvm11 * ( 1./2.*m11*m44 + 1./2.*m11*
         m44*nfm - m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm14*bbvm11 * (  - m11*m44 );

      PQM12p6x11VR +=  + pow(abvm14,2) * (  - 1./16.*m11*pow(nf,2) );

      PQM12p6x11VR +=  + abvm11*cbvm11 * ( (m44 - m11)*pow(m11,2)*nfm
          - 2*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm11*bbvm11 * ( 3./4.*pow(m11,2) + 1./2.*(
         m44 - 4*m11)*m11*nfm - (m44 - 3*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm11*abvm14 * (  - 1./4.*m11*nf );

      PQM12p6x11VR +=  + pow(abvm11,2) * ( 3./8.*m11 - 1./2.*m11*nfm + 
         1./2.*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*cbvm11 * (  - 2*(m44 - m11)*(m44 - m11)*
         pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*bbvm11 * (  - (m44 - m11)*pow(m11,2)*nfm
          - 2*(m44 - m11)*(m44 - 3*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm44 * (  - 1./2.*m11*m44 - 1./2.*m11*
         m44*nfm + m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm14 * ( 1./48.*(48*m44 + 5*m44*nf + 4*
         m44*pow(nf,2) + 3*m11*nf)*m11 );

      PQM12p6x11VR +=  + pi16*abvm11 * (  - 3./4.*pow(m11,2) + 2*(m44
          - 2*m11)*m11*pow(nfm,2) - 1./2.*(2*m44 - 5*m11)*m11*nfm );

      PQM12p6x11VR +=  + Ab(m11,mu2)*cbvm11 * ( (m44 - m11)*pow(m11,2)*
         nfm + 2*(m44 - m11)*(m44 - 2*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 3./4.*pow(m11,2) + 3./2.
         *(m44 - 2*m11)*m11*nfm + (pow(m44,2) - 7*m11*m44 + 8*pow(
         m11,2))*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( 1./2.*m44 + 1./2.*m44*
         nfm - m44*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm14 * (  - 1./4.*(4*m44 + m11*
         nf) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 3./2.*m11 + 1./2.*(m44
          - 6*m11)*nfm - (m44 - 4*m11)*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - m11*m44 );

      PQM12p6x11VR +=  + Ab(m14,mu2)*abvm14 * (  - 1./8.*m11*pow(nf,2)
          );

      PQM12p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./4.*m11*nf );

      PQM12p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( 1./2.*m11*m44 + 1./2.*
         m11*m44*nfm - m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * ( 1./3.*pow(
         m11,2) - pow(m11,2)*nfm + 2*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * ( 1./16.*(2*m44
          - 8*m11 + 5*m11*nf)*m11 );

      PQM12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * ( 1./16.*( - 1
          + nf)*(2 + nf)*m11*m44 );

      PQM12p6x11VR +=  + hhV(2,m11,m11,m11,m11,L,mu2) * ( (m44 - m11)*
         pow(m11,2)*nfm - 4*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * ( 5./8.*(m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hhV(5,m11,m11,m11,m11,L,mu2) * ( 2*(m44 - m11)
         *(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hh1V(1,m11,m14,m14,m11,L,mu2) * (  - 1./2.*
         pow(m11,2)*nf );

      PQM12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * (  - (m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m11,m14,m14,m11,L,mu2) * ( 3./16.*(2
          + nf)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m44,m14,m14,m11,L,mu2) * ( 3./16.*( - 
         1 + nf)*(2 + nf)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(2,m11,m14,m14,m11,L,mu2) * ( 3./8.*(m44
          - m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * (  - 3./16.*(
         2 + nf)*m11 );

      PQM12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * (  - 3./16.*(
          - 1 + nf)*(2 + nf)*m11 );

      PQM12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * (  - 3./8.*(
         m44 - m11)*m11 );

  return  PQM12p6x11VR/pow(mass.getf0(),4);
}
////////////////////////////////////////////////////////////////////////
// decay SON

double fnfPQSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double abvm14 = AbV(m14,L);
  //double bbvm11 = BbV(m11,L);
  double  PQF12p4x11VR =
       + abvm14 * ( 1./4.*nf );

  return  PQF12p4x11VR/pow(mass.getf0(),2);
}

double fnfPQSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfPQSONp6LV(nf,mass,Liin,L)+fnfPQSONp6RV(nf,mass,L);
}

double fnfPQSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double b23bvm11 = B23bV(m11,L);
  //double cbvm11 = CbV(m11,L);
  double   PQF12p6x11VL =
       + b23bvm11 * ( 12*(m44 - m11)*L3r*nfm + 12*(m44 - m11)*L0r*nfm )
         ;

      PQF12p6x11VL +=  + bbvm14 * ( 2*(m44 + m11)*m44*L6r*pow(nf,2) - (
         m44 + m11)*m44*L4r*pow(nf,2) + (m44 + m11)*(m44 + m11)*L8r*nf
          - 1./2.*(m44 + m11)*(m44 + m11)*L5r*nf );

      PQF12p6x11VL +=  + bbvm11 * ( 4*(m44 - m11)*m11*L5r*nfm - 12*(m44
          - m11)*m11*L3r*nfm - 12*(m44 - m11)*m11*L0r*nfm );

      PQF12p6x11VL +=  + a23bvm44 * ( 3*( - 1 + nf)*(2 + nf)*L2r );

      PQF12p6x11VL +=  + a23bvm14 * ( 3*L3r*nf + 6*L0r*nf );

      PQF12p6x11VL +=  + a23bvm11 * ( 6*L3r - 12*L3r*nfm + 6*L2r + 12*
         L1r + 6*L0r - 12*L0r*nfm );

      PQF12p6x11VL +=  + abvm44 * ( 2*( - 1 + nf)*(2 + nf)*m44*L4r - (
          - 1 + nf)*(2 + nf)*m44*L2r - 4*( - 1 + nf)*(2 + nf)*m44*L1r )
         ;

      PQF12p6x11VL +=  + abvm14 * (  - m44*L4r*pow(nf,2) + m11*L5r*nf
          - 5./2.*(m44 + m11)*L3r*nf - (m44 + m11)*L0r*nf );

      PQF12p6x11VL +=  + abvm11 * ( 2*m11*L5r - 4*m11*L5r*nfm - 6*m11*
         L3r - 10*m11*L2r - 4*m11*L1r - 6*m11*L0r - 12*(m44 - 2*m11)*
         L3r*nfm - 12*(m44 - 2*m11)*L0r*nfm );

  return PQF12p6x11VL/pow(mass.getf0(),4);
}

double fnfPQSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  //double abvm44 = AbV(m44,L);
  //double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double cbvm11 = CbV(m11,L);
  double   PQF12p6x11VR =
       + abvm14*bbvm14 * ( 1./16.*(m44 + m11)*nf );

      PQF12p6x11VR +=  + pow(abvm14,2) * ( 1./16.*nf );

      PQF12p6x11VR +=  + abvm11*bbvm14 * (  - 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + pi16*abvm14 * (  - 1./96.*(11*m44 + 4*m44*nf
          + 11*m11 + m11*nf)*nf );

      PQF12p6x11VR +=  + pi16*abvm11 * ( 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m14,mu2)*bbvm14 * ( 1./16.*(m44 + m11)*nf )
         ;

      PQF12p6x11VR +=  + Ab(m14,mu2)*abvm14 * ( 1./8.*nf + 1./16.*(m44
          + m11)*pow(m14,-1)*nf );

      PQF12p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./8.*(m44 + m11)*
         pow(m14,-1) );

      PQF12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * (  - 1./32.*(2*
         m44 + m11*nf) );

      PQF12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * (  - 1./32.*(
          - 1 + nf)*(2 + nf)*m44 );

      PQF12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * (  - 5./72.*(
         m44 - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * ( 1./144.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(3,m14,m11,m14,m11,L,mu2) * ( 1./72.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * ( 3./32.*(2
          + nf) );

      PQF12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * ( 3./32.*( - 
         1 + nf)*(2 + nf) );

      PQF12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * ( 3./16.*(m44
          - m11) );

      PQF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * ( 1./6.*pow(
         m11,2) - 1./2.*pow(m11,2)*nfm + pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hhdV(1,m11,m14,m14,m11,L,mu2) * ( 1./32.*(2*
         m44 - 8*m11 + 5*m11*nf)*m11 );

      PQF12p6x11VR +=  + hhdV(1,m44,m14,m14,m11,L,mu2) * ( 1./32.*( - 1
          + nf)*(2 + nf)*m11*m44 );

      PQF12p6x11VR +=  + hhdV(2,m11,m11,m11,m11,L,mu2) * ( 1./2.*(m44
          - m11)*pow(m11,2)*nfm - 2*(m44 - m11)*pow(m11,2)*pow(nfm,2) )
         ;

      PQF12p6x11VR +=  + hhdV(2,m11,m14,m14,m11,L,mu2) * ( 5./16.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hhdV(5,m11,m11,m11,m11,L,mu2) * ( (m44 - m11)*
         (m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hh1dV(1,m11,m14,m14,m11,L,mu2) * (  - 1./4.*
         pow(m11,2)*nf );

      PQF12p6x11VR +=  + hh1dV(2,m11,m14,m14,m11,L,mu2) * (  - 1./2.*(
         m44 - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m11,m14,m14,m11,L,mu2) * ( 3./32.*(2
          + nf)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m44,m14,m14,m11,L,mu2) * ( 3./32.*(
          - 1 + nf)*(2 + nf)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(2,m11,m14,m14,m11,L,mu2) * ( 3./16.*(
         m44 - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh27dV(1,m11,m14,m14,m11,L,mu2) * (  - 3./32.*
         (2 + nf)*m11 );

      PQF12p6x11VR +=  + hh27dV(1,m44,m14,m14,m11,L,mu2) * (  - 3./32.*
         ( - 1 + nf)*(2 + nf)*m11 );

      PQF12p6x11VR +=  + hh27dV(2,m11,m14,m14,m11,L,mu2) * (  - 3./16.*
         (m44 - m11)*m11 );

  return  PQF12p6x11VR/pow(mass.getf0(),4);
}

////////////////////////////////////////////////////////////////////////
// vev SON

double qnfPQSONp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double bbvm11 = BbV(m11,L);
  double   PQQ11p4x11VR =
       + bbvm11 * ( (m44 - m11)*nfm );

      PQQ11p4x11VR +=  + abvm14 * ( 1./2.*nf );

      PQQ11p4x11VR +=  + abvm11 * ( 1./2. - nfm );

  return  PQQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfPQSONp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfPQSONp6LV(nf,mass,Liin,L)+qnfPQSONp6RV(nf,mass,L);
}

double qnfPQSONp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  //double a23bvm11 = A23bV(m11,L);
  //double a23bvm14 = A23bV(m14,L);
  //double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQQ11p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*m11*m44*L6r - 16*(m44 - m11)*m11*m44
         *L4r + 32*(m44 - m11)*pow(m11,2)*L8r*nfm - 16*(m44 - m11)*pow(
         m11,2)*L5r*nfm );

      PQQ11p6x11VL +=  + bbvm14 * ( 4*(m44 + m11)*m44*L6r*pow(nf,2) - 2
         *(m44 + m11)*m44*L4r*pow(nf,2) + 2*(m44 + m11)*(m44 + m11)*L8r
         *nf - (m44 + m11)*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + bbvm11 * ( 8*pow(m11,2)*L8r - 4*pow(m11,2)*L5r
          + 16*(m44 - m11)*(m44 - m11)*L7r + 8*(4*m44 - 6*m11 + m11*nf)
         *m44*L6r - 4*(4*m44 - 6*m11 + m11*nf)*m44*L4r - 8*(4*m44 - 5*
         m11)*m11*L5r*nfm + 16*(pow(m44,2) + 2*m11*m44 - 4*pow(m11,2))*
         L8r*nfm );

      PQQ11p6x11VL +=  + abvm44 * ( 8*( - 1 + nf)*(2 + nf)*m44*L6r - 4*
         ( - 1 + nf)*(2 + nf)*m44*L4r );

      PQQ11p6x11VL +=  + abvm14 * ( 8*m44*L6r*pow(nf,2) - 4*m44*L4r*
         pow(nf,2) + 8*(m44 + m11)*L8r*nf - 4*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + abvm11 * ( 16*m11*L8r - 32*m11*L8r*nfm - 8*m11
         *L5r + 8*( - 2 + nf)*m44*L6r - 4*( - 2 + nf)*m44*L4r - 16*(m44
          - 2*m11)*L5r*nfm - 32*(m44 - m11)*L7r );
 
  return PQQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfPQSONp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double cbvm11 = CbV(m11,L);
  double  PQQ11p6x11VR =
       + bbvm11*cbvm11 * ( 2*(m44 - m11)*(m44 - m11)*m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + pow(bbvm11,2) * ( 1./2.*(m44 - m11)*m11*nfm + 
         1./2.*(m44 - m11)*(m44 - 5*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm44*bbvm11 * ( 1./2.*m44 + 1./2.*m44*nfm - 
         m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm14*bbvm14 * ( 1./8.*(m44 + m11)*nf );

      PQQ11p6x11VR +=  + abvm14*bbvm11 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + pow(abvm14,2) * ( 1./8.*nf );

      PQQ11p6x11VR +=  + abvm11*cbvm11 * ( (m44 - m11)*m11*nfm - 2*(m44
          - m11)*m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm11*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQQ11p6x11VR +=  + abvm11*bbvm11 * ( 3./4.*m11 + 1./2.*(m44 - 4*
         m11)*nfm - (m44 - 3*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm11*abvm14 * (  - 1./2. );

      PQQ11p6x11VR +=  + pow(abvm11,2) * ( 3./8. - 1./2.*nfm + 1./2.*
         pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*cbvm11 * (  - 2*(m44 - m11)*(m44 - m11)*
         m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*bbvm11 * (  - (m44 - m11)*m11*nfm - 2*(
         m44 - m11)*(m44 - 3*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm44 * (  - 1./2.*m44 - 1./2.*m44*nfm
          + m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm14 * (  - 1./8.*( - 4 + nf)*(m44 + 
         m11) );

      PQQ11p6x11VR +=  + pi16*abvm11 * ( 1./4.*(m44 - 2*m11) + 2*(m44
          - 2*m11)*pow(nfm,2) - 1./2.*(2*m44 - 5*m11)*nfm );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*cbvm11 * ( (m44 - m11)*m11*nfm + 2
         *(m44 - m11)*(m44 - 2*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( pow(m11,-1)*pow(m44,2)*
         pow(nfm,2) + 3./4.*m11 + 3./2.*(m44 - 2*m11)*nfm - (7*m44 - 8*
         m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( 1./2.*pow(m11,-1)*m44
          + 1./2.*pow(m11,-1)*m44*nfm - pow(m11,-1)*m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm14 * (  - 1 - 1./2.*pow(
         m11,-1)*m44 );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 3./2. + 1./2.*pow(
         m11,-1)*m44*nfm - pow(m11,-1)*m44*pow(nfm,2) - 3*nfm + 4*pow(
         nfm,2) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*bbvm14 * ( 1./8.*(m44 + m11)*nf );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*abvm14 * ( 1./4.*nf + 1./8.*(m44
          + m11)*pow(m14,-1)*nf );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./2. - 1./4.*(m44
          + m11)*pow(m14,-1) );

      PQQ11p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( 1./2.*m44 + 1./2.*m44*
         nfm - m44*pow(nfm,2) );

  return  PQQ11p6x11VR/pow(mass.getf0(),4);
}


///////////////////// SPN //////////////////////////////////////////////////
// mass SPN
double mnfPQSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double bbvm11 = BbV(m11,L);
  double   PQM12p4x11VR =
       + bbvm11 * ( 1./2.*(m44 - m11)*m11*nfm );

      PQM12p4x11VR +=  + abvm11 * (  - 1./2.*m11 - 1./2.*m11*nfm );

  return  PQM12p4x11VR/pow(mass.getf0(),2);
}

double mnfPQSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  mnfPQSPNp6LV(nf,mass,Liin,L)+mnfPQSPNp6RV(nf,mass,L);
}

double mnfPQSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQM12p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*pow(m11,2)*m44*L6r - 16*(m44 - m11)*
         pow(m11,2)*m44*L4r + 16*(m44 - m11)*pow(m11,3)*L8r*nfm - 8*(
         m44 - m11)*pow(m11,3)*L5r*nfm );

      PQM12p6x11VL +=  + b23bvm11 * (  - 12*(m44 - m11)*m11*L3r*nfm - 
         12*(m44 - m11)*m11*L0r*nfm );

      PQM12p6x11VL +=  + bbvm11 * (  - 8*pow(m11,3)*L8r + 4*pow(m11,3)*
         L5r + 8*( - 3*m44 + 4*m11 + m11*nf)*m11*m44*L4r - 16*( - 2*m44
          + 3*m11 + m11*nf)*m11*m44*L6r + 12*(m44 - m11)*pow(m11,2)*L3r
         *nfm + 12*(m44 - m11)*pow(m11,2)*L0r*nfm + 16*(m44 - m11)*(m44
          - m11)*m11*L7r - 4*(6*m44 - 7*m11)*pow(m11,2)*L5r*nfm + 8*(
         pow(m44,2) + 4*m11*m44 - 6*pow(m11,2))*m11*L8r*nfm );

      PQM12p6x11VL +=  + a23bvm44 * (  - 12*( - 1 + nf)*(1 + 2*nf)*m11*
         L2r );

      PQM12p6x11VL +=  + a23bvm14 * (  - 12*m11*L3r*nf - 24*m11*L0r*nf
          );

      PQM12p6x11VL +=  + a23bvm11 * ( 12*m11*L3r + 12*m11*L3r*nfm - 12*
         m11*L2r - 24*m11*L1r + 12*m11*L0r + 12*m11*L0r*nfm );

      PQM12p6x11VL +=  + abvm44 * ( 16*( - 1 + nf)*(1 + 2*nf)*m11*m44*
         L6r - 16*( - 1 + nf)*(1 + 2*nf)*m11*m44*L4r + 4*( - 1 + nf)*(1
          + 2*nf)*m11*m44*L2r + 16*( - 1 + nf)*(1 + 2*nf)*m11*m44*L1r )
         ;

      PQM12p6x11VL +=  + abvm14 * ( 16*(m44 + m11)*m11*L8r*nf - 8*(m44
          + m11)*m11*L5r*nf + 10*(m44 + m11)*m11*L3r*nf + 4*(m44 + m11)
         *m11*L0r*nf );

      PQM12p6x11VL +=  + abvm11 * (  - 32*pow(m11,2)*L8r - 32*pow(
         m11,2)*L8r*nfm + 16*pow(m11,2)*L5r - 12*pow(m11,2)*L3r + 20*
         pow(m11,2)*L2r + 8*pow(m11,2)*L1r - 12*pow(m11,2)*L0r + 16*(1
          + nf)*m11*m44*L4r - 16*(m44 + m44*nf - 2*m11)*m11*L6r - 8*(
         m44 - 3*m11)*m11*L5r*nfm + 12*(m44 - 2*m11)*m11*L3r*nfm + 12*(
         m44 - 2*m11)*m11*L0r*nfm - 32*(m44 - m11)*m11*L7r );
 
  return PQM12p6x11VL/pow(mass.getf0(),4);
}

double mnfPQSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQM12p6x11VR =
       + bbvm11*cbvm11 * ( 1./2.*(m44 - m11)*(m44 - m11)*pow(m11,2)*
         pow(nfm,2) );

      PQM12p6x11VR +=  + pow(bbvm11,2) * (  - 1./4.*(m44 - m11)*pow(
         m11,2)*nfm + 1./8.*(m44 - m11)*(m44 - 5*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm44*bbvm11 * ( 1./2.*m11*m44 - 1./4.*m11*
         m44*nfm - 1./4.*m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm14*bbvm11 * (  - m11*m44 );

      PQM12p6x11VR +=  + pow(abvm14,2) * (  - 1./4.*m11*pow(nf,2) );

      PQM12p6x11VR +=  + abvm11*cbvm11 * (  - 1./2.*(m44 - m11)*pow(
         m11,2)*nfm - 1./2.*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm11*bbvm11 * ( 3./4.*pow(m11,2) - 1./4.*(
         m44 - 4*m11)*m11*nfm - 1./4.*(m44 - 3*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + abvm11*abvm14 * ( 1./2.*m11*nf );

      PQM12p6x11VR +=  + pow(abvm11,2) * ( 3./8.*m11 + 1./4.*m11*nfm + 
         1./8.*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*cbvm11 * (  - 1./2.*(m44 - m11)*(m44 - 
         m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*bbvm11 * ( 1./2.*(m44 - m11)*pow(m11,2)*
         nfm - 1./2.*(m44 - m11)*(m44 - 3*m11)*m11*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm44 * (  - 1./2.*m11*m44 + 1./4.*m11*
         m44*nfm + 1./4.*m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + pi16*abvm14 * ( 1./24.*(24*m44 - 5*m44*nf + 8*
         m44*pow(nf,2) - 3*m11*nf)*m11 );

      PQM12p6x11VR +=  + pi16*abvm11 * (  - 3./4.*pow(m11,2) + 1./2.*(
         m44 - 2*m11)*m11*pow(nfm,2) + 1./4.*(2*m44 - 5*m11)*m11*nfm );

      PQM12p6x11VR +=  + Ab(m11,mu2)*cbvm11 * (  - 1./2.*(m44 - m11)*
         pow(m11,2)*nfm + 1./2.*(m44 - m11)*(m44 - 2*m11)*m11*pow(
         nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 3./4.*pow(m11,2) - 3./4.
         *(m44 - 2*m11)*m11*nfm + 1./4.*(pow(m44,2) - 7*m11*m44 + 8*
         pow(m11,2))*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( 1./2.*m44 - 1./4.*m44*
         nfm - 1./4.*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm14 * ( 1./2.*( - 2*m44 + m11*
         nf) );

      PQM12p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 3./2.*m11 - 1./4.*(m44
          - 6*m11)*nfm - 1./4.*(m44 - 4*m11)*pow(nfm,2) );

      PQM12p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - m11*m44 );

      PQM12p6x11VR +=  + Ab(m14,mu2)*abvm14 * (  - 1./2.*m11*pow(nf,2)
          );

      PQM12p6x11VR +=  + Ab(m14,mu2)*abvm11 * ( 1./2.*m11*nf );

      PQM12p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( 1./2.*m11*m44 - 1./4.*
         m11*m44*nfm - 1./4.*m11*m44*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m11,m11,m11,L,mu2) * ( 1./3.*pow(
         m11,2) + 1./2.*pow(m11,2)*nfm + 1./2.*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * (  - 1./8.*( - 
         m44 + 4*m11 + 5*m11*nf)*m11 );

      PQM12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * ( 1./8.*( - 1
          + nf)*(1 + 2*nf)*m11*m44 );

      PQM12p6x11VR +=  + hhV(2,m11,m11,m11,m11,L,mu2) * (  - 1./2.*(m44
          - m11)*pow(m11,2)*nfm - (m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * ( 5./8.*(m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hhV(5,m11,m11,m11,m11,L,mu2) * ( 1./2.*(m44 - 
         m11)*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11VR +=  + hh1V(1,m11,m14,m14,m11,L,mu2) * ( pow(m11,2)*
         nf );

      PQM12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * (  - (m44 - 
         m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m11,m14,m14,m11,L,mu2) * (  - 3./8.*(
          - 1 + nf)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(1,m44,m14,m14,m11,L,mu2) * ( 3./8.*( - 1
          + nf)*(1 + 2*nf)*pow(m11,2) );

      PQM12p6x11VR +=  + hh21V(2,m11,m14,m14,m11,L,mu2) * ( 3./8.*(m44
          - m11)*pow(m11,2) );

      PQM12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * ( 3./8.*( - 1
          + nf)*m11 );

      PQM12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * (  - 3./8.*(
          - 1 + nf)*(1 + 2*nf)*m11 );

      PQM12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * (  - 3./8.*(
         m44 - m11)*m11 );

  return  PQM12p6x11VR/pow(mass.getf0(),4);
}
////////////////////////////////////////////////////////////////////////
// decay SPN

double fnfPQSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double abvm14 = AbV(m14,L);
  //double bbvm11 = BbV(m11,L);
  double   PQF12p4x11VR =
       + abvm14 * ( 1./2.*nf );

  return  PQF12p4x11VR/pow(mass.getf0(),2);
}

double fnfPQSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  fnfPQSPNp6LV(nf,mass,Liin,L)+fnfPQSPNp6RV(nf,mass,L);
}

double fnfPQSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double a23bvm11 = A23bV(m11,L);
  double a23bvm14 = A23bV(m14,L);
  double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double b23bvm11 = B23bV(m11,L);
  //double cbvm11 = CbV(m11,L);
  double    PQF12p6x11VL =
       + b23bvm11 * ( 6*(m44 - m11)*L3r*nfm + 6*(m44 - m11)*L0r*nfm );

      PQF12p6x11VL +=  + bbvm14 * ( 8*(m44 + m11)*m44*L6r*pow(nf,2) - 4
         *(m44 + m11)*m44*L4r*pow(nf,2) + 2*(m44 + m11)*(m44 + m11)*L8r
         *nf - (m44 + m11)*(m44 + m11)*L5r*nf );

      PQF12p6x11VL +=  + bbvm11 * ( 2*(m44 - m11)*m11*L5r*nfm - 6*(m44
          - m11)*m11*L3r*nfm - 6*(m44 - m11)*m11*L0r*nfm );

      PQF12p6x11VL +=  + a23bvm44 * ( 6*( - 1 + nf)*(1 + 2*nf)*L2r );

      PQF12p6x11VL +=  + a23bvm14 * ( 6*L3r*nf + 12*L0r*nf );

      PQF12p6x11VL +=  + a23bvm11 * (  - 6*L3r - 6*L3r*nfm + 6*L2r + 12
         *L1r - 6*L0r - 6*L0r*nfm );

      PQF12p6x11VL +=  + abvm44 * ( 4*( - 1 + nf)*(1 + 2*nf)*m44*L4r - 
         2*( - 1 + nf)*(1 + 2*nf)*m44*L2r - 8*( - 1 + nf)*(1 + 2*nf)*
         m44*L1r );

      PQF12p6x11VL +=  + abvm14 * (  - 4*m44*L4r*pow(nf,2) + 2*m11*L5r*
         nf - 5*(m44 + m11)*L3r*nf - 2*(m44 + m11)*L0r*nf );

      PQF12p6x11VL +=  + abvm11 * (  - 2*m11*L5r - 2*m11*L5r*nfm + 6*
         m11*L3r - 10*m11*L2r - 4*m11*L1r + 6*m11*L0r - 6*(m44 - 2*m11)
         *L3r*nfm - 6*(m44 - 2*m11)*L0r*nfm );

  return PQF12p6x11VL/pow(mass.getf0(),4);
}

double fnfPQSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  //double abvm44 = AbV(m44,L);
  //double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double cbvm11 = CbV(m11,L);
  double   PQF12p6x11VR =
       + abvm14*bbvm14 * (  - 1./8.*(m44 + m11)*nf );

      PQF12p6x11VR +=  + pow(abvm14,2) * (  - 1./8.*nf );

      PQF12p6x11VR +=  + abvm11*bbvm14 * (  - 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + pi16*abvm14 * (  - 1./48.*( - 11*m44 + 8*m44*
         nf - 11*m11 + 2*m11*nf)*nf );

      PQF12p6x11VR +=  + pi16*abvm11 * ( 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./8.*(m44 + m11) );

      PQF12p6x11VR +=  + Ab(m14,mu2)*bbvm14 * (  - 1./8.*(m44 + m11)*nf
          );

      PQF12p6x11VR +=  + Ab(m14,mu2)*abvm14 * (  - 1./4.*nf - 1./8.*(
         m44 + m11)*pow(m14,-1)*nf );

      PQF12p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./8.*(m44 + m11)*
         pow(m14,-1) );

      PQF12p6x11VR +=  + hhV(1,m11,m14,m14,m11,L,mu2) * ( 1./16.*( - 
         m44 + m11*nf) );

      PQF12p6x11VR +=  + hhV(1,m44,m14,m14,m11,L,mu2) * (  - 1./16.*(
          - 1 + nf)*(1 + 2*nf)*m44 );

      PQF12p6x11VR +=  + hhV(2,m11,m14,m14,m11,L,mu2) * (  - 5./72.*(
         m44 - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(2,m11,m14,m14,m11,L,mu2) * ( 1./144.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh1V(3,m14,m11,m14,m11,L,mu2) * ( 1./72.*(m44
          - m11)*m11 );

      PQF12p6x11VR +=  + hh27V(1,m11,m14,m14,m11,L,mu2) * (  - 3./16.*(
          - 1 + nf) );

      PQF12p6x11VR +=  + hh27V(1,m44,m14,m14,m11,L,mu2) * ( 3./16.*( - 
         1 + nf)*(1 + 2*nf) );

      PQF12p6x11VR +=  + hh27V(2,m11,m14,m14,m11,L,mu2) * ( 3./16.*(m44
          - m11) );

      PQF12p6x11VR +=  + hhdV(1,m11,m11,m11,m11,L,mu2) * ( 1./6.*pow(
         m11,2) + 1./4.*pow(m11,2)*nfm + 1./4.*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hhdV(1,m11,m14,m14,m11,L,mu2) * (  - 1./16.*(
          - m44 + 4*m11 + 5*m11*nf)*m11 );

      PQF12p6x11VR +=  + hhdV(1,m44,m14,m14,m11,L,mu2) * ( 1./16.*( - 1
          + nf)*(1 + 2*nf)*m11*m44 );

      PQF12p6x11VR +=  + hhdV(2,m11,m11,m11,m11,L,mu2) * (  - 1./4.*(
         m44 - m11)*pow(m11,2)*nfm - 1./2.*(m44 - m11)*pow(m11,2)*pow(
         nfm,2) );

      PQF12p6x11VR +=  + hhdV(2,m11,m14,m14,m11,L,mu2) * ( 5./16.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hhdV(5,m11,m11,m11,m11,L,mu2) * ( 1./4.*(m44
          - m11)*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11VR +=  + hh1dV(1,m11,m14,m14,m11,L,mu2) * ( 1./2.*pow(
         m11,2)*nf );

      PQF12p6x11VR +=  + hh1dV(2,m11,m14,m14,m11,L,mu2) * (  - 1./2.*(
         m44 - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m11,m14,m14,m11,L,mu2) * (  - 3./16.*
         ( - 1 + nf)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(1,m44,m14,m14,m11,L,mu2) * ( 3./16.*(
          - 1 + nf)*(1 + 2*nf)*pow(m11,2) );

      PQF12p6x11VR +=  + hh21dV(2,m11,m14,m14,m11,L,mu2) * ( 3./16.*(
         m44 - m11)*pow(m11,2) );

      PQF12p6x11VR +=  + hh27dV(1,m11,m14,m14,m11,L,mu2) * ( 3./16.*(
          - 1 + nf)*m11 );

      PQF12p6x11VR +=  + hh27dV(1,m44,m14,m14,m11,L,mu2) * (  - 3./16.*
         ( - 1 + nf)*(1 + 2*nf)*m11 );

      PQF12p6x11VR +=  + hh27dV(2,m11,m14,m14,m11,L,mu2) * (  - 3./16.*
         (m44 - m11)*m11 );

  return  PQF12p6x11VR/pow(mass.getf0(),4);
}

////////////////////////////////////////////////////////////////////////
// vev SPN

double qnfPQSPNp4V(const int nf, const quarkmassnf mass, const double L){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp4RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double bbvm11 = BbV(m11,L);
  double   PQQ11p4x11VR =
       + bbvm11 * ( 1./2.*(m44 - m11)*nfm );

      PQQ11p4x11VR +=  + abvm14 * ( nf );

      PQQ11p4x11VR +=  + abvm11 * (  - 1./2. - 1./2.*nfm );

  return  PQQ11p4x11VR/pow(mass.getf0(),2);
}

double qnfPQSPNp6V(const int nf, const quarkmassnf mass, const Linf Liin,
		const double L){
  return  qnfPQSPNp6LV(nf,mass,Liin,L)+qnfPQSPNp6RV(nf,mass,L);
}

double qnfPQSPNp6LV(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp6LV probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  //double a23bvm11 = A23bV(m11,L);
  //double a23bvm14 = A23bV(m14,L);
  //double a23bvm44 = A23bV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  //double b23bvm11 = B23bV(m11,L);
  double cbvm11 = CbV(m11,L);
  double   PQQ11p6x11VL =
       + cbvm11 * ( 32*(m44 - m11)*m11*m44*L6r - 16*(m44 - m11)*m11*m44
         *L4r + 16*(m44 - m11)*pow(m11,2)*L8r*nfm - 8*(m44 - m11)*pow(
         m11,2)*L5r*nfm );

      PQQ11p6x11VL +=  + bbvm14 * ( 16*(m44 + m11)*m44*L6r*pow(nf,2) - 
         8*(m44 + m11)*m44*L4r*pow(nf,2) + 4*(m44 + m11)*(m44 + m11)*
         L8r*nf - 2*(m44 + m11)*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + bbvm11 * (  - 8*pow(m11,2)*L8r + 4*pow(m11,2)*
         L5r - 16*( - 2*m44 + 3*m11 + m11*nf)*m44*L6r + 8*( - 2*m44 + 3
         *m11 + m11*nf)*m44*L4r + 16*(m44 - m11)*(m44 - m11)*L7r - 4*(4
         *m44 - 5*m11)*m11*L5r*nfm + 8*(pow(m44,2) + 2*m11*m44 - 4*pow(
         m11,2))*L8r*nfm );

      PQQ11p6x11VL +=  + abvm44 * ( 16*( - 1 + nf)*(1 + 2*nf)*m44*L6r
          - 8*( - 1 + nf)*(1 + 2*nf)*m44*L4r );

      PQQ11p6x11VL +=  + abvm14 * ( 32*m44*L6r*pow(nf,2) - 16*m44*L4r*
         pow(nf,2) + 16*(m44 + m11)*L8r*nf - 8*(m44 + m11)*L5r*nf );

      PQQ11p6x11VL +=  + abvm11 * (  - 16*m11*L8r - 16*m11*L8r*nfm + 8*
         m11*L5r - 16*(1 + nf)*m44*L6r + 8*(1 + nf)*m44*L4r - 8*(m44 - 
         2*m11)*L5r*nfm - 32*(m44 - m11)*L7r );
 
  return PQQ11p6x11VL/pow(mass.getf0(),4);
}

double qnfPQSPNp6RV(const int nf,const quarkmassnf mass, const double L){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp6RV probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double abvm11 = AbV(m11,L);
  double abvm14 = AbV(m14,L);
  double abvm44 = AbV(m44,L);
  double bbvm11 = BbV(m11,L);
  double bbvm14 = BbV(m14,L);
  double cbvm11 = CbV(m11,L);
  double  PQQ11p6x11VR =
       + bbvm11*cbvm11 * ( 1./2.*(m44 - m11)*(m44 - m11)*m11*pow(nfm,2)
          );

      PQQ11p6x11VR +=  + pow(bbvm11,2) * (  - 1./4.*(m44 - m11)*m11*nfm
          + 1./8.*(m44 - m11)*(m44 - 5*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm44*bbvm11 * ( 1./2.*m44 - 1./4.*m44*nfm - 
         1./4.*m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm14*bbvm14 * (  - 1./4.*(m44 + m11)*nf );

      PQQ11p6x11VR +=  + abvm14*bbvm11 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + pow(abvm14,2) * (  - 1./4.*nf );

      PQQ11p6x11VR +=  + abvm11*cbvm11 * (  - 1./2.*(m44 - m11)*m11*nfm
          - 1./2.*(m44 - m11)*m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm11*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQQ11p6x11VR +=  + abvm11*bbvm11 * ( 3./4.*m11 - 1./4.*(m44 - 4*
         m11)*nfm - 1./4.*(m44 - 3*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + abvm11*abvm14 * (  - 1./2. );

      PQQ11p6x11VR +=  + pow(abvm11,2) * ( 3./8. + 1./4.*nfm + 1./8.*
         pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*cbvm11 * (  - 1./2.*(m44 - m11)*(m44 - 
         m11)*m11*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*bbvm11 * ( 1./2.*(m44 - m11)*m11*nfm - 1./
         2.*(m44 - m11)*(m44 - 3*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm44 * (  - 1./2.*m44 + 1./4.*m44*nfm
          + 1./4.*m44*pow(nfm,2) );

      PQQ11p6x11VR +=  + pi16*abvm14 * ( 1./4.*(2 + nf)*(m44 + m11) );

      PQQ11p6x11VR +=  + pi16*abvm11 * ( 1./4.*(m44 - 2*m11) + 1./2.*(
         m44 - 2*m11)*pow(nfm,2) + 1./4.*(2*m44 - 5*m11)*nfm );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*cbvm11 * (  - 1./2.*(m44 - m11)*
         m11*nfm + 1./2.*(m44 - m11)*(m44 - 2*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm14 * (  - 1./4.*(m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*bbvm11 * ( 1./4.*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 3./4.*m11 - 3./4.*(m44 - 2*m11)*nfm - 1./4.
         *(7*m44 - 8*m11)*pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm44 * ( 1./2.*pow(m11,-1)*m44
          - 1./4.*pow(m11,-1)*m44*nfm - 1./4.*pow(m11,-1)*m44*pow(
         nfm,2) );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm14 * (  - 1 - 1./2.*pow(
         m11,-1)*m44 );

      PQQ11p6x11VR +=  + Ab(m11,mu2)*abvm11 * ( 3./2. - 1./4.*pow(
         m11,-1)*m44*nfm - 1./4.*pow(m11,-1)*m44*pow(nfm,2) + 3./2.*nfm
          + pow(nfm,2) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*bbvm14 * (  - 1./4.*(m44 + m11)*nf
          );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*bbvm11 * (  - 1./2.*(m44 + m11) );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*abvm14 * (  - 1./2.*nf - 1./4.*(
         m44 + m11)*pow(m14,-1)*nf );

      PQQ11p6x11VR +=  + Ab(m14,mu2)*abvm11 * (  - 1./2. - 1./4.*(m44
          + m11)*pow(m14,-1) );

      PQQ11p6x11VR +=  + Ab(m44,mu2)*bbvm11 * ( 1./2.*m44 - 1./4.*m44*
         nfm - 1./4.*m44*pow(nfm,2) );

  return  PQQ11p6x11VR/pow(mass.getf0(),4);
}

