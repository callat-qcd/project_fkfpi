// massdecayvevloV.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.1
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains the isospin limit expressions for the finite volume
// corrections for masses and decay constants
// in terms of the infinite volume lowest order masses and F0 

// by setting the compile flag, it either produces the
// Besselfunction (CHIRONBESSEL) or the Theta function (CHIRONTHETA) output
// typically: large mL Bessel better
//            intermediate or small mL theta function version better
//   however jbdadmul as used in theta version has quite often problems
#include <iostream>
#include <cmath>

#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "oneloopintegrals.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "massdecayvevloV.h"

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
#define mpi4loV mpi4loVt
#define mpi6LloV mpi6LloVt
#define mpi6RloV mpi6RloVt
#define mpi6loV mpi6loVt
#define mk4loV mk4loVt
#define mk6LloV mk6LloVt
#define mk6RloV mk6RloVt
#define mk6loV mk6loVt
#define meta4loV meta4loVt
#define meta6LloV meta6LloVt
#define meta6RloV meta6RloVt
#define meta6loV meta6loVt
#define fpi4loV fpi4loVt
#define fpi6LloV fpi6LloVt
#define fpi6RloV fpi6RloVt
#define fpi6loV fpi6loVt
#define fk4loV fk4loVt
#define fk6LloV fk6LloVt
#define fk6RloV fk6RloVt
#define fk6loV fk6loVt
#define feta4loV feta4loVt
#define feta6LloV feta6LloVt
#define feta6RloV feta6RloVt
#define feta6loV feta6loVt
#define qqup4loV   qqup4loVt 
#define qqup6loV   qqup6loVt 
#define qqup6LloV  qqup6LloVt
#define qqup6RloV  qqup6RloVt
#define qqstrange4loV   qqstrange4loVt 
#define qqstrange6loV   qqstrange6loVt 
#define qqstrange6LloV  qqstrange6LloVt
#define qqstrange6RloV  qqstrange6RloVt
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
#define mpi4loV mpi4loVb
#define mpi6LloV mpi6LloVb
#define mpi6RloV mpi6RloVb
#define mpi6loV mpi6loVb
#define mk4loV mk4loVb
#define mk6LloV mk6LloVb
#define mk6RloV mk6RloVb
#define mk6loV mk6loVb
#define meta4loV meta4loVb
#define meta6LloV meta6LloVb
#define meta6RloV meta6RloVb
#define meta6loV meta6loVb
#define fpi4loV fpi4loVb
#define fpi6LloV fpi6LloVb
#define fpi6RloV fpi6RloVb
#define fpi6loV fpi6loVb
#define fk4loV fk4loVb
#define fk6LloV fk6LloVb
#define fk6RloV fk6RloVb
#define fk6loV fk6loVb
#define feta4loV feta4loVb
#define feta6LloV feta6LloVb
#define feta6RloV feta6RloVb
#define feta6loV feta6loVb
#define qqup4loV   qqup4loVb 
#define qqup6loV   qqup6loVb 
#define qqup6LloV  qqup6LloVb
#define qqup6RloV  qqup6RloVb
#define qqstrange4loV   qqstrange4loVb 
#define qqstrange6loV   qqstrange6loVb 
#define qqstrange6LloV  qqstrange6LloVb
#define qqstrange6RloV  qqstrange6RloVb
#else
// just some garbage to produce an error
x = 1./0.0.;
#endif
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double mpi4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double M4V = (  - 1./2.*AbV(mp2,xl)*mp2 + 1./6.*AbV(me2,xl)*mp2 ); 
  return M4V/pow(mass.getf0(),2);
}


double mpi6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LloV\n";}

  double massp6LV =
       + AbV(mp2,xl) * ( 32*L8r*pow(mp2,2) - 16*L6r*mp2*mk2 + 72*L6r*
         pow(mp2,2) - 16*L5r*pow(mp2,2) + 16*L4r*mp2*mk2 - 40*L4r*pow(
         mp2,2) + 28*L3r*pow(mp2,2) + 32*L2r*pow(mp2,2) + 56*L1r*pow(
         mp2,2) );

      massp6LV +=  + AbV(mk2,xl) * ( 32*L8r*mp2*mk2 + 64*L6r*mp2*mk2 - 
         16*L5r*mp2*mk2 - 64*L4r*mp2*mk2 + 20*L3r*mp2*mk2 + 16*L2r*mp2*
         mk2 + 64*L1r*mp2*mk2 );

      massp6LV +=  + AbV(me2,xl) * ( 32./3.*L8r*pow(mp2,2) - 64./3.*L7r
         *mp2*mk2 + 64./3.*L7r*pow(mp2,2) + 80./3.*L6r*mp2*mk2 - 8./3.*
         L6r*pow(mp2,2) - 32./9.*L5r*mp2*mk2 - 16./9.*L5r*pow(mp2,2) - 
         80./3.*L4r*mp2*mk2 + 8./3.*L4r*pow(mp2,2) + 16./3.*L3r*mp2*mk2
          - 4./3.*L3r*pow(mp2,2) + 16./3.*L2r*mp2*mk2 - 4./3.*L2r*pow(
         mp2,2) + 64./3.*L1r*mp2*mk2 - 16./3.*L1r*pow(mp2,2) );

      massp6LV +=  + A23bV(mp2,xl) * (  - 12*L3r*mp2 - 48*L2r*mp2 - 24*
         L1r*mp2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 12*L3r*mp2 - 48*L2r*mp2 );

      massp6LV +=  + A23bV(me2,xl) * (  - 4*L3r*mp2 - 12*L2r*mp2 );

      massp6LV +=  + BbV(mp2,xl) * (  - 8*L8r*pow(mp2,3) - 16*L6r*pow(
         mp2,2)*mk2 - 8*L6r*pow(mp2,3) + 4*L5r*pow(mp2,3) + 8*L4r*pow(
         mp2,2)*mk2 + 4*L4r*pow(mp2,3) );

      massp6LV +=  + BbV(me2,xl) * ( 64./9.*L8r*mp2*pow(mk2,2) - 64./9.
         *L8r*pow(mp2,2)*mk2 + 8./3.*L8r*pow(mp2,3) + 64./9.*L7r*mp2*
         pow(mk2,2) - 128./9.*L7r*pow(mp2,2)*mk2 + 64./9.*L7r*pow(
         mp2,3) + 64./9.*L6r*mp2*pow(mk2,2) + 16./9.*L6r*pow(mp2,2)*mk2
          - 8./9.*L6r*pow(mp2,3) - 64./27.*L5r*mp2*pow(mk2,2) + 32./27.
         *L5r*pow(mp2,2)*mk2 - 4./27.*L5r*pow(mp2,3) - 32./9.*L4r*mp2*
         pow(mk2,2) - 8./9.*L4r*pow(mp2,2)*mk2 + 4./9.*L4r*pow(mp2,3) )
         ;

      return massp6LV/pow(mass.getf0(),4);
}

double mpi6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 3./8.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * (  - 1./2.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * (  - 1./12.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * ( 1./4.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 1./12.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl) * ( pi16*mp2*mk2 + pi16*pow(mp2,2) - 3.
         /4.*Ab(mp2,mu2)*mp2 - 1./2.*Ab(mk2,mu2)*mp2 - 1./12.*Ab(me2,
         mu2)*mp2 + 1./4.*Bb(mp2,mu2)*pow(mp2,2) + 1./12.*Bb(me2,mu2)*
         pow(mp2,2) );

      massp6RV +=  + pow(AbV(mk2,xl),2) * (  - 1./4.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 1./2.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 2./9.*mp2*mk2 );

      massp6RV +=  + AbV(mk2,xl) * ( pi16*mp2*mk2 - 1./2.*Ab(mp2,mu2)*
         mp2 - 1./2.*Ab(mk2,mu2)*mp2 - 1./2.*Ab(me2,mu2)*mp2 - 2./9.*
         Bb(me2,mu2)*mp2*mk2 );

      massp6RV +=  + pow(AbV(me2,xl),2) * ( 1./72.*mp2 );

      massp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * (  - 1./12.*pow(mp2,2) )
         ;

      massp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 4./27.*mp2*mk2 - 7./
         108.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl) * (  - 1./12.*Ab(mp2,mu2)*mp2 - 1./2.*
         Ab(mk2,mu2)*mp2 + 1./36.*Ab(me2,mu2)*mp2 - 1./12.*Bb(mp2,mu2)*
         pow(mp2,2) + 4./27.*Bb(me2,mu2)*mp2*mk2 - 7./108.*Bb(me2,mu2)*
         pow(mp2,2) );

      massp6RV +=  + BbV(mp2,xl) * ( 1./4.*Ab(mp2,mu2)*pow(mp2,2) - 1./
         12.*Ab(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + BbV(me2,xl) * ( 1./12.*Ab(mp2,mu2)*pow(mp2,2) - 2./
         9.*Ab(mk2,mu2)*mp2*mk2 + 4./27.*Ab(me2,mu2)*mp2*mk2 - 7./108.*
         Ab(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * ( 5./6.*pow(mp2,2) )
         ;

      massp6RV +=  + hhV(mp2,mk2,mk2,mp2,xl,mu2) * ( mp2*mk2 - 5./8.*
         pow(mp2,2) );

      massp6RV +=  + hhV(mp2,me2,me2,mp2,xl,mu2) * ( 1./18.*pow(mp2,2)
          );

      massp6RV +=  + hhV(mk2,mk2,me2,mp2,xl,mu2) * ( 1./2.*mp2*mk2 + 1./
         24.*pow(mp2,2) );

      massp6RV +=  + hh1V(mp2,mk2,mk2,mp2,xl,mu2) * ( pow(mp2,2) );

      massp6RV +=  + hh1V(me2,mk2,mk2,mp2,xl,mu2) * (  - pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3*pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mk2,mk2,mp2,xl,mu2) * (  - 3./8.*pow(
         mp2,2) );

      massp6RV +=  + hh21V(mk2,mp2,mk2,mp2,xl,mu2) * ( 3*pow(mp2,2) );

      massp6RV +=  + hh21V(me2,mk2,mk2,mp2,xl,mu2) * ( 9./8.*pow(mp2,2)
          );

      massp6RV +=  + hh27V(mp2,mp2,mp2,mp2,xl,mu2) * (  - 3*mp2 );

      massp6RV +=  + hh27V(mp2,mk2,mk2,mp2,xl,mu2) * ( 3./8.*mp2 );

      massp6RV +=  + hh27V(mk2,mp2,mk2,mp2,xl,mu2) * (  - 3*mp2 );

      massp6RV +=  + hh27V(me2,mk2,mk2,mp2,xl,mu2) * (  - 9./8.*mp2 );

      return massp6RV/pow(mass.getf0(),4);
}


//***********************************************************************

double mk4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double M4V = (  - 1./4.*me2 - 1./12.*mp2 )*AbV(me2,xl);
  return M4V/pow(mass.getf0(),2);
}

double mk6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LloV\n";}

  double massp6LV =
       + AbV(mp2,xl) * ( 24*L8r*mp2*mk2 + 48*L6r*mp2*mk2 - 12*L5r*mp2*
         mk2 - 48*L4r*mp2*mk2 + 15*L3r*mp2*mk2 + 12*L2r*mp2*mk2 + 48*
         L1r*mp2*mk2 );

      massp6LV +=  + AbV(mk2,xl) * ( 48*L8r*pow(mk2,2) + 96*L6r*pow(
         mk2,2) - 24*L5r*pow(mk2,2) - 64*L4r*pow(mk2,2) + 30*L3r*pow(
         mk2,2) + 36*L2r*pow(mk2,2) + 72*L1r*pow(mk2,2) );

      massp6LV +=  + AbV(me2,xl) * ( 32./3.*L8r*pow(mk2,2) - 8*L8r*mp2*
         mk2 + 32./3.*L7r*pow(mk2,2) - 32./3.*L7r*mp2*mk2 + 32./3.*L6r*
         pow(mk2,2) - 32./3.*L6r*mp2*mk2 - 32./9.*L5r*pow(mk2,2) + 20./
         9.*L5r*mp2*mk2 - 32./3.*L4r*pow(mk2,2) + 32./3.*L4r*mp2*mk2 + 
         28./3.*L3r*pow(mk2,2) - 7./3.*L3r*mp2*mk2 + 16./3.*L2r*pow(
         mk2,2) - 4./3.*L2r*mp2*mk2 + 64./3.*L1r*pow(mk2,2) - 16./3.*
         L1r*mp2*mk2 );

      massp6LV +=  + A23bV(mp2,xl) * (  - 9*L3r*mk2 - 36*L2r*mk2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 18*L3r*mk2 - 60*L2r*mk2 - 24*
         L1r*mk2 );

      massp6LV +=  + A23bV(me2,xl) * (  - L3r*mk2 - 12*L2r*mk2 );

      massp6LV +=  + BbV(me2,xl) * (  - 128./9.*L8r*pow(mk2,3) + 128./9.
         *L8r*mp2*pow(mk2,2) - 16./3.*L8r*pow(mp2,2)*mk2 - 128./9.*L7r*
         pow(mk2,3) + 256./9.*L7r*mp2*pow(mk2,2) - 128./9.*L7r*pow(
         mp2,2)*mk2 - 128./9.*L6r*pow(mk2,3) - 32./9.*L6r*mp2*pow(
         mk2,2) + 16./9.*L6r*pow(mp2,2)*mk2 + 128./27.*L5r*pow(mk2,3)
          - 64./27.*L5r*mp2*pow(mk2,2) + 8./27.*L5r*pow(mp2,2)*mk2 + 64.
         /9.*L4r*pow(mk2,3) + 16./9.*L4r*mp2*pow(mk2,2) - 8./9.*L4r*
         pow(mp2,2)*mk2 );
  
      return massp6LV/pow(mass.getf0(),4);
}

double mk6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 3./32.*mk2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * (  - 3./4.*mk2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * (  - 3./16.*mk2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * (  - 1./6.*mp2*mk2 );

      massp6RV +=  + AbV(mp2,xl) * ( 3./4.*pi16*pow(mk2,2) - 3./16.*Ab(
         mp2,mu2)*mk2 - 3./4.*Ab(mk2,mu2)*mk2 - 3./16.*Ab(me2,mu2)*mk2
          - 1./6.*Bb(me2,mu2)*mp2*mk2 );

      massp6RV +=  + pow(AbV(mk2,xl),2) * (  - 3./4.*mk2 );

      massp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * ( 4./9.*pow(mk2,2) );

      massp6RV +=  + AbV(mk2,xl) * ( 3./4.*pi16*pow(mk2,2) + 3./4.*pi16
         *mp2*mk2 - 3./4.*Ab(mp2,mu2)*mk2 - 3./2.*Ab(mk2,mu2)*mk2 + 4./
         9.*Bb(me2,mu2)*pow(mk2,2) );

      massp6RV +=  + pow(AbV(me2,xl),2) * ( 25./288.*mk2 );

      massp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * (  - 8./27.*pow(mk2,2)
          + 7./54.*mp2*mk2 );

      massp6RV +=  + AbV(me2,xl) * ( 1./2.*pi16*pow(mk2,2) + 1./4.*pi16
         *mp2*mk2 - 3./16.*Ab(mp2,mu2)*mk2 + 25./144.*Ab(me2,mu2)*mk2
          - 8./27.*Bb(me2,mu2)*pow(mk2,2) + 7./54.*Bb(me2,mu2)*mp2*mk2
          );

      massp6RV +=  + BbV(me2,xl) * (  - 1./6.*Ab(mp2,mu2)*mp2*mk2 + 4./
         9.*Ab(mk2,mu2)*pow(mk2,2) - 8./27.*Ab(me2,mu2)*pow(mk2,2) + 7./
         54.*Ab(me2,mu2)*mp2*mk2 );

      massp6RV +=  + hhV(mp2,mp2,mk2,mk2,xl,mu2) * (  - 15./32.*pow(
         mk2,2) + 3./4.*mp2*mk2 );

      massp6RV +=  + hhV(mp2,mk2,me2,mk2,xl,mu2) * ( 13./16.*pow(mk2,2)
          );

      massp6RV +=  + hhV(mk2,mk2,mk2,mk2,xl,mu2) * ( 3./4.*pow(mk2,2) )
         ;

      massp6RV +=  + hhV(mk2,me2,me2,mk2,xl,mu2) * ( 181./288.*pow(
         mk2,2) );

      massp6RV +=  + hh1V(mk2,mp2,mp2,mk2,xl,mu2) * ( 3./4.*pow(mk2,2)
          );

      massp6RV +=  + hh1V(mk2,mp2,me2,mk2,xl,mu2) * (  - 3./2.*pow(
         mk2,2) );

      massp6RV +=  + hh1V(mk2,me2,me2,mk2,xl,mu2) * (  - 5./4.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mp2,mp2,mk2,mk2,xl,mu2) * ( 9./4.*pow(mk2,2)
          );

      massp6RV +=  + hh21V(mk2,mp2,mp2,mk2,xl,mu2) * (  - 9./32.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,mp2,me2,mk2,xl,mu2) * ( 27./16.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,mk2,mk2,mk2,xl,mu2) * ( 9./4.*pow(mk2,2)
          );

      massp6RV +=  + hh21V(mk2,me2,me2,mk2,xl,mu2) * ( 27./32.*pow(
         mk2,2) );

      massp6RV +=  + hh27V(mp2,mp2,mk2,mk2,xl,mu2) * (  - 9./4.*mk2 );

      massp6RV +=  + hh27V(mk2,mp2,mp2,mk2,xl,mu2) * ( 9./32.*mk2 );

      massp6RV +=  + hh27V(mk2,mp2,me2,mk2,xl,mu2) * (  - 27./16.*mk2 )
         ;

      massp6RV +=  + hh27V(mk2,mk2,mk2,mk2,xl,mu2) * (  - 9./4.*mk2 );

      massp6RV +=  + hh27V(mk2,me2,me2,mk2,xl,mu2) * (  - 27./32.*mk2 )
         ;

      return massp6RV/pow(mass.getf0(),4);
}

//***********************************************************************

double meta4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double M4V =   1./2.*AbV(mp2,xl)*mp2 + AbV(mk2,xl)*( - me2 - 1./3.*mp2)
    + AbV(me2,xl)*(8./9.*mk2 - 7./18.*mp2 );
  return M4V/pow(mass.getf0(),2);
}

double meta6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;

  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LloV\n";}

  double massp6LV =
       + AbV(mp2,xl) * ( 32*L8r*pow(mp2,2) - 64*L7r*mp2*mk2 + 64*L7r*
         pow(mp2,2) + 80*L6r*mp2*mk2 - 8*L6r*pow(mp2,2) - 32./3.*L5r*
         mp2*mk2 - 16./3.*L5r*pow(mp2,2) - 80*L4r*mp2*mk2 + 8*L4r*pow(
         mp2,2) + 16*L3r*mp2*mk2 - 4*L3r*pow(mp2,2) + 16*L2r*mp2*mk2 - 
         4*L2r*pow(mp2,2) + 64*L1r*mp2*mk2 - 16*L1r*pow(mp2,2) );

      massp6LV +=  + AbV(mk2,xl) * ( 128./3.*L8r*pow(mk2,2) - 32*L8r*
         mp2*mk2 + 128./3.*L7r*pow(mk2,2) - 128./3.*L7r*mp2*mk2 + 128./
         3.*L6r*pow(mk2,2) - 128./3.*L6r*mp2*mk2 - 128./9.*L5r*pow(
         mk2,2) + 80./9.*L5r*mp2*mk2 - 128./3.*L4r*pow(mk2,2) + 128./3.
         *L4r*mp2*mk2 + 112./3.*L3r*pow(mk2,2) - 28./3.*L3r*mp2*mk2 + 
         64./3.*L2r*pow(mk2,2) - 16./3.*L2r*mp2*mk2 + 256./3.*L1r*pow(
         mk2,2) - 64./3.*L1r*mp2*mk2 );

      massp6LV +=  + AbV(me2,xl) * ( 1024./9.*L8r*pow(mk2,2) - 1024./9.
         *L8r*mp2*mk2 + 32*L8r*pow(mp2,2) + 1024./9.*L7r*pow(mk2,2) - 
         1664./9.*L7r*mp2*mk2 + 640./9.*L7r*pow(mp2,2) + 1024./9.*L6r*
         pow(mk2,2) - 368./9.*L6r*mp2*mk2 - 8./9.*L6r*pow(mp2,2) - 1024.
         /27.*L5r*pow(mk2,2) + 704./27.*L5r*mp2*mk2 - 112./27.*L5r*pow(
         mp2,2) - 512./9.*L4r*pow(mk2,2) + 112./9.*L4r*mp2*mk2 + 40./9.
         *L4r*pow(mp2,2) + 64./3.*L3r*pow(mk2,2) - 32./3.*L3r*mp2*mk2
          + 4./3.*L3r*pow(mp2,2) + 128./3.*L2r*pow(mk2,2) - 64./3.*L2r*
         mp2*mk2 + 8./3.*L2r*pow(mp2,2) + 128./3.*L1r*pow(mk2,2) - 64./
         3.*L1r*mp2*mk2 + 8./3.*L1r*pow(mp2,2) );

      massp6LV +=  + A23bV(mp2,xl) * (  - 16*L3r*mk2 + 4*L3r*mp2 - 48*
         L2r*mk2 + 12*L2r*mp2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 16./3.*L3r*mk2 + 4./3.*L3r*
         mp2 - 64*L2r*mk2 + 16*L2r*mp2 );

      massp6LV +=  + A23bV(me2,xl) * (  - 16*L3r*mk2 + 4*L3r*mp2 - 32*
         L2r*mk2 + 8*L2r*mp2 - 32*L1r*mk2 + 8*L1r*mp2 );

      massp6LV +=  + BbV(mp2,xl) * ( 8*L8r*pow(mp2,3) + 16*L6r*pow(
         mp2,2)*mk2 + 8*L6r*pow(mp2,3) - 4*L5r*pow(mp2,3) - 8*L4r*pow(
         mp2,2)*mk2 - 4*L4r*pow(mp2,3) );

      massp6LV +=  + BbV(mk2,xl) * (  - 64./3.*L8r*pow(mk2,3) - 128./3.
         *L6r*pow(mk2,3) - 64./3.*L6r*mp2*pow(mk2,2) + 32./3.*L5r*pow(
         mk2,3) + 64./3.*L4r*pow(mk2,3) + 32./3.*L4r*mp2*pow(mk2,2) );

      massp6LV +=  + BbV(me2,xl) * ( 1024./27.*L8r*pow(mk2,3) - 1472./
         27.*L8r*mp2*pow(mk2,2) + 832./27.*L8r*pow(mp2,2)*mk2 - 56./9.*
         L8r*pow(mp2,3) + 1024./27.*L7r*pow(mk2,3) - 832./9.*L7r*mp2*
         pow(mk2,2) + 640./9.*L7r*pow(mp2,2)*mk2 - 448./27.*L7r*pow(
         mp2,3) + 1024./27.*L6r*pow(mk2,3) - 64./9.*L6r*mp2*pow(mk2,2)
          - 80./9.*L6r*pow(mp2,2)*mk2 + 56./27.*L6r*pow(mp2,3) - 1024./
         81.*L5r*pow(mk2,3) + 320./27.*L5r*mp2*pow(mk2,2) - 32./9.*L5r*
         pow(mp2,2)*mk2 + 28./81.*L5r*pow(mp2,3) - 512./27.*L4r*pow(
         mk2,3) + 32./9.*L4r*mp2*pow(mk2,2) + 40./9.*L4r*pow(mp2,2)*mk2
          - 28./27.*L4r*pow(mp2,3) );

      return massp6LV/pow(mass.getf0(),4);
}

double meta6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 1./8.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * (  - 3./2.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * ( 1./12.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * (  - 1./4.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 4./9.*mp2*mk2 - 7./36.
         *pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl) * (  - 1./4.*Ab(mp2,mu2)*mp2 - 3./2.*
         Ab(mk2,mu2)*mp2 + 1./12.*Ab(me2,mu2)*mp2 - 1./4.*Bb(mp2,mu2)*
         pow(mp2,2) + 4./9.*Bb(me2,mu2)*mp2*mk2 - 7./36.*Bb(me2,mu2)*
         pow(mp2,2) );

      massp6RV +=  + pow(AbV(mk2,xl),2) * ( mk2 + 3./4.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 32./9.*mk2 + 3./2.*
         mp2 );

      massp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 32./27.*pow(mk2,2)
          + 14./27.*mp2*mk2 );

      massp6RV +=  + AbV(mk2,xl) * ( 8./3.*pi16*pow(mk2,2) + 2./3.*pi16
         *mp2*mk2 - 1./3.*pi16*pow(mp2,2) - 3./2.*Ab(mp2,mu2)*mp2 + 2*
         Ab(mk2,mu2)*mk2 + 3./2.*Ab(mk2,mu2)*mp2 - 32./9.*Ab(me2,mu2)*
         mk2 + 3./2.*Ab(me2,mu2)*mp2 - 32./27.*Bb(me2,mu2)*pow(mk2,2)
          + 14./27.*Bb(me2,mu2)*mp2*mk2 );

      massp6RV +=  + pow(AbV(me2,xl),2) * ( 8./27.*mk2 - 31./216.*mp2 )
         ;

      massp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./12.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * ( 4./9.*pow(mk2,2) );

      massp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 64./81.*pow(mk2,2) - 
         56./81.*mp2*mk2 + 49./324.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl) * ( 1./12.*Ab(mp2,mu2)*mp2 - 32./9.*
         Ab(mk2,mu2)*mk2 + 3./2.*Ab(mk2,mu2)*mp2 + 16./27.*Ab(me2,mu2)*
         mk2 - 31./108.*Ab(me2,mu2)*mp2 + 1./12.*Bb(mp2,mu2)*pow(mp2,2)
          + 4./9.*Bb(mk2,mu2)*pow(mk2,2) + 64./81.*Bb(me2,mu2)*pow(
         mk2,2) - 56./81.*Bb(me2,mu2)*mp2*mk2 + 49./324.*Bb(me2,mu2)*
         pow(mp2,2) );

      massp6RV +=  + BbV(mp2,xl) * (  - 1./4.*Ab(mp2,mu2)*pow(mp2,2) + 
         1./12.*Ab(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + BbV(mk2,xl) * ( 4./9.*Ab(me2,mu2)*pow(mk2,2) );

      massp6RV +=  + BbV(me2,xl) * ( 4./9.*Ab(mp2,mu2)*mp2*mk2 - 7./36.
         *Ab(mp2,mu2)*pow(mp2,2) - 32./27.*Ab(mk2,mu2)*pow(mk2,2) + 14./
         27.*Ab(mk2,mu2)*mp2*mk2 + 64./81.*Ab(me2,mu2)*pow(mk2,2) - 56./
         81.*Ab(me2,mu2)*mp2*mk2 + 49./324.*Ab(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + hhV(mp2,mp2,me2,me2,xl,mu2) * ( 1./6.*pow(mp2,2) )
         ;

      massp6RV +=  + hhV(mp2,mk2,mk2,me2,xl,mu2) * ( 3./2.*mp2*mk2 + 1./
         8.*pow(mp2,2) );

      massp6RV +=  + hhV(mk2,mk2,me2,me2,xl,mu2) * ( 50./9.*pow(mk2,2)
          - 11./3.*mp2*mk2 + 5./8.*pow(mp2,2) );

      massp6RV +=  + hhV(me2,me2,me2,me2,xl,mu2) * ( 128./243.*pow(
         mk2,2) - 112./243.*mp2*mk2 + 49./486.*pow(mp2,2) );

      massp6RV +=  + hh1V(mp2,mk2,mk2,me2,xl,mu2) * (  - 4*mp2*mk2 + 
         pow(mp2,2) );

      massp6RV +=  + hh1V(me2,mk2,mk2,me2,xl,mu2) * (  - 32./3.*pow(
         mk2,2) + 20./3.*mp2*mk2 - pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mk2,mk2,me2,xl,mu2) * ( 6*pow(mk2,2) - 3
         *mp2*mk2 + 3./8.*pow(mp2,2) );

      massp6RV +=  + hh21V(me2,mk2,mk2,me2,xl,mu2) * ( 6*pow(mk2,2) - 3
         *mp2*mk2 + 3./8.*pow(mp2,2) );

      massp6RV +=  + hh27V(mp2,mk2,mk2,me2,xl,mu2) * (  - 9./2.*mk2 + 9.
         /8.*mp2 );

      massp6RV +=  + hh27V(me2,mk2,mk2,me2,xl,mu2) * (  - 9./2.*mk2 + 9.
         /8.*mp2 );

      return massp6RV/pow(mass.getf0(),4);
}

//**********************************************************************
double fpi4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double F4V = (  AbV(mp2,xl) + 1./2.*AbV(mk2,xl) ); 
  return F4V/pow(mass.getf0(),2);
}


double fpi6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6LloV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 2*L5r*mp2 - 8*L4r*mk2 + 8*L4r*mp2 - 14*L3r*mp2
          - 16*L2r*mp2 - 28*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 2*L5r*mp2 + 12*L4r*mk2 - 2*L4r*
         mp2 - 10*L3r*mk2 - 8*L2r*mk2 - 32*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 2./3.*L5r*mp2 + 16./3.*L4r*mk2 - 
         4./3.*L4r*mp2 - 8./3.*L3r*mk2 + 2./3.*L3r*mp2 - 8./3.*L2r*mk2
          + 2./3.*L2r*mp2 - 32./3.*L1r*mk2 + 8./3.*L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 6*L3r + 24*L2r + 12*L1r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 6*L3r + 24*L2r );

      decayp6LV +=  + A23bV(me2,xl) * ( 2*L3r + 6*L2r );

      decayp6LV +=  + BbV(mp2,xl) * ( 16*L8r*pow(mp2,2) + 32*L6r*mp2*
         mk2 + 16*L6r*pow(mp2,2) - 8*L5r*pow(mp2,2) - 16*L4r*mp2*mk2 - 
         8*L4r*pow(mp2,2) );

      decayp6LV +=  + BbV(mk2,xl) * ( 8*L8r*pow(mk2,2) + 16*L6r*pow(
         mk2,2) + 8*L6r*mp2*mk2 - 4*L5r*pow(mk2,2) - 8*L4r*pow(mk2,2)
          - 4*L4r*mp2*mk2 );

      return decayp6LV/pow(mass.getf0(),4);
}

double fpi6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + AbV(mp2,xl)*BbV(mp2,xl) * (  - 1./2.*mp2 );

      decayp6RV +=  + AbV(mp2,xl) * (  - 1./2.*pi16*mk2 - 1./4.*pi16*
         mp2 - 1./2.*Ab(mp2,mu2) );

      decayp6RV +=  + AbV(mk2,xl) * (  - 1./2.*pi16*mk2 - 1./8.*pi16*
         mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./6.*mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./6.*mk2 );

      decayp6RV +=  + AbV(me2,xl) * ( 1./6.*pi16*mk2 - 1./6.*pi16*mp2
          + 1./6.*Ab(mp2,mu2) - 1./6.*Ab(mk2,mu2) );

      decayp6RV +=  + BbV(mp2,xl) * (  - 1./2.*Ab(mp2,mu2)*mp2 + 1./6.*
         Ab(me2,mu2)*mp2 );

      decayp6RV +=  + BbV(mk2,xl) * (  - 1./6.*Ab(me2,mu2)*mk2 );

      decayp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * (  - 1./2.*mp2 );

      decayp6RV +=  + hhV(mp2,mk2,mk2,mp2,xl,mu2) * (  - 1./2.*mk2 + 1./
         16.*mp2 );

      decayp6RV +=  + hhV(mk2,mk2,me2,mp2,xl,mu2) * (  - 1./4.*mk2 + 1./
         16.*mp2 );

      decayp6RV +=  + hh27V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3./2. );

      decayp6RV +=  + hh27V(mp2,mk2,mk2,mp2,xl,mu2) * (  - 3./16. );

      decayp6RV +=  + hh27V(mk2,mp2,mk2,mp2,xl,mu2) * ( 3./2. );

      decayp6RV +=  + hh27V(me2,mk2,mk2,mp2,xl,mu2) * ( 9./16. );

      decayp6RV +=  + hhdV(mp2,mp2,mp2,mp2,xl,mu2) * ( 5./12.*pow(
         mp2,2) );

      decayp6RV +=  + hhdV(mp2,mk2,mk2,mp2,xl,mu2) * ( 1./2.*mp2*mk2 - 
         5./16.*pow(mp2,2) );

      decayp6RV +=  + hhdV(mp2,me2,me2,mp2,xl,mu2) * ( 1./36.*pow(
         mp2,2) );

      decayp6RV +=  + hhdV(mk2,mk2,me2,mp2,xl,mu2) * ( 1./4.*mp2*mk2 + 
         1./48.*pow(mp2,2) );

      decayp6RV +=  + hh1dV(mp2,mk2,mk2,mp2,xl,mu2) * ( 1./2.*pow(
         mp2,2) );

      decayp6RV +=  + hh1dV(me2,mk2,mk2,mp2,xl,mu2) * (  - 1./2.*pow(
         mp2,2) );

      decayp6RV +=  + hh21dV(mp2,mp2,mp2,mp2,xl,mu2) * ( 3./2.*pow(
         mp2,2) );

      decayp6RV +=  + hh21dV(mp2,mk2,mk2,mp2,xl,mu2) * (  - 3./16.*pow(
         mp2,2) );

      decayp6RV +=  + hh21dV(mk2,mp2,mk2,mp2,xl,mu2) * ( 3./2.*pow(
         mp2,2) );

      decayp6RV +=  + hh21dV(me2,mk2,mk2,mp2,xl,mu2) * ( 9./16.*pow(
         mp2,2) );

      decayp6RV +=  + hh27dV(mp2,mp2,mp2,mp2,xl,mu2) * (  - 3./2.*mp2 )
         ;

      decayp6RV +=  + hh27dV(mp2,mk2,mk2,mp2,xl,mu2) * ( 3./16.*mp2 );

      decayp6RV +=  + hh27dV(mk2,mp2,mk2,mp2,xl,mu2) * (  - 3./2.*mp2 )
         ;

      decayp6RV +=  + hh27dV(me2,mk2,mk2,mp2,xl,mu2) * (  - 9./16.*mp2
          );

      return decayp6RV/pow(mass.getf0(),4);
}

//***********************************************************************

double fk4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double F4V = ( AbV(mp2,xl) * ( 3./8. )
		 + AbV(mk2,xl) * ( 3./4. )
		 + AbV(me2,xl) * ( 3./8. ));
  return F4V/pow(mass.getf0(),2);
}

double fk6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk6LloV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 3./2.*L5r*mk2 - 3*L4r*mk2 + 21./2.*L4r*mp2 - 
         15./2.*L3r*mp2 - 6*L2r*mp2 - 24*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 3*L5r*mk2 + 10*L4r*mk2 - 3*L4r*
         mp2 - 15*L3r*mk2 - 18*L2r*mk2 - 36*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 1./6.*L5r*mk2 + 7./3.*L4r*mk2 - 
         17./6.*L4r*mp2 - 14./3.*L3r*mk2 + 7./6.*L3r*mp2 - 8./3.*L2r*
         mk2 + 2./3.*L2r*mp2 - 32./3.*L1r*mk2 + 8./3.*L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 9./2.*L3r + 18*L2r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 9*L3r + 30*L2r + 12*L1r );

      decayp6LV +=  + A23bV(me2,xl) * ( 1./2.*L3r + 6*L2r );

      decayp6LV +=  + BbV(mp2,xl) * ( 6*L8r*pow(mp2,2) + 12*L6r*mp2*mk2
          + 6*L6r*pow(mp2,2) - 3*L5r*pow(mp2,2) - 6*L4r*mp2*mk2 - 3*L4r
         *pow(mp2,2) );

      decayp6LV +=  + BbV(mk2,xl) * ( 12*L8r*pow(mk2,2) + 24*L6r*pow(
         mk2,2) + 12*L6r*mp2*mk2 - 6*L5r*pow(mk2,2) - 12*L4r*pow(mk2,2)
          - 6*L4r*mp2*mk2 );

      decayp6LV +=  + BbV(me2,xl) * ( 16*L8r*pow(mk2,2) - 16*L8r*mp2*
         mk2 + 6*L8r*pow(mp2,2) + 16*L7r*pow(mk2,2) - 32*L7r*mp2*mk2 + 
         16*L7r*pow(mp2,2) + 16*L6r*pow(mk2,2) + 4*L6r*mp2*mk2 - 2*L6r*
         pow(mp2,2) - 16./3.*L5r*pow(mk2,2) + 8./3.*L5r*mp2*mk2 - 1./3.
         *L5r*pow(mp2,2) - 8*L4r*pow(mk2,2) - 2*L4r*mp2*mk2 + L4r*pow(
         mp2,2) );

      return decayp6LV/pow(mass.getf0(),4);
}

double fk6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + pow(AbV(mp2,xl),2) * (  - 15./128. );

      decayp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * ( 3./32. );

      decayp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * ( 9./64. );

      decayp6RV +=  + AbV(mp2,xl)*Bb(me2,mu2) * ( 3./16.*mp2 );

      decayp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * (  - 3./16.*mp2 );

      decayp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 3./16.*mp2 );

      decayp6RV +=  + AbV(mp2,xl) * (  - 15./32.*pi16*mk2 + 3./16.*pi16
         *mp2 - 27./64.*Ab(mp2,mu2) + 3./32.*Ab(mk2,mu2) + 9./64.*Ab(
         me2,mu2) );

      decayp6RV +=  + pow(AbV(mk2,xl),2) * ( 3./32. );

      decayp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 9./32. );

      decayp6RV +=  + AbV(mk2,xl)*Bb(me2,mu2) * (  - 1./2.*mk2 );

      decayp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 1./2.*mk2 );

      decayp6RV +=  + AbV(mk2,xl) * (  - 9./16.*pi16*mk2 - 3./8.*pi16*
         mp2 + 3./32.*Ab(mp2,mu2) + 3./16.*Ab(mk2,mu2) - 9./32.*Ab(me2,
         mu2) );

      decayp6RV +=  + pow(AbV(me2,xl),2) * ( 9./128. );

      decayp6RV +=  + AbV(me2,xl)*Bb(me2,mu2) * ( 1./3.*mk2 - 7./48.*
         mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./16.*mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./4.*mk2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 1./3.*mk2 - 7./48.*
         mp2 );

      decayp6RV +=  + AbV(me2,xl) * (  - 3./32.*pi16*mk2 - 3./16.*pi16*
         mp2 + 13./64.*Ab(mp2,mu2) - 17./32.*Ab(mk2,mu2) + 9./64.*Ab(
         me2,mu2) );

      decayp6RV +=  + BbV(mp2,xl) * (  - 3./16.*Ab(mp2,mu2)*mp2 + 1./16.
         *Ab(me2,mu2)*mp2 );

      decayp6RV +=  + BbV(mk2,xl) * (  - 1./4.*Ab(me2,mu2)*mk2 );

      decayp6RV +=  + BbV(me2,xl) * ( 3./16.*Ab(mp2,mu2)*mp2 - 1./2.*
         Ab(mk2,mu2)*mk2 + 1./3.*Ab(me2,mu2)*mk2 - 7./48.*Ab(me2,mu2)*
         mp2 );

      decayp6RV +=  + hhV(mp2,mp2,mk2,mk2,xl,mu2) * ( 3./64.*mk2 - 3./8.
         *mp2 );

      decayp6RV +=  + hhV(mp2,mk2,me2,mk2,xl,mu2) * (  - 9./32.*mk2 );

      decayp6RV +=  + hhV(mk2,mk2,mk2,mk2,xl,mu2) * (  - 3./8.*mk2 );

      decayp6RV +=  + hhV(mk2,me2,me2,mk2,xl,mu2) * (  - 9./64.*mk2 );

      decayp6RV +=  + hh27V(mp2,mp2,mk2,mk2,xl,mu2) * ( 9./8. );

      decayp6RV +=  + hh27V(mk2,mp2,mp2,mk2,xl,mu2) * (  - 9./64. );

      decayp6RV +=  + hh27V(mk2,mp2,me2,mk2,xl,mu2) * ( 27./32. );

      decayp6RV +=  + hh27V(mk2,mk2,mk2,mk2,xl,mu2) * ( 9./8. );

      decayp6RV +=  + hh27V(mk2,me2,me2,mk2,xl,mu2) * ( 27./64. );

      decayp6RV +=  + hhdV(mp2,mp2,mk2,mk2,xl,mu2) * (  - 15./64.*pow(
         mk2,2) + 3./8.*mp2*mk2 );

      decayp6RV +=  + hhdV(mp2,mk2,me2,mk2,xl,mu2) * ( 13./32.*pow(
         mk2,2) );

      decayp6RV +=  + hhdV(mk2,mk2,mk2,mk2,xl,mu2) * ( 3./8.*pow(mk2,2)
          );

      decayp6RV +=  + hhdV(mk2,me2,me2,mk2,xl,mu2) * ( 181./576.*pow(
         mk2,2) );

      decayp6RV +=  + hh1dV(mk2,mp2,mp2,mk2,xl,mu2) * ( 3./8.*pow(
         mk2,2) );

      decayp6RV +=  + hh1dV(mk2,mp2,me2,mk2,xl,mu2) * (  - 3./4.*pow(
         mk2,2) );

      decayp6RV +=  + hh1dV(mk2,me2,me2,mk2,xl,mu2) * (  - 5./8.*pow(
         mk2,2) );

      decayp6RV +=  + hh21dV(mp2,mp2,mk2,mk2,xl,mu2) * ( 9./8.*pow(
         mk2,2) );

      decayp6RV +=  + hh21dV(mk2,mp2,mp2,mk2,xl,mu2) * (  - 9./64.*pow(
         mk2,2) );

      decayp6RV +=  + hh21dV(mk2,mp2,me2,mk2,xl,mu2) * ( 27./32.*pow(
         mk2,2) );

      decayp6RV +=  + hh21dV(mk2,mk2,mk2,mk2,xl,mu2) * ( 9./8.*pow(
         mk2,2) );

      decayp6RV +=  + hh21dV(mk2,me2,me2,mk2,xl,mu2) * ( 27./64.*pow(
         mk2,2) );

      decayp6RV +=  + hh27dV(mp2,mp2,mk2,mk2,xl,mu2) * (  - 9./8.*mk2 )
         ;

      decayp6RV +=  + hh27dV(mk2,mp2,mp2,mk2,xl,mu2) * ( 9./64.*mk2 );

      decayp6RV +=  + hh27dV(mk2,mp2,me2,mk2,xl,mu2) * (  - 27./32.*mk2
          );

      decayp6RV +=  + hh27dV(mk2,mk2,mk2,mk2,xl,mu2) * (  - 9./8.*mk2 )
         ;

      decayp6RV +=  + hh27dV(mk2,me2,me2,mk2,xl,mu2) * (  - 27./64.*mk2
          );

      return decayp6RV/pow(mass.getf0(),4);
}

//***********************************************************************

double feta4loV(const lomass mass, const double xl){
  double mk2 = pow(mass.getmk0(),2);
  double F4V =  AbV(mk2,xl) * ( 3./2. );
  return F4V/pow(mass.getf0(),2);
}

double feta6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta6LloV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 2*L5r*mp2 + 12*L4r*mp2 - 6*L3r*mp2 - 6*L2r*mp2
          - 24*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 8./3.*L5r*mk2 - 2*L5r*mp2 + 4*L4r
         *mk2 - 6*L4r*mp2 - 14*L3r*mk2 - 8*L2r*mk2 - 32*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 32./9.*L5r*mk2 - 14./9.*L5r*mp2
          + 16./3.*L4r*mk2 - 4./3.*L4r*mp2 - 8*L3r*mk2 + 2*L3r*mp2 - 16
         *L2r*mk2 + 4*L2r*mp2 - 16*L1r*mk2 + 4*L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 6*L3r + 18*L2r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 2*L3r + 24*L2r );

      decayp6LV +=  + A23bV(me2,xl) * ( 6*L3r + 12*L2r + 12*L1r );

      decayp6LV +=  + BbV(mk2,xl) * ( 24*L8r*pow(mk2,2) + 48*L6r*pow(
         mk2,2) + 24*L6r*mp2*mk2 - 12*L5r*pow(mk2,2) - 24*L4r*pow(
         mk2,2) - 12*L4r*mp2*mk2 );


      return decayp6LV/pow(mass.getf0(),4);
}

double feta6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + AbV(mk2,xl) * (  - 3./2.*pi16*mk2 - 3./8.*pi16*mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./2.*mk2 );

      decayp6RV +=  + AbV(me2,xl) * ( 1./2.*pi16*mk2 - 1./2.*Ab(mk2,mu2
         ) );

      decayp6RV +=  + BbV(mk2,xl) * (  - 1./2.*Ab(me2,mu2)*mk2 );

      decayp6RV +=  + hhV(mp2,mk2,mk2,me2,xl,mu2) * (  - 9./16.*mp2 );

      decayp6RV +=  + hhV(mk2,mk2,me2,me2,xl,mu2) * (  - 3./4.*mk2 + 3./
         16.*mp2 );

      decayp6RV +=  + hh27V(mp2,mk2,mk2,me2,xl,mu2) * ( 27./16. );

      decayp6RV +=  + hh27V(me2,mk2,mk2,me2,xl,mu2) * ( 27./16. );

      decayp6RV +=  + hhdV(mp2,mp2,me2,me2,xl,mu2) * ( 1./12.*pow(
         mp2,2) );

      decayp6RV +=  + hhdV(mp2,mk2,mk2,me2,xl,mu2) * ( 3./4.*mp2*mk2 + 
         1./16.*pow(mp2,2) );

      decayp6RV +=  + hhdV(mk2,mk2,me2,me2,xl,mu2) * ( 25./9.*pow(
         mk2,2) - 11./6.*mp2*mk2 + 5./16.*pow(mp2,2) );

      decayp6RV +=  + hhdV(me2,me2,me2,me2,xl,mu2) * ( 64./243.*pow(
         mk2,2) - 56./243.*mp2*mk2 + 49./972.*pow(mp2,2) );

      decayp6RV +=  + hh1dV(mp2,mk2,mk2,me2,xl,mu2) * (  - 2*mp2*mk2 + 
         1./2.*pow(mp2,2) );

      decayp6RV +=  + hh1dV(me2,mk2,mk2,me2,xl,mu2) * (  - 16./3.*pow(
         mk2,2) + 10./3.*mp2*mk2 - 1./2.*pow(mp2,2) );

      decayp6RV +=  + hh21dV(mp2,mk2,mk2,me2,xl,mu2) * ( 3*pow(mk2,2)
          - 3./2.*mp2*mk2 + 3./16.*pow(mp2,2) );

      decayp6RV +=  + hh21dV(me2,mk2,mk2,me2,xl,mu2) * ( 3*pow(mk2,2)
          - 3./2.*mp2*mk2 + 3./16.*pow(mp2,2) );

      decayp6RV +=  + hh27dV(mp2,mk2,mk2,me2,xl,mu2) * (  - 9./4.*mk2
          + 9./16.*mp2 );

      decayp6RV +=  + hh27dV(me2,mk2,mk2,me2,xl,mu2) * (  - 9./4.*mk2
          + 9./16.*mp2 );

      return decayp6RV/pow(mass.getf0(),4);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double qqup4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = (4.*mk2-mp2)/3.;
  double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  double   P4Vx1 =
       + abvme2 * ( 1./6. );

      P4Vx1 +=  + abvmk2 * ( 1 );

      P4Vx1 +=  + abvmp2 * ( 3./2. );
  return P4Vx1/pow(mass.getf0(),2);
}


double qqup6loV(const lomass mass, const Li Liin, const double xl){
  return qqup6LloV(mass,Liin,xl)+qqup6RloV(mass,xl);
}

double qqup6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup6LloV\n";}
  double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  double bbvmp2 = BbV(mp2,xl);
  double bbvmk2 = BbV(mk2,xl);
  double bbvme2 = BbV(me2,xl);

  double   P6LVx1 =
       + bbvme2 * ( 64./9.*pow(mk2,2)*L8r + 64./9.*pow(mk2,2)*L7r + 64./
         9.*pow(mk2,2)*L6r - 64./27.*pow(mk2,2)*L5r - 32./9.*pow(mk2,2)
         *L4r - 64./9.*mp2*mk2*L8r - 128./9.*mp2*mk2*L7r + 16./9.*mp2*
         mk2*L6r + 32./27.*mp2*mk2*L5r - 8./9.*mp2*mk2*L4r + 8./3.*pow(
         mp2,2)*L8r + 64./9.*pow(mp2,2)*L7r - 8./9.*pow(mp2,2)*L6r - 4./
         27.*pow(mp2,2)*L5r + 4./9.*pow(mp2,2)*L4r );

      P6LVx1 +=  + bbvmk2 * ( 16*pow(mk2,2)*L8r + 32*pow(mk2,2)*L6r - 8
         *pow(mk2,2)*L5r - 16*pow(mk2,2)*L4r + 16*mp2*mk2*L6r - 8*mp2*
         mk2*L4r );

      P6LVx1 +=  + bbvmp2 * ( 48*mp2*mk2*L6r - 24*mp2*mk2*L4r + 24*pow(
         mp2,2)*L8r + 24*pow(mp2,2)*L6r - 12*pow(mp2,2)*L5r - 12*pow(
         mp2,2)*L4r );

      P6LVx1 +=  + abvme2 * (  - 64./3.*mk2*L7r + 80./3.*mk2*L6r - 32./
         9.*mk2*L5r - 40./3.*mk2*L4r + 16./3.*mp2*L8r + 64./3.*mp2*L7r
          - 8./3.*mp2*L6r + 8./9.*mp2*L5r + 4./3.*mp2*L4r );

      P6LVx1 +=  + abvmk2 * ( 32*mk2*L8r + 96*mk2*L6r - 16*mk2*L5r - 48
         *mk2*L4r + 16*mp2*L6r - 8*mp2*L4r );

      P6LVx1 +=  + abvmp2 * ( 48*mk2*L6r - 24*mk2*L4r + 48*mp2*L8r + 72
         *mp2*L6r - 24*mp2*L5r - 36*mp2*L4r );

      return P6LVx1/pow(mass.getf0(),4);
}

double qqup6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  double bbvmp2 = BbV(mp2,xl);
  double bbvmk2 = BbV(mk2,xl);
  double bbvme2 = BbV(me2,xl);

  double   P6RVx1 =
       + abvme2*bbvme2 * ( 4./27.*mk2 - 7./108.*mp2 );

      P6RVx1 +=  + abvme2*bbvmk2 * (  - 1./3.*mk2 );

      P6RVx1 +=  + abvme2*bbvmp2 * ( 1./4.*mp2 );

      P6RVx1 +=  + pow(abvme2,2) * ( 1./72. );

      P6RVx1 +=  + abvmk2*bbvme2 * (  - 2./9.*mk2 );

      P6RVx1 +=  + abvmk2*abvme2 * (  - 1./3. );

      P6RVx1 +=  + abvmp2*bbvme2 * ( 1./12.*mp2 );

      P6RVx1 +=  + abvmp2*bbvmp2 * (  - 3./4.*mp2 );

      P6RVx1 +=  + abvmp2*abvme2 * ( 1./4. );

      P6RVx1 +=  + pow(abvmp2,2) * (  - 3./8. );

      P6RVx1 +=  + Ab(mp2,mu2)*bbvme2 * ( 1./12.*mp2 );

      P6RVx1 +=  + Ab(mp2,mu2)*bbvmp2 * (  - 3./4.*mp2 );

      P6RVx1 +=  + Ab(mp2,mu2)*abvme2 * ( 1./4. );

      P6RVx1 +=  + Ab(mp2,mu2)*abvmp2 * (  - 3./4. );

      P6RVx1 +=  + Ab(mk2,mu2)*bbvme2 * (  - 2./9.*mk2 );

      P6RVx1 +=  + Ab(mk2,mu2)*abvme2 * (  - 1./3. );

      P6RVx1 +=  + Ab(me2,mu2)*bbvme2 * ( 4./27.*mk2 - 7./108.*mp2 );

      P6RVx1 +=  + Ab(me2,mu2)*bbvmk2 * (  - 1./3.*mk2 );

      P6RVx1 +=  + Ab(me2,mu2)*bbvmp2 * ( 1./4.*mp2 );

      P6RVx1 +=  + Ab(me2,mu2)*abvme2 * ( 1./36. );

      P6RVx1 +=  + Ab(me2,mu2)*abvmk2 * (  - 1./3. );

      P6RVx1 +=  + Ab(me2,mu2)*abvmp2 * ( 1./4. );

      P6RVx1 +=  + Bb(mp2,mu2)*abvme2 * ( 1./4.*mp2 );

      P6RVx1 +=  + Bb(mp2,mu2)*abvmp2 * (  - 3./4.*mp2 );

      P6RVx1 +=  + Bb(mk2,mu2)*abvme2 * (  - 1./3.*mk2 );

      P6RVx1 +=  + Bb(me2,mu2)*abvme2 * ( 4./27.*mk2 - 7./108.*mp2 );

      P6RVx1 +=  + Bb(me2,mu2)*abvmk2 * (  - 2./9.*mk2 );

      P6RVx1 +=  + Bb(me2,mu2)*abvmp2 * ( 1./12.*mp2 );

      return P6RVx1/pow(mass.getf0(),4);
}

double qqstrange4loV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = (4.*mk2-mp2)/3.;
  //double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  double  P4Vx3 =
       + abvme2 * ( 2./3. );

      P4Vx3 +=  + abvmk2 * ( 2 );
  return P4Vx3/pow(mass.getf0(),2);
}


double qqstrange6loV(const lomass mass, const Li Liin, const double xl){
  return qqstrange6LloV(mass,Liin,xl)+qqstrange6RloV(mass,xl);
}

double qqstrange6LloV(const lomass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange6LloV\n";}
  double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  //double bbvmp2 = BbV(mp2,xl);
  double bbvmk2 = BbV(mk2,xl);
  double bbvme2 = BbV(me2,xl);

  double  P6LVx3 =
       + bbvme2 * ( 256./9.*pow(mk2,2)*L8r + 256./9.*pow(mk2,2)*L7r + 
         256./9.*pow(mk2,2)*L6r - 256./27.*pow(mk2,2)*L5r - 128./9.*
         pow(mk2,2)*L4r - 256./9.*mp2*mk2*L8r - 512./9.*mp2*mk2*L7r + 
         64./9.*mp2*mk2*L6r + 128./27.*mp2*mk2*L5r - 32./9.*mp2*mk2*L4r
          + 32./3.*pow(mp2,2)*L8r + 256./9.*pow(mp2,2)*L7r - 32./9.*
         pow(mp2,2)*L6r - 16./27.*pow(mp2,2)*L5r + 16./9.*pow(mp2,2)*
         L4r );

      P6LVx3 +=  + bbvmk2 * ( 32*pow(mk2,2)*L8r + 64*pow(mk2,2)*L6r - 
         16*pow(mk2,2)*L5r - 32*pow(mk2,2)*L4r + 32*mp2*mk2*L6r - 16*
         mp2*mk2*L4r );

      P6LVx3 +=  + abvme2 * ( 128./3.*mk2*L8r + 128./3.*mk2*L7r + 128./
         3.*mk2*L6r - 128./9.*mk2*L5r - 64./3.*mk2*L4r - 64./3.*mp2*L8r
          - 128./3.*mp2*L7r + 16./3.*mp2*L6r + 32./9.*mp2*L5r - 8./3.*
         mp2*L4r );

      P6LVx3 +=  + abvmk2 * ( 64*mk2*L8r + 128*mk2*L6r - 32*mk2*L5r - 
         64*mk2*L4r + 32*mp2*L6r - 16*mp2*L4r );

      P6LVx3 +=  + abvmp2 * ( 48*mp2*L6r - 24*mp2*L4r );

      return P6LVx3/pow(mass.getf0(),4);
}

double qqstrange6RloV(const lomass mass, const double xl){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double abvmp2 = AbV(mp2,xl);
  double abvmk2 = AbV(mk2,xl);
  double abvme2 = AbV(me2,xl);
  //double bbvmp2 = BbV(mp2,xl);
  double bbvmk2 = BbV(mk2,xl);
  double bbvme2 = BbV(me2,xl);

  double   P6RVx3 =
       + abvme2*bbvme2 * ( 16./27.*mk2 - 7./27.*mp2 );

      P6RVx3 +=  + abvme2*bbvmk2 * (  - 2./3.*mk2 );

      P6RVx3 +=  + pow(abvme2,2) * ( 2./9. );

      P6RVx3 +=  + abvmk2*bbvme2 * (  - 8./9.*mk2 );

      P6RVx3 +=  + abvmk2*abvme2 * (  - 2./3. );

      P6RVx3 +=  + abvmp2*bbvme2 * ( 1./3.*mp2 );

      P6RVx3 +=  + Ab(mp2,mu2)*bbvme2 * ( 1./3.*mp2 );

      P6RVx3 +=  + Ab(mk2,mu2)*bbvme2 * (  - 8./9.*mk2 );

      P6RVx3 +=  + Ab(mk2,mu2)*abvme2 * (  - 2./3. );

      P6RVx3 +=  + Ab(me2,mu2)*bbvme2 * ( 16./27.*mk2 - 7./27.*mp2 );

      P6RVx3 +=  + Ab(me2,mu2)*bbvmk2 * (  - 2./3.*mk2 );

      P6RVx3 +=  + Ab(me2,mu2)*abvme2 * ( 4./9. );

      P6RVx3 +=  + Ab(me2,mu2)*abvmk2 * (  - 2./3. );

      P6RVx3 +=  + Bb(mk2,mu2)*abvme2 * (  - 2./3.*mk2 );

      P6RVx3 +=  + Bb(me2,mu2)*abvme2 * ( 16./27.*mk2 - 7./27.*mp2 );

      P6RVx3 +=  + Bb(me2,mu2)*abvmk2 * (  - 8./9.*mk2 );

      P6RVx3 +=  + Bb(me2,mu2)*abvmp2 * ( 1./3.*mp2 );

      return P6RVx3/pow(mass.getf0(),4);
}

