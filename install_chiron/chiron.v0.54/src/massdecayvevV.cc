// massdecayvevV.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains the isospin limit expressions for the finite volume
// corrections for masses and decay constants
// in terms of the infinite volume physical masses and Fpi 

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
#include "massdecayvevV.h"

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
#define mpi4V mpi4Vt
#define mpi6LV mpi6LVt
#define mpi6RV mpi6RVt
#define mpi6V mpi6Vt
#define mk4V mk4Vt
#define mk6LV mk6LVt
#define mk6RV mk6RVt
#define mk6V mk6Vt
#define meta4V meta4Vt
#define meta6LV meta6LVt
#define meta6RV meta6RVt
#define meta6V meta6Vt
#define fpi4V fpi4Vt
#define fpi6LV fpi6LVt
#define fpi6RV fpi6RVt
#define fpi6V fpi6Vt
#define fk4V fk4Vt
#define fk6LV fk6LVt
#define fk6RV fk6RVt
#define fk6V fk6Vt
#define feta4V feta4Vt
#define feta6LV feta6LVt
#define feta6RV feta6RVt
#define feta6V feta6Vt
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
#define mpi4V mpi4Vb
#define mpi6LV mpi6LVb
#define mpi6RV mpi6RVb
#define mpi6V mpi6Vb
#define mk4V mk4Vb
#define mk6LV mk6LVb
#define mk6RV mk6RVb
#define mk6V mk6Vb
#define meta4V meta4Vb
#define meta6LV meta6LVb
#define meta6RV meta6RVb
#define meta6V meta6Vb
#define fpi4V fpi4Vb
#define fpi6LV fpi6LVb
#define fpi6RV fpi6RVb
#define fpi6V fpi6Vb
#define fk4V fk4Vb
#define fk6LV fk6LVb
#define fk6RV fk6RVb
#define fk6V fk6Vb
#define feta4V feta4Vb
#define feta6LV feta6LVb
#define feta6RV feta6RVb
#define feta6V feta6Vb
#else
// just some garbage to produce an error
x = 1./0.0.;
#endif
#endif
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double mpi4V(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double me2 = pow(mass.getmeta(),2);
  double M4V = (  - 1./2.*AbV(mp2,xl)*mp2 + 1./6.*AbV(me2,xl)*mp2 ); 
  return M4V/pow(mass.getfpi(),2);
}

double mpi6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double  massp6LV =
       + AbV(mp2,xl) * ( 40*L8r*pow(mp2,2) + 80*L6r*pow(mp2,2) - 24*L5r
         *pow(mp2,2) - 48*L4r*pow(mp2,2) + 28*L3r*pow(mp2,2) + 32*L2r*
         pow(mp2,2) + 56*L1r*pow(mp2,2) );

      massp6LV +=  + AbV(mk2,xl) * ( 32*L8r*mp2*mk2 + 64*L6r*mp2*mk2 - 
         16*L5r*mp2*mk2 - 64*L4r*mp2*mk2 + 20*L3r*mp2*mk2 + 16*L2r*mp2*
         mk2 + 64*L1r*mp2*mk2 );

      massp6LV +=  + AbV(me2,xl) * ( 8*L8r*pow(mp2,2) - 64./3.*L7r*mp2*
         mk2 + 64./3.*L7r*pow(mp2,2) + 64./3.*L6r*mp2*mk2 - 16./3.*L6r*
         pow(mp2,2) - 32./9.*L5r*mp2*mk2 + 8./9.*L5r*pow(mp2,2) - 64./3.
         *L4r*mp2*mk2 + 16./3.*L4r*pow(mp2,2) + 16./3.*L3r*mp2*mk2 - 4./
         3.*L3r*pow(mp2,2) + 16./3.*L2r*mp2*mk2 - 4./3.*L2r*pow(mp2,2)
          + 64./3.*L1r*mp2*mk2 - 16./3.*L1r*pow(mp2,2) );

      massp6LV +=  + A23bV(mp2,xl) * (  - 12*L3r*mp2 - 48*L2r*mp2 - 24*
         L1r*mp2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 12*L3r*mp2 - 48*L2r*mp2 );

      massp6LV +=  + A23bV(me2,xl) * (  - 4*L3r*mp2 - 12*L2r*mp2 );

      return massp6LV/pow(mass.getfpi(),4);
}

double mpi6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 3./8.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * (  - 1./2.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * (  - 1./12.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * ( 1./4.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 1./12.*pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl) * ( pi16*mp2*mk2 + 3./4.*pi16*pow(
         mp2,2) - 7./4.*Ab(mp2,mu2)*mp2 - Ab(mk2,mu2)*mp2 + 1./12.*Bb(
         me2,mu2)*pow(mp2,2) );

      massp6RV +=  + pow(AbV(mk2,xl),2) * (  - 1./4.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 1./2.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 2./9.*mp2*mk2 );

      massp6RV +=  + AbV(mk2,xl) * ( pi16*mp2*mk2 - 1./2.*Ab(mp2,mu2)*
         mp2 - 1./2.*Ab(mk2,mu2)*mp2 - 1./2.*Ab(me2,mu2)*mp2 - 2./9.*
         Bb(me2,mu2)*mp2*mk2 );

      massp6RV +=  + pow(AbV(me2,xl),2) * ( 1./72.*mp2 );

      massp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * (  - 1./12.*pow(mp2,2)
          );

      massp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 4./27.*mp2*mk2 - 7./
         108.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl) * ( 1./12.*pi16*pow(mp2,2) + 1./4.*Ab(
         mp2,mu2)*mp2 - 1./3.*Ab(mk2,mu2)*mp2 + 4./27.*Bb(me2,mu2)*mp2
         *mk2 - 7./108.*Bb(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * ( 5./6.*pow(mp2,2)
          );

      massp6RV +=  + hhV(mp2,mk2,mk2,mp2,xl,mu2) * ( mp2*mk2 - 5./8.*
         pow(mp2,2) );

      massp6RV +=  + hhV(mp2,me2,me2,mp2,xl,mu2) * ( 1./18.*pow(mp2,2)
          );

      massp6RV +=  + hhV(mk2,mk2,me2,mp2,xl,mu2) * ( 1./2.*mp2*mk2 + 1.
         /24.*pow(mp2,2) );

      massp6RV +=  + hh1V(mp2,mk2,mk2,mp2,xl,mu2) * ( pow(mp2,2) );

      massp6RV +=  + hh1V(me2,mk2,mk2,mp2,xl,mu2) * (  - pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mp2,mp2,mp2,xl,mu2) * ( 3*pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mk2,mk2,mp2,xl,mu2) * (  - 3./8.*pow(
         mp2,2) );

      massp6RV +=  + hh21V(mk2,mp2,mk2,mp2,xl,mu2) * ( 3*pow(mp2,2) );

      massp6RV +=  + hh21V(me2,mk2,mk2,mp2,xl,mu2) * ( 9./8.*pow(
         mp2,2) );

      massp6RV +=  + hh27V(mp2,mp2,mp2,mp2,xl,mu2) * (  - 3*mp2 );

      massp6RV +=  + hh27V(mp2,mk2,mk2,mp2,xl,mu2) * ( 3./8.*mp2 );

      massp6RV +=  + hh27V(mk2,mp2,mk2,mp2,xl,mu2) * (  - 3*mp2 );

      massp6RV +=  + hh27V(me2,mk2,mk2,mp2,xl,mu2) * (  - 9./8.*mp2 );

      return massp6RV/pow(mass.getfpi(),4);
}

//***********************************************************************

double mk4V(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double me2 = pow(mass.getmeta(),2);
  double M4V = (  - 1./4.*me2 - 1./12.*mp2 )*AbV(me2,xl);
  return M4V/pow(mass.getfpi(),2);
}

double mk6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double  massp6LV =
       + AbV(mp2,xl) * ( 24*L8r*mp2*mk2 + 48*L6r*mp2*mk2 - 12*L5r*mp2*
         mk2 - 48*L4r*mp2*mk2 + 15*L3r*mp2*mk2 + 12*L2r*mp2*mk2 + 48*
         L1r*mp2*mk2 );

      massp6LV +=  + AbV(mk2,xl) * ( 48*L8r*pow(mk2,2) + 96*L6r*pow(
         mk2,2) - 24*L5r*pow(mk2,2) - 64*L4r*pow(mk2,2) + 30*L3r*pow(
         mk2,2) + 36*L2r*pow(mk2,2) + 72*L1r*pow(mk2,2) );

      massp6LV +=  + AbV(me2,xl) * ( 64./3.*L8r*pow(mk2,2) - 56./3.*L8r
         *mp2*mk2 + 16./3.*L8r*pow(mp2,2) + 64./3.*L7r*pow(mk2,2) - 32*
         L7r*mp2*mk2 + 32./3.*L7r*pow(mp2,2) + 64./3.*L6r*pow(mk2,2) - 
         16./3.*L6r*mp2*mk2 - 64./9.*L5r*pow(mk2,2) + 4./3.*L5r*mp2*mk2
          - 8./9.*L5r*pow(mp2,2) - 64./3.*L4r*pow(mk2,2) + 16./3.*L4r*
         mp2*mk2 + 28./3.*L3r*pow(mk2,2) - 7./3.*L3r*mp2*mk2 + 16./3.*
         L2r*pow(mk2,2) - 4./3.*L2r*mp2*mk2 + 64./3.*L1r*pow(mk2,2) - 
         16./3.*L1r*mp2*mk2 );

      massp6LV +=  + A23bV(mp2,xl) * (  - 9*L3r*mk2 - 36*L2r*mk2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 18*L3r*mk2 - 60*L2r*mk2 - 24*
         L1r*mk2 );

      massp6LV +=  + A23bV(me2,xl) * (  - L3r*mk2 - 12*L2r*mk2 );

      return massp6LV/pow(mass.getfpi(),4);
}

double mk6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double  massp6RV =
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
         *mp2*mk2 - 41./48.*Ab(mp2,mu2)*mk2 + 1./12.*Ab(mp2,mu2)*mp2 - 
         2./3.*Ab(mk2,mu2)*mk2 + 19./48.*Ab(me2,mu2)*mk2 - 1./12.*Ab(
         me2,mu2)*mp2 - 8./27.*Bb(me2,mu2)*pow(mk2,2) + 7./54.*Bb(me2
         ,mu2)*mp2*mk2 );

      massp6RV +=  + hhV(mp2,mp2,mk2,mk2,xl,mu2) * (  - 15./32.*pow(
         mk2,2) + 3./4.*mp2*mk2 );

      massp6RV +=  + hhV(mp2,mk2,me2,mk2,xl,mu2) * ( 13./16.*pow(
         mk2,2) );

      massp6RV +=  + hhV(mk2,mk2,mk2,mk2,xl,mu2) * ( 3./4.*pow(mk2,2)
          );

      massp6RV +=  + hhV(mk2,me2,me2,mk2,xl,mu2) * ( 181./288.*pow(
         mk2,2) );

      massp6RV +=  + hh1V(mk2,mp2,mp2,mk2,xl,mu2) * ( 3./4.*pow(mk2,2)
          );

      massp6RV +=  + hh1V(mk2,mp2,me2,mk2,xl,mu2) * (  - 3./2.*pow(
         mk2,2) );

      massp6RV +=  + hh1V(mk2,me2,me2,mk2,xl,mu2) * (  - 5./4.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mp2,mp2,mk2,mk2,xl,mu2) * ( 9./4.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,mp2,mp2,mk2,xl,mu2) * (  - 9./32.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,mp2,me2,mk2,xl,mu2) * ( 27./16.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,mk2,mk2,mk2,xl,mu2) * ( 9./4.*pow(
         mk2,2) );

      massp6RV +=  + hh21V(mk2,me2,me2,mk2,xl,mu2) * ( 27./32.*pow(
         mk2,2) );

      massp6RV +=  + hh27V(mp2,mp2,mk2,mk2,xl,mu2) * (  - 9./4.*mk2 );

      massp6RV +=  + hh27V(mk2,mp2,mp2,mk2,xl,mu2) * ( 9./32.*mk2 );

      massp6RV +=  + hh27V(mk2,mp2,me2,mk2,xl,mu2) * (  - 27./16.*mk2
          );

      massp6RV +=  + hh27V(mk2,mk2,mk2,mk2,xl,mu2) * (  - 9./4.*mk2 );

      massp6RV +=  + hh27V(mk2,me2,me2,mk2,xl,mu2) * (  - 27./32.*mk2
          );

      return massp6RV/pow(mass.getfpi(),4);
}

//***********************************************************************

double meta4V(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double M4V =   1./2.*AbV(mp2,xl)*mp2 + AbV(mk2,xl)*( - me2 - 1./3.*mp2)
    + AbV(me2,xl)*(8./9.*mk2 - 7./18.*mp2 );
  return M4V/pow(mass.getfpi(),2);
}

double meta6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double   massp6LV =
       + AbV(mp2,xl) * ( 24*L8r*pow(mp2,2) - 64*L7r*mp2*mk2 + 64*L7r*
         pow(mp2,2) + 64*L6r*mp2*mk2 - 16*L6r*pow(mp2,2) - 32./3.*L5r*
         mp2*mk2 + 8./3.*L5r*pow(mp2,2) - 64*L4r*mp2*mk2 + 16*L4r*pow(
         mp2,2) + 16*L3r*mp2*mk2 - 4*L3r*pow(mp2,2) + 16*L2r*mp2*mk2 - 
         4*L2r*pow(mp2,2) + 64*L1r*mp2*mk2 - 16*L1r*pow(mp2,2) );

      massp6LV +=  + AbV(mk2,xl) * ( 256./3.*L8r*pow(mk2,2) - 224./3.*
         L8r*mp2*mk2 + 64./3.*L8r*pow(mp2,2) + 256./3.*L7r*pow(mk2,2)
          - 128*L7r*mp2*mk2 + 128./3.*L7r*pow(mp2,2) + 256./3.*L6r*pow(
         mk2,2) - 64./3.*L6r*mp2*mk2 - 256./9.*L5r*pow(mk2,2) + 16./3.*
         L5r*mp2*mk2 - 32./9.*L5r*pow(mp2,2) - 256./3.*L4r*pow(mk2,2)
          + 64./3.*L4r*mp2*mk2 + 112./3.*L3r*pow(mk2,2) - 28./3.*L3r*
         mp2*mk2 + 64./3.*L2r*pow(mk2,2) - 16./3.*L2r*mp2*mk2 + 256./3.
         *L1r*pow(mk2,2) - 64./3.*L1r*mp2*mk2 );

      massp6LV +=  + AbV(me2,xl) * ( 896./9.*L8r*pow(mk2,2) - 1024./9.*
         L8r*mp2*mk2 + 344./9.*L8r*pow(mp2,2) + 1024./9.*L7r*pow(mk2,2)
          - 1664./9.*L7r*mp2*mk2 + 640./9.*L7r*pow(mp2,2) + 256./3.*L6r
         *pow(mk2,2) - 128./3.*L6r*mp2*mk2 + 16./3.*L6r*pow(mp2,2) - 
         832./27.*L5r*pow(mk2,2) + 896./27.*L5r*mp2*mk2 - 280./27.*L5r*
         pow(mp2,2) - 256./9.*L4r*pow(mk2,2) + 128./9.*L4r*mp2*mk2 - 16.
         /9.*L4r*pow(mp2,2) + 64./3.*L3r*pow(mk2,2) - 32./3.*L3r*mp2*
         mk2 + 4./3.*L3r*pow(mp2,2) + 128./3.*L2r*pow(mk2,2) - 64./3.*
         L2r*mp2*mk2 + 8./3.*L2r*pow(mp2,2) + 128./3.*L1r*pow(mk2,2) - 
         64./3.*L1r*mp2*mk2 + 8./3.*L1r*pow(mp2,2) );

      massp6LV +=  + A23bV(mp2,xl) * (  - 16*L3r*mk2 + 4*L3r*mp2 - 48*
         L2r*mk2 + 12*L2r*mp2 );

      massp6LV +=  + A23bV(mk2,xl) * (  - 16./3.*L3r*mk2 + 4./3.*L3r*
         mp2 - 64*L2r*mk2 + 16*L2r*mp2 );

      massp6LV +=  + A23bV(me2,xl) * (  - 16*L3r*mk2 + 4*L3r*mp2 - 32*
         L2r*mk2 + 8*L2r*mp2 - 32*L1r*mk2 + 8*L1r*mp2 );

      return massp6LV/pow(mass.getfpi(),4);
}

double meta6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double massp6RV =
       + pow(AbV(mp2,xl),2) * (  - 1./8.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * (  - 3./2.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * ( 1./12.*mp2 );

      massp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * (  - 1./4.*pow(mp2,2) )
         ;

      massp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 4./9.*mp2*mk2 - 7./36.
         *pow(mp2,2) );

      massp6RV +=  + AbV(mp2,xl) * ( 1./4.*pi16*pow(mp2,2) + 3./4.*Ab(
         mp2,mu2)*mp2 - Ab(mk2,mu2)*mp2 + 4./9.*Bb(me2,mu2)*mp2*mk2 - 
         7./36.*Bb(me2,mu2)*pow(mp2,2) );

      massp6RV +=  + pow(AbV(mk2,xl),2) * ( mk2 + 3./4.*mp2 );

      massp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 32./9.*mk2 + 3./2.*
         mp2 );

      massp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 32./27.*pow(mk2,2)
          + 14./27.*mp2*mk2 );

      massp6RV +=  + AbV(mk2,xl) * ( 8./3.*pi16*pow(mk2,2) + 2./3.*pi16
         *mp2*mk2 - 1./3.*pi16*pow(mp2,2) - 8./3.*Ab(mp2,mu2)*mk2 - 7./
         6.*Ab(mp2,mu2)*mp2 - 2./3.*Ab(mk2,mu2)*mk2 + 3./2.*Ab(mk2,mu2)
         *mp2 - 8./3.*Ab(me2,mu2)*mk2 + 7./6.*Ab(me2,mu2)*mp2 - 32./27.
         *Bb(me2,mu2)*pow(mk2,2) + 14./27.*Bb(me2,mu2)*mp2*mk2 );

      massp6RV +=  + pow(AbV(me2,xl),2) * ( 8./27.*mk2 - 31./216.*mp2 )
         ;

      massp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./12.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * ( 4./9.*pow(mk2,2) );

      massp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 64./81.*pow(mk2,2) - 
         56./81.*mp2*mk2 + 49./324.*pow(mp2,2) );

      massp6RV +=  + AbV(me2,xl) * (  - 4./9.*pi16*pow(mk2,2) - 1./12.*
         pi16*pow(mp2,2) + 16./9.*Ab(mp2,mu2)*mk2 - 29./36.*Ab(mp2,mu2)
         *mp2 - 20./9.*Ab(mk2,mu2)*mk2 + 10./9.*Ab(mk2,mu2)*mp2 + 8./9.
         *Ab(me2,mu2)*mk2 - 2./9.*Ab(me2,mu2)*mp2 + 64./81.*Bb(me2,mu2
         )*pow(mk2,2) - 56./81.*Bb(me2,mu2)*mp2*mk2 + 49./324.*Bb(me2
         ,mu2)*pow(mp2,2) );

      massp6RV +=  + hhV(mp2,mp2,me2,me2,xl,mu2) * ( 1./6.*pow(mp2,2)
          );

      massp6RV +=  + hhV(mp2,mk2,mk2,me2,xl,mu2) * ( 3./2.*mp2*mk2 + 1.
         /8.*pow(mp2,2) );

      massp6RV +=  + hhV(mk2,mk2,me2,me2,xl,mu2) * ( 50./9.*pow(mk2,2)
          - 11./3.*mp2*mk2 + 5./8.*pow(mp2,2) );

      massp6RV +=  + hhV(me2,me2,me2,me2,xl,mu2) * ( 128./243.*pow(
         mk2,2) - 112./243.*mp2*mk2 + 49./486.*pow(mp2,2) );

      massp6RV +=  + hh1V(mp2,mk2,mk2,me2,xl,mu2) * (  - 4*mp2*mk2 + 
         pow(mp2,2) );

      massp6RV +=  + hh1V(me2,mk2,mk2,me2,xl,mu2) * (  - 32./3.*pow(
         mk2,2) + 20./3.*mp2*mk2 - pow(mp2,2) );

      massp6RV +=  + hh21V(mp2,mk2,mk2,me2,xl,mu2) * ( 6*pow(mk2,2) - 
         3*mp2*mk2 + 3./8.*pow(mp2,2) );

      massp6RV +=  + hh21V(me2,mk2,mk2,me2,xl,mu2) * ( 6*pow(mk2,2) - 
         3*mp2*mk2 + 3./8.*pow(mp2,2) );

      massp6RV +=  + hh27V(mp2,mk2,mk2,me2,xl,mu2) * (  - 9./2.*mk2 + 
         9./8.*mp2 );

      massp6RV +=  + hh27V(me2,mk2,mk2,me2,xl,mu2) * (  - 9./2.*mk2 + 
         9./8.*mp2 );

      return massp6RV/pow(mass.getfpi(),4);
}

//**********************************************************************
double fpi4V(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double F4V = (  AbV(mp2,xl) + 1./2.*AbV(mk2,xl) ); 
  return F4V/pow(mass.getfpi(),1);
}

double fpi6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 6*L5r*mp2 + 12*L4r*mp2 - 14*L3r*mp2 - 16*L2r*
         mp2 - 28*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 4*L5r*mp2 + 16*L4r*mk2 - 10*L3r*
         mk2 - 8*L2r*mk2 - 32*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 2./3.*L5r*mp2 + 16./3.*L4r*mk2 - 
         4./3.*L4r*mp2 - 8./3.*L3r*mk2 + 2./3.*L3r*mp2 - 8./3.*L2r*mk2
          + 2./3.*L2r*mp2 - 32./3.*L1r*mk2 + 8./3.*L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 6*L3r + 24*L2r + 12*L1r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 6*L3r + 24*L2r );

      decayp6LV +=  + A23bV(me2,xl) * ( 2*L3r + 6*L2r );

      return decayp6LV/pow(mass.getfpi(),3);
}

double fpi6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + AbV(mp2,xl)*BbV(mp2,xl) * (  - 1./2.*mp2 );

      decayp6RV +=  + AbV(mp2,xl) * (  - 1./2.*pi16*mk2 - 1./4.*pi16*
         mp2 + 1./2.*Ab(mp2,mu2) + 1./2.*Ab(mk2,mu2) );

      decayp6RV +=  + AbV(mk2,xl) * (  - 1./2.*pi16*mk2 - 1./8.*pi16*
         mp2 + 1./2.*Ab(mp2,mu2) + 1./4.*Ab(mk2,mu2) );

      decayp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./6.*mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./6.*mk2 );

      decayp6RV +=  + AbV(me2,xl) * ( 1./6.*pi16*mk2 - 1./6.*pi16*mp2
          + 1./6.*Ab(mp2,mu2) - 1./6.*Ab(mk2,mu2) );

      decayp6RV +=  + hhV(mp2,mp2,mp2,mp2,xl,mu2) * (  - 1./2.*mp2 );

      decayp6RV +=  + hhV(mp2,mk2,mk2,mp2,xl,mu2) * (  - 1./2.*mk2 + 1.
         /16.*mp2 );

      decayp6RV +=  + hhV(mk2,mk2,me2,mp2,xl,mu2) * (  - 1./4.*mk2 + 1.
         /16.*mp2 );

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

      return decayp6RV/pow(mass.getfpi(),3);
}

//***********************************************************************

double fk4V(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double F4V = ( AbV(mp2,xl) * ( 3./8. )
		 + AbV(mk2,xl) * ( 3./4. )
		 + AbV(me2,xl) * ( 3./8. ));
  return F4V/pow(mass.getfpi(),1);
}

double fk6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 3./2.*L5r*mk2 + 3./2.*L5r*mp2 + 12*L4r*mp2 - 
         15./2.*L3r*mp2 - 6*L2r*mp2 - 24*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 3*L5r*mk2 + 3*L5r*mp2 + 16*L4r*
         mk2 - 15*L3r*mk2 - 18*L2r*mk2 - 36*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 1./6.*L5r*mk2 + 3./2.*L5r*mp2 + 
         16./3.*L4r*mk2 - 4./3.*L4r*mp2 - 14./3.*L3r*mk2 + 7./6.*L3r*
         mp2 - 8./3.*L2r*mk2 + 2./3.*L2r*mp2 - 32./3.*L1r*mk2 + 8./3.*
         L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 9./2.*L3r + 18*L2r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 9*L3r + 30*L2r + 12*L1r );

      decayp6LV +=  + A23bV(me2,xl) * ( 1./2.*L3r + 6*L2r );

      return decayp6LV/pow(mass.getfpi(),3);
}

double fk6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + pow(AbV(mp2,xl),2) * (  - 15./128. );

      decayp6RV +=  + AbV(mp2,xl)*AbV(mk2,xl) * ( 3./32. );

      decayp6RV +=  + AbV(mp2,xl)*AbV(me2,xl) * ( 9./64. );

      decayp6RV +=  + AbV(mp2,xl)*BbV(mp2,xl) * (  - 3./16.*mp2 );

      decayp6RV +=  + AbV(mp2,xl)*BbV(me2,xl) * ( 3./16.*mp2 );

      decayp6RV +=  + AbV(mp2,xl) * (  - 15./32.*pi16*mk2 + 3./16.*pi16
         *mp2 - 3./64.*Ab(mp2,mu2) + 9./32.*Ab(mk2,mu2) + 9./64.*Ab(me2
         ,mu2) + 3./16.*Bb(me2,mu2)*mp2 );

      decayp6RV +=  + pow(AbV(mk2,xl),2) * ( 3./32. );

      decayp6RV +=  + AbV(mk2,xl)*AbV(me2,xl) * (  - 9./32. );

      decayp6RV +=  + AbV(mk2,xl)*BbV(me2,xl) * (  - 1./2.*mk2 );

      decayp6RV +=  + AbV(mk2,xl) * (  - 9./16.*pi16*mk2 - 3./8.*pi16*
         mp2 + 27./32.*Ab(mp2,mu2) + 9./16.*Ab(mk2,mu2) - 9./32.*Ab(me2
         ,mu2) - 1./2.*Bb(me2,mu2)*mk2 );

      decayp6RV +=  + pow(AbV(me2,xl),2) * ( 9./128. );

      decayp6RV +=  + AbV(me2,xl)*BbV(mp2,xl) * ( 1./16.*mp2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./4.*mk2 );

      decayp6RV +=  + AbV(me2,xl)*BbV(me2,xl) * ( 1./3.*mk2 - 7./48.*
         mp2 );

      decayp6RV +=  + AbV(me2,xl) * (  - 3./32.*pi16*mk2 - 3./16.*pi16*
         mp2 + 37./64.*Ab(mp2,mu2) - 11./32.*Ab(mk2,mu2) + 9./64.*Ab(
         me2,mu2) + 1./3.*Bb(me2,mu2)*mk2 - 7./48.*Bb(me2,mu2)*mp2 );

      decayp6RV +=  + hhV(mp2,mp2,mk2,mk2,xl,mu2) * ( 3./64.*mk2 - 3./
         8.*mp2 );

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

      return decayp6RV/pow(mass.getfpi(),3);
}

//***********************************************************************

double feta4V(const physmass mass, const double xl){
  double mk2 = pow(mass.getmk(),2);
  double F4V =  AbV(mk2,xl) * ( 3./2. );
  return F4V/pow(mass.getfpi(),1);
}

double feta6LV(const physmass mass, const Li Liin, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  //double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6LV\n";}

  double decayp6LV =
       + AbV(mp2,xl) * ( 2*L5r*mp2 + 12*L4r*mp2 - 6*L3r*mp2 - 6*L2r*mp2
          - 24*L1r*mp2 );

      decayp6LV +=  + AbV(mk2,xl) * ( 8./3.*L5r*mk2 + 4*L5r*mp2 + 16*
         L4r*mk2 - 14*L3r*mk2 - 8*L2r*mk2 - 32*L1r*mk2 );

      decayp6LV +=  + AbV(me2,xl) * ( 32./9.*L5r*mk2 - 14./9.*L5r*mp2
          + 16./3.*L4r*mk2 - 4./3.*L4r*mp2 - 8*L3r*mk2 + 2*L3r*mp2 - 16
         *L2r*mk2 + 4*L2r*mp2 - 16*L1r*mk2 + 4*L1r*mp2 );

      decayp6LV +=  + A23bV(mp2,xl) * ( 6*L3r + 18*L2r );

      decayp6LV +=  + A23bV(mk2,xl) * ( 2*L3r + 24*L2r );

      decayp6LV +=  + A23bV(me2,xl) * ( 6*L3r + 12*L2r + 12*L1r );


      return decayp6LV/pow(mass.getfpi(),3);
}

double feta6RV(const physmass mass, const double xl){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double decayp6RV =
       + AbV(mk2,xl) * (  - 3./2.*pi16*mk2 - 3./8.*pi16*mp2 + 3./2.*Ab(
         mp2,mu2) + 3./4.*Ab(mk2,mu2) );

      decayp6RV +=  + AbV(me2,xl)*BbV(mk2,xl) * (  - 1./2.*mk2 );

      decayp6RV +=  + AbV(me2,xl) * ( 1./2.*pi16*mk2 - 1./2.*Ab(mk2,mu2
         ) );

      decayp6RV +=  + hhV(mp2,mk2,mk2,me2,xl,mu2) * (  - 9./16.*mp2 );

      decayp6RV +=  + hhV(mk2,mk2,me2,me2,xl,mu2) * (  - 3./4.*mk2 + 3.
         /16.*mp2 );

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

      return decayp6RV/pow(mass.getfpi(),3);
}
