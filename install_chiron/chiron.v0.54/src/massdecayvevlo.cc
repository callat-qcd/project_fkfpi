// massdecayvevlo.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results derived in
//  G.~Amoros, J.~Bijnens and P.~Talavera,
//  Two point functions at two loops in three flavor chiral perturbation theory,
//  Nucl.\ Phys.\ B {\bf 568} (2000) 319 [hep-ph/9907264]
//
// contains the isospin limit expressions for masses, decay
// constants and vacuum expectation values 
// in terms of the lowest order masses and decay constnats

#include <iostream>
#include <iomanip>
#include <string>

#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "oneloopintegrals.h"
#include "sunsetintegrals.h"
#include "massdecayvevlo.h"

const double pi = M_PI;
const double pi16 = 1./(16.*pi*pi);

// mpi ******************************************************************

double mpi4lo(const lomass mass, const Li Liin){
  return mpi4Llo(mass,Liin)+mpi4Rlo(mass);
}

double mpi4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mpi4L =
       + 32*mp2*mk2*L6r - 16*mp2*mk2*L4r + 16*pow(mp2,2)*L8r + 16*pow(
         mp2,2)*L6r - 8*pow(mp2,2)*L5r - 8*pow(mp2,2)*L4r;
  return mpi4L/pow(mass.getf0(),2);
}

double mpi4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double mpi4R =
       + Ab(mp2,mu2) * (  - 1./2.*mp2 );

      mpi4R +=  + Ab(me2,mu2) * ( 1./6.*mp2 );

  return mpi4R/pow(mass.getf0(),2);
}

double mpi6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return mpi6Rlo(mass)+mpi6Llo(mass,Liin)+mpi6Clo(mass,Ciin);
}

double mpi6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return mpi6Llo(mass,Liin)+mpi6Clo(mass,Ciin);
}

double mpi6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6Llo\n";}

  double    mpi6L =
       + Ab(mp2,mu2) * (  - 32*mp2*mk2*L6r + 24*mp2*mk2*L4r + 24*pow(
         mp2,2)*L8r + 64*pow(mp2,2)*L6r - 12*pow(mp2,2)*L5r - 36*pow(
         mp2,2)*L4r + 28*pow(mp2,2)*L3r + 32*pow(mp2,2)*L2r + 56*pow(
         mp2,2)*L1r );

      mpi6L +=  + Ab(mk2,mu2) * ( 32*mp2*mk2*L8r + 64*mp2*mk2*L6r - 16*
         mp2*mk2*L5r - 64*mp2*mk2*L4r + 20*mp2*mk2*L3r + 16*mp2*mk2*L2r
          + 64*mp2*mk2*L1r );

      mpi6L +=  + Ab(me2,mu2) * ( 16./3.*mp2*mk2*L8r - 16*mp2*mk2*L7r
          + 32*mp2*mk2*L6r - 16./3.*mp2*mk2*L5r - 88./3.*mp2*mk2*L4r + 
         16./3.*mp2*mk2*L3r + 16./3.*mp2*mk2*L2r + 64./3.*mp2*mk2*L1r
          + 20./3.*pow(mp2,2)*L8r + 12*pow(mp2,2)*L7r - 4./3.*pow(
         mp2,2)*L5r + 4./3.*pow(mp2,2)*L4r - 4./3.*pow(mp2,2)*L3r - 4./
         3.*pow(mp2,2)*L2r - 16./3.*pow(mp2,2)*L1r + 4./3.*pow(mp2,3)*
         pow(me2,-1)*L8r + 4*pow(mp2,3)*pow(me2,-1)*L7r );

      mpi6L +=  - 512*mp2*pow(mk2,2)*L4r*L6r + 256*mp2*pow(mk2,2)*pow(
         L4r,2) - 64./9.*mp2*pow(mk2,2)*pi16*L8r - 64./9.*mp2*pow(
         mk2,2)*pi16*L7r - 64./9.*mp2*pow(mk2,2)*pi16*L6r + 64./27.*mp2
         *pow(mk2,2)*pi16*L5r + 32./9.*mp2*pow(mk2,2)*pi16*L4r + 86./27.
         *mp2*pow(mk2,2)*pi16*L3r + 104./9.*mp2*pow(mk2,2)*pi16*L2r - 
         256*pow(mp2,2)*mk2*L5r*L6r - 256*pow(mp2,2)*mk2*L4r*L8r - 512*
         pow(mp2,2)*mk2*L4r*L6r + 256*pow(mp2,2)*mk2*L4r*L5r + 256*pow(
         mp2,2)*mk2*pow(L4r,2) + 64./9.*pow(mp2,2)*mk2*pi16*L8r + 128./
         9.*pow(mp2,2)*mk2*pi16*L7r + 128./9.*pow(mp2,2)*mk2*pi16*L6r
          - 32./27.*pow(mp2,2)*mk2*pi16*L5r - 64./9.*pow(mp2,2)*mk2*
         pi16*L4r - 16./27.*pow(mp2,2)*mk2*pi16*L3r - 16./9.*pow(mp2,2)
         *mk2*pi16*L2r - 128*pow(mp2,3)*L5r*L8r - 128*pow(mp2,3)*L5r*
         L6r + 64*pow(mp2,3)*pow(L5r,2) - 128*pow(mp2,3)*L4r*L8r - 128*
         pow(mp2,3)*L4r*L6r + 128*pow(mp2,3)*L4r*L5r + 64*pow(mp2,3)*
         pow(L4r,2) + 16./3.*pow(mp2,3)*pi16*L8r - 64./9.*pow(mp2,3)*
         pi16*L7r;
      mpi6L +=  + 80./9.*pow(mp2,3)*pi16*L6r - 104./27.*pow(mp2,3)*pi16
         *L5r - 40./9.*pow(mp2,3)*pi16*L4r + 56./27.*pow(mp2,3)*pi16*
         L3r + 74./9.*pow(mp2,3)*pi16*L2r + 4*pow(mp2,3)*pi16*L1r;

      return mpi6L/pow(mass.getf0(),4);
}

double mpi6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6Clo\n";}


  double    mpi6C =
       - 32*CC[12]*pow(mp2,3) - 64*CC[13]*pow(mp2,2)*mk2 - 32*CC[13]*
         pow(mp2,3) - 16*CC[14]*pow(mp2,3) - 32*CC[15]*pow(mp2,2)*mk2
          - 16*CC[15]*pow(mp2,3) - 64*CC[16]*mp2*pow(mk2,2) + 64*CC[16]
         *pow(mp2,2)*mk2 - 48*CC[16]*pow(mp2,3) - 16*CC[17]*pow(mp2,3)
          + 48*CC[19]*pow(mp2,3) + 64*CC[20]*mp2*pow(mk2,2) + 80*CC[20]
         *pow(mp2,3) + 192*CC[21]*mp2*pow(mk2,2) + 192*CC[21]*pow(
         mp2,2)*mk2 + 48*CC[21]*pow(mp2,3) + 32*CC[31]*pow(mp2,3) + 64*
         CC[32]*pow(mp2,2)*mk2 + 32*CC[32]*pow(mp2,3);

  return mpi6C/pow(mass.getf0(),4);
}

double mpi6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double mpi6R =
       + pow(Ab(mp2,mu2),2) * (  - 1./2.*mk2 - 373./144.*mp2 );

      mpi6R +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - 1./2.*mp2 );

      mpi6R +=  + Ab(mp2,mu2)*Ab(me2,mu2) * (  - 1./6.*mp2 + 1./12.*
         pow(mp2,2)*pow(me2,-1) );

      mpi6R +=  + Ab(mp2,mu2) * ( mp2*mk2*pi16 + 2./3.*pow(mp2,2)*pi16
          );

      mpi6R +=  + pow(Ab(mk2,mu2),2) * (  - 7./4.*mp2 - 1./6.*pow(
         mp2,2)*pow(mk2,-1) );

      mpi6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 2./3.*mp2 - 1./18.*
         pow(mp2,2)*pow(me2,-1) );

      mpi6R +=  + Ab(mk2,mu2) * ( 11./9.*mp2*mk2*pi16 );

      mpi6R +=  + pow(Ab(me2,mu2),2) * (  - 1./16.*mp2 - 1./6.*pow(
         mp2,2)*pow(me2,-1) );

      mpi6R +=  + Ab(me2,mu2) * (  - 4./27.*mp2*mk2*pi16 + 4./27.*pow(
         mp2,2)*pi16 );

      mpi6R +=  + hh(mp2,mp2,mp2,mp2,mu2) * ( 5./6.*pow(mp2,2) );

      mpi6R +=  + hh(mp2,mk2,mk2,mp2,mu2) * (  - 5./8.*pow(mp2,2) );

      mpi6R +=  + hh(mp2,me2,me2,mp2,mu2) * ( 1./18.*pow(mp2,2) );

      mpi6R +=  + hh(mk2,mp2,mk2,mp2,mu2) * ( mp2*mk2 );

      mpi6R +=  + hh(mk2,mk2,me2,mp2,mu2) * (  - 5./6.*pow(mp2,2) );

      mpi6R +=  + hh(me2,mk2,mk2,mp2,mu2) * ( 1./2.*mp2*mk2 - 1./8.*
         pow(mp2,2) );

      mpi6R +=  + hh1(mp2,mk2,mk2,mp2,mu2) * ( pow(mp2,2) );

      mpi6R +=  + hh1(mk2,mk2,me2,mp2,mu2) * ( 2*pow(mp2,2) );

      mpi6R +=  + hh21(mp2,mp2,mp2,mp2,mu2) * ( 3*pow(mp2,2) );

      mpi6R +=  + hh21(mp2,mk2,mk2,mp2,mu2) * (  - 3./8.*pow(mp2,2) );

      mpi6R +=  + hh21(mk2,mp2,mk2,mp2,mu2) * ( 3*pow(mp2,2) );

      mpi6R +=  + hh21(me2,mk2,mk2,mp2,mu2) * ( 9./8.*pow(mp2,2) );

      mpi6R +=  - 15./16.*mp2*pow(mk2,2)*pow(pi16,2) - 11./36.*mp2*pow(
         mk2,2)*pow(pi16,2)*pow(pi,2) - 139./216.*pow(mp2,2)*mk2*pow(
         pi16,2) - 37./324.*pow(mp2,2)*mk2*pow(pi16,2)*pow(pi,2) - 3217.
         /1728.*pow(mp2,3)*pow(pi16,2) - 527./1296.*pow(mp2,3)*pow(
         pi16,2)*pow(pi,2);

      return mpi6R/pow(mass.getf0(),4);
}

// *************************************************************************

double mk4lo(const lomass mass, const Li Liin){
  return mk4Llo(mass,Liin)+mk4Rlo(mass);
}

double mk4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mk4L =
       + 16*pow(mk2,2)*L8r + 32*pow(mk2,2)*L6r - 8*pow(mk2,2)*L5r - 16*
         pow(mk2,2)*L4r + 16*mp2*mk2*L6r - 8*mp2*mk2*L4r;

  return mk4L/pow(mass.getf0(),2);
}

double mk4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double mk4R =
       + Ab(me2,mu2) * (  - 1./3.*mk2 );

  return mk4R/pow(mass.getf0(),2);
}

double mk6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return mk6Rlo(mass)+mk6Llo(mass,Liin)+mk6Clo(mass,Ciin);
}

double mk6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return mk6Llo(mass,Liin)+mk6Clo(mass,Ciin);
}

double mk6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk6Llo\n";}

  double mk6L =
       + Ab(mp2,mu2) * ( 24*mp2*mk2*L8r + 48*mp2*mk2*L6r - 12*mp2*mk2*
         L5r - 48*mp2*mk2*L4r + 15*mp2*mk2*L3r + 12*mp2*mk2*L2r + 48*
         mp2*mk2*L1r );

      mk6L +=  + Ab(mk2,mu2) * ( 48*pow(mk2,2)*L8r + 96*pow(mk2,2)*L6r
          - 24*pow(mk2,2)*L5r - 64*pow(mk2,2)*L4r + 30*pow(mk2,2)*L3r
          + 36*pow(mk2,2)*L2r + 72*pow(mk2,2)*L1r );

      mk6L +=  + Ab(me2,mu2) * (  - 16./3.*pow(mk2,2)*L4r + 28./3.*pow(
         mk2,2)*L3r + 16./3.*pow(mk2,2)*L2r + 64./3.*pow(mk2,2)*L1r + 8
         *mp2*mk2*L7r - 16*mp2*mk2*L6r + 4./3.*mp2*mk2*L5r + 40./3.*mp2
         *mk2*L4r - 7./3.*mp2*mk2*L3r - 4./3.*mp2*mk2*L2r - 16./3.*mp2*
         mk2*L1r - 2*pow(mp2,2)*L8r - 6*pow(mp2,2)*L7r - 2./3.*pow(
         mp2,3)*pow(me2,-1)*L8r - 2*pow(mp2,3)*pow(me2,-1)*L7r );

      mk6L +=  - 128*pow(mk2,3)*L5r*L8r - 256*pow(mk2,3)*L5r*L6r + 64*
         pow(mk2,3)*pow(L5r,2) - 256*pow(mk2,3)*L4r*L8r - 512*pow(
         mk2,3)*L4r*L6r + 256*pow(mk2,3)*L4r*L5r + 256*pow(mk2,3)*pow(
         L4r,2) + 128./9.*pow(mk2,3)*pi16*L8r + 128./9.*pow(mk2,3)*pi16
         *L7r + 128./9.*pow(mk2,3)*pi16*L6r - 128./27.*pow(mk2,3)*pi16*
         L5r - 64./9.*pow(mk2,3)*pi16*L4r + 89./27.*pow(mk2,3)*pi16*L3r
          + 122./9.*pow(mk2,3)*pi16*L2r + 4*pow(mk2,3)*pi16*L1r - 128*
         mp2*pow(mk2,2)*L5r*L6r - 128*mp2*pow(mk2,2)*L4r*L8r - 512*mp2*
         pow(mk2,2)*L4r*L6r + 128*mp2*pow(mk2,2)*L4r*L5r + 256*mp2*pow(
         mk2,2)*pow(L4r,2) - 128./9.*mp2*pow(mk2,2)*pi16*L8r - 256./9.*
         mp2*pow(mk2,2)*pi16*L7r + 32./9.*mp2*pow(mk2,2)*pi16*L6r + 64./
         27.*mp2*pow(mk2,2)*pi16*L5r - 16./9.*mp2*pow(mk2,2)*pi16*L4r
          - 4./27.*mp2*pow(mk2,2)*pi16*L3r - 16./9.*mp2*pow(mk2,2)*pi16
         *L2r - 128*pow(mp2,2)*mk2*L4r*L6r + 64*pow(mp2,2)*mk2*pow(
         L4r,2) + 16./3.*pow(mp2,2)*mk2*pi16*L8r + 128./9.*pow(mp2,2)*
         mk2*pi16*L7r;
      mk6L +=  - 16./9.*pow(mp2,2)*mk2*pi16*L6r - 8./27.*pow(mp2,2)*mk2
         *pi16*L5r + 8./9.*pow(mp2,2)*mk2*pi16*L4r + 41./27.*pow(mp2,2)
         *mk2*pi16*L3r + 56./9.*pow(mp2,2)*mk2*pi16*L2r;

      return mk6L/pow(mass.getf0(),4);
}

double mk6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk6Clo\n";}

  double mk6C =
       - 32*CC[12]*pow(mk2,3) - 64*CC[13]*pow(mk2,3) - 32*CC[13]*mp2*
         pow(mk2,2) - 32*CC[14]*pow(mk2,3) + 32*CC[14]*mp2*pow(mk2,2)
          - 16*CC[14]*pow(mp2,2)*mk2 - 32*CC[15]*pow(mk2,3) - 16*CC[15]
         *mp2*pow(mk2,2) - 64*CC[16]*pow(mk2,3) + 64*CC[16]*mp2*pow(
         mk2,2) - 48*CC[16]*pow(mp2,2)*mk2 - 32*CC[17]*mp2*pow(mk2,2)
          + 16*CC[17]*pow(mp2,2)*mk2 + 96*CC[19]*pow(mk2,3) - 96*CC[19]
         *mp2*pow(mk2,2) + 48*CC[19]*pow(mp2,2)*mk2 + 128*CC[20]*pow(
         mk2,3) - 32*CC[20]*mp2*pow(mk2,2) + 48*CC[20]*pow(mp2,2)*mk2
          + 192*CC[21]*pow(mk2,3) + 192*CC[21]*mp2*pow(mk2,2) + 48*CC[
         21]*pow(mp2,2)*mk2 + 32*CC[31]*pow(mk2,3) + 64*CC[32]*pow(
         mk2,3) + 32*CC[32]*mp2*pow(mk2,2);

  return mk6C/pow(mass.getf0(),4);
}

double mk6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double mk6R =
       + pow(Ab(mp2,mu2),2) * (  - 1./2.*pow(mp2,-1)*pow(mk2,2) - 27./
         32.*mk2 );

      mk6R +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - 3./4.*mk2 );

      mk6R +=  + Ab(mp2,mu2)*Ab(me2,mu2) * (  - 3./16.*mk2 - 1./8.*mp2
          - 1./24.*pow(mp2,2)*pow(me2,-1) );

      mk6R +=  + Ab(mp2,mu2) * ( 3./4.*pow(mk2,2)*pi16 + 1./6.*mp2*mk2*
         pi16 );

      mk6R +=  + pow(Ab(mk2,mu2),2) * (  - 251./72.*mk2 - 3./8.*mp2 );

      mk6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * ( 1./3.*mk2 + 1./12.*mp2 + 1.
         /36.*pow(mp2,2)*pow(me2,-1) );

      mk6R +=  + Ab(mk2,mu2) * ( 11./36.*pow(mk2,2)*pi16 + 3./4.*mp2*
         mk2*pi16 );

      mk6R +=  + pow(Ab(me2,mu2),2) * (  - 7./12.*mk2 - 9./128.*mp2 - 3.
         /128.*pow(mp2,2)*pow(me2,-1) );

      mk6R +=  + Ab(me2,mu2) * ( 43./54.*pow(mk2,2)*pi16 + 13./108.*mp2
         *mk2*pi16 );

      mk6R +=  + hh(mp2,mp2,mk2,mk2,mu2) * ( 3./8.*pow(mk2,2) + 3./4.*
         mp2*mk2 );

      mk6R +=  + hh(mp2,mk2,me2,mk2,mu2) * ( 1./4.*pow(mk2,2) );

      mk6R +=  + hh(mk2,mp2,mp2,mk2,mu2) * (  - 3./32.*pow(mk2,2) );

      mk6R +=  + hh(mk2,mp2,me2,mk2,mu2) * ( 9./16.*pow(mk2,2) );

      mk6R +=  + hh(mk2,mk2,mk2,mk2,mu2) * ( 3./4.*pow(mk2,2) );

      mk6R +=  + hh(mk2,me2,me2,mk2,mu2) * ( 181./288.*pow(mk2,2) );

      mk6R +=  + hh1(mp2,mp2,mk2,mk2,mu2) * (  - 3./2.*pow(mk2,2) );

      mk6R +=  + hh1(mk2,mp2,me2,mk2,mu2) * (  - 3./2.*pow(mk2,2) );

      mk6R +=  + hh1(mk2,me2,me2,mk2,mu2) * (  - 5./4.*pow(mk2,2) );

      mk6R +=  + hh21(mp2,mp2,mk2,mk2,mu2) * ( 9./4.*pow(mk2,2) );

      mk6R +=  + hh21(mk2,mp2,mp2,mk2,mu2) * (  - 9./32.*pow(mk2,2) );

      mk6R +=  + hh21(mk2,mp2,me2,mk2,mu2) * ( 27./16.*pow(mk2,2) );

      mk6R +=  + hh21(mk2,mk2,mk2,mk2,mu2) * ( 9./4.*pow(mk2,2) );

      mk6R +=  + hh21(mk2,me2,me2,mk2,mu2) * ( 27./32.*pow(mk2,2) );

      mk6R +=  - 4709./1728.*pow(mk2,3)*pow(pi16,2) - 763./1296.*pow(
         mk2,3)*pow(pi16,2)*pow(pi,2) - 19./108.*mp2*pow(mk2,2)*pow(
         pi16,2) - 73./648.*mp2*pow(mk2,2)*pow(pi16,2)*pow(pi,2) - 13./
         24.*pow(mp2,2)*mk2*pow(pi16,2) - 1./8.*pow(mp2,2)*mk2*pow(
         pi16,2)*pow(pi,2);

      return mk6R/pow(mass.getf0(),4);
}


// *************************************************************************

double meta4lo(const lomass mass, const Li Liin){
  return meta4Llo(mass,Liin)+meta4Rlo(mass);
}

double meta4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  //double me2 = 4./3.*mk2-1./3.*mp2;

  double  meta4L =
       + 128./3.*pow(mk2,2)*L8r + 128./3.*pow(mk2,2)*L7r + 128./3.*pow(
         mk2,2)*L6r - 128./9.*pow(mk2,2)*L5r - 64./3.*pow(mk2,2)*L4r - 
         128./3.*mp2*mk2*L8r - 256./3.*mp2*mk2*L7r + 32./3.*mp2*mk2*L6r
          + 64./9.*mp2*mk2*L5r - 16./3.*mp2*mk2*L4r + 16*pow(mp2,2)*L8r
          + 128./3.*pow(mp2,2)*L7r - 16./3.*pow(mp2,2)*L6r - 8./9.*pow(
         mp2,2)*L5r + 8./3.*pow(mp2,2)*L4r;

  return meta4L/pow(mass.getf0(),2);
}

double meta4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double meta4R =
       + Ab(mp2,mu2) * ( 1./2.*mp2 );

      meta4R +=  + Ab(mk2,mu2) * (  - 4./3.*mk2 );

      meta4R +=  + Ab(me2,mu2) * ( 8./9.*mk2 - 7./18.*mp2 );

  return meta4R/pow(mass.getf0(),2);
}

double meta6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return meta6Rlo(mass)+meta6Llo(mass,Liin)+meta6Clo(mass,Ciin);
}

double meta6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return meta6Llo(mass,Liin)+meta6Clo(mass,Ciin);
}

double meta6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta6Llo\n";}

  double meta6L =
       + Ab(mp2,mu2) * (  - 64*mp2*mk2*L7r + 96*mp2*mk2*L6r - 32./3.*
         mp2*mk2*L5r - 88*mp2*mk2*L4r + 16*mp2*mk2*L3r + 16*mp2*mk2*L2r
          + 64*mp2*mk2*L1r + 40*pow(mp2,2)*L8r + 64*pow(mp2,2)*L7r - 28.
         /3.*pow(mp2,2)*L5r + 4*pow(mp2,2)*L4r - 4*pow(mp2,2)*L3r - 4*
         pow(mp2,2)*L2r - 16*pow(mp2,2)*L1r );

      meta6L +=  + Ab(mk2,mu2) * ( 64./3.*pow(mk2,2)*L8r + 128./3.*pow(
         mk2,2)*L7r - 32./9.*pow(mk2,2)*L5r - 64./3.*pow(mk2,2)*L4r + 
         112./3.*pow(mk2,2)*L3r + 64./3.*pow(mk2,2)*L2r + 256./3.*pow(
         mk2,2)*L1r - 32*mp2*mk2*L8r - 128./3.*mp2*mk2*L7r - 64*mp2*mk2
         *L6r + 80./9.*mp2*mk2*L5r + 160./3.*mp2*mk2*L4r - 28./3.*mp2*
         mk2*L3r - 16./3.*mp2*mk2*L2r - 64./3.*mp2*mk2*L1r );

      meta6L +=  + Ab(me2,mu2) * ( 1280./9.*pow(mk2,2)*L8r + 1280./9.*
         pow(mk2,2)*L7r + 1280./9.*pow(mk2,2)*L6r - 1280./27.*pow(
         mk2,2)*L5r - 640./9.*pow(mk2,2)*L4r + 64./3.*pow(mk2,2)*L3r + 
         128./3.*pow(mk2,2)*L2r + 128./3.*pow(mk2,2)*L1r - 1328./9.*mp2
         *mk2*L8r - 2224./9.*mp2*mk2*L7r - 352./9.*mp2*mk2*L6r + 880./
         27.*mp2*mk2*L5r + 104./9.*mp2*mk2*L4r - 32./3.*mp2*mk2*L3r - 
         64./3.*mp2*mk2*L2r - 64./3.*mp2*mk2*L1r + 404./9.*pow(mp2,2)*
         L8r + 964./9.*pow(mp2,2)*L7r - 80./9.*pow(mp2,2)*L6r - 140./27.
         *pow(mp2,2)*L5r + 68./9.*pow(mp2,2)*L4r + 4./3.*pow(mp2,2)*L3r
          + 8./3.*pow(mp2,2)*L2r + 8./3.*pow(mp2,2)*L1r + 64./27.*pow(
         mp2,2)*mk2*pow(me2,-1)*L8r + 64./27.*pow(mp2,2)*mk2*pow(
         me2,-1)*L7r + 64./27.*pow(mp2,2)*mk2*pow(me2,-1)*L6r - 52./27.
         *pow(mp2,3)*pow(me2,-1)*L8r - 124./27.*pow(mp2,3)*pow(me2,-1)*
         L7r - 16./27.*pow(mp2,3)*pow(me2,-1)*L6r );

      meta6L +=  - 4096./9.*pow(mk2,3)*L5r*L8r - 4096./9.*pow(mk2,3)*
         L5r*L7r - 4096./9.*pow(mk2,3)*L5r*L6r + 4096./27.*pow(mk2,3)*
         pow(L5r,2) - 2048./3.*pow(mk2,3)*L4r*L8r - 2048./3.*pow(mk2,3)
         *L4r*L7r - 2048./3.*pow(mk2,3)*L4r*L6r + 4096./9.*pow(mk2,3)*
         L4r*L5r + 1024./3.*pow(mk2,3)*pow(L4r,2) - 448./27.*pow(mk2,3)
         *pi16*L8r - 1024./27.*pow(mk2,3)*pi16*L7r + 128./27.*pow(
         mk2,3)*pi16*L6r + 160./81.*pow(mk2,3)*pi16*L5r - 64./27.*pow(
         mk2,3)*pi16*L4r + 152./27.*pow(mk2,3)*pi16*L3r + 544./27.*pow(
         mk2,3)*pi16*L2r + 256./27.*pow(mk2,3)*pi16*L1r + 5120./9.*mp2*
         pow(mk2,2)*L5r*L8r + 1024*mp2*pow(mk2,2)*L5r*L7r - 1024./9.*
         mp2*pow(mk2,2)*pow(L5r,2) + 1024./3.*mp2*pow(mk2,2)*L4r*L8r + 
         1024*mp2*pow(mk2,2)*L4r*L7r - 512*mp2*pow(mk2,2)*L4r*L6r + 256
         *mp2*pow(mk2,2)*pow(L4r,2) + 1472./27.*mp2*pow(mk2,2)*pi16*L8r
          + 832./9.*mp2*pow(mk2,2)*pi16*L7r + 256./9.*mp2*pow(mk2,2)*
         pi16*L6r - 320./27.*mp2*pow(mk2,2)*pi16*L5r - 128./9.*mp2*pow(
         mk2,2)*pi16*L4r;
      meta6L +=  - 34./9.*mp2*pow(mk2,2)*pi16*L3r - 88./9.*mp2*pow(
         mk2,2)*pi16*L2r - 64./9.*mp2*pow(mk2,2)*pi16*L1r - 2560./9.*
         pow(mp2,2)*mk2*L5r*L8r - 2048./3.*pow(mp2,2)*mk2*L5r*L7r + 256.
         /3.*pow(mp2,2)*mk2*L5r*L6r + 256./9.*pow(mp2,2)*mk2*pow(L5r,2)
          + 256./3.*pow(mp2,2)*mk2*L4r*L8r - 256./3.*pow(mp2,2)*mk2*L4r
         *L5r - 832./27.*pow(mp2,2)*mk2*pi16*L8r - 640./9.*pow(mp2,2)*
         mk2*pi16*L7r - 64./9.*pow(mp2,2)*mk2*pi16*L6r + 32./9.*pow(
         mp2,2)*mk2*pi16*L5r + 32./9.*pow(mp2,2)*mk2*pi16*L4r + 32./9.*
         pow(mp2,2)*mk2*pi16*L3r + 88./9.*pow(mp2,2)*mk2*pi16*L2r + 16./
         9.*pow(mp2,2)*mk2*pi16*L1r + 128./3.*pow(mp2,3)*L5r*L8r + 1024.
         /9.*pow(mp2,3)*L5r*L7r - 128./9.*pow(mp2,3)*L5r*L6r - 64./27.*
         pow(mp2,3)*pow(L5r,2) - 128*pow(mp2,3)*L4r*L8r - 1024./3.*pow(
         mp2,3)*L4r*L7r + 128./3.*pow(mp2,3)*L4r*L6r + 128./9.*pow(
         mp2,3)*L4r*L5r - 64./3.*pow(mp2,3)*pow(L4r,2) - 16./9.*pow(
         mp2,3)*pi16*L8r + 448./27.*pow(mp2,3)*pi16*L7r - 272./27.*pow(
         mp2,3)*pi16*L6r;
      meta6L +=  + 296./81.*pow(mp2,3)*pi16*L5r + 136./27.*pow(mp2,3)*
         pi16*L4r - 20./27.*pow(mp2,3)*pi16*L3r - 58./27.*pow(mp2,3)*
         pi16*L2r - 4./27.*pow(mp2,3)*pi16*L1r;

      return meta6L/pow(mass.getf0(),4);
}

double meta6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta6Clo\n";}

  double meta6C =
       - 2048./27.*CC[12]*pow(mk2,3) + 512./9.*CC[12]*mp2*pow(mk2,2) - 
         128./9.*CC[12]*pow(mp2,2)*mk2 + 32./27.*CC[12]*pow(mp2,3) -  
         1024./9.*CC[13]*pow(mk2,3) + 64./3.*CC[13]*pow(mp2,2)*mk2 - 32.
         /9.*CC[13]*pow(mp2,3) - 512./9.*CC[14]*pow(mk2,3) + 640./9.*
         CC[14]*mp2*pow(mk2,2) - 320./9.*CC[14]*pow(mp2,2)*mk2 + 16./3.
         *CC[14]*pow(mp2,3) - 512./9.*CC[15]*pow(mk2,3) + 32./3.*CC[15]
         *pow(mp2,2)*mk2 - 16./9.*CC[15]*pow(mp2,3) - 256./3.*CC[16]*
         pow(mk2,3) + 320./3.*CC[16]*mp2*pow(mk2,2) - 256./3.*CC[16]*
         pow(mp2,2)*mk2 + 16*CC[16]*pow(mp2,3) - 512./9.*CC[17]*pow(
         mk2,3) + 640./9.*CC[17]*mp2*pow(mk2,2) - 320./9.*CC[17]*pow(
         mp2,2)*mk2 + 16./3.*CC[17]*pow(mp2,3) - 512./9.*CC[18]*pow(
         mk2,3) + 128*CC[18]*mp2*pow(mk2,2) - 256./3.*CC[18]*pow(mp2,2)
         *mk2 + 128./9.*CC[18]*pow(mp2,3) + 256*CC[19]*pow(mk2,3) - 384
         *CC[19]*mp2*pow(mk2,2) + 192*CC[19]*pow(mp2,2)*mk2 - 16*CC[19]
         *pow(mp2,3);
      meta6C +=  + 256*CC[20]*pow(mk2,3) - 192*CC[20]*mp2*pow(mk2,2) + 
	 64*CC[20]*pow(mp2,2)*mk2 + 16*CC[20]*pow(mp2,3) + 256*CC[21]*
         pow(mk2,3) + 192*CC[21]*mp2*pow(mk2,2) - 16*CC[21]*pow(mp2,3)
          + 512./3.*CC[31]*pow(mk2,3) - 256*CC[31]*mp2*pow(mk2,2) + 128
         *CC[31]*pow(mp2,2)*mk2 - 32./3.*CC[31]*pow(mp2,3) + 512./3.*
         CC[32]*pow(mk2,3) - 256./3.*CC[32]*mp2*pow(mk2,2) - 64./3.*CC[
         32]*pow(mp2,2)*mk2 + 32*CC[32]*pow(mp2,3) + 512./3.*CC[33]*
         pow(mk2,3) - 1024./3.*CC[33]*mp2*pow(mk2,2) + 512./3.*CC[33]*
         pow(mp2,2)*mk2;

  return meta6C/pow(mass.getf0(),4);
}

double meta6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double meta6R =
       + pow(Ab(mp2,mu2),2) * (  - 3./4.*mk2 - 29./48.*mp2 );

      meta6R +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - 3./2.*mp2 );

      meta6R +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 1./2.*mp2 - 1./12.*pow(
         mp2,2)*pow(me2,-1) );

      meta6R +=  + Ab(mp2,mu2) * (  - 4./9.*mp2*mk2*pi16 + 4./9.*pow(
         mp2,2)*pi16 );

      meta6R +=  + pow(Ab(mk2,mu2),2) * (  - 29./9.*mk2 + 43./12.*mp2
          - pow(mp2,2)*pow(mk2,-1) );

      meta6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 4*mk2 + 5./3.*mp2 + 1.
         /18.*pow(mp2,2)*pow(me2,-1) );

      meta6R +=  + Ab(mk2,mu2) * ( 104./27.*pow(mk2,2)*pi16 + 4./27.*
         mp2*mk2*pi16 - 1./3.*pow(mp2,2)*pi16 );

      meta6R +=  + pow(Ab(me2,mu2),2) * (  - 193./108.*mk2 + 307./432.*
         mp2 - 1./9.*pow(mp2,2)*pow(me2,-1) );

      meta6R +=  + Ab(me2,mu2) * (  - 100./81.*pow(mk2,2)*pi16 + 56./81.
         *mp2*mk2*pi16 - 32./81.*pow(mp2,2)*pi16 + 52./243.*pow(mp2,2)*
         mk2*pow(me2,-1)*pi16 - 13./243.*pow(mp2,3)*pow(me2,-1)*pi16 );

      meta6R +=  + hh(mp2,mp2,me2,me2,mu2) * ( 1./6.*pow(mp2,2) );

      meta6R +=  + hh(mp2,mk2,mk2,me2,mu2) * ( 3./2.*mp2*mk2 + 1./8.*
         pow(mp2,2) );

      meta6R +=  + hh(mk2,mk2,me2,me2,mu2) * (  - 64./9.*pow(mk2,2) + 4
         *mp2*mk2 - 1./2.*pow(mp2,2) );

      meta6R +=  + hh(me2,mk2,mk2,me2,mu2) * ( 2*pow(mk2,2) - mp2*mk2
          + 1./8.*pow(mp2,2) );

      meta6R +=  + hh(me2,me2,me2,me2,mu2) * ( 128./243.*pow(mk2,2) - 
         112./243.*mp2*mk2 + 49./486.*pow(mp2,2) );

      meta6R +=  + hh1(mp2,mk2,mk2,me2,mu2) * (  - 4*mp2*mk2 + pow(
         mp2,2) );

      meta6R +=  + hh1(mk2,mk2,me2,me2,mu2) * ( 64./3.*pow(mk2,2) - 40./
         3.*mp2*mk2 + 2*pow(mp2,2) );

      meta6R +=  + hh21(mp2,mk2,mk2,me2,mu2) * ( 6*pow(mk2,2) - 3*mp2*
         mk2 + 3./8.*pow(mp2,2) );

      meta6R +=  + hh21(me2,mk2,mk2,me2,mu2) * ( 6*pow(mk2,2) - 3*mp2*
         mk2 + 3./8.*pow(mp2,2) );

      meta6R +=  - 7567./972.*pow(mk2,3)*pow(pi16,2) - 1091./729.*pow(
         mk2,3)*pow(pi16,2)*pow(pi,2) + 367./48.*mp2*pow(mk2,2)*pow(
         pi16,2) + 133./108.*mp2*pow(mk2,2)*pow(pi16,2)*pow(pi,2) - 
         4133./1296.*pow(mp2,2)*mk2*pow(pi16,2) - 269./486.*pow(mp2,2)*
         mk2*pow(pi16,2)*pow(pi,2) - 1781./15552.*pow(mp2,3)*pow(
         pi16,2) - 91./11664.*pow(mp2,3)*pow(pi16,2)*pow(pi,2);

      return meta6R/pow(mass.getf0(),4);
}

// *************************************************************************

double fpi4lo(const lomass mass, const Li Liin){
  return fpi4Llo(mass,Liin)+fpi4Rlo(mass);
}

double fpi4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);


  double fpi4L =
       + 8*mk2*L4r + 4*mp2*L5r + 4*mp2*L4r;

  return fpi4L/pow(mass.getf0(),2);
}

double fpi4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double mu2 = pow(mass.getmu(),2);

  double  fpi4R =
       + Ab(mp2,mu2) * ( 1 );

      fpi4R +=  + Ab(mk2,mu2) * ( 1./2. );

  return fpi4R/pow(mass.getf0(),2);
}

double fpi6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return fpi6Rlo(mass)+fpi6Llo(mass,Liin)+fpi6Clo(mass,Ciin);
}

double fpi6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return fpi6Llo(mass,Liin)+fpi6Clo(mass,Ciin);
}

double fpi6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6Llo\n";}

  double fpi6L =
       + Ab(mp2,mu2) * ( 32*mk2*L6r - 24*mk2*L4r + 16*mp2*L8r + 16*mp2*
         L6r - 6*mp2*L5r - 14*mp2*L3r - 16*mp2*L2r - 28*mp2*L1r );

      fpi6L +=  + Ab(mk2,mu2) * ( 8*mk2*L8r + 16*mk2*L6r - 4*mk2*L5r + 
         4*mk2*L4r - 10*mk2*L3r - 8*mk2*L2r - 32*mk2*L1r + 8*mp2*L6r + 
         2*mp2*L5r - 6*mp2*L4r );

      fpi6L +=  + Ab(me2,mu2) * ( 16./3.*mk2*L4r - 8./3.*mk2*L3r - 8./3.
         *mk2*L2r - 32./3.*mk2*L1r + 2./3.*mp2*L5r - 4./3.*mp2*L4r + 2./
         3.*mp2*L3r + 2./3.*mp2*L2r + 8./3.*mp2*L1r );

      fpi6L +=  - 32*pow(mk2,2)*pow(L4r,2) - 8*pow(mk2,2)*pi16*L8r - 16
         *pow(mk2,2)*pi16*L6r + 4*pow(mk2,2)*pi16*L5r + 8*pow(mk2,2)*
         pi16*L4r - 43./27.*pow(mk2,2)*pi16*L3r - 52./9.*pow(mk2,2)*
         pi16*L2r - 32*mp2*mk2*L4r*L5r - 32*mp2*mk2*pow(L4r,2) - 40*mp2
         *mk2*pi16*L6r + 20*mp2*mk2*pi16*L4r + 8./27.*mp2*mk2*pi16*L3r
          + 8./9.*mp2*mk2*pi16*L2r - 8*pow(mp2,2)*pow(L5r,2) - 16*pow(
         mp2,2)*L4r*L5r - 8*pow(mp2,2)*pow(L4r,2) - 16*pow(mp2,2)*pi16*
         L8r - 16*pow(mp2,2)*pi16*L6r + 8*pow(mp2,2)*pi16*L5r + 8*pow(
         mp2,2)*pi16*L4r - 28./27.*pow(mp2,2)*pi16*L3r - 37./9.*pow(
         mp2,2)*pi16*L2r - 2*pow(mp2,2)*pi16*L1r;

  return fpi6L/pow(mass.getf0(),4);
}

double fpi6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6Clo\n";}

  double fpi6C =
       + 8*CC[14]*pow(mp2,2) + 16*CC[15]*mp2*mk2 + 8*CC[15]*pow(mp2,2)
          + 32*CC[16]*pow(mk2,2) - 32*CC[16]*mp2*mk2 + 24*CC[16]*pow(
         mp2,2) + 8*CC[17]*pow(mp2,2);

  return fpi6C/pow(mass.getf0(),4);
}

double fpi6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double fpi6R =
       + pow(Ab(mp2,mu2),2) * ( 7./32. + 1./4.*pow(mp2,-1)*mk2 );

      fpi6R +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 1./6. );

      fpi6R +=  + Ab(mp2,mu2) * (  - 1./2.*mk2*pi16 - 1./4.*mp2*pi16 );

      fpi6R +=  + pow(Ab(mk2,mu2),2) * ( 3./4. - 1./8.*mp2*pow(mk2,-1)
          );

      fpi6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 1./6. );

      fpi6R +=  + Ab(mk2,mu2) * (  - 1./2.*mk2*pi16 - 1./8.*mp2*pi16 );

      fpi6R +=  + pow(Ab(me2,mu2),2) * ( 3./32. );

      fpi6R +=  + Ab(me2,mu2) * ( 1./6.*mk2*pi16 - 1./6.*mp2*pi16 );

      fpi6R +=  + hh(mp2,mp2,mp2,mp2,mu2) * (  - 1./2.*mp2 );

      fpi6R +=  + hh(mp2,mk2,mk2,mp2,mu2) * ( 1./16.*mp2 );

      fpi6R +=  + hh(mk2,mp2,mk2,mp2,mu2) * (  - 1./2.*mk2 );

      fpi6R +=  + hh(me2,mk2,mk2,mp2,mu2) * (  - 1./4.*mk2 + 1./16.*mp2
          );

      fpi6R +=  + hhd(mp2,mp2,mp2,mp2,mu2) * ( 5./12.*pow(mp2,2) );

      fpi6R +=  + hhd(mp2,mk2,mk2,mp2,mu2) * (  - 5./16.*pow(mp2,2) );

      fpi6R +=  + hhd(mp2,me2,me2,mp2,mu2) * ( 1./36.*pow(mp2,2) );

      fpi6R +=  + hhd(mk2,mp2,mk2,mp2,mu2) * ( 1./2.*mp2*mk2 );

      fpi6R +=  + hhd(mk2,mk2,me2,mp2,mu2) * (  - 5./12.*pow(mp2,2) );

      fpi6R +=  + hhd(me2,mk2,mk2,mp2,mu2) * ( 1./4.*mp2*mk2 - 1./16.*
         pow(mp2,2) );

      fpi6R +=  + hh1d(mp2,mk2,mk2,mp2,mu2) * ( 1./2.*pow(mp2,2) );

      fpi6R +=  + hh1d(mk2,mk2,me2,mp2,mu2) * ( pow(mp2,2) );

      fpi6R +=  + hh21d(mp2,mp2,mp2,mp2,mu2) * ( 3./2.*pow(mp2,2) );

      fpi6R +=  + hh21d(mp2,mk2,mk2,mp2,mu2) * (  - 3./16.*pow(mp2,2) )
         ;

      fpi6R +=  + hh21d(mk2,mp2,mk2,mp2,mu2) * ( 3./2.*pow(mp2,2) );

      fpi6R +=  + hh21d(me2,mk2,mk2,mp2,mu2) * ( 9./16.*pow(mp2,2) );

      fpi6R +=  + 15./32.*pow(mk2,2)*pow(pi16,2) + 11./72.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 5./32.*mp2*mk2*pow(pi16,2) + 1./144.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 41./128.*pow(mp2,2)*pow(
         pi16,2) + 35./288.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return fpi6R/pow(mass.getf0(),4);
}

// *************************************************************************

double fk4lo(const lomass mass, const Li Liin){
  return fk4Llo(mass,Liin)+fk4Rlo(mass);
}

double fk4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);

  double  fk4L =
       + 4*mk2*L5r + 8*mk2*L4r + 4*mp2*L4r;

  return fk4L/pow(mass.getf0(),2);
}

double fk4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double fk4R =
       + Ab(mp2,mu2) * ( 3./8. );

      fk4R +=  + Ab(mk2,mu2) * ( 3./4. );

      fk4R +=  + Ab(me2,mu2) * ( 3./8. );

  return fk4R/pow(mass.getf0(),2);
}

double fk6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return fk6Rlo(mass)+fk6Llo(mass,Liin)+fk6Clo(mass,Ciin);
}

double fk6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return fk6Llo(mass,Liin)+fk6Clo(mass,Ciin);
}

double fk6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk6Llo\n";}

  double fk6L =
       + Ab(mp2,mu2) * ( 12*mk2*L6r + 3./2.*mk2*L5r - 9*mk2*L4r + 6*mp2
         *L8r + 6*mp2*L6r - 3*mp2*L5r + 15./2.*mp2*L4r - 15./2.*mp2*L3r
          - 6*mp2*L2r - 24*mp2*L1r );

      fk6L +=  + Ab(mk2,mu2) * ( 12*mk2*L8r + 24*mk2*L6r - 3*mk2*L5r - 
         2*mk2*L4r - 15*mk2*L3r - 18*mk2*L2r - 36*mk2*L1r + 12*mp2*L6r
          - 9*mp2*L4r );

      fk6L +=  + Ab(me2,mu2) * ( 12*mk2*L8r + 12*mk2*L7r + 12*mk2*L6r
          - 23./6.*mk2*L5r - 11./3.*mk2*L4r - 14./3.*mk2*L3r - 8./3.*
         mk2*L2r - 32./3.*mk2*L1r - 9*mp2*L8r - 21*mp2*L7r + 6*mp2*L6r
          + mp2*L5r - 35./6.*mp2*L4r + 7./6.*mp2*L3r + 2./3.*mp2*L2r + 
         8./3.*mp2*L1r + 3*pow(mp2,2)*pow(me2,-1)*L8r + 9*pow(mp2,2)*
         pow(me2,-1)*L7r );

      fk6L +=  - 8*pow(mk2,2)*pow(L5r,2) - 32*pow(mk2,2)*L4r*L5r - 32*
         pow(mk2,2)*pow(L4r,2) - 28*pow(mk2,2)*pi16*L8r - 16*pow(mk2,2)
         *pi16*L7r - 40*pow(mk2,2)*pi16*L6r + 34./3.*pow(mk2,2)*pi16*
         L5r + 20*pow(mk2,2)*pi16*L4r - 89./54.*pow(mk2,2)*pi16*L3r - 
         61./9.*pow(mk2,2)*pi16*L2r - 2*pow(mk2,2)*pi16*L1r - 16*mp2*
         mk2*L4r*L5r - 32*mp2*mk2*pow(L4r,2) + 16*mp2*mk2*pi16*L8r + 32
         *mp2*mk2*pi16*L7r - 28*mp2*mk2*pi16*L6r - 8./3.*mp2*mk2*pi16*
         L5r + 14*mp2*mk2*pi16*L4r + 2./27.*mp2*mk2*pi16*L3r + 8./9.*
         mp2*mk2*pi16*L2r - 8*pow(mp2,2)*pow(L4r,2) - 12*pow(mp2,2)*
         pi16*L8r - 16*pow(mp2,2)*pi16*L7r - 4*pow(mp2,2)*pi16*L6r + 10.
         /3.*pow(mp2,2)*pi16*L5r + 2*pow(mp2,2)*pi16*L4r - 41./54.*pow(
         mp2,2)*pi16*L3r - 28./9.*pow(mp2,2)*pi16*L2r;

      return fk6L/pow(mass.getf0(),4);
}

double fk6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk6Clo\n";}

  double fk6C =
       + 16*CC[14]*pow(mk2,2) - 16*CC[14]*mp2*mk2 + 8*CC[14]*pow(mp2,2)
          + 16*CC[15]*pow(mk2,2) + 8*CC[15]*mp2*mk2 + 32*CC[16]*pow(
         mk2,2) - 32*CC[16]*mp2*mk2 + 24*CC[16]*pow(mp2,2) + 16*CC[17]*
         mp2*mk2 - 8*CC[17]*pow(mp2,2);

  return fk6C/pow(mass.getf0(),4);
}

double fk6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double fk6R =
       + pow(Ab(mp2,mu2),2) * ( 9./128. + 3./32.*pow(mp2,-1)*mk2 );

      fk6R +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 3./32. );

      fk6R +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 13./64. + 3./16.*mp2*pow(
         me2,-1) );

      fk6R +=  + Ab(mp2,mu2) * (  - 15./32.*mk2*pi16 );

      fk6R +=  + pow(Ab(mk2,mu2),2) * ( 27./32. + 3./16.*mp2*pow(
         mk2,-1) );

      fk6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 29./32. - 1./8.*mp2*
         pow(me2,-1) );

      fk6R +=  + Ab(mk2,mu2) * (  - 1./16.*mk2*pi16 - 3./8.*mp2*pi16 );

      fk6R +=  + pow(Ab(me2,mu2),2) * ( 17./32. + 1./128.*mp2*pow(
         me2,-1) );

      fk6R +=  + Ab(me2,mu2) * (  - 41./96.*mk2*pi16 - 1./24.*mp2*pi16
          );

      fk6R +=  + hh(mp2,mp2,mk2,mk2,mu2) * (  - 3./8.*mp2 );

      fk6R +=  + hh(mk2,mp2,mp2,mk2,mu2) * ( 3./64.*mk2 );

      fk6R +=  + hh(mk2,mp2,me2,mk2,mu2) * (  - 9./32.*mk2 );

      fk6R +=  + hh(mk2,mk2,mk2,mk2,mu2) * (  - 3./8.*mk2 );

      fk6R +=  + hh(mk2,me2,me2,mk2,mu2) * (  - 9./64.*mk2 );

      fk6R +=  + hhd(mp2,mp2,mk2,mk2,mu2) * ( 3./16.*pow(mk2,2) + 3./8.
         *mp2*mk2 );

      fk6R +=  + hhd(mp2,mk2,me2,mk2,mu2) * ( 1./8.*pow(mk2,2) );

      fk6R +=  + hhd(mk2,mp2,mp2,mk2,mu2) * (  - 3./64.*pow(mk2,2) );

      fk6R +=  + hhd(mk2,mp2,me2,mk2,mu2) * ( 9./32.*pow(mk2,2) );

      fk6R +=  + hhd(mk2,mk2,mk2,mk2,mu2) * ( 3./8.*pow(mk2,2) );

      fk6R +=  + hhd(mk2,me2,me2,mk2,mu2) * ( 181./576.*pow(mk2,2) );

      fk6R +=  + hh1d(mp2,mp2,mk2,mk2,mu2) * (  - 3./4.*pow(mk2,2) );

      fk6R +=  + hh1d(mk2,mp2,me2,mk2,mu2) * (  - 3./4.*pow(mk2,2) );

      fk6R +=  + hh1d(mk2,me2,me2,mk2,mu2) * (  - 5./8.*pow(mk2,2) );

      fk6R +=  + hh21d(mp2,mp2,mk2,mk2,mu2) * ( 9./8.*pow(mk2,2) );

      fk6R +=  + hh21d(mk2,mp2,mp2,mk2,mu2) * (  - 9./64.*pow(mk2,2) );

      fk6R +=  + hh21d(mk2,mp2,me2,mk2,mu2) * ( 27./32.*pow(mk2,2) );

      fk6R +=  + hh21d(mk2,mk2,mk2,mk2,mu2) * ( 9./8.*pow(mk2,2) );

      fk6R +=  + hh21d(mk2,me2,me2,mk2,mu2) * ( 27./64.*pow(mk2,2) );

      fk6R +=  + 197./384.*pow(mk2,2)*pow(pi16,2) + 3./16.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 29./192.*mp2*mk2*pow(pi16,2) + 1./32.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 13./48.*pow(mp2,2)*pow(pi16,2)
          + 1./16.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return fk6R/pow(mass.getf0(),4);
}

// *************************************************************************

double feta4lo(const lomass mass, const Li Liin){
  return feta4Llo(mass,Liin)+feta4Rlo(mass);
}

double feta4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta4Llo\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);

  double feta4L =
       + 16./3.*mk2*L5r + 8*mk2*L4r - 4./3.*mp2*L5r + 4*mp2*L4r;

  return feta4L/pow(mass.getf0(),2);
}

double feta4Rlo(const lomass mass){
  double mk2 = pow(mass.getmk0(),2);
  double mu2 = pow(mass.getmu(),2);

  double feta4R =
       + Ab(mk2,mu2) * ( 3./2. );

  return feta4R/pow(mass.getf0(),2);
}

double feta6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return feta6Rlo(mass)+feta6Llo(mass,Liin)+feta6Clo(mass,Ciin);
}

double feta6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return feta6Llo(mass,Liin)+feta6Clo(mass,Ciin);
}

double feta6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta6Llo\n";}

  double feta6L =
       + Ab(mp2,mu2) * ( 2*mp2*L5r + 12*mp2*L4r - 6*mp2*L3r - 6*mp2*L2r
          - 24*mp2*L1r );

      feta6L +=  + Ab(mk2,mu2) * ( 24*mk2*L8r + 48*mk2*L6r - 28./3.*mk2
         *L5r - 20*mk2*L4r - 14*mk2*L3r - 8*mk2*L2r - 32*mk2*L1r + 24*
         mp2*L6r - 2*mp2*L5r - 18*mp2*L4r );

      feta6L +=  + Ab(me2,mu2) * ( 32./9.*mk2*L5r + 16./3.*mk2*L4r - 8*
         mk2*L3r - 16*mk2*L2r - 16*mk2*L1r - 14./9.*mp2*L5r - 4./3.*mp2
         *L4r + 2*mp2*L3r + 4*mp2*L2r + 4*mp2*L1r );

      feta6L +=  - 128./9.*pow(mk2,2)*pow(L5r,2) - 128./3.*pow(mk2,2)*
         L4r*L5r - 32*pow(mk2,2)*pow(L4r,2) - 24*pow(mk2,2)*pi16*L8r - 
         48*pow(mk2,2)*pi16*L6r + 12*pow(mk2,2)*pi16*L5r + 24*pow(
         mk2,2)*pi16*L4r - 19./9.*pow(mk2,2)*pi16*L3r - 68./9.*pow(
         mk2,2)*pi16*L2r - 32./9.*pow(mk2,2)*pi16*L1r + 64./9.*mp2*mk2*
         pow(L5r,2) - 32./3.*mp2*mk2*L4r*L5r - 32*mp2*mk2*pow(L4r,2) - 
         24*mp2*mk2*pi16*L6r + 12*mp2*mk2*pi16*L4r + 8./9.*mp2*mk2*pi16
         *L3r + 16./9.*mp2*mk2*pi16*L2r + 16./9.*mp2*mk2*pi16*L1r - 8./
         9.*pow(mp2,2)*pow(L5r,2) + 16./3.*pow(mp2,2)*L4r*L5r - 8*pow(
         mp2,2)*pow(L4r,2) - 10./9.*pow(mp2,2)*pi16*L3r - 29./9.*pow(
         mp2,2)*pi16*L2r - 2./9.*pow(mp2,2)*pi16*L1r;

      return feta6L/pow(mass.getf0(),4);
}

double feta6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta6Clo\n";}

  double feta6C =
       + 64./3.*CC[14]*pow(mk2,2) - 64./3.*CC[14]*mp2*mk2 + 8*CC[14]*
         pow(mp2,2) + 64./3.*CC[15]*pow(mk2,2) + 16./3.*CC[15]*mp2*mk2
          - 8./3.*CC[15]*pow(mp2,2) + 32*CC[16]*pow(mk2,2) - 32*CC[16]*
         mp2*mk2 + 24*CC[16]*pow(mp2,2) + 64./3.*CC[17]*pow(mk2,2) - 64.
         /3.*CC[17]*mp2*mk2 + 8*CC[17]*pow(mp2,2) + 64./3.*CC[18]*pow(
         mk2,2) - 128./3.*CC[18]*mp2*mk2 + 64./3.*CC[18]*pow(mp2,2);

  return feta6C/pow(mass.getf0(),4);
}

double feta6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double feta6R =
       + pow(Ab(mp2,mu2),2) * ( 9./32. );

      feta6R +=  + pow(Ab(mk2,mu2),2) * ( 3./4. + 3./8.*mp2*pow(mk2,-1)
          );

      feta6R +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 1./2. );

      feta6R +=  + Ab(mk2,mu2) * (  - 3./2.*mk2*pi16 - 3./8.*mp2*pi16 )
         ;

      feta6R +=  + pow(Ab(me2,mu2),2) * ( 9./32. );

      feta6R +=  + Ab(me2,mu2) * ( 1./2.*mk2*pi16 );

      feta6R +=  + hh(mp2,mk2,mk2,me2,mu2) * (  - 9./16.*mp2 );

      feta6R +=  + hh(me2,mk2,mk2,me2,mu2) * (  - 3./4.*mk2 + 3./16.*
         mp2 );

      feta6R +=  + hhd(mp2,mp2,me2,me2,mu2) * ( 1./12.*pow(mp2,2) );

      feta6R +=  + hhd(mp2,mk2,mk2,me2,mu2) * ( 3./4.*mp2*mk2 + 1./16.*
         pow(mp2,2) );

      feta6R +=  + hhd(mk2,mk2,me2,me2,mu2) * (  - 32./9.*pow(mk2,2) + 
         2*mp2*mk2 - 1./4.*pow(mp2,2) );

      feta6R +=  + hhd(me2,mk2,mk2,me2,mu2) * ( pow(mk2,2) - 1./2.*mp2*
         mk2 + 1./16.*pow(mp2,2) );

      feta6R +=  + hhd(me2,me2,me2,me2,mu2) * ( 64./243.*pow(mk2,2) - 
         56./243.*mp2*mk2 + 49./972.*pow(mp2,2) );

      feta6R +=  + hh1d(mp2,mk2,mk2,me2,mu2) * (  - 2*mp2*mk2 + 1./2.*
         pow(mp2,2) );

      feta6R +=  + hh1d(mk2,mk2,me2,me2,mu2) * ( 32./3.*pow(mk2,2) - 20.
         /3.*mp2*mk2 + pow(mp2,2) );

      feta6R +=  + hh21d(mp2,mk2,mk2,me2,mu2) * ( 3*pow(mk2,2) - 3./2.*
         mp2*mk2 + 3./16.*pow(mp2,2) );

      feta6R +=  + hh21d(me2,mk2,mk2,me2,mu2) * ( 3*pow(mk2,2) - 3./2.*
         mp2*mk2 + 3./16.*pow(mp2,2) );

      feta6R +=  + 49./96.*pow(mk2,2)*pow(pi16,2) + 5./24.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 11./96.*mp2*mk2*pow(pi16,2) + 1./48.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 91./384.*pow(mp2,2)*pow(
         pi16,2) + 5./96.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return feta6R/pow(mass.getf0(),4);
}

// *************************************************************************

double qqup4lo(const lomass mass, const Li Liin){
  return qqup4Llo(mass,Liin)+qqup4Rlo(mass);
}

double qqup4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup4L\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);

  double Qnr4L1 =
       + 32*mk2*L6r + 4*mp2*H2r + 8*mp2*L8r + 16*mp2*L6r;

  return Qnr4L1/pow(mass.getf0(),2);
}

double qqup4Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double Qnr4R1 =
       + Ab(mp2,mu2) * ( 3./2. );

      Qnr4R1 +=  + Ab(mk2,mu2) * ( 1 );

      Qnr4R1 +=  + Ab(me2,mu2) * ( 1./6. );

  return Qnr4R1/pow(mass.getf0(),2);
}

double qqup6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return qqup6Rlo(mass)+qqup6Llo(mass,Liin)+qqup6Clo(mass,Ciin);
}

double qqup6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return qqup6Llo(mass,Liin)+qqup6Clo(mass,Ciin);
}

double qqup6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup6L\n";}


  double Qnr6L1 =
       + Ab(mp2,mu2) * ( 96*mk2*L6r - 48*mk2*L4r + 72*mp2*L8r + 96*mp2*
         L6r - 36*mp2*L5r - 48*mp2*L4r );

      Qnr6L1 +=  + Ab(mk2,mu2) * ( 48*mk2*L8r + 128*mk2*L6r - 24*mk2*
         L5r - 64*mk2*L4r + 32*mp2*L6r - 16*mp2*L4r );

      Qnr6L1 +=  + Ab(me2,mu2) * ( 16./3.*mk2*L8r - 16*mk2*L7r + 32*mk2
         *L6r - 16./3.*mk2*L5r - 16*mk2*L4r + 4./3.*mp2*L8r + 12*mp2*
         L7r + 4./3.*mp2*L5r + 4./3.*pow(mp2,2)*pow(me2,-1)*L8r + 4*
         pow(mp2,2)*pow(me2,-1)*L7r );

      Qnr6L1 +=  - 208./9.*pow(mk2,2)*pi16*L8r - 64./9.*pow(mk2,2)*pi16
         *L7r - 352./9.*pow(mk2,2)*pi16*L6r + 280./27.*pow(mk2,2)*pi16*
         L5r + 176./9.*pow(mk2,2)*pi16*L4r + 64./9.*mp2*mk2*pi16*L8r + 
         128./9.*mp2*mk2*pi16*L7r - 592./9.*mp2*mk2*pi16*L6r - 32./27.*
         mp2*mk2*pi16*L5r + 296./9.*mp2*mk2*pi16*L4r - 80./3.*pow(
         mp2,2)*pi16*L8r - 64./9.*pow(mp2,2)*pi16*L7r - 208./9.*pow(
         mp2,2)*pi16*L6r + 328./27.*pow(mp2,2)*pi16*L5r + 104./9.*pow(
         mp2,2)*pi16*L4r;

      return Qnr6L1/pow(mass.getf0(),4);
}

double qqup6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup6C\n";}

  double Qnr6C1 =
       + 192*pow(mk2,2)*CC[21] + 64*pow(mk2,2)*CC[20] + 8*mp2*mk2*CC[94] + 
         192*mp2*mk2*CC[21] - 4*pow(mp2,2)*CC[94] + 48*pow(mp2,2)*CC[21] + 80
         *pow(mp2,2)*CC[20] + 48*pow(mp2,2)*CC[19];

  return Qnr6C1/pow(mass.getf0(),4);
}

double qqup6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double Qnr6R1 =
       + pow(Ab(mp2,mu2),2) * (  - 9./8. );

      Qnr6R1 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 1./2. + 1./12.*mp2*pow(
         me2,-1) );

      Qnr6R1 +=  + Ab(mp2,mu2) * ( 2./3.*mp2*pi16 );

      Qnr6R1 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 5./6. - 1./18.*mp2*
         pow(me2,-1) );

      Qnr6R1 +=  + Ab(mk2,mu2) * ( 2./9.*mk2*pi16 );

      Qnr6R1 +=  + pow(Ab(me2,mu2),2) * ( 1./8. - 1./36.*mp2*pow(
         me2,-1) );

      Qnr6R1 +=  + Ab(me2,mu2) * ( 5./27.*mk2*pi16 - 5./27.*mp2*pi16 );

      return Qnr6R1/pow(mass.getf0(),4);
}

// *************************************************************************

double qqstrange4lo(const lomass mass, const Li Liin){
  return qqstrange4Llo(mass,Liin)+qqstrange4Rlo(mass);
}

double qqstrange4Llo(const lomass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange4L\n";}
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);

  double Qnr4L3 =
       + 8*mk2*H2r + 16*mk2*L8r + 32*mk2*L6r - 4*mp2*H2r - 8*mp2*L8r + 
         16*mp2*L6r;

  return Qnr4L3/pow(mass.getf0(),2);
}

double qqstrange4Rlo(const lomass mass){
  double mk2 = pow(mass.getmk0(),2);
  double mp2 = pow(mass.getmp0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);

  double  Qnr4R3 =
       + Ab(mk2,mu2) * ( 2 );

      Qnr4R3 +=  + Ab(me2,mu2) * ( 2./3. );

  return Qnr4R3/pow(mass.getf0(),2);
}

double qqstrange6lo(const lomass mass, const Li Liin, const Ci Ciin){
  return qqstrange6Rlo(mass)+qqstrange6Llo(mass,Liin)+qqstrange6Clo(mass,Ciin);
}

double qqstrange6LClo(const lomass mass, const Li Liin, const Ci Ciin){
  return qqstrange6Llo(mass,Liin)+qqstrange6Clo(mass,Ciin);
}

double qqstrange6Llo(const lomass mass, const Li Liin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange6L\n";}

  double Qnr6L3 =
       + Ab(mp2,mu2) * ( 48*mp2*L6r - 24*mp2*L4r );

      Qnr6L3 +=  + Ab(mk2,mu2) * ( 96*mk2*L8r + 192*mk2*L6r - 48*mk2*
         L5r - 96*mk2*L4r + 64*mp2*L6r - 32*mp2*L4r );

      Qnr6L3 +=  + Ab(me2,mu2) * ( 64*mk2*L8r + 64*mk2*L7r + 64*mk2*L6r
          - 64./3.*mk2*L5r - 32*mk2*L4r - 112./3.*mp2*L8r - 80*mp2*L7r
          + 16*mp2*L6r + 16./3.*mp2*L5r - 8*mp2*L4r + 16./3.*pow(mp2,2)
         *pow(me2,-1)*L8r + 16*pow(mp2,2)*pow(me2,-1)*L7r );

      Qnr6L3 +=  - 544./9.*pow(mk2,2)*pi16*L8r - 256./9.*pow(mk2,2)*
         pi16*L7r - 832./9.*pow(mk2,2)*pi16*L6r + 688./27.*pow(mk2,2)*
         pi16*L5r + 416./9.*pow(mk2,2)*pi16*L4r + 256./9.*mp2*mk2*pi16*
         L8r + 512./9.*mp2*mk2*pi16*L7r - 352./9.*mp2*mk2*pi16*L6r - 
         128./27.*mp2*mk2*pi16*L5r + 176./9.*mp2*mk2*pi16*L4r - 32./3.*
         pow(mp2,2)*pi16*L8r - 256./9.*pow(mp2,2)*pi16*L7r + 32./9.*
         pow(mp2,2)*pi16*L6r + 16./27.*pow(mp2,2)*pi16*L5r - 16./9.*
         pow(mp2,2)*pi16*L4r;

      return Qnr6L3/pow(mass.getf0(),4);
}

double qqstrange6Clo(const lomass mass, const Ci Ciin){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange6C\n";}


  double Qnr6C3 =
       + 192*pow(mk2,2)*CC[21] + 192*pow(mk2,2)*CC[20] + 192*pow(mk2,2)*
         CC[19] + 192*mp2*mk2*CC[21] - 64*mp2*mk2*CC[20] - 192*mp2*mk2*CC[19]
          + 4*pow(mp2,2)*CC[94] + 48*pow(mp2,2)*CC[21] + 16*pow(mp2,2)*CC[20]
          + 48*pow(mp2,2)*CC[19];

  return Qnr6C3/pow(mass.getf0(),4);
}

double qqstrange6Rlo(const lomass mass){
  double mp2 = pow(mass.getmp0(),2);
  double mk2 = pow(mass.getmk0(),2);
  double me2 = 4./3.*mk2-1./3.*mp2;
  double mu2 = pow(mass.getmu(),2);


  double Qnr6R3 =
       + Ab(mp2,mu2)*Ab(me2,mu2) * ( 1./3.*mp2*pow(me2,-1) );

      Qnr6R3 +=  + Ab(mp2,mu2) * (  - 1./3.*mp2*pi16 );

      Qnr6R3 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 2 - 2./9.*mp2*pow(
         me2,-1) );

      Qnr6R3 +=  + Ab(mk2,mu2) * ( 8./9.*mk2*pi16 );

      Qnr6R3 +=  + pow(Ab(me2,mu2),2) * ( 2./3. - 1./9.*mp2*pow(me2,-1)
          );

      Qnr6R3 +=  + Ab(me2,mu2) * ( 2./27.*mk2*pi16 + 7./27.*mp2*pi16 );

      return Qnr6R3/pow(mass.getf0(),4);
}
