// massdecayvev.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results derived in
//  G.~Amoros, J.~Bijnens and P.~Talavera,
//  Two point functions at two loops in three flavor chiral perturbation theory,
//  Nucl.\ Phys.\ B {\bf 568} (2000) 319 [hep-ph/9907264]
//
// G.~Amoros, J.~Bijnens and P.~Talavera,
// ``K(lepton 4) form-factors and pi pi scattering,''
// Nucl.\ Phys.\ B {\bf 585} (2000) 293
// [Erratum-ibid.\ B {\bf 598} (2001) 665] [hep-ph/0003258].

// contains the isospin limit expressions for masses, decay
// constants and vacuum expectation values 
// in terms of physical masses and Fpi

#include <iostream>
#include <iomanip>
#include <string>

#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "oneloopintegrals.h"
#include "sunsetintegrals.h"
#include "massdecayvev.h"

const double pi = M_PI;
const double pi16 = 1./(16.*pi*pi);

// mpi ******************************************************************

double mpi4(const physmass mass, const Li Liin){
  return mpi4L(mass,Liin)+mpi4R(mass);
}

double mpi4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double   M4L11 =
    + 32*mp2*mk2*L6r - 16*mp2*mk2*L4r + 16*pow(mp2,2)*L8r + 16*pow(
         mp2,2)*L6r - 8*pow(mp2,2)*L5r - 8*pow(mp2,2)*L4r;
  return M4L11/pow(mass.getfpi(),2);
}

double mpi4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double M4R11 =
    + Ab(mp2,mu2) * (  - 1./2.*mp2 );
  M4R11 +=  + Ab(me2,mu2) * ( 1./6.*mp2 );
  return M4R11/pow(mass.getfpi(),2);
}

double mpi6(const physmass mass, const Li Liin, const Ci Ciin){
  return mpi6R(mass)+mpi6L(mass,Liin)+mpi6C(mass,Ciin);
}

double mpi6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return mpi6L(mass,Liin)+mpi6C(mass,Ciin);
}

double mpi6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6L\n";}

  double M6L11 =
       + Ab(mp2,mu2) * ( 80*mp2*mk2*L6r - 40*mp2*mk2*L4r + 88*pow(
         mp2,2)*L8r + 128*pow(mp2,2)*L6r - 48*pow(mp2,2)*L5r - 72*pow(
         mp2,2)*L4r + 28*pow(mp2,2)*L3r + 32*pow(mp2,2)*L2r + 56*pow(
         mp2,2)*L1r );

      M6L11 +=  + Ab(mk2,mu2) * ( 32*mp2*mk2*L8r + 96*mp2*mk2*L6r -
         16*mp2*mk2*L5r - 80*mp2*mk2*L4r + 20*mp2*mk2*L3r + 16*mp2*mk2*
         L2r + 64*mp2*mk2*L1r + 16*pow(mp2,2)*L8r + 16*pow(mp2,2)*L6r
          - 8*pow(mp2,2)*L5r - 8*pow(mp2,2)*L4r );

      M6L11 +=  + Ab(me2,mu2) * (  - 64./3.*mp2*mk2*L7r + 80./3.*mp2
         *mk2*L6r - 32./9.*mp2*mk2*L5r - 24*mp2*mk2*L4r + 16./3.*mp2*
         mk2*L3r + 16./3.*mp2*mk2*L2r + 64./3.*mp2*mk2*L1r + 8./3.*pow(
         mp2,2)*L8r + 64./3.*pow(mp2,2)*L7r - 32./3.*pow(mp2,2)*L6r +
         32./9.*pow(mp2,2)*L5r + 8*pow(mp2,2)*L4r - 4./3.*pow(mp2,2)*
         L3r - 4./3.*pow(mp2,2)*L2r - 16./3.*pow(mp2,2)*L1r );

      M6L11 +=  - 512*mp2*pow(mk2,2)*L6r*L8r - 2048*mp2*pow(mk2,2)*pow(
         L6r,2) + 256*mp2*pow(mk2,2)*L5r*L6r + 256*mp2*pow(mk2,2)*L4r*
         L8r + 2048*mp2*pow(mk2,2)*L4r*L6r - 128*mp2*pow(mk2,2)*L4r*L5r
          - 512*mp2*pow(mk2,2)*pow(L4r,2) + 86./27.*mp2*pow(mk2,2)*pi16
         *L3r + 104./9.*mp2*pow(mk2,2)*pi16*L2r - 1536*pow(mp2,2)*mk2*
         L6r*L8r - 2048*pow(mp2,2)*mk2*pow(L6r,2) + 768*pow(mp2,2)*mk2*
         L5r*L6r + 768*pow(mp2,2)*mk2*L4r*L8r + 2048*pow(mp2,2)*mk2*L4r
         *L6r - 384*pow(mp2,2)*mk2*L4r*L5r - 512*pow(mp2,2)*mk2*pow(
         L4r,2) - 16./27.*pow(mp2,2)*mk2*pi16*L3r - 16./9.*pow(mp2,2)*
         mk2*pi16*L2r - 512*pow(mp2,3)*pow(L8r,2) - 1024*pow(mp2,3)*L6r
         *L8r - 512*pow(mp2,3)*pow(L6r,2) + 512*pow(mp2,3)*L5r*L8r +
         512*pow(mp2,3)*L5r*L6r - 128*pow(mp2,3)*pow(L5r,2) + 512*pow(
         mp2,3)*L4r*L8r + 512*pow(mp2,3)*L4r*L6r - 256*pow(mp2,3)*L4r*
         L5r - 128*pow(mp2,3)*pow(L4r,2) + 56./27.*pow(mp2,3)*pi16*L3r
          + 74./9.*pow(mp2,3)*pi16*L2r + 4*pow(mp2,3)*pi16*L1r;
      return M6L11/pow(mass.getfpi(),4);
}

double mpi6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mpi6C\n";}


  double M6C11 =
      + 192*mp2*pow(mk2,2)*CC[21] + 64*mp2*pow(mk2,2)*CC[20] - 64*
         mp2*pow(mk2,2)*CC[16] + 64*pow(mp2,2)*mk2*CC[32] + 192*
         pow(mp2,2)*mk2*CC[21] + 64*pow(mp2,2)*mk2*CC[16] - 32*pow(
         mp2,2)*mk2*CC[15] - 64*pow(mp2,2)*mk2*CC[13] + 32*pow(
         mp2,3)*CC[32] + 32*pow(mp2,3)*CC[31] + 48*pow(mp2,3)*
         CC[21] + 80*pow(mp2,3)*CC[20] + 48*pow(mp2,3)*CC[19] -
         16*pow(mp2,3)*CC[17] - 48*pow(mp2,3)*CC[16] - 16*pow(
         mp2,3)*CC[15] - 16*pow(mp2,3)*CC[14] - 32*pow(mp2,3)*
         CC[13] - 32*pow(mp2,3)*CC[12];

  return M6C11/pow(mass.getfpi(),4);
}

double mpi6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double M6R11 =
       + pow(Ab(mp2,mu2),2) * (  - 1./2.*mk2 - 589./144.*mp2 );

      M6R11 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - mp2 );

      M6R11 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 5./12.*mp2 );

      M6R11 +=  + Ab(mp2,mu2) * ( mp2*mk2*pi16 + pow(mp2,2)*pi16 );

      M6R11 +=  + pow(Ab(mk2,mu2),2) * (  - 7./4.*mp2 - 1./6.*pow(
         mp2,2)*pow(mk2,-1) );

      M6R11 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 1./3.*mp2 );

      M6R11 +=  + Ab(mk2,mu2) * ( mp2*mk2*pi16 );

      M6R11 +=  + pow(Ab(me2,mu2),2) * (  - 29./144.*mp2 - 5./36.*
         pow(mp2,2)*pow(me2,-1) );

      M6R11 +=  + hh(mp2,mp2,mp2,mp2,mu2) * ( 5./6.*pow(mp2,2) );

      M6R11 +=  + hh(mp2,mk2,mk2,mp2,mu2) * (  - 5./8.*pow(mp2,2) );

      M6R11 +=  + hh(mp2,me2,me2,mp2,mu2) * ( 1./18.*pow(mp2,2) );

      M6R11 +=  + hh(mk2,mp2,mk2,mp2,mu2) * ( mp2*mk2 );

      M6R11 +=  + hh(mk2,mk2,me2,mp2,mu2) * (  - 5./6.*pow(mp2,2) );

      M6R11 +=  + hh(me2,mk2,mk2,mp2,mu2) * ( 1./2.*mp2*mk2 - 1./8.*
         pow(mp2,2) );

      M6R11 +=  + hh1(mp2,mk2,mk2,mp2,mu2) * ( pow(mp2,2) );

      M6R11 +=  + hh1(mk2,mk2,me2,mp2,mu2) * ( 2*pow(mp2,2) );

      M6R11 +=  + hh21(mp2,mp2,mp2,mp2,mu2) * ( 3*pow(mp2,2) );

      M6R11 +=  + hh21(mp2,mk2,mk2,mp2,mu2) * (  - 3./8.*pow(mp2,2) );

      M6R11 +=  + hh21(mk2,mp2,mk2,mp2,mu2) * ( 3*pow(mp2,2) );

      M6R11 +=  + hh21(me2,mk2,mk2,mp2,mu2) * ( 9./8.*pow(mp2,2) );

      M6R11 +=  - 15./16.*mp2*pow(mk2,2)*pow(pi16,2) - 11./36.*mp2*pow(
         mk2,2)*pow(pi16,2)*pow(pi,2) - 139./216.*pow(mp2,2)*mk2*pow(
         pi16,2) - 37./324.*pow(mp2,2)*mk2*pow(pi16,2)*pow(pi,2) - 3217.
         /1728.*pow(mp2,3)*pow(pi16,2) - 527./1296.*pow(mp2,3)*pow(
         pi16,2)*pow(pi,2);

      return M6R11/pow(mass.getfpi(),4);
}

// *************************************************************************

double mk4(const physmass mass, const Li Liin){
  return mk4L(mass,Liin)+mk4R(mass);
}

double mk4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double   M4L44 =
       + 16*pow(mk2,2)*L8r + 32*pow(mk2,2)*L6r - 8*pow(mk2,2)*L5r - 16*
         pow(mk2,2)*L4r + 16*mp2*mk2*L6r - 8*mp2*mk2*L4r;

  return M4L44/pow(mass.getfpi(),2);
}

double mk4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

double   M4R44 =
       + Ab(me2,mu2) * (  - 1./4.*me2 - 1./12.*mp2 );

  return M4R44/pow(mass.getfpi(),2);
}

double mk6(const physmass mass, const Li Liin, const Ci Ciin){
  return mk6R(mass)+mk6L(mass,Liin)+mk6C(mass,Ciin);
}

double mk6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return mk6L(mass,Liin)+mk6C(mass,Ciin);
}

double mk6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk6L\n";}

  double   M6L44 =
       + Ab(mp2,mu2) * ( 32*pow(mk2,2)*L8r + 64*pow(mk2,2)*L6r - 16*
         pow(mk2,2)*L5r - 32*pow(mk2,2)*L4r + 24*mp2*mk2*L8r + 88*mp2*
         mk2*L6r - 12*mp2*mk2*L5r - 68*mp2*mk2*L4r + 15*mp2*mk2*L3r + 
         12*mp2*mk2*L2r + 48*mp2*mk2*L1r );

      M6L44 +=  + Ab(mk2,mu2) * ( 64*pow(mk2,2)*L8r + 128*pow(mk2,2)
         *L6r - 32*pow(mk2,2)*L5r - 80*pow(mk2,2)*L4r + 30*pow(mk2,2)*
         L3r + 36*pow(mk2,2)*L2r + 72*pow(mk2,2)*L1r + 16*mp2*mk2*L6r
          - 8*mp2*mk2*L4r );

      M6L44 +=  + Ab(me2,mu2) * ( 32*pow(mk2,2)*L8r + 64./3.*pow(
         mk2,2)*L7r + 128./3.*pow(mk2,2)*L6r - 112./9.*pow(mk2,2)*L5r
          - 32*pow(mk2,2)*L4r + 28./3.*pow(mk2,2)*L3r + 16./3.*pow(
         mk2,2)*L2r + 64./3.*pow(mk2,2)*L1r - 56./3.*mp2*mk2*L8r - 32*
         mp2*mk2*L7r - 8./3.*mp2*mk2*L6r + 4./3.*mp2*mk2*L5r + 4*mp2*
         mk2*L4r - 7./3.*mp2*mk2*L3r - 4./3.*mp2*mk2*L2r - 16./3.*mp2*
         mk2*L1r + 16./3.*pow(mp2,2)*L8r + 32./3.*pow(mp2,2)*L7r - 8./9.
         *pow(mp2,2)*L5r );

      M6L44 +=  - 512*pow(mk2,3)*pow(L8r,2) - 2048*pow(mk2,3)*L6r*L8r
          - 2048*pow(mk2,3)*pow(L6r,2) + 384*pow(mk2,3)*L5r*L8r + 768*
         pow(mk2,3)*L5r*L6r - 64*pow(mk2,3)*pow(L5r,2) + 1024*pow(
         mk2,3)*L4r*L8r + 2048*pow(mk2,3)*L4r*L6r - 384*pow(mk2,3)*L4r*
         L5r - 512*pow(mk2,3)*pow(L4r,2) + 89./27.*pow(mk2,3)*pi16*L3r
          + 122./9.*pow(mk2,3)*pi16*L2r + 4*pow(mk2,3)*pi16*L1r - 768*
         mp2*pow(mk2,2)*L6r*L8r - 2048*mp2*pow(mk2,2)*pow(L6r,2) + 128*
         mp2*pow(mk2,2)*L5r*L8r + 512*mp2*pow(mk2,2)*L5r*L6r - 64*mp2*
         pow(mk2,2)*pow(L5r,2) + 384*mp2*pow(mk2,2)*L4r*L8r + 2048*mp2*
         pow(mk2,2)*L4r*L6r - 256*mp2*pow(mk2,2)*L4r*L5r - 512*mp2*pow(
         mk2,2)*pow(L4r,2) - 4./27.*mp2*pow(mk2,2)*pi16*L3r - 16./9.*
         mp2*pow(mk2,2)*pi16*L2r - 256*pow(mp2,2)*mk2*L6r*L8r - 512*
         pow(mp2,2)*mk2*pow(L6r,2) + 256*pow(mp2,2)*mk2*L5r*L6r + 128*
         pow(mp2,2)*mk2*L4r*L8r + 512*pow(mp2,2)*mk2*L4r*L6r - 128*pow(
         mp2,2)*mk2*L4r*L5r - 128*pow(mp2,2)*mk2*pow(L4r,2) + 41./27.*
         pow(mp2,2)*mk2*pi16*L3r;
      M6L44 +=  + 56./9.*pow(mp2,2)*mk2*pi16*L2r;

      return M6L44/pow(mass.getfpi(),4);
}

double mk6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in mk6C\n";}

  double   M6C44 =
       + 64*pow(mk2,3)*CC[32] + 32*pow(mk2,3)*CC[31] + 192*pow(
         mk2,3)*CC[21] + 128*pow(mk2,3)*CC[20] + 96*pow(mk2,3)*
         CC[19] - 64*pow(mk2,3)*CC[16] - 32*pow(mk2,3)*CC[15] - 
         32*pow(mk2,3)*CC[14] - 64*pow(mk2,3)*CC[13] - 32*pow(
         mk2,3)*CC[12] + 32*mp2*pow(mk2,2)*CC[32] + 192*mp2*pow(
         mk2,2)*CC[21] - 32*mp2*pow(mk2,2)*CC[20] - 96*mp2*pow(
         mk2,2)*CC[19] - 32*mp2*pow(mk2,2)*CC[17] + 64*mp2*pow(
         mk2,2)*CC[16] - 16*mp2*pow(mk2,2)*CC[15] + 32*mp2*pow(
         mk2,2)*CC[14] - 32*mp2*pow(mk2,2)*CC[13] + 48*pow(mp2,2)*
         mk2*CC[21] + 48*pow(mp2,2)*mk2*CC[20] + 48*pow(mp2,2)*mk2*
         CC[19] + 16*pow(mp2,2)*mk2*CC[17] - 48*pow(mp2,2)*mk2*
         CC[16] - 16*pow(mp2,2)*mk2*CC[14];

  return M6C44/pow(mass.getfpi(),4);
}

double mk6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   M6R44 =
       + pow(Ab(mp2,mu2),2) * (  - 1./2.*pow(mp2,-1)*pow(mk2,2) - 27.
         /32.*mk2 );

      M6R44 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - 3./4.*mk2 );

      M6R44 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * (  - 41./48.*mk2 + 1./
         12.*mp2 );

      M6R44 +=  + Ab(mp2,mu2) * ( 3./4.*pow(mk2,2)*pi16 );

      M6R44 +=  + pow(Ab(mk2,mu2),2) * (  - 251./72.*mk2 - 3./8.*mp2
          );

      M6R44 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 2./3.*mk2 );

      M6R44 +=  + Ab(mk2,mu2) * ( 3./4.*pow(mk2,2)*pi16 + 3./4.*mp2*
         mk2*pi16 );

      M6R44 +=  + pow(Ab(me2,mu2),2) * (  - 5./36.*mk2 - 25./128.*
         mp2 - 43./1152.*pow(mp2,2)*pow(me2,-1) );

      M6R44 +=  + Ab(me2,mu2) * ( 1./2.*pow(mk2,2)*pi16 + 1./4.*mp2*
         mk2*pi16 );

      M6R44 +=  + hh(mp2,mp2,mk2,mk2,mu2) * ( 3./8.*pow(mk2,2) + 3./4.*
         mp2*mk2 );

      M6R44 +=  + hh(mp2,mk2,me2,mk2,mu2) * ( 1./4.*pow(mk2,2) );

      M6R44 +=  + hh(mk2,mp2,mp2,mk2,mu2) * (  - 3./32.*pow(mk2,2) );

      M6R44 +=  + hh(mk2,mp2,me2,mk2,mu2) * ( 9./16.*pow(mk2,2) );

      M6R44 +=  + hh(mk2,mk2,mk2,mk2,mu2) * ( 3./4.*pow(mk2,2) );

      M6R44 +=  + hh(mk2,me2,me2,mk2,mu2) * ( 181./288.*pow(mk2,2) );

      M6R44 +=  + hh1(mp2,mp2,mk2,mk2,mu2) * (  - 3./2.*pow(mk2,2) );

      M6R44 +=  + hh1(mk2,mp2,me2,mk2,mu2) * (  - 3./2.*pow(mk2,2) );

      M6R44 +=  + hh1(mk2,me2,me2,mk2,mu2) * (  - 5./4.*pow(mk2,2) );

      M6R44 +=  + hh21(mp2,mp2,mk2,mk2,mu2) * ( 9./4.*pow(mk2,2) );

      M6R44 +=  + hh21(mk2,mp2,mp2,mk2,mu2) * (  - 9./32.*pow(mk2,2) );

      M6R44 +=  + hh21(mk2,mp2,me2,mk2,mu2) * ( 27./16.*pow(mk2,2) );

      M6R44 +=  + hh21(mk2,mk2,mk2,mk2,mu2) * ( 9./4.*pow(mk2,2) );

      M6R44 +=  + hh21(mk2,me2,me2,mk2,mu2) * ( 27./32.*pow(mk2,2) );

      M6R44 +=  - 4709./1728.*pow(mk2,3)*pow(pi16,2) - 763./1296.*pow(
         mk2,3)*pow(pi16,2)*pow(pi,2) - 19./108.*mp2*pow(mk2,2)*pow(
         pi16,2) - 73./648.*mp2*pow(mk2,2)*pow(pi16,2)*pow(pi,2) - 13./
         24.*pow(mp2,2)*mk2*pow(pi16,2) - 1./8.*pow(mp2,2)*mk2*pow(
         pi16,2)*pow(pi,2);

      return M6R44/pow(mass.getfpi(),4);
}


// *************************************************************************

double meta4(const physmass mass, const Li Liin){
  return meta4L(mass,Liin)+meta4R(mass);
}

double meta4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);

  double   M4L88 =
       - 32./3.*mk2*me2*L5r - 16*mk2*me2*L4r + 128./3.*pow(mk2,2)*L8r
          + 128./3.*pow(mk2,2)*L7r + 128./3.*pow(mk2,2)*L6r + 8./3.*mp2
         *me2*L5r - 8*mp2*me2*L4r - 128./3.*mp2*mk2*L8r - 256./3.*mp2*
         mk2*L7r + 32./3.*mp2*mk2*L6r + 16*pow(mp2,2)*L8r + 128./3.*
         pow(mp2,2)*L7r - 16./3.*pow(mp2,2)*L6r;

  return M4L88/pow(mass.getfpi(),2);
}

double meta4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

double   M4R88 =
       + Ab(mp2,mu2) * ( 1./2.*mp2 );

      M4R88 +=  + Ab(mk2,mu2) * (  - me2 - 1./3.*mp2 );

      M4R88 +=  + Ab(me2,mu2) * ( 8./9.*mk2 - 7./18.*mp2 );

  return M4R88/pow(mass.getfpi(),2);
}

double meta6(const physmass mass, const Li Liin, const Ci Ciin){
  return meta6R(mass)+meta6L(mass,Liin)+meta6C(mass,Ciin);
}

double meta6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return meta6L(mass,Liin)+meta6C(mass,Ciin);
}

double meta6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta6L\n";}

  double M6L88 =
       + Ab(mp2,mu2) * ( 256./3.*pow(mk2,2)*L8r + 256./3.*pow(mk2,2)
         *L7r + 256./3.*pow(mk2,2)*L6r - 256./9.*pow(mk2,2)*L5r - 128./
         3.*pow(mk2,2)*L4r - 320./3.*mp2*mk2*L8r - 832./3.*mp2*mk2*L7r
          + 272./3.*mp2*mk2*L6r + 32./3.*mp2*mk2*L5r - 72*mp2*mk2*L4r
          + 16*mp2*mk2*L3r + 16*mp2*mk2*L2r + 64*mp2*mk2*L1r + 72*pow(
         mp2,2)*L8r + 192*pow(mp2,2)*L7r - 32*pow(mp2,2)*L6r - 8./9.*
         pow(mp2,2)*L5r + 80./3.*pow(mp2,2)*L4r - 4*pow(mp2,2)*L3r - 4*
         pow(mp2,2)*L2r - 16*pow(mp2,2)*L1r );

      M6L88 +=  + Ab(mk2,mu2) * ( 128*pow(mk2,2)*L8r + 128*pow(
         mk2,2)*L7r + 128*pow(mk2,2)*L6r - 512./9.*pow(mk2,2)*L5r - 128
         *pow(mk2,2)*L4r + 112./3.*pow(mk2,2)*L3r + 64./3.*pow(mk2,2)*
         L2r + 256./3.*pow(mk2,2)*L1r - 352./3.*mp2*mk2*L8r - 640./3.*
         mp2*mk2*L7r - 32./3.*mp2*mk2*L6r + 16*mp2*mk2*L5r + 16./3.*mp2
         *mk2*L4r - 28./3.*mp2*mk2*L3r - 16./3.*mp2*mk2*L2r - 64./3.*
         mp2*mk2*L1r + 112./3.*pow(mp2,2)*L8r + 256./3.*pow(mp2,2)*L7r
          - 16./3.*pow(mp2,2)*L6r - 40./9.*pow(mp2,2)*L5r + 8./3.*pow(
         mp2,2)*L4r );

      M6L88 +=  + Ab(me2,mu2) * ( 128*pow(mk2,2)*L8r + 1280./9.*pow(
         mk2,2)*L7r + 1024./9.*pow(mk2,2)*L6r - 704./27.*pow(mk2,2)*L5r
          - 64./3.*pow(mk2,2)*L4r + 64./3.*pow(mk2,2)*L3r + 128./3.*
         pow(mk2,2)*L2r + 128./3.*pow(mk2,2)*L1r - 1088./9.*mp2*mk2*L8r
          - 1792./9.*mp2*mk2*L7r - 368./9.*mp2*mk2*L6r + 736./27.*mp2*
         mk2*L5r + 56./3.*mp2*mk2*L4r - 32./3.*mp2*mk2*L3r - 64./3.*mp2
         *mk2*L2r - 64./3.*mp2*mk2*L1r + 296./9.*pow(mp2,2)*L8r + 512./
         9.*pow(mp2,2)*L7r + 64./9.*pow(mp2,2)*L6r - 248./27.*pow(
         mp2,2)*L5r - 16./3.*pow(mp2,2)*L4r + 4./3.*pow(mp2,2)*L3r + 8./
         3.*pow(mp2,2)*L2r + 8./3.*pow(mp2,2)*L1r );

      M6L88 +=  - 4096./3.*pow(mk2,3)*pow(L8r,2) - 4096./3.*pow(mk2,3)*
         L7r*L8r - 4096*pow(mk2,3)*L6r*L8r - 8192./3.*pow(mk2,3)*L6r*
         L7r - 8192./3.*pow(mk2,3)*pow(L6r,2) + 8192./9.*pow(mk2,3)*L5r
         *L8r + 2048./3.*pow(mk2,3)*L5r*L7r + 10240./9.*pow(mk2,3)*L5r*
         L6r - 1024./9.*pow(mk2,3)*pow(L5r,2) + 7168./3.*pow(mk2,3)*L4r
         *L8r + 2048*pow(mk2,3)*L4r*L7r + 8192./3.*pow(mk2,3)*L4r*L6r
          - 5632./9.*pow(mk2,3)*L4r*L5r - 2048./3.*pow(mk2,3)*pow(
         L4r,2) + 152./27.*pow(mk2,3)*pi16*L3r + 544./27.*pow(mk2,3)*
         pi16*L2r + 256./27.*pow(mk2,3)*pi16*L1r + 2048./3.*mp2*pow(
         mk2,2)*pow(L8r,2) + 4096./3.*mp2*pow(mk2,2)*L7r*L8r + 3584./3.
         *mp2*pow(mk2,2)*L6r*L8r + 4096*mp2*pow(mk2,2)*L6r*L7r - 2048*
         mp2*pow(mk2,2)*pow(L6r,2) - 512./9.*mp2*pow(mk2,2)*L5r*L8r - 
         1024./3.*mp2*pow(mk2,2)*L5r*L7r + 1280./3.*mp2*pow(mk2,2)*L5r*
         L6r - 256./3.*mp2*pow(mk2,2)*pow(L5r,2) - 3328./3.*mp2*pow(
         mk2,2)*L4r*L8r - 3072*mp2*pow(mk2,2)*L4r*L7r + 2048*mp2*pow(
         mk2,2)*L4r*L6r;
      M6L88 +=  - 128*mp2*pow(mk2,2)*L4r*L5r - 512*mp2*pow(mk2,2)*pow(
         L4r,2) - 34./9.*mp2*pow(mk2,2)*pi16*L3r - 88./9.*mp2*pow(
         mk2,2)*pi16*L2r - 64./9.*mp2*pow(mk2,2)*pi16*L1r + 2048./3.*
         pow(mp2,2)*mk2*pow(L8r,2) + 4096./3.*pow(mp2,2)*mk2*L7r*L8r + 
         512./3.*pow(mp2,2)*mk2*L6r*L8r - 6656./9.*pow(mp2,2)*mk2*L5r*
         L8r - 4096./3.*pow(mp2,2)*mk2*L5r*L7r + 256./3.*pow(mp2,2)*mk2
         *L5r*L6r + 256./3.*pow(mp2,2)*mk2*pow(L5r,2) - 256./3.*pow(
         mp2,2)*mk2*L4r*L8r - 128./3.*pow(mp2,2)*mk2*L4r*L5r + 32./9.*
         pow(mp2,2)*mk2*pi16*L3r + 88./9.*pow(mp2,2)*mk2*pi16*L2r + 16./
         9.*pow(mp2,2)*mk2*pi16*L1r - 512*pow(mp2,3)*pow(L8r,2) - 4096./
         3.*pow(mp2,3)*L7r*L8r - 1024./3.*pow(mp2,3)*L6r*L8r - 4096./3.
         *pow(mp2,3)*L6r*L7r + 512./3.*pow(mp2,3)*pow(L6r,2) + 3584./9.
         *pow(mp2,3)*L5r*L8r + 1024*pow(mp2,3)*L5r*L7r - 1024./9.*pow(
         mp2,3)*L5r*L6r - 128./9.*pow(mp2,3)*pow(L5r,2) + 1024./3.*pow(
         mp2,3)*L4r*L8r + 1024*pow(mp2,3)*L4r*L7r - 512./3.*pow(mp2,3)*
         L4r*L6r;
      M6L88 +=  + 256./9.*pow(mp2,3)*L4r*L5r + 128./3.*pow(mp2,3)*pow(
         L4r,2) - 20./27.*pow(mp2,3)*pi16*L3r - 58./27.*pow(mp2,3)*pi16
         *L2r - 4./27.*pow(mp2,3)*pi16*L1r;

      return M6L88/pow(mass.getfpi(),4);
}

double meta6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in meta6C\n";}

  double   M6C88 =
       + 512./3.*pow(mk2,3)*CC[33] + 512./3.*pow(mk2,3)*CC[32] + 
         512./3.*pow(mk2,3)*CC[31] + 256*pow(mk2,3)*CC[21] + 256*
         pow(mk2,3)*CC[20] + 256*pow(mk2,3)*CC[19] - 512./9.*pow(
         mk2,3)*CC[18] - 512./9.*pow(mk2,3)*CC[17] - 256./3.*pow(
         mk2,3)*CC[16] - 512./9.*pow(mk2,3)*CC[15] - 512./9.*pow(
         mk2,3)*CC[14] - 1024./9.*pow(mk2,3)*CC[13] - 2048./27.*
         pow(mk2,3)*CC[12] - 1024./3.*mp2*pow(mk2,2)*CC[33] - 256./
         3.*mp2*pow(mk2,2)*CC[32] - 256*mp2*pow(mk2,2)*CC[31] + 192
         *mp2*pow(mk2,2)*CC[21] - 192*mp2*pow(mk2,2)*CC[20] - 384*
         mp2*pow(mk2,2)*CC[19] + 128*mp2*pow(mk2,2)*CC[18] + 640./9.
         *mp2*pow(mk2,2)*CC[17] + 320./3.*mp2*pow(mk2,2)*CC[16] + 
         640./9.*mp2*pow(mk2,2)*CC[14] + 512./9.*mp2*pow(mk2,2)*
         CC[12] + 512./3.*pow(mp2,2)*mk2*CC[33] - 64./3.*pow(mp2,2)
         *mk2*CC[32] + 128*pow(mp2,2)*mk2*CC[31] + 64*pow(mp2,2)*
         mk2*CC[20];
      M6C88 +=  + 192*pow(mp2,2)*mk2*CC[19] - 256./3.*pow(mp2,2)*mk2*
         CC[18] - 320./9.*pow(mp2,2)*mk2*CC[17] - 256./3.*pow(
         mp2,2)*mk2*CC[16] + 32./3.*pow(mp2,2)*mk2*CC[15] - 320./9.
         *pow(mp2,2)*mk2*CC[14] + 64./3.*pow(mp2,2)*mk2*CC[13] - 
         128./9.*pow(mp2,2)*mk2*CC[12] + 32*pow(mp2,3)*CC[32] - 32./
         3.*pow(mp2,3)*CC[31] - 16*pow(mp2,3)*CC[21] + 16*pow(
         mp2,3)*CC[20] - 16*pow(mp2,3)*CC[19] + 128./9.*pow(mp2,3)*
         CC[18] + 16./3.*pow(mp2,3)*CC[17] + 16*pow(mp2,3)*CC[16]
          - 16./9.*pow(mp2,3)*CC[15] + 16./3.*pow(mp2,3)*CC[14] - 
         32./9.*pow(mp2,3)*CC[13] + 32./27.*pow(mp2,3)*CC[12];

  return M6C88/pow(mass.getfpi(),4);
}

double meta6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   M6R88 =
       + pow(Ab(mp2,mu2),2) * (  - 3./4.*mk2 + 43./48.*mp2 );

      M6R88 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * (  - 8./3.*mk2 - 2./3.
         *mp2 );

      M6R88 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 16./9.*mk2 - 35./36.
         *mp2 );

      M6R88 +=  + pow(Ab(mk2,mu2),2) * (  - 53./9.*mk2 + 43./12.*mp2
          - pow(mp2,2)*pow(mk2,-1) );

      M6R88 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 16./9.*mk2 + 7./
         9.*mp2 );

      M6R88 +=  + Ab(mk2,mu2) * ( 8./3.*pow(mk2,2)*pi16 + 2./3.*mp2*
         mk2*pi16 - 1./3.*pow(mp2,2)*pi16 );

      M6R88 +=  + pow(Ab(me2,mu2),2) * (  - 25./12.*mk2 + 55./48.*
         mp2 - 5./36.*pow(mp2,2)*pow(me2,-1) );

      M6R88 +=  + hh(mp2,mp2,me2,me2,mu2) * ( 1./6.*pow(mp2,2) );

      M6R88 +=  + hh(mp2,mk2,mk2,me2,mu2) * ( 3./2.*mp2*mk2 + 1./8.*
         pow(mp2,2) );

      M6R88 +=  + hh(mk2,mk2,me2,me2,mu2) * (  - 64./9.*pow(mk2,2) + 4*
         mp2*mk2 - 1./2.*pow(mp2,2) );

      M6R88 +=  + hh(me2,mk2,mk2,me2,mu2) * ( 2*pow(mk2,2) - mp2*mk2 + 
         1./8.*pow(mp2,2) );

      M6R88 +=  + hh(me2,me2,me2,me2,mu2) * ( 128./243.*pow(mk2,2) - 
         112./243.*mp2*mk2 + 49./486.*pow(mp2,2) );

      M6R88 +=  + hh1(mp2,mk2,mk2,me2,mu2) * (  - 4*mp2*mk2 + pow(
         mp2,2) );

      M6R88 +=  + hh1(mk2,mk2,me2,me2,mu2) * ( 64./3.*pow(mk2,2) - 40./
         3.*mp2*mk2 + 2*pow(mp2,2) );

      M6R88 +=  + hh21(mp2,mk2,mk2,me2,mu2) * ( 6*pow(mk2,2) - 3*mp2*
         mk2 + 3./8.*pow(mp2,2) );

      M6R88 +=  + hh21(me2,mk2,mk2,me2,mu2) * ( 6*pow(mk2,2) - 3*mp2*
         mk2 + 3./8.*pow(mp2,2) );

      M6R88 +=  - 7567./972.*pow(mk2,3)*pow(pi16,2) - 1091./729.*pow(
         mk2,3)*pow(pi16,2)*pow(pi,2) + 367./48.*mp2*pow(mk2,2)*pow(
         pi16,2) + 133./108.*mp2*pow(mk2,2)*pow(pi16,2)*pow(pi,2) - 
         4133./1296.*pow(mp2,2)*mk2*pow(pi16,2) - 269./486.*pow(mp2,2)*
         mk2*pow(pi16,2)*pow(pi,2) - 1781./15552.*pow(mp2,3)*pow(
         pi16,2) - 91./11664.*pow(mp2,3)*pow(pi16,2)*pow(pi,2);

      return M6R88/pow(mass.getfpi(),4);
}

// *************************************************************************

double fpi4(const physmass mass, const Li Liin){
  return fpi4L(mass,Liin)+fpi4R(mass);
}

double fpi4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);


  double   F4L11 =
       + 8*mk2*L4r + 4*mp2*L5r + 4*mp2*L4r;

  return F4L11/pow(mass.getfpi(),2);
}

double fpi4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double mu2 = pow(mass.getmu(),2);

  double   F4R11 =
       + Ab(mp2,mu2) * ( 1 );

  F4R11 +=  + Ab(mk2,mu2) * ( 1./2. );

  return F4R11/pow(mass.getfpi(),2);
}

double fpi6(const physmass mass, const Li Liin, const Ci Ciin){
  return fpi6R(mass)+fpi6L(mass,Liin)+fpi6C(mass,Ciin);
}

double fpi6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return fpi6L(mass,Liin)+fpi6C(mass,Ciin);
}

double fpi6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6L\n";}

  double   F6L11 =
       + Ab(mp2,mu2) * ( 24*mk2*L4r + 20*mp2*L5r + 26*mp2*L4r - 14*
         mp2*L3r - 16*mp2*L2r - 28*mp2*L1r );

  F6L11 +=  + Ab(mk2,mu2) * ( 28*mk2*L4r - 10*mk2*L3r - 8*mk2*
         L2r - 32*mk2*L1r + 10*mp2*L5r + 6*mp2*L4r );

  F6L11 +=  + Ab(me2,mu2) * ( 8*mk2*L4r - 8./3.*mk2*L3r - 8./3.*
         mk2*L2r - 32./3.*mk2*L1r - 2*mp2*L4r + 2./3.*mp2*L3r + 2./3.*
         mp2*L2r + 8./3.*mp2*L1r );

  F6L11 +=  - 128*pow(mk2,2)*L4r*L8r - 256*pow(mk2,2)*L4r*L6r + 64*
         pow(mk2,2)*L4r*L5r + 224*pow(mk2,2)*pow(L4r,2) - 43./27.*pow(
         mk2,2)*pi16*L3r - 52./9.*pow(mk2,2)*pi16*L2r - 128*mp2*mk2*L5r
         *L6r - 256*mp2*mk2*L4r*L6r + 160*mp2*mk2*L4r*L5r + 224*mp2*mk2
         *pow(L4r,2) + 8./27.*mp2*mk2*pi16*L3r + 8./9.*mp2*mk2*pi16*L2r
          - 64*pow(mp2,2)*L5r*L8r - 64*pow(mp2,2)*L5r*L6r + 56*pow(
         mp2,2)*pow(L5r,2) - 64*pow(mp2,2)*L4r*L8r - 64*pow(mp2,2)*L4r*
         L6r + 112*pow(mp2,2)*L4r*L5r + 56*pow(mp2,2)*pow(L4r,2) - 28./
         27.*pow(mp2,2)*pi16*L3r - 37./9.*pow(mp2,2)*pi16*L2r - 2*pow(
         mp2,2)*pi16*L1r;

  return F6L11/pow(mass.getfpi(),4);
}

double fpi6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fpi6C\n";}

  double   F6C11 =
       + 32*pow(mk2,2)*CC[16] - 32*mp2*mk2*CC[16] + 16*mp2*mk2*
         CC[15] + 8*pow(mp2,2)*CC[17] + 24*pow(mp2,2)*CC[16] + 8*
         pow(mp2,2)*CC[15] + 8*pow(mp2,2)*CC[14];

  return F6C11/pow(mass.getfpi(),4);
}

double fpi6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   F6R11 =
       + pow(Ab(mp2,mu2),2) * ( 87./32. + 1./4.*pow(mp2,-1)*mk2 );

      F6R11 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 2 );

      F6R11 +=  + Ab(mp2,mu2) * (  - 1./2.*mk2*pi16 - 3./4.*mp2*pi16
          );

      F6R11 +=  + pow(Ab(mk2,mu2),2) * ( 5./4. - 1./8.*mp2*pow(
         mk2,-1) );

      F6R11 +=  + Ab(mk2,mu2) * (  - 1./2.*mk2*pi16 - 1./8.*mp2*pi16
          );

      F6R11 +=  + pow(Ab(me2,mu2),2) * ( 3./32. );

      F6R11 +=  + hh(mp2,mp2,mp2,mp2,mu2) * (  - 1./2.*mp2 );

      F6R11 +=  + hh(mp2,mk2,mk2,mp2,mu2) * ( 1./16.*mp2 );

      F6R11 +=  + hh(mk2,mp2,mk2,mp2,mu2) * (  - 1./2.*mk2 );

      F6R11 +=  + hh(me2,mk2,mk2,mp2,mu2) * (  - 1./4.*mk2 + 1./16.*mp2
          );

      F6R11 +=  + hhd(mp2,mp2,mp2,mp2,mu2) * ( 5./12.*pow(mp2,2) );

      F6R11 +=  + hhd(mp2,mk2,mk2,mp2,mu2) * (  - 5./16.*pow(mp2,2) );

      F6R11 +=  + hhd(mp2,me2,me2,mp2,mu2) * ( 1./36.*pow(mp2,2) );

      F6R11 +=  + hhd(mk2,mp2,mk2,mp2,mu2) * ( 1./2.*mp2*mk2 );

      F6R11 +=  + hhd(mk2,mk2,me2,mp2,mu2) * (  - 5./12.*pow(mp2,2) );

      F6R11 +=  + hhd(me2,mk2,mk2,mp2,mu2) * ( 1./4.*mp2*mk2 - 1./16.*
         pow(mp2,2) );

      F6R11 +=  + hh1d(mp2,mk2,mk2,mp2,mu2) * ( 1./2.*pow(mp2,2) );

      F6R11 +=  + hh1d(mk2,mk2,me2,mp2,mu2) * ( pow(mp2,2) );

      F6R11 +=  + hh21d(mp2,mp2,mp2,mp2,mu2) * ( 3./2.*pow(mp2,2) );

      F6R11 +=  + hh21d(mp2,mk2,mk2,mp2,mu2) * (  - 3./16.*pow(mp2,2) )
         ;

      F6R11 +=  + hh21d(mk2,mp2,mk2,mp2,mu2) * ( 3./2.*pow(mp2,2) );

      F6R11 +=  + hh21d(me2,mk2,mk2,mp2,mu2) * ( 9./16.*pow(mp2,2) );

      F6R11 +=  + 15./32.*pow(mk2,2)*pow(pi16,2) + 11./72.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 5./32.*mp2*mk2*pow(pi16,2) + 1./144.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 41./128.*pow(mp2,2)*pow(
         pi16,2) + 35./288.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return F6R11/pow(mass.getfpi(),4);
}

// *************************************************************************

double fk4(const physmass mass, const Li Liin){
  return fk4L(mass,Liin)+fk4R(mass);
}

double fk4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);

  double   F4L44 =
       + 4*mk2*L5r + 8*mk2*L4r + 4*mp2*L4r;

  return F4L44/pow(mass.getfpi(),2);
}

double fk4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double   F4R44 =
       + Ab(mp2,mu2) * ( 3./8. );

      F4R44 +=  + Ab(mk2,mu2) * ( 3./4. );

      F4R44 +=  + Ab(me2,mu2) * ( 3./8. );

  return F4R44/pow(mass.getfpi(),2);
}

double fk6(const physmass mass, const Li Liin, const Ci Ciin){
  return fk6R(mass)+fk6L(mass,Liin)+fk6C(mass,Ciin);
}

double fk6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return fk6L(mass,Liin)+fk6C(mass,Ciin);
}

double fk6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk6L\n";}

  double   F6L44 =
       + Ab(mp2,mu2) * ( 19./2.*mk2*L5r + 19*mk2*L4r + 3*mp2*L5r + 
         47./2.*mp2*L4r - 15./2.*mp2*L3r - 6*mp2*L2r - 24*mp2*L1r );

      F6L44 +=  + Ab(mk2,mu2) * ( 7*mk2*L5r + 30*mk2*L4r - 15*mk2*
         L3r - 18*mk2*L2r - 36*mk2*L1r + 6*mp2*L5r + 7*mp2*L4r );

      F6L44 +=  + Ab(me2,mu2) * ( 3./2.*mk2*L5r + 11*mk2*L4r - 14./3.
         *mk2*L3r - 8./3.*mk2*L2r - 32./3.*mk2*L1r + 3*mp2*L5r - 1./2.*
         mp2*L4r + 7./6.*mp2*L3r + 2./3.*mp2*L2r + 8./3.*mp2*L1r );

      F6L44 +=  - 64*pow(mk2,2)*L5r*L8r - 128*pow(mk2,2)*L5r*L6r + 24*
         pow(mk2,2)*pow(L5r,2) - 128*pow(mk2,2)*L4r*L8r - 256*pow(
         mk2,2)*L4r*L6r + 160*pow(mk2,2)*L4r*L5r + 224*pow(mk2,2)*pow(
         L4r,2) - 89./54.*pow(mk2,2)*pi16*L3r - 61./9.*pow(mk2,2)*pi16*
         L2r - 2*pow(mk2,2)*pi16*L1r - 64*mp2*mk2*L5r*L6r + 32*mp2*mk2*
         pow(L5r,2) - 256*mp2*mk2*L4r*L6r + 112*mp2*mk2*L4r*L5r + 224*
         mp2*mk2*pow(L4r,2) + 2./27.*mp2*mk2*pi16*L3r + 8./9.*mp2*mk2*
         pi16*L2r - 64*pow(mp2,2)*L4r*L8r - 64*pow(mp2,2)*L4r*L6r + 64*
         pow(mp2,2)*L4r*L5r + 56*pow(mp2,2)*pow(L4r,2) - 41./54.*pow(
         mp2,2)*pi16*L3r - 28./9.*pow(mp2,2)*pi16*L2r;

      return F6L44/pow(mass.getfpi(),4);
}

double fk6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in fk6C\n";}

  double   F6C44 =
       + 32*pow(mk2,2)*CC[16] + 16*pow(mk2,2)*CC[15] + 16*pow(
         mk2,2)*CC[14] + 16*mp2*mk2*CC[17] - 32*mp2*mk2*CC[16] + 
         8*mp2*mk2*CC[15] - 16*mp2*mk2*CC[14] - 8*pow(mp2,2)*
         CC[17] + 24*pow(mp2,2)*CC[16] + 8*pow(mp2,2)*CC[14];

  return F6C44/pow(mass.getfpi(),4);
}

double fk6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double   F6R44 =
       + pow(Ab(mp2,mu2),2) * ( 129./128. + 3./32.*pow(mp2,-1)*mk2 )
         ;

      F6R44 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 63./32. );

      F6R44 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 57./64. );

      F6R44 +=  + Ab(mp2,mu2) * (  - 15./32.*mk2*pi16 );

      F6R44 +=  + pow(Ab(mk2,mu2),2) * ( 51./32. + 3./16.*mp2*pow(
         mk2,-1) );

      F6R44 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * ( 3./32. );

      F6R44 +=  + Ab(mk2,mu2) * (  - 9./16.*mk2*pi16 - 3./8.*mp2*
         pi16 );

      F6R44 +=  + pow(Ab(me2,mu2),2) * ( 9./32. + 9./128.*mp2*pow(
         me2,-1) );

      F6R44 +=  + Ab(me2,mu2) * (  - 11./32.*mk2*pi16 - 1./8.*mp2*
         pi16 );

      F6R44 +=  + hh(mp2,mp2,mk2,mk2,mu2) * (  - 3./8.*mp2 );

      F6R44 +=  + hh(mk2,mp2,mp2,mk2,mu2) * ( 3./64.*mk2 );

      F6R44 +=  + hh(mk2,mp2,me2,mk2,mu2) * (  - 9./32.*mk2 );

      F6R44 +=  + hh(mk2,mk2,mk2,mk2,mu2) * (  - 3./8.*mk2 );

      F6R44 +=  + hh(mk2,me2,me2,mk2,mu2) * (  - 9./64.*mk2 );

      F6R44 +=  + hhd(mp2,mp2,mk2,mk2,mu2) * ( 3./16.*pow(mk2,2) + 3./8.
         *mp2*mk2 );

      F6R44 +=  + hhd(mp2,mk2,me2,mk2,mu2) * ( 1./8.*pow(mk2,2) );

      F6R44 +=  + hhd(mk2,mp2,mp2,mk2,mu2) * (  - 3./64.*pow(mk2,2) );

      F6R44 +=  + hhd(mk2,mp2,me2,mk2,mu2) * ( 9./32.*pow(mk2,2) );

      F6R44 +=  + hhd(mk2,mk2,mk2,mk2,mu2) * ( 3./8.*pow(mk2,2) );

      F6R44 +=  + hhd(mk2,me2,me2,mk2,mu2) * ( 181./576.*pow(mk2,2) );

      F6R44 +=  + hh1d(mp2,mp2,mk2,mk2,mu2) * (  - 3./4.*pow(mk2,2) );

      F6R44 +=  + hh1d(mk2,mp2,me2,mk2,mu2) * (  - 3./4.*pow(mk2,2) );

      F6R44 +=  + hh1d(mk2,me2,me2,mk2,mu2) * (  - 5./8.*pow(mk2,2) );

      F6R44 +=  + hh21d(mp2,mp2,mk2,mk2,mu2) * ( 9./8.*pow(mk2,2) );

      F6R44 +=  + hh21d(mk2,mp2,mp2,mk2,mu2) * (  - 9./64.*pow(mk2,2) )
         ;

      F6R44 +=  + hh21d(mk2,mp2,me2,mk2,mu2) * ( 27./32.*pow(mk2,2) );

      F6R44 +=  + hh21d(mk2,mk2,mk2,mk2,mu2) * ( 9./8.*pow(mk2,2) );

      F6R44 +=  + hh21d(mk2,me2,me2,mk2,mu2) * ( 27./64.*pow(mk2,2) );

      F6R44 +=  + 197./384.*pow(mk2,2)*pow(pi16,2) + 3./16.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 29./192.*mp2*mk2*pow(pi16,2) + 1./32.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 13./48.*pow(mp2,2)*pow(pi16,2)
          + 1./16.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return F6R44/pow(mass.getfpi(),4);
}

// *************************************************************************

double feta4(const physmass mass, const Li Liin){
  return feta4L(mass,Liin)+feta4R(mass);
}

double feta4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);

  double   F4L88 =
       + 16./3.*mk2*L5r + 8*mk2*L4r - 4./3.*mp2*L5r + 4*mp2*L4r;

  return F4L88/pow(mass.getfpi(),2);
}

double feta4R(const physmass mass){
  double mk2 = pow(mass.getmk(),2);
  double mu2 = pow(mass.getmu(),2);

  double  F4R88 =
       + Ab(mk2,mu2) * ( 3./2. );

  return F4R88/pow(mass.getfpi(),2);
}

double feta6(const physmass mass, const Li Liin, const Ci Ciin){
  return feta6R(mass)+feta6L(mass,Liin)+feta6C(mass,Ciin);
}

double feta6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return feta6L(mass,Liin)+feta6C(mass,Ciin);
}

double feta6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta6L\n";}

  double   F6L88 =
       + Ab(mp2,mu2) * ( 32./3.*mk2*L5r + 16*mk2*L4r - 4./3.*mp2*L5r
          + 22*mp2*L4r - 6*mp2*L3r - 6*mp2*L2r - 24*mp2*L1r );

      F6L88 +=  + Ab(mk2,mu2) * ( 8*mk2*L5r + 36*mk2*L4r - 14*mk2*
         L3r - 8*mk2*L2r - 32*mk2*L1r + 26./3.*mp2*L5r + 10*mp2*L4r );

      F6L88 +=  + Ab(me2,mu2) * ( 16./3.*mk2*L5r + 8*mk2*L4r - 8*mk2
         *L3r - 16*mk2*L2r - 16*mk2*L1r - 4./3.*mp2*L5r - 2*mp2*L4r + 2
         *mp2*L3r + 4*mp2*L2r + 4*mp2*L1r );

      F6L88 +=  - 256./3.*pow(mk2,2)*L5r*L8r - 512./3.*pow(mk2,2)*L5r*
         L6r + 256./9.*pow(mk2,2)*pow(L5r,2) - 128*pow(mk2,2)*L4r*L8r
          - 256*pow(mk2,2)*L4r*L6r + 192*pow(mk2,2)*L4r*L5r + 224*pow(
         mk2,2)*pow(L4r,2) - 19./9.*pow(mk2,2)*pi16*L3r - 68./9.*pow(
         mk2,2)*pi16*L2r - 32./9.*pow(mk2,2)*pi16*L1r - 128./3.*mp2*mk2
         *L5r*L6r + 448./9.*mp2*mk2*pow(L5r,2) - 256*mp2*mk2*L4r*L6r + 
         96*mp2*mk2*L4r*L5r + 224*mp2*mk2*pow(L4r,2) + 8./9.*mp2*mk2*
         pi16*L3r + 16./9.*mp2*mk2*pi16*L2r + 16./9.*mp2*mk2*pi16*L1r
          + 64./3.*pow(mp2,2)*L5r*L8r + 64./3.*pow(mp2,2)*L5r*L6r - 200.
         /9.*pow(mp2,2)*pow(L5r,2) - 64*pow(mp2,2)*L4r*L8r - 64*pow(
         mp2,2)*L4r*L6r + 48*pow(mp2,2)*L4r*L5r + 56*pow(mp2,2)*pow(
         L4r,2) - 10./9.*pow(mp2,2)*pi16*L3r - 29./9.*pow(mp2,2)*pi16*
         L2r - 2./9.*pow(mp2,2)*pi16*L1r;

      return F6L88/pow(mass.getfpi(),4);
}

double feta6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in feta6C\n";}

  double   F6C88 =
       + 64./3.*pow(mk2,2)*CC[18] + 64./3.*pow(mk2,2)*CC[17] + 32*
         pow(mk2,2)*CC[16] + 64./3.*pow(mk2,2)*CC[15] + 64./3.*pow(
         mk2,2)*CC[14] - 128./3.*mp2*mk2*CC[18] - 64./3.*mp2*mk2*
         CC[17] - 32*mp2*mk2*CC[16] + 16./3.*mp2*mk2*CC[15] - 64./
         3.*mp2*mk2*CC[14] + 64./3.*pow(mp2,2)*CC[18] + 8*pow(
         mp2,2)*CC[17] + 24*pow(mp2,2)*CC[16] - 8./3.*pow(mp2,2)*
         CC[15] + 8*pow(mp2,2)*CC[14];

  return F6C88/pow(mass.getfpi(),4);
}

double feta6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   F6R88 =
       + pow(Ab(mp2,mu2),2) * ( 9./32. );

      F6R88 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 3 );

      F6R88 +=  + pow(Ab(mk2,mu2),2) * ( 9./4. + 3./8.*mp2*pow(
         mk2,-1) );

      F6R88 +=  + Ab(mk2,mu2) * (  - 3./2.*mk2*pi16 - 3./8.*mp2*pi16
          );

      F6R88 +=  + pow(Ab(me2,mu2),2) * ( 9./32. );

      F6R88 +=  + hh(mp2,mk2,mk2,me2,mu2) * (  - 9./16.*mp2 );

      F6R88 +=  + hh(me2,mk2,mk2,me2,mu2) * (  - 3./4.*mk2 + 3./16.*mp2
          );

      F6R88 +=  + hhd(mp2,mp2,me2,me2,mu2) * ( 1./12.*pow(mp2,2) );

      F6R88 +=  + hhd(mp2,mk2,mk2,me2,mu2) * ( 3./4.*mp2*mk2 + 1./16.*
         pow(mp2,2) );

      F6R88 +=  + hhd(mk2,mk2,me2,me2,mu2) * (  - 32./9.*pow(mk2,2) + 2
         *mp2*mk2 - 1./4.*pow(mp2,2) );

      F6R88 +=  + hhd(me2,mk2,mk2,me2,mu2) * ( pow(mk2,2) - 1./2.*mp2*
         mk2 + 1./16.*pow(mp2,2) );

      F6R88 +=  + hhd(me2,me2,me2,me2,mu2) * ( 64./243.*pow(mk2,2) - 56.
         /243.*mp2*mk2 + 49./972.*pow(mp2,2) );

      F6R88 +=  + hh1d(mp2,mk2,mk2,me2,mu2) * (  - 2*mp2*mk2 + 1./2.*
         pow(mp2,2) );

      F6R88 +=  + hh1d(mk2,mk2,me2,me2,mu2) * ( 32./3.*pow(mk2,2) - 20./
         3.*mp2*mk2 + pow(mp2,2) );

      F6R88 +=  + hh21d(mp2,mk2,mk2,me2,mu2) * ( 3*pow(mk2,2) - 3./2.*
         mp2*mk2 + 3./16.*pow(mp2,2) );

      F6R88 +=  + hh21d(me2,mk2,mk2,me2,mu2) * ( 3*pow(mk2,2) - 3./2.*
         mp2*mk2 + 3./16.*pow(mp2,2) );

      F6R88 +=  + 49./96.*pow(mk2,2)*pow(pi16,2) + 5./24.*pow(mk2,2)*
         pow(pi16,2)*pow(pi,2) - 11./96.*mp2*mk2*pow(pi16,2) + 1./48.*
         mp2*mk2*pow(pi16,2)*pow(pi,2) + 91./384.*pow(mp2,2)*pow(
         pi16,2) + 5./96.*pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return F6R88/pow(mass.getfpi(),4);
}

// *************************************************************************

double qqup4(const physmass mass, const Li Liin){
  return qqup4L(mass,Liin)+qqup4R(mass);
}

double qqup4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);

  double Q4L1 =
       + 32*mk2*L6r + 4*mp2*H2r + 8*mp2*L8r + 16*mp2*L6r;

  return Q4L1/pow(mass.getfpi(),2);
}

double qqup4R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double  Q4R1 =
       + Ab(mp2,mu2) * ( 3./2. );

      Q4R1 +=  + Ab(mk2,mu2) * ( 1 );

      Q4R1 +=  + Ab(me2,mu2) * ( 1./6. );


  return Q4R1/pow(mass.getfpi(),2);
}

double qqup6(const physmass mass, const Li Liin, const Ci Ciin){
  return qqup6R(mass)+qqup6L(mass,Liin)+qqup6C(mass,Ciin);
}

double qqup6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return qqup6L(mass,Liin)+qqup6C(mass,Ciin);
}

double qqup6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup6L\n";}


  double   Q6L1 =
       + Ab(mp2,mu2) * ( 112*mk2*L6r + 10*mp2*H2r + 68*mp2*L8r + 112
         *mp2*L6r - 12*mp2*L5r - 24*mp2*L4r );

      Q6L1 +=  + Ab(mk2,mu2) * ( 32*mk2*L8r + 128*mk2*L6r - 16*mk2*
         L5r - 32*mk2*L4r + 4*mp2*H2r + 8*mp2*L8r + 32*mp2*L6r + 8*mp2*
         L5r );

      Q6L1 +=  + Ab(me2,mu2) * (  - 8*me2*L4r - 64./3.*mk2*L7r + 112.
         /3.*mk2*L6r - 32./9.*mk2*L5r - 2./3.*mp2*H2r + 4*mp2*L8r + 64./
         3.*mp2*L7r - 16./3.*mp2*L6r + 20./9.*mp2*L5r );

      Q6L1 +=  - 512*pow(mk2,2)*L6r*L8r - 1024*pow(mk2,2)*pow(L6r,2) + 
         256*pow(mk2,2)*L5r*L6r + 1024*pow(mk2,2)*L4r*L6r - 128*mp2*mk2
         *L6r*H2r - 256*mp2*mk2*L6r*L8r - 1024*mp2*mk2*pow(L6r,2) + 256
         *mp2*mk2*L5r*L6r + 128*mp2*mk2*L4r*H2r + 256*mp2*mk2*L4r*L8r
          + 1024*mp2*mk2*L4r*L6r - 64*pow(mp2,2)*L8r*H2r - 128*pow(
         mp2,2)*pow(L8r,2) - 64*pow(mp2,2)*L6r*H2r - 384*pow(mp2,2)*L6r
         *L8r - 256*pow(mp2,2)*pow(L6r,2) + 64*pow(mp2,2)*L5r*H2r + 128
         *pow(mp2,2)*L5r*L8r + 256*pow(mp2,2)*L5r*L6r + 64*pow(mp2,2)*
         L4r*H2r + 128*pow(mp2,2)*L4r*L8r + 256*pow(mp2,2)*L4r*L6r;

      return Q6L1/pow(mass.getfpi(),4);
}

double qqup6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqup6C\n";}


  double   Q6C1 =
       + 192*pow(mk2,2)*CC[21] + 64*pow(mk2,2)*CC[20] + 8*mp2*mk2*
         CC[94] + 192*mp2*mk2*CC[21] - 4*pow(mp2,2)*CC[94] + 48*
         pow(mp2,2)*CC[21] + 80*pow(mp2,2)*CC[20] + 48*pow(mp2,2)*
         CC[19];


  return Q6C1/pow(mass.getfpi(),4);
}

double qqup6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   Q6R1 =
       + pow(Ab(mp2,mu2),2) * ( 21./8. );

      Q6R1 +=  + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 7./2. );

      Q6R1 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 7./12. );

      Q6R1 +=  + pow(Ab(mk2,mu2),2) * ( 1 );

      Q6R1 +=  + Ab(mk2,mu2)*Ab(me2,mu2) * (  - 1./6. );

      Q6R1 +=  + pow(Ab(me2,mu2),2) * ( 1./72. );

      Q6R1 +=  + 5./36.*pow(me2,2)*pow(pi16,2) + 5./216.*pow(me2,2)*
         pow(pi16,2)*pow(pi,2) - 13./36.*mk2*me2*pow(pi16,2) - 13./216.
         *mk2*me2*pow(pi16,2)*pow(pi,2) + 19./81.*pow(mk2,2)*pow(
         pi16,2) + 19./486.*pow(mk2,2)*pow(pi16,2)*pow(pi,2) + 19./72.*
         mp2*me2*pow(pi16,2) + 19./432.*mp2*me2*pow(pi16,2)*pow(pi,2)
          - 113./324.*mp2*mk2*pow(pi16,2) - 113./1944.*mp2*mk2*pow(
         pi16,2)*pow(pi,2) + 47./648.*pow(mp2,2)*pow(pi16,2) + 47./3888.
         *pow(mp2,2)*pow(pi16,2)*pow(pi,2);

      return Q6R1/pow(mass.getfpi(),4);
}

// *************************************************************************

double qqstrange4(const physmass mass, const Li Liin){
  return qqstrange4L(mass,Liin)+qqstrange4R(mass);
}

double qqstrange4L(const physmass mass, const Li Liin){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange4L\n";}
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);

  double   Q4L3 =
       + 8*mk2*H2r + 16*mk2*L8r + 32*mk2*L6r - 4*mp2*H2r - 8*mp2*L8r + 
         16*mp2*L6r;

  return Q4L3/pow(mass.getfpi(),2);
}

double qqstrange4R(const physmass mass){
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);

  double     Q4R3 =
       + Ab(mk2,mu2) * ( 2 );

      Q4R3 +=  + Ab(me2,mu2) * ( 2./3. );

  return Q4R3/pow(mass.getfpi(),2);
}

double qqstrange6(const physmass mass, const Li Liin, const Ci Ciin){
  return qqstrange6R(mass)+qqstrange6L(mass,Liin)+qqstrange6C(mass,Ciin);
}

double qqstrange6LC(const physmass mass, const Li Liin, const Ci Ciin){
  return qqstrange6L(mass,Liin)+qqstrange6C(mass,Ciin);
}

double qqstrange6L(const physmass mass, const Li Liin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu;
  Liin.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange6L\n";}

  double   Q6L3 =
       + Ab(mp2,mu2) * ( 16*mk2*H2r + 32*mk2*L8r + 64*mk2*L6r - 10*
         mp2*H2r - 20*mp2*L8r + 88*mp2*L6r - 24*mp2*L4r );

      Q6L3 +=  + Ab(mk2,mu2) * ( 8*mk2*H2r + 80*mk2*L8r + 160*mk2*
         L6r - 32*mk2*L5r - 32*mk2*L4r - 4*mp2*H2r - 8*mp2*L8r + 48*mp2
         *L6r + 16*mp2*L5r );

      Q6L3 +=  + Ab(me2,mu2) * (  - 8*me2*L4r + 8./3.*mk2*H2r + 48*
         mk2*L8r + 128./3.*mk2*L7r + 160./3.*mk2*L6r - 128./9.*mk2*L5r
          + 2./3.*mp2*H2r - 20*mp2*L8r - 128./3.*mp2*L7r + 8./3.*mp2*
         L6r + 80./9.*mp2*L5r );

      Q6L3 +=  - 128*pow(mk2,2)*L8r*H2r - 256*pow(mk2,2)*pow(L8r,2) - 
         256*pow(mk2,2)*L6r*H2r - 1024*pow(mk2,2)*L6r*L8r - 1024*pow(
         mk2,2)*pow(L6r,2) + 64*pow(mk2,2)*L5r*H2r + 128*pow(mk2,2)*L5r
         *L8r + 256*pow(mk2,2)*L5r*L6r + 256*pow(mk2,2)*L4r*H2r + 512*
         pow(mk2,2)*L4r*L8r + 1024*pow(mk2,2)*L4r*L6r - 1024*mp2*mk2*
         pow(L6r,2) + 64*mp2*mk2*L5r*H2r + 128*mp2*mk2*L5r*L8r + 256*
         mp2*mk2*L5r*L6r + 1024*mp2*mk2*L4r*L6r + 64*pow(mp2,2)*L8r*H2r
          + 128*pow(mp2,2)*pow(L8r,2) + 64*pow(mp2,2)*L6r*H2r - 128*
         pow(mp2,2)*L6r*L8r - 256*pow(mp2,2)*pow(L6r,2) - 64*pow(mp2,2)
         *L5r*H2r - 128*pow(mp2,2)*L5r*L8r + 256*pow(mp2,2)*L5r*L6r - 
         64*pow(mp2,2)*L4r*H2r - 128*pow(mp2,2)*L4r*L8r + 256*pow(
         mp2,2)*L4r*L6r;

      return Q6L3/pow(mass.getfpi(),4);
}

double qqstrange6C(const physmass mass, const Ci Ciin){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double CC[95],mu;
  Ciin.out(CC,mu);
  if (mass.getmu() != mu){
    cout << "ERROR: differing mu's in qqstrange6C\n";}


  double   Q6C3 =
       + 192*pow(mk2,2)*CC[21] + 192*pow(mk2,2)*CC[20] + 192*pow(
         mk2,2)*CC[19] + 192*mp2*mk2*CC[21] - 64*mp2*mk2*CC[20]
          - 192*mp2*mk2*CC[19] + 4*pow(mp2,2)*CC[94] + 48*pow(
         mp2,2)*CC[21] + 16*pow(mp2,2)*CC[20] + 48*pow(mp2,2)*
         CC[19];

  return Q6C3/pow(mass.getfpi(),4);
}

double qqstrange6R(const physmass mass){
  double mp2 = pow(mass.getmpi(),2);
  double mk2 = pow(mass.getmk(),2);
  double me2 = pow(mass.getmeta(),2);
  double mu2 = pow(mass.getmu(),2);


  double   Q6R3 =
       + Ab(mp2,mu2)*Ab(mk2,mu2) * ( 4 );

      Q6R3 +=  + Ab(mp2,mu2)*Ab(me2,mu2) * ( 4./3. );

      Q6R3 +=  + pow(Ab(mk2,mu2),2) * ( 2 );

      Q6R3 +=  + pow(Ab(me2,mu2),2) * ( 2./9. );

      Q6R3 +=  + 17./36.*pow(me2,2)*pow(pi16,2) + 17./216.*pow(me2,2)*
         pow(pi16,2)*pow(pi,2) - 2./3.*mk2*me2*pow(pi16,2) - 1./9.*mk2*
         me2*pow(pi16,2)*pow(pi,2) + 4./81.*pow(mk2,2)*pow(pi16,2) + 2./
         243.*pow(mk2,2)*pow(pi16,2)*pow(pi,2) + 1./36.*mp2*me2*pow(
         pi16,2) + 1./216.*mp2*me2*pow(pi16,2)*pow(pi,2) + 13./81.*mp2*
         mk2*pow(pi16,2) + 13./486.*mp2*mk2*pow(pi16,2)*pow(pi,2) - 7./
         162.*pow(mp2,2)*pow(pi16,2) - 7./972.*pow(mp2,2)*pow(pi16,2)*
         pow(pi,2);

      return Q6R3/pow(mass.getfpi(),4);
}
