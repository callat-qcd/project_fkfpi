// getfpimeta.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// getting consistent valeus of masses and decay constants
#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "massdecayvev.h"
#include "getfpimeta.h"

// mpi and mk as input p^6 calculation
physmass getfpimeta6(const double mpiin, const double mkin,
		const physmass massin, const Li liin, const Ci ciin){
  //calculate F0 from the standard input
  double mu = massin.getmu();
  double F0 = massin.getfpi()/(1.+fpi4(massin,liin)+fpi6(massin,liin,ciin));
  double metain = sqrt(4./3.*mkin*mkin-1./3.*mpiin*mpiin);
  physmass mass0(mpiin,mkin,metain,massin.getfpi(),mu);
  double mp2lo = pow(mpiin,2)-mpi4(mass0,liin)-mpi6(mass0,liin,ciin);
  double mk2lo = pow(mkin,2)-mk4(mass0,liin)-mk6(mass0,liin,ciin);
  double me2lo = 4./3.*mk2lo-1./3.*mp2lo;
  double me0 = sqrt(me2lo+meta4(mass0,liin)+meta6(mass0,liin,ciin));
  double fpi0 = F0*(1.+fpi4(mass0,liin)+fpi6(mass0,liin,ciin));
  physmass mass1;
  double me1,fpi1;
  for(int i=1;i<=1000;i++){
    mass1 = physmass(mpiin,mkin,me0,fpi0,mu);
    mp2lo = pow(mpiin,2)-mpi4(mass1,liin)-mpi6(mass1,liin,ciin);
    mk2lo = pow(mkin,2)-mk4(mass1,liin)-mk6(mass1,liin,ciin);
    me2lo = 4./3.*mk2lo-1./3.*mp2lo;
    me1 = sqrt(me2lo+meta4(mass1,liin)+meta6(mass1,liin,ciin));
    double fpi40 = fpi4(mass1,liin);
    double fpi60 = fpi6(mass1,liin,ciin);
    fpi1 = F0*(1.+fpi40+fpi60);
    if((fabs((me1-me0)/me1)+fabs((fpi1-fpi0)/fpi1)) < 1e-6){
      return physmass(mpiin,mkin,me1,fpi1,mu);
    }
    me0 = me1; fpi0=fpi1;
  }
  std::cout << "precision in getfpimeta6 not reached in 1000 steps\n";
  return physmass(mpiin,mkin,me1,fpi1,mu);
}
// mpi and mk as input p^4 calculation
physmass getfpimeta4(const double mpiin, const double mkin,
		     const physmass massin, const Li liin){
  //calculate F0 from the standard input
  double mu = massin.getmu();
  double F0 = massin.getfpi()/(1.+fpi4(massin,liin));
  double metain = sqrt(4./3.*mkin*mkin-1./3.*mpiin*mpiin);
  physmass mass0(mpiin,mkin,metain,massin.getfpi(),mu);
  double mp2lo = pow(mpiin,2)-mpi4(mass0,liin);
  double mk2lo = pow(mkin,2)-mk4(mass0,liin);
  double me2lo = 4./3.*mk2lo-1./3.*mp2lo;
  double me0 = sqrt(me2lo+meta4(mass0,liin));
  double fpi0 = F0*(1.+fpi4(mass0,liin));
  physmass mass1;
  double me1,fpi1;
  for(int i=1;i<=1000;i++){
    mass1 = physmass(mpiin,mkin,me0,fpi0,mu);
    mp2lo = pow(mpiin,2)-mpi4(mass1,liin);
    mk2lo = pow(mkin,2)-mk4(mass1,liin);
    me2lo = 4./3.*mk2lo-1./3.*mp2lo;
    me1 = sqrt(me2lo+meta4(mass1,liin));
    double fpi40 = fpi4(mass1,liin);
    fpi1 = F0*(1.+fpi40);
    if((fabs((me1-me0)/me1)+fabs((fpi1-fpi0)/fpi1)) < 1e-6){
      return physmass(mpiin,mkin,me1,fpi1,mu);
    }
    me0 = me1; fpi0=fpi1;
  }
  std::cout << "precision in getfpimeta4 not reached in 1000 steps\n";
  return physmass(mpiin,mkin,me1,fpi1,mu);
}
