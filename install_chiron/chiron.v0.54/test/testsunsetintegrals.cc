// testsunsetintegrals.cc is part of 
// the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// prints out the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]


// an example program to print out some values of the sunsetintegrals

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
using namespace std;

typedef std::complex<double> dcomplex;
#include "sunsetintegrals.h"


void printall(double m1sq,double m2sq,double m3sq,double qsq,double mu2){
  cout << "output for m1sq,m2sq,m3sq,qsq,musq: ";
  cout <<  m1sq<<' '<<m2sq<<' '<<m3sq<<' '<<qsq<<' '<<mu2<<endl;
  cout << "hh    "<< hh(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "zhh   "<< zhh(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hh1   "<< hh1(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "zhh1  "<< zhh1(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hh21  "<< hh21(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "zhh21 "<< zhh21(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hh31  "<< hh31(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "zhh31 "<< zhh31(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hhd   "<< hhd(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hh1d  "<< hh1d(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
  cout << "hh21d "<< hh21d(m1sq,m2sq,m3sq,qsq,mu2)<<endl;
}

int main(void){
  const double  pi = M_PI;
  const double pi16 = 1./(16.*pi*pi);
  const double pi162 = pi16*pi16;
  cout << "note that above threshold only the zh... should be used\n"
       << "below theshold the others are usually faster\n\n";
  setprecisionsunsetintegrals(1e-11*pi162);
  double  mpi = 0.1395;
  double  mk = 0.495;
  double  meta = 0.548;
  double  mu = 0.77;

  double mpi2 = mpi*mpi;
  double mk2 = mk*mk;
  double me2 = meta*meta;
  double mu2 = mu*mu;
  cout << "precision sunsetintegrals: " << getprecisionsunsetintegrals()
       <<"\n\n";
  
  cout.precision(10);

  cout <<"##### mpi2,mpi2,mpi2,mpi2 ###################################\n";
       printall(mpi2,mpi2,mpi2,mpi2,mu2);
  cout <<"##### mpi2,mk2,mk2,mpi2 ###################################\n";
       printall(mpi2,mk2,mk2,mpi2,mu2);
  cout <<"##### me2,mk2,mk2,mpi2 ###################################\n";
       printall(me2,mk2,mk2,mpi2,mu2);
  cout <<"##### mpi2,mpi2,mk2,mk2 ###################################\n";
       printall(mpi2,mpi2,mk2,mk2,mu2);
  cout <<"##### me2,mpi2,mk2,mk2 ###################################\n";
       printall(me2,mpi2,mk2,mk2,mu2);
  cout <<"##### me2,me2,mk2,mk2 ###################################\n";
       printall(me2,me2,mk2,mk2,mu2);
  cout <<"##### mk2,mk2,mk2,mk2 ###################################\n";
       printall(mk2,mk2,mk2,mk2,mu2);

  return 0;
}
