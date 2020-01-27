// testoneloopintegrals.cc is part of 
// the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// prints out the one-loop functions and some tests of them

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
using namespace std;

typedef std::complex<double> dcomplex;
#include "oneloopintegrals.h"

const double pi16 = 1./pow(4.*M_PI,2);

int main(void){
  setprecisiononeloopintegrals(1e-12*pi16);
  double  mpi = 0.1395;
  double  mk = 0.495;
  double  meta = 0.548;
  double  mu = 0.77;

  double mpi2 = mpi*mpi;
  double mk2 = mk*mk;
  double me2 = meta*meta;
  double mu2 = mu*mu;
  cout << "\nprecision oneloopintegrals: " << getprecisiononeloopintegrals()
       <<"\n";
  
  cout << "\ntadpoles\n"
       << Ab(mpi2,mu2) <<' '<<Ab(mk2,mu2)<<'\n'
       << Bb(mpi2,mu2) <<' '<<Bb(mk2,mu2)<<'\n'
       << Cb(mpi2,mu2) <<' '<<Cb(mk2,mu2)<<'\n';
  cout << "  different interface\n"
       << Ab(1,mpi2,mu2) <<' '<<Ab(1,mk2,mu2)<<'\n'
       << Ab(2,mpi2,mu2) <<' '<<Ab(2,mk2,mu2)<<'\n'
       << Ab(3,mpi2,mu2) <<' '<<Ab(3,mk2,mu2)<<'\n';

  cout << "\nBubble integrals Bb\n"
       << "  different masses,analytical and numerical\n"
       << Bb(mpi2,mk2,-me2,mu2) <<' '<<Bbnum(mpi2,mk2,-me2,mu2)<<'\n'
       << Bb(mpi2,mk2,2.*mk2,mu2) <<' '<<Bbnum(mpi2,mk2,2.*mk2,mu2)<<'\n';
  cout << "  equal masses two analytical, one numerical\n"
       << Bb(mpi2,-me2,mu2) <<' '<<Bb(mpi2,mpi2,-me2,mu2)
       <<' '<<Bbnum(mpi2,mk2,-me2,mu2)<<'\n'
       << Bb(mpi2,me2,mu2) <<' '<<Bb(mpi2,mpi2,me2,mu2)
       <<' '<<Bbnum(mpi2,mpi2,me2,mu2)<<'\n';
  cout << "\nBubble integrals B1b\n"
       << B1b(mpi2,mk2,-me2,mu2) <<' '<<B1bnum(mpi2,mk2,-me2,mu2)<<'\n'
       << B1b(mpi2,mk2,1.5*me2,mu2) <<' '<<B1bnum(mpi2,mk2,1.5*me2,mu2)<<'\n'
       << B1b(mpi2,mk2,0.,mu2) <<' '<<B1bnum(mpi2,mk2,0.,mu2)
       << (B1b(mpi2,mk2,1e-4,mu2)+B1b(mpi2,mk2,-1e-4,mu2))/2.
       <<'\n';
  cout << "\nBubble integrals B21b\n"
       << B21b(mpi2,mk2,-me2,mu2) <<' '<<B21bnum(mpi2,mk2,-me2,mu2)<<'\n'
       << B21b(mpi2,mk2,1.5*me2,mu2) <<' '<<B21bnum(mpi2,mk2,1.5*me2,mu2)<<'\n';
  cout << "\nBubble integrals B22b\n"
       << B22b(mpi2,mk2,-me2,mu2) <<' '<<B22bnum(mpi2,mk2,-me2,mu2)<<'\n'
       << B22b(mpi2,mk2,1.5*me2,mu2) <<' '<<B22bnum(mpi2,mk2,1.5*me2,mu2)<<'\n'
       << B22b(mpi2,-me2,mu2) <<' '<<B22bnum(mpi2,mpi2,-me2,mu2)
       << B22b(mpi2,me2,mu2) <<' '<<B22bnum(mpi2,mpi2,me2,mu2)
       <<'\n';
  cout << "\nBubble integrals B21b, B32b only done via integration over x\n"
       << "This tests the consistency relation with B1\n"
       << -me2*B31bnum(mpi2,mk2,-me2,mu2)+6.*B32bnum(mpi2,mk2,-me2,mu2) <<' '
       << mpi2*B1b(mpi2,mk2,-me2,mu2)+Ab(mk2,mu2)+pi16*(mpi2/6.+mk2/3.+me2/12.)
       << '\n'
       <<1.5*me2*B31bnum(mpi2,mk2,1.5*me2,mu2)+6.*B32bnum(mpi2,mk2,1.5*me2,mu2)
       <<' '
       << mpi2*B1b(mpi2,mk2,1.5*me2,mu2)+Ab(mk2,mu2)
            +pi16*(mpi2/6.+mk2/3.-1.5*+me2/12.)
       << '\n';
  return 0;
}
