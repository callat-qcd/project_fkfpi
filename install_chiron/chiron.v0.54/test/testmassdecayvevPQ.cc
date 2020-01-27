// testmassdecayvevPQ.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"
#include "quenchedsunsetintegrals.h"
#include "massdecayvevPQ.h"

const double pi162 = 1./pow(4.*M_PI,4);

using namespace std;
int main(void){
  setprecisionquenchedsunsetintegrals(1e-11*pi162);
  cout << setprecision(10);
  cout << "precision set : "<< getprecisionquenchedsunsetintegrals() <<'\n';

  double f0 = 0.0877;
  double mu = 0.77;
  ifstream fin("test/LinfKiRandom.dat");
  string temps;
  for(int i=0; i < 5;i++){
    getline(fin,temps);
  }
  Linf lirandom,li0;
  Ki kirandom,ki0;
  fin >> lirandom;
  fin >> kirandom;
  fin.close();
  //lirandom = li0;  kirandom = ki0;
  const vector<double> Bmq11={0.055,0.065};
  quarkmassnf B0mq11(Bmq11,f0,mu);
  cout << lirandom;
  cout << kirandom;
  cout << "MASSES for the various cases, three sea flavours\n";
  cout << B0mq11;
  cout << "mv1s1nf3p4     : " << mv1s1nf3p4(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p6     : " << mv1s1nf3p6(B0mq11,lirandom,kirandom) << '\n';
  cout << "mv1s1nf3p4L    : " << mv1s1nf3p4L(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p4R    : " << mv1s1nf3p4R(B0mq11) << '\n';
  cout << "mv1s1nf3p6L    : " << mv1s1nf3p6L(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p6K    : " << mv1s1nf3p6K(B0mq11,kirandom) << '\n';
  cout << "mv1s1nf3p6R    : " << mv1s1nf3p6R(B0mq11) << '\n';
  const vector<double> Bmq21={0.055,0.085,0.065};
  quarkmassnf B0mq21(Bmq21,f0,mu);
  cout << B0mq21;
  cout << "mv2s1nf3p4     : " << mv2s1nf3p4(B0mq21,lirandom) << '\n';
  cout << "mv2s1nf3p6     : " << mv2s1nf3p6(B0mq21,lirandom,kirandom) << '\n';
  cout << "mv2s1nf3p4L    : " << mv2s1nf3p4L(B0mq21,lirandom) << '\n';
  cout << "mv2s1nf3p4R    : " << mv2s1nf3p4R(B0mq21) << '\n';
  cout << "mv2s1nf3p6L    : " << mv2s1nf3p6L(B0mq21,lirandom) << '\n';
  cout << "mv2s1nf3p6K    : " << mv2s1nf3p6K(B0mq21,kirandom) << '\n';
  cout << "mv2s1nf3p6R    : " << mv2s1nf3p6R(B0mq21) << '\n';
  const vector<double> Bmq12={0.055,0.065,0.145};
  quarkmassnf B0mq12(Bmq12,f0,mu);
  cout << B0mq12;
  cout << "mv1s2nf3p4     : " << mv1s2nf3p4(B0mq12,lirandom) << '\n';
  cout << "mv1s2nf3p6     : " << mv1s2nf3p6(B0mq12,lirandom,kirandom) << '\n';
  cout << "mv1s2nf3p4L    : " << mv1s2nf3p4L(B0mq12,lirandom) << '\n';
  cout << "mv1s2nf3p4R    : " << mv1s2nf3p4R(B0mq12) << '\n';
  cout << "mv1s2nf3p6L    : " << mv1s2nf3p6L(B0mq12,lirandom) << '\n';
  cout << "mv1s2nf3p6K    : " << mv1s2nf3p6K(B0mq12,kirandom) << '\n';
  cout << "mv1s2nf3p6R    : " << mv1s2nf3p6R(B0mq12) << '\n';
  const vector<double> Bmq22={0.055,0.085,0.065,0.145};
  quarkmassnf B0mq22(Bmq22,f0,mu);
  cout << B0mq22;
  cout << "mv2s2nf3p4     : " << mv2s2nf3p4(B0mq22,lirandom) << '\n';
  cout << "mv2s2nf3p6     : " << mv2s2nf3p6(B0mq22,lirandom,kirandom) << '\n';
  cout << "mv2s2nf3p4L    : " << mv2s2nf3p4L(B0mq22,lirandom) << '\n';
  cout << "mv2s2nf3p4R    : " << mv2s2nf3p4R(B0mq22) << '\n';
  cout << "mv2s2nf3p6L    : " << mv2s2nf3p6L(B0mq22,lirandom) << '\n';
  cout << "mv2s2nf3p6K    : " << mv2s2nf3p6K(B0mq22,kirandom) << '\n';
  cout << "mv2s2nf3p6R    : " << mv2s2nf3p6R(B0mq22) << '\n';
  const vector<double> Bmq13={0.055,0.065,0.105,0.145};
  quarkmassnf B0mq13(Bmq13,f0,mu);
  cout << B0mq13;
  cout << "mv1s3nf3p4     : " << mv1s3nf3p4(B0mq13,lirandom) << '\n';
  cout << "mv1s3nf3p6     : " << mv1s3nf3p6(B0mq13,lirandom,kirandom) << '\n';
  cout << "mv1s3nf3p4L    : " << mv1s3nf3p4L(B0mq13,lirandom) << '\n';
  cout << "mv1s3nf3p4R    : " << mv1s3nf3p4R(B0mq13) << '\n';
  cout << "mv1s3nf3p6L    : " << mv1s3nf3p6L(B0mq13,lirandom) << '\n';
  cout << "mv1s3nf3p6K    : " << mv1s3nf3p6K(B0mq13,kirandom) << '\n';
  cout << "mv1s3nf3p6R    : " << mv1s3nf3p6R(B0mq13) << '\n';
  const vector<double> Bmq23={0.055,0.085,0.065,0.105,0.145};
  quarkmassnf B0mq23(Bmq23,f0,mu);
  cout << B0mq23;
  cout << "mv2s3nf3p4     : " << mv2s3nf3p4(B0mq23,lirandom) << '\n';
  cout << "mv2s3nf3p6     : " << mv2s3nf3p6(B0mq23,lirandom,kirandom) << '\n';
  cout << "mv2s3nf3p4L    : " << mv2s3nf3p4L(B0mq23,lirandom) << '\n';
  cout << "mv2s3nf3p4R    : " << mv2s3nf3p4R(B0mq23) << '\n';
  cout << "mv2s3nf3p6L    : " << mv2s3nf3p6L(B0mq23,lirandom) << '\n';
  cout << "mv2s3nf3p6K    : " << mv2s3nf3p6K(B0mq23,kirandom) << '\n';
  cout << "mv2s3nf3p6R    : " << mv2s3nf3p6R(B0mq23) << '\n';

  cout << "DECAY CONSTANTS for the various cases, three sea flavours\n";
  cout << B0mq11;
  cout << "fv1s1nf3p4     : " << fv1s1nf3p4(B0mq11,lirandom) << '\n';
  cout << "fv1s1nf3p6     : " << fv1s1nf3p6(B0mq11,lirandom,kirandom) << '\n';
  cout << "fv1s1nf3p4L    : " << fv1s1nf3p4L(B0mq11,lirandom) << '\n';
  cout << "fv1s1nf3p4R    : " << fv1s1nf3p4R(B0mq11) << '\n';
  cout << "fv1s1nf3p6L    : " << fv1s1nf3p6L(B0mq11,lirandom) << '\n';
  cout << "fv1s1nf3p6K    : " << fv1s1nf3p6K(B0mq11,kirandom) << '\n';
  cout << "fv1s1nf3p6R    : " << fv1s1nf3p6R(B0mq11) << '\n';
  cout << B0mq21;
  cout << "fv2s1nf3p4     : " << fv2s1nf3p4(B0mq21,lirandom) << '\n';
  cout << "fv2s1nf3p6     : " << fv2s1nf3p6(B0mq21,lirandom,kirandom) << '\n';
  cout << "fv2s1nf3p4L    : " << fv2s1nf3p4L(B0mq21,lirandom) << '\n';
  cout << "fv2s1nf3p4R    : " << fv2s1nf3p4R(B0mq21) << '\n';
  cout << "fv2s1nf3p6L    : " << fv2s1nf3p6L(B0mq21,lirandom) << '\n';
  cout << "fv2s1nf3p6K    : " << fv2s1nf3p6K(B0mq21,kirandom) << '\n';
  cout << "fv2s1nf3p6R    : " << fv2s1nf3p6R(B0mq21) << '\n';
  cout << B0mq12;
  cout << "fv1s2nf3p4     : " << fv1s2nf3p4(B0mq12,lirandom) << '\n';
  cout << "fv1s2nf3p6     : " << fv1s2nf3p6(B0mq12,lirandom,kirandom) << '\n';
  cout << "fv1s2nf3p4L    : " << fv1s2nf3p4L(B0mq12,lirandom) << '\n';
  cout << "fv1s2nf3p4R    : " << fv1s2nf3p4R(B0mq12) << '\n';
  cout << "fv1s2nf3p6L    : " << fv1s2nf3p6L(B0mq12,lirandom) << '\n';
  cout << "fv1s2nf3p6K    : " << fv1s2nf3p6K(B0mq12,kirandom) << '\n';
  cout << "fv1s2nf3p6R    : " << fv1s2nf3p6R(B0mq12) << '\n';
  cout << B0mq22;
  cout << "fv2s2nf3p4     : " << fv2s2nf3p4(B0mq22,lirandom) << '\n';
  cout << "fv2s2nf3p6     : " << fv2s2nf3p6(B0mq22,lirandom,kirandom) << '\n';
  cout << "fv2s2nf3p4L    : " << fv2s2nf3p4L(B0mq22,lirandom) << '\n';
  cout << "fv2s2nf3p4R    : " << fv2s2nf3p4R(B0mq22) << '\n';
  cout << "fv2s2nf3p6L    : " << fv2s2nf3p6L(B0mq22,lirandom) << '\n';
  cout << "fv2s2nf3p6K    : " << fv2s2nf3p6K(B0mq22,kirandom) << '\n';
  cout << "fv2s2nf3p6R    : " << fv2s2nf3p6R(B0mq22) << '\n';
  cout << B0mq13;
  cout << "fv1s3nf3p4     : " << fv1s3nf3p4(B0mq13,lirandom) << '\n';
  cout << "fv1s3nf3p6     : " << fv1s3nf3p6(B0mq13,lirandom,kirandom) << '\n';
  cout << "fv1s3nf3p4L    : " << fv1s3nf3p4L(B0mq13,lirandom) << '\n';
  cout << "fv1s3nf3p4R    : " << fv1s3nf3p4R(B0mq13) << '\n';
  cout << "fv1s3nf3p6L    : " << fv1s3nf3p6L(B0mq13,lirandom) << '\n';
  cout << "fv1s3nf3p6K    : " << fv1s3nf3p6K(B0mq13,kirandom) << '\n';
  cout << "fv1s3nf3p6R    : " << fv1s3nf3p6R(B0mq13) << '\n';
  cout << B0mq23;
  cout << "fv2s3nf3p4     : " << fv2s3nf3p4(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p6     : " << fv2s3nf3p6(B0mq23,lirandom,kirandom) << '\n';
  cout << "fv2s3nf3p4L    : " << fv2s3nf3p4L(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p4R    : " << fv2s3nf3p4R(B0mq23) << '\n';
  cout << "fv2s3nf3p6L    : " << fv2s3nf3p6L(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p6K    : " << fv2s3nf3p6K(B0mq23,kirandom) << '\n';
  cout << "fv2s3nf3p6R    : " << fv2s3nf3p6R(B0mq23) << '\n';


  // checking scalechanges
  cout << "same but now with the scale set to 0.95 and the LECs"
       << "changed accordingly\n";
  kirandom.changescale(0.95,lirandom);
  B0mq11.setmu(0.95);
  B0mq23.setmu(0.95);
  cout << B0mq11;
  cout << "mv1s1nf3p4     : " << mv1s1nf3p4(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p6     : " << mv1s1nf3p6(B0mq11,lirandom,kirandom) << '\n';
  cout << "mv1s1nf3p4L    : " << mv1s1nf3p4L(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p4R    : " << mv1s1nf3p4R(B0mq11) << '\n';
  cout << "mv1s1nf3p6L    : " << mv1s1nf3p6L(B0mq11,lirandom) << '\n';
  cout << "mv1s1nf3p6K    : " << mv1s1nf3p6K(B0mq11,kirandom) << '\n';
  cout << "mv1s1nf3p6R    : " << mv1s1nf3p6R(B0mq11) << '\n';
  cout << B0mq23;
  cout << "fv2s3nf3p4     : " << fv2s3nf3p4(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p6     : " << fv2s3nf3p6(B0mq23,lirandom,kirandom) << '\n';
  cout << "fv2s3nf3p4L    : " << fv2s3nf3p4L(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p4R    : " << fv2s3nf3p4R(B0mq23) << '\n';
  cout << "fv2s3nf3p6L    : " << fv2s3nf3p6L(B0mq23,lirandom) << '\n';
  cout << "fv2s3nf3p6K    : " << fv2s3nf3p6K(B0mq23,kirandom) << '\n';
  cout << "fv2s3nf3p6R    : " << fv2s3nf3p6R(B0mq23) << '\n';
  return 0;
}
