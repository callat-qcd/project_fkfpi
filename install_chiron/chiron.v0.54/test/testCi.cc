// testCi.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>  // for srand


using namespace std;
#include "Ci.h"

int main(void){
  Li lifit10(0.43173E-03,
	    0.73457E-03,
	    -0.23469E-02,
	    0.0000,
	    0.97205E-03,
	    0.0000,
	    -0.30833E-03,
	    0.60467E-03,
	    0.69000E-02,0.,0.,0.,0.77,"fit 10");
  Ci Ciset1;
  Ciset1.setmu(0.77);
  cout << lifit10;
  Ciset1.changescale(lifit10,1.2);
  cout << lifit10;
  cout << Ciset1;
  Li lifit10scale = lifit10;
  Ci Ciset2 = Ciset1;
  Ciset2.changescale(lifit10scale,1.0);
  cout << "changed scale also for Li";
  cout << lifit10scale;
  cout << Ciset2;
  ofstream fout("temp.dat");
  fout << Ciset2;
  fout.close();
  ifstream fin("temp.dat");
  Ci Ciset3;
  fin >> Ciset3;
  fin.close();
  cout << "now written to and read back in from a file\n";
  cout << Ciset3;
  // examples of the out functions
  double Cr[95];
  const double MV = 0.77;
  const double fV = 0.20;
  const double gV = 0.09;
  const double fChi = -0.025;
  const double aV =  -0.014;
  const double MS = 0.98;
  const double cd = 0.032;
  const double cm = 0.042;
  const double l3ss = 0.;
  const double MEP = 0.958;
  const double dmt = 0.020;
  const double MP = 1.4;
  const double dm = 0.;
  const double fpi = 0.0922;
  const double mureso = 0.77;
  Ci Cireso = Ci(fpi,mureso,MV,fV,gV,fChi,aV,MS,cd,cm,l3ss,MEP,dmt,MP,dm);
  cout << "The Ci from resonances\n";
  cout << Cireso;
  // sum and difference
  Ciset3.setmu(0.77);
  cout << "sum\n";
  cout << Ciset3+Cireso;
  cout << "difference\n";
  cout << Ciset3-Cireso;
  cout << "multiplying by a number from either side\n";
  cout << 3.*Cireso;
  cout << Cireso*5;
  cout << "theoutput from the random functions\n";
  srand(time (0));// initializes the random number generator
  cout << Cirandom();
  cout << CirandomlargeNc();
  cout << CirandomlargeNc2();

  return 0;


}
