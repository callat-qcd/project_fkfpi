// testKi.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>  // for srand


using namespace std;
#include "Ki.h"

int main(void){
  Linf linffit10(0.1e-3,0.43173E-03,
	    0.73457E-03,
	    -0.23469E-02,
	    0.0000,
	    0.97205E-03,
	    0.0000,
	    -0.30833E-03,
	    0.60467E-03,
		 0.69000E-02,0.,0.1e-4,0.,0.,0.77,"almost fit 10",5);
  Ki Kiset1(0.77,"just a try",5);
  Kiset1.setmu(0.77);
  cout << linffit10;
  cout << Kiset1;
  Kiset1.changescale(linffit10,1.2);
  cout << "changed scale also for Li";
  cout << linffit10;
  cout << Kiset1;
  Linf linffit10scale = linffit10;
  Ki Kiset2 = Kiset1;
  cout << "Setting the scale back\n";
  Kiset2.changescale(linffit10scale,0.77);
  cout << linffit10scale;
  cout << Kiset2;
  ofstream fout("temp.dat");
  fout << Kiset1;
  fout.close();
  ifstream fin("temp.dat");
  Ki Kiset3;
  fin >> Kiset3;
  fin.close();
  cout << "now written to and read back in from a file\n";
  cout << Kiset3;
  // examples of the out functions
  double Kr[116];
  // sum and difference
  Kiset3.setmu(0.77);
  cout << "sum\n";
  cout << Kiset3+Kiset3;
  cout << "difference\n";
  cout << Kiset3-Kiset2;
  cout << "multiplying by a number from either side\n";
  cout << 3.*Kiset3;
  cout << Kiset3*5.;
  cout << "output from the random functions\n";
  srand(time (0));// initializes the random number generator
  cout << Kirandom();

  return 0;


}
