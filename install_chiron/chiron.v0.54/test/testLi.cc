// testLi.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// test the Li class implementations

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>  // for srand


using namespace std;
#include "Li.h"

int main(void){
  Li liset1;
  cout << "A default set if Li inputs:\n";
  cout << liset1;
  // an example with setting all numerical values and functions with the scale
  // and name
  Li liset2(1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,10e-3,11e-3,12e-3);
  liset2.setmu(0.9);
  liset2.setname("Set 2");
  cout << "Another simple example:\n";
  cout << liset2;
  // an example setting the name and some individual Li
  Li liset3;
  liset3.setli(3,6e-3);
  liset3.setli(5e-3,4);
  liset3.setname("set 3");
  cout << "The output of some inputs\n";
  cout << liset3;
  // summing some, note the first should complain
  Li liset4 = liset2+liset3;
  liset2.setmu(0.77);
  Li liset5 = liset2+liset3;
  liset5.setname("sum");
  cout << liset5;
  Li liset5b = liset2-liset3;
  liset5b.setname("difference");
  cout << liset5b;
  // multiplying a set with a number from both sides
  Li liset6 = 2.*liset3;
  Li liset7 = liset3*2.;
  liset6.setname("2 times set 3");
  liset7.setname("set 3 times 2");
  cout << liset6;
  cout << liset7;
  // changin the scale
  Li liset8 = liset3;
  Li liset9 = liset3;
  liset8.changescale(0.5);
  liset9.changescale(1.0);
  cout << liset8;
  cout << liset9;
  // writing to a stream and reading it back in, here with a filestream
  ofstream ouf("temp.dat");
  ouf << liset9;
  ouf.close();
  ifstream inf("temp.dat");
  Li liset10;
  inf >> liset10;
  inf.close();
  cout << "Written to and read back from file:\n";
  cout << liset10;
  // testing the random sets, some examples
  srand(time (0));// initializes the random number generator
  cout << "Output for some random Li:\n"
       << "NOTE: these can (should) differ from the outputs in testLi.dat\n";
  cout << Lirandom();
  cout << LirandomlargeNc();
  cout << LirandomlargeNc2();
  return 0;
}
