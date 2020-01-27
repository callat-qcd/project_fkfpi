// testmassdecayvev.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include "sunsetintegrals.h"
#include "massdecayvevlo.h"

const double pi162 = 1./pow(4.*M_PI,4);

const double mp = 0.135;
const double mk = 0.450;
//const double metaGMO = sqrt((4.*mk*mk-mpi*mpi)/3.);

const lomass stdmass(mp,mk,0.09,0.77);

using namespace std;
int main(void){
  setprecisionsunsetintegrals(1e-9*pi162);
  // for comparing with fortran, otherwise less precision and digits OK
  setprecisionsunsetintegrals(1e-11*pi162);
  cout << setprecision(14);
  ifstream fin("test/LiCiBE14.dat");
  string temps;
  for(int i=0; i < 5;i++){
    getline(fin,temps);
  }
  Li liBE14;
  Ci ciBE14;
  fin >> liBE14;
  fin >> ciBE14;
  fin.close();
  cout << "Overall outputs:\n";
  cout << "mpi4       : " <<mpi4lo(stdmass,liBE14) << '\n';
  cout << "mpi6       : " <<mpi6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "mk4        : " <<mk4lo(stdmass,liBE14) << '\n';
  cout << "mk6        : " <<mk6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "meta4      : " <<meta4lo(stdmass,liBE14) << '\n';
  cout << "meta6      : " <<meta6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "fpi4       : " <<fpi4lo(stdmass,liBE14) << '\n';
  cout << "fpi6       : " <<fpi6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "fk4        : " <<fk4lo(stdmass,liBE14) << '\n';
  cout << "fk6        : " <<fk6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "feta4      : " <<feta4lo(stdmass,liBE14) << '\n';
  cout << "feta6      : " <<feta6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "qqup4      : " <<qqup4lo(stdmass,liBE14) << '\n';
  cout << "qqup6      : " <<qqup6lo(stdmass,liBE14,ciBE14) << '\n';
  cout << "qqstrange4 : " <<qqstrange4lo(stdmass,liBE14) << '\n';
  cout << "qqstrange6 : " <<qqstrange6lo(stdmass,liBE14,ciBE14) << '\n';

  // same outputs but split in Li dependent, Ci dependent, pure loops
  cout << "The Li,Ci and pure loop parts\n";
  cout << "mpi4L       : " <<mpi4Llo(stdmass,liBE14) << '\n';
  cout << "mpi4R       : " <<mpi4Rlo(stdmass) << '\n';
  cout << "mpi6L       : " <<mpi6Llo(stdmass,liBE14) << '\n';
  cout << "mpi6C       : " <<mpi6Clo(stdmass,ciBE14) << '\n';
  cout << "mpi6R       : " <<mpi6Rlo(stdmass) << '\n';
  cout << "mk4L        : " <<mk4Llo(stdmass,liBE14) << '\n';
  cout << "mk4R        : " <<mk4Rlo(stdmass) << '\n';
  cout << "mk6L        : " <<mk6Llo(stdmass,liBE14) << '\n';
  cout << "mk6C        : " <<mk6Clo(stdmass,ciBE14) << '\n';
  cout << "mk6R        : " <<mk6Rlo(stdmass) << '\n';
  cout << "meta4L      : " <<meta4Llo(stdmass,liBE14) << '\n';
  cout << "meta4R      : " <<meta4Rlo(stdmass) << '\n';
  cout << "meta6L      : " <<meta6Llo(stdmass,liBE14) << '\n';
  cout << "meta6C      : " <<meta6Clo(stdmass,ciBE14) << '\n';
  cout << "meta6R      : " <<meta6Rlo(stdmass) << '\n';
  cout << "fpi4L       : " <<fpi4Llo(stdmass,liBE14) << '\n';
  cout << "fpi4R       : " <<fpi4Rlo(stdmass) << '\n';
  cout << "fpi6L       : " <<fpi6Llo(stdmass,liBE14) << '\n';
  cout << "fpi6C       : " <<fpi6Clo(stdmass,ciBE14) << '\n';
  cout << "fpi6R       : " <<fpi6Rlo(stdmass) << '\n';
  cout << "fk4L        : " <<fk4Llo(stdmass,liBE14) << '\n';
  cout << "fk4R        : " <<fk4Rlo(stdmass) << '\n';
  cout << "fk6L        : " <<fk6Llo(stdmass,liBE14) << '\n';
  cout << "fk6C        : " <<fk6Clo(stdmass,ciBE14) << '\n';
  cout << "fk6R        : " <<fk6Rlo(stdmass) << '\n';
  cout << "feta4L      : " <<feta4Llo(stdmass,liBE14) << '\n';
  cout << "feta4R      : " <<feta4Rlo(stdmass) << '\n';
  cout << "feta6L      : " <<feta6Llo(stdmass,liBE14) << '\n';
  cout << "feta6C      : " <<feta6Clo(stdmass,ciBE14) << '\n';
  cout << "feta6R      : " <<feta6Rlo(stdmass) << '\n';
  cout << "qqup4L      : " <<qqup4Llo(stdmass,liBE14) << '\n';
  cout << "qqup4R      : " <<qqup4Rlo(stdmass) << '\n';
  cout << "qqup6L      : " <<qqup6Llo(stdmass,liBE14) << '\n';
  cout << "qqup6C      : " <<qqup6Clo(stdmass,ciBE14) << '\n';
  cout << "qqup6R      : " <<qqup6Rlo(stdmass) << '\n';
  cout << "qqstrange4L : " <<qqstrange4Llo(stdmass,liBE14) << '\n';
  cout << "qqstrange4R : " <<qqstrange4Rlo(stdmass) << '\n';
  cout << "qqstrange6L : " <<qqstrange6Llo(stdmass,liBE14) << '\n';
  cout << "qqstrange6C : " <<qqstrange6Clo(stdmass,ciBE14) << '\n';
  cout << "qqstrange6R : " <<qqstrange6Rlo(stdmass) << '\n';
  return 0;
}
