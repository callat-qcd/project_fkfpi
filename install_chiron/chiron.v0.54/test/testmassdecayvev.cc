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
#include "massdecayvev.h"

const double pi162 = 1./pow(4.*M_PI,4);

const double mpp = 0.13956995;
const double mpo = 0.1349764;
const double mkp = 0.493677;

const double mko = 0.497672;
const double mka = 0.49453; // average mK with em removed as well
const double meta = 0.54730;

const double mpi = mpo;
const double mk = mka;
//const double metaGMO = sqrt((4.*mk*mk-mpi*mpi)/3.);

const physmass stdmass(mpo,mka,meta,0.0922,0.77);

using namespace std;
int main(void){
  setprecisionsunsetintegrals(1e-9*pi162);
  // for comparing with fortran, otherwise less precision and digits OK
  //setprecisionsunsetintegrals(1e-11*pi162);
  //cout << setprecision(14);
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
  cout << "mpi4       : " <<mpi4(stdmass,liBE14) << '\n';
  cout << "mpi6       : " <<mpi6(stdmass,liBE14,ciBE14) << '\n';
  cout << "mk4        : " <<mk4(stdmass,liBE14) << '\n';
  cout << "mk6        : " <<mk6(stdmass,liBE14,ciBE14) << '\n';
  cout << "meta4      : " <<meta4(stdmass,liBE14) << '\n';
  cout << "meta6      : " <<meta6(stdmass,liBE14,ciBE14) << '\n';
  cout << "fpi4       : " <<fpi4(stdmass,liBE14) << '\n';
  cout << "fpi6       : " <<fpi6(stdmass,liBE14,ciBE14) << '\n';
  cout << "fk4        : " <<fk4(stdmass,liBE14) << '\n';
  cout << "fk6        : " <<fk6(stdmass,liBE14,ciBE14) << '\n';
  cout << "feta4      : " <<feta4(stdmass,liBE14) << '\n';
  cout << "feta6      : " <<feta6(stdmass,liBE14,ciBE14) << '\n';
  cout << "qqup4      : " <<qqup4(stdmass,liBE14) << '\n';
  cout << "qqup6      : " <<qqup6(stdmass,liBE14,ciBE14) << '\n';
  cout << "qqstrange4 : " <<qqstrange4(stdmass,liBE14) << '\n';
  cout << "qqstrange6 : " <<qqstrange6(stdmass,liBE14,ciBE14) << '\n';

  // same outputs but split in Li dependent, Ci dependent, pure loops
  cout << "The Li,Ci and pure loop parts\n";
  cout << "mpi4L       : " <<mpi4L(stdmass,liBE14) << '\n';
  cout << "mpi4R       : " <<mpi4R(stdmass) << '\n';
  cout << "mpi6L       : " <<mpi6L(stdmass,liBE14) << '\n';
  cout << "mpi6C       : " <<mpi6C(stdmass,ciBE14) << '\n';
  cout << "mpi6R       : " <<mpi6R(stdmass) << '\n';
  cout << "mk4L        : " <<mk4L(stdmass,liBE14) << '\n';
  cout << "mk4R        : " <<mk4R(stdmass) << '\n';
  cout << "mk6L        : " <<mk6L(stdmass,liBE14) << '\n';
  cout << "mk6C        : " <<mk6C(stdmass,ciBE14) << '\n';
  cout << "mk6R        : " <<mk6R(stdmass) << '\n';
  cout << "meta4L      : " <<meta4L(stdmass,liBE14) << '\n';
  cout << "meta4R      : " <<meta4R(stdmass) << '\n';
  cout << "meta6L      : " <<meta6L(stdmass,liBE14) << '\n';
  cout << "meta6C      : " <<meta6C(stdmass,ciBE14) << '\n';
  cout << "meta6R      : " <<meta6R(stdmass) << '\n';
  cout << "fpi4L       : " <<fpi4L(stdmass,liBE14) << '\n';
  cout << "fpi4R       : " <<fpi4R(stdmass) << '\n';
  cout << "fpi6L       : " <<fpi6L(stdmass,liBE14) << '\n';
  cout << "fpi6C       : " <<fpi6C(stdmass,ciBE14) << '\n';
  cout << "fpi6R       : " <<fpi6R(stdmass) << '\n';
  cout << "fk4L        : " <<fk4L(stdmass,liBE14) << '\n';
  cout << "fk4R        : " <<fk4R(stdmass) << '\n';
  cout << "fk6L        : " <<fk6L(stdmass,liBE14) << '\n';
  cout << "fk6C        : " <<fk6C(stdmass,ciBE14) << '\n';
  cout << "fk6R        : " <<fk6R(stdmass) << '\n';
  cout << "feta4L      : " <<feta4L(stdmass,liBE14) << '\n';
  cout << "feta4R      : " <<feta4R(stdmass) << '\n';
  cout << "feta6L      : " <<feta6L(stdmass,liBE14) << '\n';
  cout << "feta6C      : " <<feta6C(stdmass,ciBE14) << '\n';
  cout << "feta6R      : " <<feta6R(stdmass) << '\n';
  cout << "qqup4L      : " <<qqup4L(stdmass,liBE14) << '\n';
  cout << "qqup4R      : " <<qqup4R(stdmass) << '\n';
  cout << "qqup6L      : " <<qqup6L(stdmass,liBE14) << '\n';
  cout << "qqup6C      : " <<qqup6C(stdmass,ciBE14) << '\n';
  cout << "qqup6R      : " <<qqup6R(stdmass) << '\n';
  cout << "qqstrange4L : " <<qqstrange4L(stdmass,liBE14) << '\n';
  cout << "qqstrange4R : " <<qqstrange4R(stdmass) << '\n';
  cout << "qqstrange6L : " <<qqstrange6L(stdmass,liBE14) << '\n';
  cout << "qqstrange6C : " <<qqstrange6C(stdmass,ciBE14) << '\n';
  cout << "qqstrange6R : " <<qqstrange6R(stdmass) << '\n';
  return 0;
}
