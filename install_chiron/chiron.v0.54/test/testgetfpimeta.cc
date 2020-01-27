// testgetfpimeta.cc is part of the CHIRON ChPT at two loops program
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
#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "getfpimeta.h"
#include "sunsetintegrals.h"

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
  setprecisionsunsetintegrals(1e-7*pi162);
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
  double mpiin[7]={mpo,0.1,0.3,0.1,0.1,0.3,0.495};
  double mkin[7]={mka,0.48714,0.5496,0.4,0.495,0.495,0.495};
  cout <<"#  mpi         mk        meta(p^6)  fpi(p^6)  meta(p^4)  fpi(p^4)\n";
  for(int i=0; i< 7;i++){
    physmass mass6 = getfpimeta6(mpiin[i],mkin[i],stdmass,liBE14,ciBE14);
    physmass mass4 = getfpimeta4(mpiin[i],mkin[i],stdmass,liBE14);
    cout<< setw(10)<< mpiin[i] <<' '<< setw(10)<<mkin[i]
	<<' '<< setw(10)<<mass6.getmeta()<<' '<< setw(10)<<mass6.getfpi()
	<<' '<< setw(10)<<mass4.getmeta()<<' '<< setw(10)<<mass4.getfpi()
	<< '\n';
  }

  return 0;
}
