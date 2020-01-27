// testfinitevolumesunsetintegrals.cc is part of 
// the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include<iomanip>
#include<cmath>
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"

const double pi16 = 1./pow(4.*M_PI,2);
const double pi162 = 1./pow(4.*M_PI,4);

int main(void){
  using namespace std;
  const double mpp = 0.13956995;
  const double mpo = 0.1349764;
  const double mkp = 0.493677;
  const double mko = 0.497672;
  const double mka = 0.49453; // average mK with em removed as well
  const double meta = 0.54730;

  double m1sq = mpo*mpo;
  double m2sq = mka*mka;
  double m3sq = meta*meta;
  double qsq = m1sq;
  double mu2 = pow(0.77,2);
  double xl = 3./mpo;
  double msq = mpo*mpo;
  cout << "A  number of results for mpi L=3, theta and Bessel\n";
  //cout << "AbVacc maxsum AbVt AbVb A22bVt A22bVb A23bVt A23bVb BbVt BbVb\n"; 
  double epst[8] = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};
  int nb[8] =      {1,2,5,10,20,40,60,100};
  int nb2[8] =     {1,2,5,10,20,40,40,40};
  // sunset actually contain some calls to oneloop integrals AbV,A22bV and A23bV
  setprecisionfinitevolumeoneloopt(1e-6*pi162,1e-5,false);
  setprecisionfinitevolumeoneloopb(60,1e-5,false);

  cout <<"standard cases\n";
  cout <<"eps_Theta sumBessel  hhVt      hhVb        hh1Vt          hh1Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hhVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1Vt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1Vb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh21Vt    hh21Vb       hh22Vt         hh22Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh21Vt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21Vb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22Vt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22Vb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh27Vt    hh27Vb       hhdVt         hhdVb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh27Vt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27Vb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhdVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhdVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh1dVt    hh1dVb       hh21dVt       hh21dVb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh1dVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1dVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21dVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21dVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh22dVt   hh22dVb      hh27dVt       hhd27Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh22dVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22dVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27dVt(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27dVb(m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout << "############################################################\n";
  cout << "cases with different propagators\n";

  for(int nprop=1;nprop<9;nprop++){
    cout << "case propagator = "<<nprop<<'\n';

  cout <<"eps_Theta sumBessel  hhVt      hhVb        hh1Vt          hh1Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hhVt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhVb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1Vt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1Vb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh21Vt    hh21Vb       hh22Vt         hh22Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh21Vt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21Vb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22Vt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22Vb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh27Vt    hh27Vb       hhdVt         hhdVb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh27Vt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27Vb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhdVt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hhdVb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh1dVt    hh1dVb       hh21dVt       hh21dVb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh1Vt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh1Vb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21dVt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh21dVb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  cout <<"eps_Theta sumBessel  hh22dVt   hh22dVb      hh27dVt       hhd27Vb\n";
  for(int i=0; i<6; i++){
    setprecisionfinitevolumesunsett(epst[i],epst[i],false);
    setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst[i],false);
    cout << setw(8)<<epst[i]<<' '<<setw(3)<<nb2[i]<<' '
	 <<setw(13)<<hh22dVt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh22dVb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27dVt(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<setw(13)<<hh27dVb(nprop,m1sq,m2sq,m3sq,m1sq,xl,mu2)<<' '
	 <<'\n';
  }
  } // nprop
  return 0;
}
