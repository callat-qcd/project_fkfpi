// testfinitevolumeoneloopintegrals.cc is part of 
// the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#include<iostream>
#include<iomanip>
#include<cmath>
#include "finitevolumeoneloopintegrals.h"

const double pi16 = 1./pow(4.*M_PI,2);

int main(void){
  using namespace std;
  double mpi = 0.135;
  double xl = 2./mpi;
  double msq = mpi*mpi;
  cout <<"A  number of results for mpi L=2, theta and Bessel\n";
  cout << "AbVacc maxsum AbVt AbVb A22bVt A22bVb A23bVt A23bVb\n"; 
  double epst[8] = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};
  int nb[8] = {1,2,5,10,20,40,60,100};
  for(int i=0; i<8; i++){
    setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
    setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
    cout << epst[i]<<' '<<nb[i]<<' '
	 <<AbVt(msq,xl) <<' '<<AbVb(msq,xl)<<' '
	 <<A22bVt(msq,xl) <<' '<<A22bVb(msq,xl)<<' '
	 <<A23bVt(msq,xl) <<' '<<A23bVb(msq,xl)<<'\n';
  }
  cout << "AbVacc maxsum BbVt BbVb B22bVt B22bVb B23bVt B23bVb\n"; 
  for(int i=0; i<8; i++){
    setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
    setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
    cout << epst[i]<<' '<<nb[i]<<' '
	 <<BbVt(msq,xl) <<' '<<BbVb(msq,xl)<<' '
	 <<B22bVt(msq,xl) <<' '<<B22bVb(msq,xl)<<' '
	 <<B23bVt(msq,xl) <<' '<<B23bVb(msq,xl)<<'\n';
  }
  cout << "AbVacc maxsum CbVt CbVb C22bVt C22bVb C23bVt C23bVb\n"; 
  for(int i=0; i<8; i++){
    setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
    setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
    cout << epst[i]<<' '<<nb[i]<<' '
	 <<CbVt(msq,xl) <<' '<<CbVb(msq,xl)<<' '
	 <<C22bVt(msq,xl) <<' '<<C22bVb(msq,xl)<<' '
	 <<C23bVt(msq,xl) <<' '<<C23bVb(msq,xl)<<'\n';
  }
  cout << "AbVacc maxsum DbVt DbVb D22bVt D22bVb D23bVt D23bVb\n"; 
  for(int i=0; i<8; i++){
    setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
    setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
    cout << epst[i]<<' '<<nb[i]<<' '
	 <<DbVt(msq,xl) <<' '<<DbVb(msq,xl)<<' '
	 <<D22bVt(msq,xl) <<' '<<D22bVb(msq,xl)<<' '
	 <<D23bVt(msq,xl) <<' '<<D23bVb(msq,xl)<<'\n';
  }

  cout <<"checking the relation 4*A22bV+3*A23bV=msq*AbV\n";
  cout << 4.*A22bVt(msq,xl)+3.*A23bVt(msq,xl) << ' '<<msq*AbVt(msq,xl)<<'\n';
  cout << 4.*A22bVb(msq,xl)+3.*A23bVb(msq,xl) << ' '<<msq*AbVb(msq,xl)<<'\n';
  cout <<"checking the relation 4*B22bV+3*B23bV=AbV+msq*BbV\n";
  cout << 4.*B22bVt(msq,xl)+3.*B23bVt(msq,xl) << ' '<<AbVt(msq,xl)+msq*BbVt(msq,xl)<<'\n';
  cout << 4.*B22bVb(msq,xl)+3.*B23bVb(msq,xl) << ' '<<AbVb(msq,xl)+msq*BbVb(msq,xl)<<'\n';
  cout <<"checking the relation 4*C22bV+3*C23bV=BbV+msq*CbV\n";
  cout << 4.*C22bVt(msq,xl)+3.*C23bVt(msq,xl) << ' '<<BbVt(msq,xl)+msq*CbVt(msq,xl)<<'\n';
  cout << 4.*C22bVb(msq,xl)+3.*C23bVb(msq,xl) << ' '<<BbVb(msq,xl)+msq*CbVb(msq,xl)<<'\n';
  cout <<"checking the relation 4*D22bV+3*D23bV=CbV+msq*DbV\n";
  cout << 4.*D22bVt(msq,xl)+3.*D23bVt(msq,xl) << ' '<<CbVt(msq,xl)+msq*DbVt(msq,xl)<<'\n';
  cout << 4.*D22bVb(msq,xl)+3.*D23bVb(msq,xl) << ' '<<CbVb(msq,xl)+msq*DbVb(msq,xl)<<'\n';
  cout <<"checking B22bV=1/2 AbV\n";
  cout << B22bVb(msq,xl) <<' '<<0.5*AbVb(msq,xl)<<'\n';
  cout <<"checking C22bV=1/4 BbV\n";
  cout << C22bVb(msq,xl) <<' '<<0.25*BbVb(msq,xl)<<'\n';
  cout <<"checking D22bV=1/6 AbV\n";
  cout << D22bVb(msq,xl) <<' '<<1./6.*CbVb(msq,xl)<<'\n';
  cout <<"checking B23bV=1/3 (m^2 BbV-AbV)\n";
  cout << B23bVt(msq,xl) <<' '<<1./3.*(msq*BbVt(msq,xl)-AbVt(msq,xl))<<'\n';
  cout <<"checking C23bV=1/3 (m^2 CbV)\n";
  cout << C23bVt(msq,xl) <<' '<<1./3.*(msq*CbVt(msq,xl))<<'\n';
  cout <<"checking D23bV=1/3 (m^2 DbV+1/3 CbV)\n";
  cout << D23bVt(msq,xl) <<' '<<1./3.*(msq*DbVt(msq,xl)+CbVt(msq,xl)/3.)<<'\n';
  cout <<"checking A23bV=1/3 (m^2 AbV)\n";
  cout << A23bVt(msq,xl) <<' '<<1./3.*(msq*AbVt(msq,xl))<<'\n';

  return 0;
}
