// testmassdecayvevPQV.cc is part of the CHIRON ChPT at two loops program
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
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "massdecayvevPQV.h"

const double pi162 = 1./pow(4.*M_PI,4);

using namespace std;
int main(void){

  double f0 = 0.0877;
  double mu = 0.77;
  ifstream fin("test/LinfKiRandom.dat");
  string temps;
  for(int i=0; i < 5;i++){
    getline(fin,temps);
  }
  Linf lirandom,li0;
  fin >> lirandom;
  fin.close();

  // masses of 150, 250,450 MeV for valence mesons (lowest order)
  double m11 = pow(0.150,2), mq1 = m11/2.;
  double m22 = pow(0.250,2), mq2 = m22/2.;
  double m33 = pow(0.450,2), mq3 = m33/2.;
  double mq4 = 1.5*mq1;
  double mq5 = 1.3*mq2;
  double mq6 = 1.1*mq3;

  quarkmassnf B0mq11({mq1,mq4},f0,mu);
  quarkmassnf B0mq21({mq1,mq3,mq4},f0,mu);
  quarkmassnf B0mq12({mq1,mq4,mq6},f0,mu);
  quarkmassnf B0mq22({mq1,mq3,mq4,mq6},f0,mu);
  quarkmassnf B0mq13({mq1,mq4,mq5,mq6},f0,mu);
  quarkmassnf B0mq23({mq1,mq3,mq4,mq5,mq6},f0,mu);

  // setting all needed precisions
  double epst[9]  = {10.,1.,1e-1,1e-2,1e-3,1e-4,1e-5,1e-5,1e-6};
  double epst2[9] = {10.,1.,1e-1,1e-2,1e-3,1e-4,1e-4,1e-5,1e-6};
  int nb[9] =  {1,2,5,10,20,40,80,80,80};
  int nb2[9] = {1,2,5,10,20,40,40,40,40};
  int maxprint = 9;


  cout << "#Partially quenched finite volume corrections three flavours\n";
  cout << lirandom;

  cout << "\n#Mass one valence one sea three sea quarks\n";
  cout << B0mq11;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv1s1nf3p4Vt       mv1s1nf3p4Vb        "
	 << "mv1s1nf3p6LVt       mv1s1nf3p6LVb        mv1s1nf3p6RVt       mv1s1nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv1s1nf3p4Vt(B0mq11,xl) << flush
	  <<' '<<setw(13)<< mv1s1nf3p4Vb(B0mq11,xl) << flush
	  <<' '<<setw(13)<< mv1s1nf3p6LVt(B0mq11,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s1nf3p6LVb(B0mq11,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s1nf3p6RVt(B0mq11,xl) << flush
	  <<' '<<setw(13)<< mv1s1nf3p6RVb(B0mq11,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Mass two valence one sea three sea quarks\n";
  cout << B0mq21;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv2s1nf3p4Vt       mv2s1nf3p4Vb        "
	 << "mv2s1nf3p6LVt       mv2s1nf3p6LVb        mv2s1nf3p6RVt       mv2s1nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv2s1nf3p4Vt(B0mq21,xl) << flush
	  <<' '<<setw(13)<< mv2s1nf3p4Vb(B0mq21,xl) << flush
	  <<' '<<setw(13)<< mv2s1nf3p6LVt(B0mq21,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s1nf3p6LVb(B0mq21,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s1nf3p6RVt(B0mq21,xl) << flush
	  <<' '<<setw(13)<< mv2s1nf3p6RVb(B0mq21,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Mass one valence two sea three sea quarks\n";
  cout << B0mq12;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv1s2nf3p4Vt       mv1s2nf3p4Vb        "
	 << "mv1s2nf3p6LVt       mv1s2nf3p6LVb        mv1s2nf3p6RVt       mv1s2nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv1s2nf3p4Vt(B0mq12,xl) << flush
	  <<' '<<setw(13)<< mv1s2nf3p4Vb(B0mq12,xl) << flush
	  <<' '<<setw(13)<< mv1s2nf3p6LVt(B0mq12,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s2nf3p6LVb(B0mq12,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s2nf3p6RVt(B0mq12,xl) << flush
	  <<' '<<setw(13)<< mv1s2nf3p6RVb(B0mq12,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Mass two valence two sea three sea quarks\n";
  cout << B0mq22;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv2s2nf3p4Vt       mv2s2nf3p4Vb        "
	 << "mv2s2nf3p6LVt       mv2s2nf3p6LVb        mv2s2nf3p6RVt       mv2s2nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv2s2nf3p4Vt(B0mq22,xl) << flush
	  <<' '<<setw(13)<< mv2s2nf3p4Vb(B0mq22,xl) << flush
	  <<' '<<setw(13)<< mv2s2nf3p6LVt(B0mq22,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s2nf3p6LVb(B0mq22,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s2nf3p6RVt(B0mq22,xl) << flush
	  <<' '<<setw(13)<< mv2s2nf3p6RVb(B0mq22,xl) << flush
	  <<'\n';
    }
  }
  cout << "\n#Mass one valence three sea three sea quarks\n";
  cout << B0mq13;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv1s3nf3p4Vt       mv1s3nf3p4Vb        "
	 << "mv1s3nf3p6LVt       mv1s3nf3p6LVb        mv1s3nf3p6RVt       mv1s3nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv1s3nf3p4Vt(B0mq13,xl) << flush
	  <<' '<<setw(13)<< mv1s3nf3p4Vb(B0mq13,xl) << flush
	  <<' '<<setw(13)<< mv1s3nf3p6LVt(B0mq13,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s3nf3p6LVb(B0mq13,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv1s3nf3p6RVt(B0mq13,xl) << flush
	  <<' '<<setw(13)<< mv1s3nf3p6RVb(B0mq13,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Mass two valence three sea three sea quarks\n";
  cout << B0mq23;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  mv2s3nf3p4Vt       mv2s3nf3p4Vb        "
	 << "mv2s3nf3p6LVt       mv2s3nf3p6LVb        mv2s3nf3p6RVt       mv2s3nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mv2s3nf3p4Vt(B0mq23,xl) << flush
	  <<' '<<setw(13)<< mv2s3nf3p4Vb(B0mq23,xl) << flush
	  <<' '<<setw(13)<< mv2s3nf3p6LVt(B0mq23,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s3nf3p6LVb(B0mq23,lirandom,xl) << flush
	  <<' '<<setw(13)<< mv2s3nf3p6RVt(B0mq23,xl) << flush
	  <<' '<<setw(13)<< mv2s3nf3p6RVb(B0mq23,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n Decay constants\n";

  cout << "\n#Decay constant one valence one sea three sea quarks\n";
  cout << B0mq11;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv1s1nf3p4Vt       fv1s1nf3p4Vb        "
	 << "fv1s1nf3p6LVt       fv1s1nf3p6LVb        fv1s1nf3p6RVt       fv1s1nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv1s1nf3p4Vt(B0mq11,xl) << flush
	  <<' '<<setw(13)<< fv1s1nf3p4Vb(B0mq11,xl) << flush
	  <<' '<<setw(13)<< fv1s1nf3p6LVt(B0mq11,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s1nf3p6LVb(B0mq11,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s1nf3p6RVt(B0mq11,xl) << flush
	  <<' '<<setw(13)<< fv1s1nf3p6RVb(B0mq11,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Decay constant two valence one sea three sea quarks\n";
  cout << B0mq21;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv2s1nf3p4Vt       fv2s1nf3p4Vb        "
	 << "fv2s1nf3p6LVt       fv2s1nf3p6LVb        fv2s1nf3p6RVt       fv2s1nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv2s1nf3p4Vt(B0mq21,xl) << flush
	  <<' '<<setw(13)<< fv2s1nf3p4Vb(B0mq21,xl) << flush
	  <<' '<<setw(13)<< fv2s1nf3p6LVt(B0mq21,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s1nf3p6LVb(B0mq21,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s1nf3p6RVt(B0mq21,xl) << flush
	  <<' '<<setw(13)<< fv2s1nf3p6RVb(B0mq21,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Decay constant one valence two sea three sea quarks\n";
  cout << B0mq12;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv1s2nf3p4Vt       fv1s2nf3p4Vb        "
	 << "fv1s2nf3p6LVt       fv1s2nf3p6LVb        fv1s2nf3p6RVt       fv1s2nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv1s2nf3p4Vt(B0mq12,xl) << flush
	  <<' '<<setw(13)<< fv1s2nf3p4Vb(B0mq12,xl) << flush
	  <<' '<<setw(13)<< fv1s2nf3p6LVt(B0mq12,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s2nf3p6LVb(B0mq12,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s2nf3p6RVt(B0mq12,xl) << flush
	  <<' '<<setw(13)<< fv1s2nf3p6RVb(B0mq12,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Decay constant two valence two sea three sea quarks\n";
  cout << B0mq22;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv2s2nf3p4Vt       fv2s2nf3p4Vb        "
	 << "fv2s2nf3p6LVt       fv2s2nf3p6LVb        fv2s2nf3p6RVt       fv2s2nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv2s2nf3p4Vt(B0mq22,xl) << flush
	  <<' '<<setw(13)<< fv2s2nf3p4Vb(B0mq22,xl) << flush
	  <<' '<<setw(13)<< fv2s2nf3p6LVt(B0mq22,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s2nf3p6LVb(B0mq22,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s2nf3p6RVt(B0mq22,xl) << flush
	  <<' '<<setw(13)<< fv2s2nf3p6RVb(B0mq22,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Decay constant one valence three sea three sea quarks\n";
  cout << B0mq13;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv1s3nf3p4Vt       fv1s3nf3p4Vb        "
	 << "fv1s3nf3p6LVt       fv1s3nf3p6LVb        fv1s3nf3p6RVt       fv1s3nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv1s3nf3p4Vt(B0mq13,xl) << flush
	  <<' '<<setw(13)<< fv1s3nf3p4Vb(B0mq13,xl) << flush
	  <<' '<<setw(13)<< fv1s3nf3p6LVt(B0mq13,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s3nf3p6LVb(B0mq13,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv1s3nf3p6RVt(B0mq13,xl) << flush
	  <<' '<<setw(13)<< fv1s3nf3p6RVb(B0mq13,xl) << flush
	  <<'\n';
    }
  }

  cout << "\n#Decay constant two valence three sea three sea quarks\n";
  cout << B0mq23;
  for(int il = 2; il < 5; il++){
    double xl = double(il)/sqrt(m11);
    cout << "output for m11 L = "<<il<<'\n';
    cout << " eps    eps2  maxonesum maxtwosum  fv2s3nf3p4Vt       fv2s3nf3p4Vb        "
	 << "fv2s3nf3p6LVt       fv2s3nf3p6LVb        fv2s3nf3p6RVt       fv2s3nf3p6RVb"  	 <<'\n';
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<setw(8)<<epst2[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fv2s3nf3p4Vt(B0mq23,xl) << flush
	  <<' '<<setw(13)<< fv2s3nf3p4Vb(B0mq23,xl) << flush
	  <<' '<<setw(13)<< fv2s3nf3p6LVt(B0mq23,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s3nf3p6LVb(B0mq23,lirandom,xl) << flush
	  <<' '<<setw(13)<< fv2s3nf3p6RVt(B0mq23,xl) << flush
	  <<' '<<setw(13)<< fv2s3nf3p6RVb(B0mq23,xl) << flush
	  <<'\n';
    }
  }

  return 0;
}
