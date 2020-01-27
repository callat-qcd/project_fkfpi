// testmassdecayvevloV.cc is part of the CHIRON ChPT at two loops program
// collection
// Copyright (C) 2015 Johan Bijnens, v1.0
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
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "massdecayvevloV.h"

const double pi16 = 1./pow(4.*M_PI,2);
const double pi162 = 1./pow(4.*M_PI,4);

const double mpo = 0.130;
const double mk = 0.450;

const lomass stdmass(mpo,mk,0.0877,0.77);

using namespace std;
int main(void){
  // reading in an example set of LECS
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

  cout << "input used:\n";
  cout << stdmass << liBE14 << setprecision(6);

  // setting all needed precision one-loop
  double epst[8]  = {10.,1.,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6};
  double epst2[8] = {10.,1.,1e-1,1e-2,1e-3,1e-4,1e-4,1e-5};
  int nb[8] =  {1,2,5,10,20,40,80,160};
  int nb2[8] = {1,2,5,10,20,40,40,40};
  int maxprint = 7;

  for(int il = 2; il < 5; il++){
    double xl = double(il)/mpo;
    cout << "output for m_pi L = "<<il<<'\n';
    cout << " eps maxonesum maxtwosum  mpi4loVt    mpi4loVb    mpi6LloVt "
	 << "    mpi6LloVb      mpi6RloVt      mpi6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mpi4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<mpi4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<< mpi6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<mpi6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< mpi6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<mpi6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum   mk4loVt     mk4loVb     mk6LloVt "
	 << "     mk6LloVb       mk6RloVt       mk6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<  mk4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< mk4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<  mk6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< mk6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<  mk6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< mk6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum meta4loVt   meta4loVb   meta6LloVt "
	 << "   meta6LloVb     meta6RloVt     meta6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<meta4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<meta6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<meta6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta6RloVb(stdmass,xl)
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum  fpi4loVt    fpi4loVb    fpi6LloVt "
	 << "    fpi6LloVb      fpi6RloVt      fpi6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fpi4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<fpi4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<< fpi6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<fpi6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< fpi6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<fpi6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum   fk4loVt     fk4loVb     fk6LloVt "
	 << "     fk6LloVb       fk6RloVt       fk6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<  fk4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< fk4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<  fk6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< fk6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<  fk6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< fk6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum feta4loVt   feta4loVb   feta6LloVt "
	 << "    feta6LloVb     feta6RloVt     feta6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<feta4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<feta6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<feta6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum qqup4loVt   qqup4loVb   qqup6LloVt "
	 << "    qqup6LloVb     qqup6RloVt     qqup6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<qqup4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqup4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqup6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<qqup6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<qqup6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqup6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum qqstrange4loVt qqstrange4loVb qqstrange6LloVt"
	 << " qqstrange6LloVb qqstrange6RloVt qqstrange6RloVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<qqstrange4loVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqstrange4loVb(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqstrange6LloVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<qqstrange6LloVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<qqstrange6RloVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<qqstrange6RloVb(stdmass,xl) << flush
	  <<'\n';
    }
  }

  return 0;
}
