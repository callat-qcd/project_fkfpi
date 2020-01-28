// testmassdecayvevV.cc is part of the CHIRON ChPT at two loops program
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
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
#include "massdecayvevV.h"

const double pi16 = 1./pow(4.*M_PI,2);
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
    cout << " eps maxonesum maxtwosum  mpi4Vt      mpi4Vb      mpi6LVt "
	 << "      mpi6LVb        mpi6RVt        mpi6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< mpi4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<<mpi4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<< mpi6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<mpi6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< mpi6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<mpi6RVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum   mk4Vt       mk4Vb       mk6LVt "
	 << "       mk6LVb         mk6RVt         mk6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<  mk4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<< mk4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<<  mk6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< mk6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<  mk6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< mk6RVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum meta4Vt     meta4Vb     meta6LVt "
	 << "     meta6LVb       meta6RVt       meta6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<meta4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<meta6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<meta6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<meta6RVb(stdmass,xl)
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum  fpi4Vt      fpi4Vb      fpi6LVt "
	 << "      fpi6LVb        fpi6RVt        fpi6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<< fpi4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<<fpi4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<< fpi6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<fpi6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< fpi6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<fpi6RVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum   fk4Vt       fk4Vb       fk6LVt "
	 << "       fk6LVb         fk6RVt         fk6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<  fk4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<< fk4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<<  fk6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<< fk6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<  fk6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<< fk6RVb(stdmass,xl) << flush
	  <<'\n';
    }
    cout << " eps maxonesum maxtwosum feta4Vt     feta4Vb     feta6LVt "
	 << "     feta6LVb       feta6RVt       feta6RVb\n";
    for(int i=0;i<maxprint;i++){
      setprecisionfinitevolumeoneloopt(epst[i],1e-5,false);
      setprecisionfinitevolumeoneloopb(nb[i],1e-5,false);
      setprecisionfinitevolumesunsett(epst[i],epst2[i],false);
      setprecisionfinitevolumesunsetb(nb[i],nb2[i],epst[i],epst2[i],false);
      cout<<setw(8)<< epst[i]<<' '<<setw(4)<<nb[i]<<' '<<setw(4)<<nb2[i]
	  <<' '<<setw(13)<<feta4Vt(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta4Vb(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta6LVt(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<feta6LVb(stdmass,liBE14,xl) << flush
	  <<' '<<setw(13)<<feta6RVt(stdmass,xl) << flush
	  <<' '<<setw(13)<<feta6RVb(stdmass,xl) << flush
	  <<'\n';
    }
  }

  return 0;
}
