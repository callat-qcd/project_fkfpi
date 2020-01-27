// testquenchedsunsetintegrals.cc is part of the CHIRON ChPT
// at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains all the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]
// has the extra ones needed for the partiallyquenched case with the
// extra first integer index indicating the power of the propagators added
// Derived in:
//
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 73} (2006) 074509
//  [hep-lat/0602003].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 72} (2005) 074502
//  [hep-lat/0506004].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 71} (2005) 094502
//  [hep-lat/0501014].
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 70} (2004) 111503
//  [hep-lat/0406017].

// g++ testquenchedsunsetintegrals.cc quenchedsunsetintegrals.o -ljbnumlib

#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>

#include "quenchedsunsetintegrals.h"


using namespace std;
int main(void){
  const double pi162 = 1./pow(4.*M_PI,4);
  setprecisionquenchedsunsetintegrals(1e-10*pi162);
  double mpi = 0.135;
  double mk = 0.495;
  double meta = 0.548;
  double mu = 0.77;
  double mu2 = mu*mu;
  double m1sq = mpi*mpi;
  double m2sq = mk*mk;
  double m3sq = meta*meta;
  double qsq = -m3sq;
  ifstream fortrandata("test/comparequenchedsunsetfortran.dat");
  string temps;
  getline(fortrandata,temps);
  cout << "precision set to :"<<getprecisionquenchedsunsetintegrals()<<'\n';
  int iin;
  double fortranin;
  // for comparing with the old fortran program
  getline(fortrandata,temps);
  cout <<"m1sq m2sq m3sq qsq mu2: "<<m1sq<<' '<<m2sq<<' '<<m3sq<<' '<<mu2<<'\n';
  getline(fortrandata,temps);
  cout << "hh, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hh(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout <<"here "<<temps<<'\n';
  cout << "hh1, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hhd, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hhd(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1d(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21d(i,m1sq,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  // lm = 0 case 1
  cout << "lm = 0, case 1\n";
  double m1s =pow(sqrt(m2sq)+sqrt(m3sq),2);
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout <<"m1sq m2sq m3sq qsq mu2: "<<m1s<<' '<<m2sq<<' '<<m3sq<<' '<<mu2<<'\n';
  getline(fortrandata,temps);
  cout << "hh, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hh(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hhd, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hhd(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1d(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21d(i,m1s,m2sq,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  // lm = 0 case 2
  cout << "lm = 0, case 2\n";
  double m2s =pow(sqrt(m1sq)+sqrt(m3sq),2);
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout <<"m1sq m2sq m3sq qsq mu2: "<<m1sq<<' '<<m2s<<' '<<m3sq<<' '<<mu2<<'\n';
  getline(fortrandata,temps);
  cout << "hh, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hh(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hhd, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hhd(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1d(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21d(i,m1sq,m2s,m3sq,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  // lm = 0 case 3
  cout << "lm = 0, case 3\n";
  double m3s =pow(sqrt(m1sq)+sqrt(m2sq),2);
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout <<"m1sq m2sq m3sq qsq mu2: "<<m1sq<<' '<<m2sq<<' '<<m3s<<' '<<mu2<<'\n';
  getline(fortrandata,temps);
  cout << "hh, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hh(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hhd, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i <<' '
	 <<setprecision(10) << hhd(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh1d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh1d(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }
  getline(fortrandata,temps);// blanks after the last number
  getline(fortrandata,temps);
  cout << "hh21d, iprop=1,8\n";
  for(int i=1;i<9;i++){
    fortrandata >> iin >> fortranin;
    cout << setw(2) << i 
	 <<' '<<setprecision(10) << hh21d(i,m1sq,m2sq,m3s,qsq,mu2)
	 <<' '<<setprecision(10)<<fortranin
	 << '\n';
  }


  return 0;
}



