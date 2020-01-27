// testintegralsreal.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.
#include<ctime>
#include<cmath>
#include<iostream>
#include<iomanip>
#include "jbnumlib.h"

/////////////////////////////////////////////////////////////////////
// real integral
double f1(const double x){
  return sqrt(x);
}
// with wiggles
double f2(const double x){
  return 1.+sin(10.*x);
}
// one that return zero
double f3(const double x){
  return 1.-5.*x*x*x*x;
}

double f4(const double x){
  return 1.-sqrt(x);
}

// integrable sqrt divergence at 1
double f5(const double x){
  return 1./sqrt(1.-x*x);
}
double f6(const double x){
  return x/sqrt(1.-x*x);
}
// integral log divergence at 0
double f7(const double x){
  return log(x);
}

int main(){
  using namespace std;
  double accu = 0.;
  double result;
  clock_t t,tp;
  int imax;

  ///////real integrals no singularities//////////////////////////////
  // used for testing how fast they are
  cout << 
    "printed are obtained values and relative difference with exact result\n";
  cout << "int_0^4 sqrt(x) dx\n";
  accu = 1e-14;
  imax = 300000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f1,0.,4.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3./16.*result-1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f1,0.,4.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3./16.*result-1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f1,0.,4.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3./16.*result-1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f1,0.,4.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3./16.*result-1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^pi (1+sin(10x)) dx\n";
  accu = 1e-14;
  imax = 1000000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f2,0.,M_PI,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/M_PI-1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f2,0.,M_PI,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/M_PI-1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f2,0.,M_PI,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/M_PI-1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f2,0.,M_PI,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/M_PI-1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1(1-5x^4)dx\n";
  accu = 1e-15;
  imax = 10000000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f3,0.,1.,accu);
  t = clock();
  cout << result << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f3,0.,1.,accu);
  t = clock();
  cout << result << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f3,0.,1.,accu);
  t = clock();
  cout << result << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f3,0.,1.,accu);
  t = clock();
  cout << result << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 (1-sqrt(x)) dx\n";
  accu = 1e-14;
  imax = 300000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f4,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3.*result-1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f4,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3.*result-1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f4,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3.*result-1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f4,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (3.*result-1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 (1/sqrt(1-x^2)) dx\n";
  accu = 1e-7;
  imax = 100000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f5,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (2./M_PI*result-1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f5,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (2./M_PI*result-1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f5,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (2./M_PI*result-1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f5,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (2./M_PI*result-1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 (x/sqrt(1-x^2)) dx\n";
  accu = 1e-7;
  imax = 100000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f6,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result-1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f6,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result-1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f6,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result-1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f6,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result-1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 log(x) dx\n";
  accu = 1e-10;
  imax = 100000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdgauss(f7,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result+1.) 
       << " via jbdgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdgauss2(f7,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result+1.) 
       << " via jbdgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad15(f7,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result+1.) 
       << " via jbdquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdquad21(f7,0.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result+1.) 
       << " via jbdquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  return 0;
}
