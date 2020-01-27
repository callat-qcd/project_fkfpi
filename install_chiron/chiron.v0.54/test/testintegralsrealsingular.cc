// g++ testintegrals.cc -ljbnumlib
#include<ctime>
#include<cmath>
#include<iostream>
#include<iomanip>
#include "jbnumlib.h"

/////////////////////////////////////////////////////////////////////////
// real integral with singularity

double f8(const double x){
  return 1./(1.-x*x);
}

double f9(const double x){
  return x/(1.-x*x);
}

double f10(const double x){
  return sqrt(x)/(x-0.3);
}

double f11(const double x){
  return 2.*x*x/(x*x-0.3);
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

// real integrals singularities /////////////////////////////////
  cout << "int_0^3 1/(1-x^2) dx\n";
  accu = 1e-13;
  imax = 1000000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdcauch(f8,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(0.5*log(2.))-1.) 
       << " via jbdcauch in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdcauch2(f8,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(0.5*log(2.))-1.) 
       << " via jbdcauch2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing15(f8,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(0.5*log(2.))-1.) 
       << " via jbdsing15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing21(f8,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(0.5*log(2.))-1.) 
       << " via jbdsing21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;


  cout << "int_0^3 x/(1-x^2) dx\n";
  accu = 1e-13;
  imax = 1000000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdcauch(f9,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(-0.5*log(8.))-1.) 
       << " via jbdcauch in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdcauch2(f9,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(-0.5*log(8.))-1.) 
       << " via jbdcauch2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing15(f9,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(-0.5*log(8.))-1.) 
       << " via jbdsing15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing21(f9,0.,3.,1.,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/(-0.5*log(8.))-1.) 
       << " via jbdsing21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 sqrt(x)/(x-0.3) dx\n";
  accu = 1e-13;
  imax = 100000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  double a = 0.3;
  double sqrta = sqrt(a);
  double resulta = 2.-sqrta*log((1.+sqrta)/(1-sqrta));
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdcauch(f10,0.,1.,a,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdcauch in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdcauch2(f10,0.,1.,a,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdcauch2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing15(f10,0.,1.,a,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdsing15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing21(f10,0.,1.,a,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdsing21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^1 2x^2/(x^2-0.3) dx\n";
  accu = 1e-13;
  imax = 5000000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  a = 0.3;
  sqrta = sqrt(a);
  resulta = 2.-sqrta*log((1.+sqrta)/(1-sqrta));
  tp = clock();
  for(int i=0; i<imax; i++) result = jbdcauch(f11,0.,1.,sqrta,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdcauch in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdcauch2(f11,0.,1.,sqrta,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdcauch2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing15(f11,0.,1.,sqrta,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdsing15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbdsing21(f11,0.,1.,sqrta,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< (result/resulta-1.)
       << " via jbdsing21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  return 0;
}
