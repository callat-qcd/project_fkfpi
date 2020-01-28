// g++ testintegrals.cc -ljbnumlib
#include<ctime>
#include<cmath>
#include<complex>
#include<iostream>
#include<iomanip>
#include "jbnumlib.h"
typedef std::complex<double> dcomplex;


// complex line integral

dcomplex alpha = dcomplex(6.,5.),beta = dcomplex(-25.,0.5);
dcomplex zero = dcomplex(0.,0.);

dcomplex myexp(const dcomplex x){
  return exp(beta*x);
}

dcomplex intmyexp(){
  return exp(beta*alpha)/beta-dcomplex(1.,0.)/beta;
}

// beware to not cross the cut
dcomplex f12(const dcomplex x){
  return sqrt(x);
}

int main(){
  using namespace std;
  double accu = 0.;
  dcomplex result,resultint;
  clock_t t,tp;
  int imax;

  // used for testing how fast they are
  cout << 
    "printed are obtained values and relative difference with exact result\n";
  cout << "int_0^alpha exp(beta*x) dx\n";
  accu = 1e-15;
  imax = 10000;
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbwgauss(myexp,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/intmyexp()-1.) << " via jbwgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwgauss2(myexp,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/intmyexp()-1.) << " via jbwgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad15(myexp,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/intmyexp()-1.) << " via jbwquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad21(myexp,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/intmyexp()-1.) << " via jbwquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_beta^alpha sqrt(x) dx\n";
  accu = 1e-15;
  imax = 100000;
  resultint = 2./3.*(pow(alpha,1.5)-pow(beta,1.5));
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbwgauss(f12,beta,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  tp = clock();
  for(int i=0; i<imax; i++) result = jbwgauss2(f12,beta,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad15(f12,beta,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad21(f12,beta,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;

  cout << "int_0^alpha sqrt(x) dx\n";
  accu = 1e-15;
  imax = 10000;
  resultint = 2./3.*(pow(alpha,1.5));
  cout << "accuracy asked " << accu << " number of times evaluated "
       << imax <<'\n';
  tp = clock();
  for(int i=0; i<imax; i++) result = jbwgauss(f12,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwgauss in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  tp = clock();
  for(int i=0; i<imax; i++) result = jbwgauss2(f12,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwgauss2 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad15(f12,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwquad15 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  for(int i=0; i<imax; i++) result = jbwquad21(f12,zero,alpha,accu);
  t = clock();
  cout << setprecision(15) <<result <<' '<< setprecision(3)
       << (result/resultint-1.) << " via jbwquad21 in "
       <<double(t-tp)/double(CLOCKS_PER_SEC)<<" sec\n";
  tp = t;
  return 0;
}
