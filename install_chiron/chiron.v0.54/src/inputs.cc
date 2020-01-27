// inputs.cc is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.


// implementation of input classes
// and setting them to particular values

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "inputs.h"

//++++++++ physmass ++++++++++++++++++++++++++++++++++++++++++++++++++++++

physmass::~physmass(void){}

physmass::physmass(const double mpiin, const double mkin,
		   const double metain,  const double fpiin,
		   const double muin){
  mpi = mpiin;
  mk = mkin;
  meta = metain;
  fpi = fpiin;
  mu = muin;
}

void physmass::setmpi(const double mpiin){
  mpi = mpiin;
}
void physmass::setmk(const double mkin){
  mk = mkin;
}
void physmass::setmeta(const double metain){
  meta = metain;
}
void physmass::setfpi(const double fpiin){
  fpi = fpiin;
}
void physmass::setmu(const double muin){
  mu = muin;
}

void physmass::out(double &mpiout, double &mkout, double &metaout,
		   double &fpiout, double &muout) const{
  mpiout = mpi;
  mkout = mk;
  metaout = meta;
  fpiout = fpi;
  muout = mu;
}

double physmass::getmpi(void ) const{
  return mpi;
}
double physmass::getmk(void) const{
  return mk;
}
double physmass::getmeta(void) const{
  return meta;
}
double physmass::getfpi(void) const{
  return fpi;
}
double physmass::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const physmass & bb){
  os   << setprecision(15) << fixed
       << "# mpi :  "<< bb.mpi << "\n"
       << "# mk  :  "<<bb.mk <<"\n"
       << "# meta:  "<<bb.meta <<"\n"
       << "# fpi :  "<<bb.fpi <<"\n"
       << "# mu  :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, physmass & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double mpi,mk,meta,fpi,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> mpi;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# mpi ") std::cout << "trouble reading in mass mpi "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mk;
  std::getline(is,temps);
  if (temps2 != "# mk  ") std::cout << "trouble reading in mass mk "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> meta;
  std::getline(is,temps);
  if (temps2 != "# meta") std::cout <<"trouble reading in mass meta "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> fpi;
  std::getline(is,temps);
  if (temps2 != "# fpi ") std::cout <<"trouble reading in mass fpi "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu  ") std::cout << "trouble reading in mass mu "<< temps2<<'\n';
  mass = physmass(mpi,mk,meta,fpi,mu);
  return is;
}

bool operator==(const physmass mass1,const physmass mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.mpi-mass2.mpi)/mass1.mpi) < massprecision) &&
       (abs((mass1.mk-mass2.mk)/mass1.mk) < massprecision) &&
       (abs((mass1.meta-mass2.meta)/mass1.meta) < massprecision) &&
       (abs((mass1.fpi-mass2.fpi)/mass1.fpi) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

//+++++++++++lomass+++++++++++++++++++++++++++++++++++++++++++++++++++++++

lomass::~lomass(void){}

lomass::lomass(const double mp0in, const double mk0in,
	       const double f0in, const double muin){
  mp0 = mp0in;
  mk0 = mk0in;
  f0 = f0in;
  mu = muin;
}
lomass::lomass(const quarkmass mass){
  mp0 = sqrt(2.*mass.getB0mhat());
  mk0 = sqrt(mass.getB0mhat()+mass.getB0ms());
  f0 = mass.getf0();
  mu = mass.getmu();
}
void lomass::setmp0(const double mp0in){
  mp0 = mp0in;
}
void lomass::setmk0(const double mk0in){
  mk0 = mk0in;
}
void lomass::setf0(const double f0in){
  f0 = f0in;
}
void lomass::setmu(const double muin){
  mu = muin;
}
void lomass::out(double &mp0out, double &mk0out,
		   double &f0out, double &muout) const{
  mp0out = mp0;
  mk0out = mk0;
  f0out = f0;
  muout = mu;
}

double lomass::getmp0(void ) const{
  return mp0;
}
double lomass::getmk0(void) const{
  return mk0;
}
double lomass::getf0(void) const{
  return f0;
}
double lomass::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const lomass & bb){
  os   << setprecision(15) << fixed
       << "# mp0 :  "<< bb.mp0 << "\n"
       << "# mk0 :  "<<bb.mk0 <<"\n"
       << "# f0  :  "<<bb.f0 <<"\n"
       << "# mu  :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, lomass & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double mp0,mk0,f0,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> mp0;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# mp0 ") std::cout << "trouble reading in lomass mp0 "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mk0;
  std::getline(is,temps);
  if (temps2 != "# mk0 ") std::cout << "trouble reading in lomass mk0 "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> f0;
  std::getline(is,temps);
  if (temps2 != "# f0  ") std::cout <<"trouble reading in lomass f0 "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu  ") std::cout << "trouble reading in lomass mu "<< temps2<<'\n';
  mass = lomass(mp0,mk0,f0,mu);
  return is;
}

bool operator==(const lomass mass1,const lomass mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.mp0-mass2.mp0)/mass1.mp0) < massprecision) &&
       (abs((mass1.mk0-mass2.mk0)/mass1.mk0) < massprecision) &&
       (abs((mass1.f0-mass2.f0)/mass1.f0) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

//+++++++++++quarkmass+++++++++++++++++++++++++++++++++++++++++++++++++++++++

quarkmass::~quarkmass(void){}

quarkmass::quarkmass(const double B0mhatin, const double B0msin,
	       const double f0in, const double muin){
  B0mhat = B0mhatin;
  B0ms = B0msin;
  f0 = f0in;
  mu = muin;
}
quarkmass::quarkmass(const lomass mass){
  B0mhat = pow(mass.getmp0(),2)/2.;
  B0ms = pow(mass.getmk0(),2)-B0mhat;
  f0 = mass.getf0();
  mu = mass.getmu();
}
void quarkmass::setB0mhat(const double B0mhatin){
  B0mhat = B0mhatin;
}
void quarkmass::setB0ms(const double B0msin){
  B0ms = B0msin;
}
void quarkmass::setf0(const double f0in){
  f0 = f0in;
}
void quarkmass::setmu(const double muin){
  mu = muin;
}
void quarkmass::out(double &B0mhatout, double &B0msout,
		   double &f0out, double &muout) const{
  B0mhatout = B0mhat;
  B0msout = B0ms;
  f0out = f0;
  muout = mu;
}

double quarkmass::getB0mhat(void ) const{
  return B0mhat;
}
double quarkmass::getB0ms(void) const{
  return B0ms;
}
double quarkmass::getf0(void) const{
  return f0;
}
double quarkmass::getmu(void) const{
  return mu;
}

std::ostream & operator<<(std::ostream & os, const quarkmass & bb){
  os   << setprecision(15) << fixed
       << "# B0mhat :  "<< bb.B0mhat << "\n"
       << "# B0ms   :  "<<bb.B0ms <<"\n"
       << "# f0     :  "<<bb.f0 <<"\n"
       << "# mu     :  "<<bb.mu <<"\n";
  os.unsetf(std::ios_base::floatfield);
  return os;
}

std::istream & operator>>(std::istream & is, quarkmass & mass){
  string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  string temps2;
  double B0mhat,B0ms,f0,mu;
  std::getline(is,temps2,':');// reads in characters up to ':'
  is >> B0mhat;
  std::getline(is,temps); // to remove end of line
  if (temps2 != "# B0mhat ") std::cout << "trouble reading in quarkmass B0mhat "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> B0ms;
  std::getline(is,temps);
  if (temps2 != "# B0ms   ") std::cout << "trouble reading in quarkmass B0ms "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> f0;
  std::getline(is,temps);
  if (temps2 != "# f0     ") std::cout <<"trouble reading in quarkmass f0 "<< temps2<<'\n';
  std::getline(is,temps2,':');
  is >> mu;
  std::getline(is,temps);
  if (temps2 != "# mu     ") std::cout << "trouble reading in quarkmass mu "<< temps2<<'\n';
  mass = quarkmass(B0mhat,B0ms,f0,mu);
  return is;
}

bool operator==(const quarkmass mass1,const quarkmass mass2){
  const double massprecision = 1e-7;
  if ( (abs((mass1.B0mhat-mass2.B0mhat)/mass1.B0mhat) < massprecision) &&
       (abs((mass1.B0ms-mass2.B0ms)/mass1.B0ms) < massprecision) &&
       (abs((mass1.f0-mass2.f0)/mass1.f0) < massprecision) &&
       (abs((mass1.mu-mass2.mu)/mass1.mu) < massprecision) ) return true;
  return false;
}

