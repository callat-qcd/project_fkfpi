// Li.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// implementation of the Li class

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>  // for rand


#include "Li.h"

const double pi16 = 1./(16.*M_PI*M_PI);

Li::Li(const double l1r,const double l2r,const double l3r,
       const double l4r,const double l5r,const double l6r,
       const double l7r,const double l8r,const double l9r,
       const double l10r,const double h1r,const double h2r,
       const double muin,
       const std::string Name){
    L1r = l1r;
    L2r = l2r; 
    L3r = l3r;
    L4r = l4r;
    L5r = l5r;
    L6r = l6r;
    L7r = l7r;
    L8r = l8r;
    L9r = l9r;
    L10r = l10r;
    H1r = h1r;
    H2r = h2r;
    mu = muin;
    name = Name;
}

Li::~Li(void){}

std::ostream & operator<<(std::ostream & os, const Li & bb){
  os   << std::setprecision(15)
       << "# Li set : " << bb.name << "\n"
       << "#mu : "<<bb.mu <<"\n"
       << "#L1r  L2r  L3r :" <<bb.L1r <<' '<<bb.L2r<<' '<<bb.L3r<<"\n"
       << "#L4r  L5r  L6r :" <<bb.L4r <<' '<<bb.L5r<<' '<<bb.L6r<<"\n"
       << "#L7r  L8r  L9r :" <<bb.L7r <<' '<<bb.L8r<<' '<<bb.L9r<<"\n"
       << "#L10r H1r  H2r :" <<bb.L10r<<' '<<bb.H1r<<' '<<bb.H2r<<"\n";
  os.unsetf(std::ios_base::floatfield);
  os << std::setprecision(6);
  return os;
}

std::istream & operator>>(std::istream & is, Li & Liout){
  std::string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  std::string temps2;
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,h1r,h2r,mu;
  std::string name;
  std::getline(is,temps2,':');// reads in characters up to ':'
  std::getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# Li set ") std::cout << "trouble reading in Li\n";
  std::getline(is,temps2,':');   is >> mu;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l1r >> l2r >> l3r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l4r >> l5r >> l6r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l7r >> l8r >> l9r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l10r>> h1r >> h2r;  getline(is,temps);
  Liout = Li(l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,h1r,h2r,mu,name);
  return is;
}


void Li::setli(const int n, const double lin){
  switch(n){
  case 1: L1r = lin; break;
  case 2: L2r = lin; break;  
  case 3: L3r = lin; break;  
  case 4: L4r = lin; break;  
  case 5: L5r = lin; break;  
  case 6: L6r = lin; break;  
  case 7: L7r = lin; break;  
  case 8: L8r = lin; break;  
  case 9: L9r = lin; break;  
  case 10: L10r = lin; break;  
  case 11: H1r = lin; break;  
  case 12: H2r = lin; break;  
  default: std::cout << "Trying to set Li number "<<n<<" in set "<<name<<'\n';
  }
}

void Li::setli(const double lin,const int n){
  setli(n,lin);
}

void Li::setmu(const double muin){
  mu = muin;
}

void Li::setname(const std::string inputname){
  name = inputname;
}

// summing two sets of Li
Li Li::operator+(const Li & bb) const{
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,l9r,l8r,l10r,h1r,h2r,muout;
  l1r = L1r+bb.L1r;
  l2r = L2r+bb.L2r;
  l3r = L3r+bb.L3r;
  l4r = L4r+bb.L4r;
  l5r = L5r+bb.L5r;
  l6r = L6r+bb.L6r;
  l7r = L7r+bb.L7r;
  l8r = L8r+bb.L8r;
  l9r = L9r+bb.L9r;
  l10r = L10r+bb.L10r;
  h1r = H1r+bb.H1r;
  h2r = H2r+bb.H2r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: summing Li \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return Li(l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,h1r,h2r,muout);
}
// difference of two sets of Li
Li Li::operator-(const Li & bb) const{
  double l1r,l2r,l3r,l4r,l5r,l6r,l7r,l9r,l8r,l10r,h1r,h2r,muout;
  l1r = L1r-bb.L1r;
  l2r = L2r-bb.L2r;
  l3r = L3r-bb.L3r;
  l4r = L4r-bb.L4r;
  l5r = L5r-bb.L5r;
  l6r = L6r-bb.L6r;
  l7r = L7r-bb.L7r;
  l8r = L8r-bb.L8r;
  l9r = L9r-bb.L9r;
  l10r = L10r-bb.L10r;
  h1r = H1r-bb.H1r;
  h2r = H2r-bb.H2r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: difference of Li \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return Li(l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,h1r,h2r,muout);
}

// multiplying an Li with a number

Li Li::operator*(const double & aa) const{
  return Li(aa*L1r,aa*L2r,aa*L3r,aa*L4r,aa*L5r,aa*L6r,aa*L7r,aa*L8r,
	    aa*L9r,aa*L10r,aa*H1r,aa*H2r,mu);
}

Li operator*(const double & aa,const Li & bb){
  return Li(aa*bb.L1r,aa*bb.L2r,aa*bb.L3r,aa*bb.L4r,aa*bb.L5r,aa*bb.L6r,
	    aa*bb.L7r,aa*bb.L8r,aa*bb.L9r,aa*bb.L10r,aa*bb.H1r,aa*bb.H2r,
	    bb.mu);
}

// changing scale

void Li::changescale(const double newmu){
  const double g1 = 3./32.;
  const double g2 = 3./16.;
  const double g3 = 0.;
  const double g4 = 1./8.;
  const double g5 = 3./8.;
  const double g6 = 11./144.;
  const double g7 = 0.;
  const double g8 = 5./48.;
  const double g9 = 1./4.;
  const double g10=-1./4.; 
  const double g11=-1./8.; 
  const double g12= 5./24.; 

  double logmm = log(mu/newmu);

  L1r =  L1r+pi16*g1*logmm;
  L2r =  L2r+pi16*g2*logmm;
  L3r =  L3r+pi16*g3*logmm;
  L4r =  L4r+pi16*g4*logmm;
  L5r =  L5r+pi16*g5*logmm;
  L6r =  L6r+pi16*g6*logmm;
  L7r =  L7r+pi16*g7*logmm;
  L8r =  L8r+pi16*g8*logmm;
  L9r =  L9r+pi16*g9*logmm;
  L10r=L10r+pi16*g10*logmm;
  H1r = H1r+pi16*g11*logmm;
  H2r = H2r+pi16*g12*logmm;
  mu = newmu;
}

void Li::out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & H1t,double & H2t,
	     double & mut, std::string nameout) const{
  L1t = L1r;
  L2t = L2r;
  L3t = L3r;
  L4t = L4r;
  L5t = L5r;
  L6t = L6r;
  L7t = L7r;
  L8t = L8r;
  L9t = L9r;
  L10t= L10r;
  H1t = H1r;
  H2t = H2r;
  mut =mu;
  nameout = name;
}

void Li::out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & H1t,double & H2t,
	     double & mut) const{
  L1t = L1r;
  L2t = L2r;
  L3t = L3r;
  L4t = L4r;
  L5t = L5r;
  L6t = L6r;
  L7t = L7r;
  L8t = L8r;
  L9t = L9r;
  L10t= L10r;
  H1t = H1r;
  H2t = H2r;
  mut =mu;
}

void Li::out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	       double & L10t,double & H1t,double & H2t) const{
  L1t = L1r;
  L2t = L2r;
  L3t = L3r;
  L4t = L4r;
  L5t = L5r;
  L6t = L6r;
  L7t = L7r;
  L8t = L8r;
  L9t = L9r;
  L10t= L10r;
  H1t = H1r;
  H2t = H2r;
}

void Li::out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	       double & L10t) const{
  L1t = L1r;
  L2t = L2r;
  L3t = L3r;
  L4t = L4r;
  L5t = L5r;
  L6t = L6r;
  L7t = L7r;
  L8t = L8r;
  L9t = L9r;
  L10t= L10r;
}

void Li::out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	       double & L10t, double & mut) const{
  L1t = L1r;
  L2t = L2r;
  L3t = L3r;
  L4t = L4r;
  L5t = L5r;
  L6t = L6r;
  L7t = L7r;
  L8t = L8r;
  L9t = L9r;
  L10t= L10r;
  mut =mu;
}

double Li::out(const int n) const{
  switch(n){
  case 1:  return L1r ;
  case 2:  return L2r ;
  case 3:  return L3r ;
  case 4:  return L4r ;
  case 5:  return L5r ;
  case 6:  return L6r ;
  case 7:  return L7r ;
  case 8:  return L8r ;
  case 9:  return L9r ;
  case 10: return L10r;
  case 11: return  H1r;
  case 12: return  H2r;  
  default:
      std::cout << "Trying to output Li number "<<n<<" in set "<<name<<'\n';
  }
  return 0.;
}

Li Lirandom(void){
  double Lir[13];
  for(int i=1;i<=12;i++){
     Lir[i] = pi16*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Li(Lir[1],Lir[2],Lir[3],Lir[4],Lir[5],Lir[6],Lir[7],Lir[8],
	   Lir[9],Lir[10],Lir[11],Lir[12],mu,"Li_random");
}

Li LirandomlargeNc(void){
  double Lir[13];
  const int largeNc[10]={1,2,3,5,8,9,10,11,12};
  for(int i=0;i<13;++i){
    Lir[i] = 0;
  }
  for(int i=0;i<9;++i){
     Lir[largeNc[i]] = pi16*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Li(Lir[1],Lir[2],Lir[3],Lir[4],Lir[5],Lir[6],Lir[7],Lir[8],
	   Lir[9],Lir[10],Lir[11],Lir[12],mu,"Li_random large Nc");
}

Li LirandomlargeNc2(void){
  double Lir[13];
  const int largeNc[10]={1,2,3,5,8,9,10,11,12};
  for(int i=1;i<13;++i){
    Lir[i] = 1./3.*pi16*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
  }
  for(int i=0;i<9;++i){
    Lir[largeNc[i]] *= 3.;
 }
 const double mu = 0.77;
 return Li(Lir[1],Lir[2],Lir[3],Lir[4],Lir[5],Lir[6],Lir[7],Lir[8],
	   Lir[9],Lir[10],Lir[11],Lir[12],mu,"Li_random large Nc 2");
}

