// Linf.cc is part of the CHIRON ChPT at two loops program collection
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


#include "Linf.h"

const double pi16 = 1./(16.*M_PI*M_PI);

Linf::Linf(const double l0r,const double l1r,const double l2r,
	   const double l3r,const double l4r,const double l5r,
	   const double l6r,const double l7r,const double l8r,
	   const double l9r,const double l10r,const double l11r,
	   const double h1r,const double h2r,
	   const double muin,
	   const std::string Name,
	   const int nfin){
    L0r = l0r;
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
    L11r = l11r;
    H1r = h1r;
    H2r = h2r;
    mu = muin;
    name = Name;
    nf = nfin;
}

Linf::~Linf(void){}

std::ostream & operator<<(std::ostream & os, const Linf & bb){
  os   << std::setprecision(15)
       << "# Linf set : " << bb.name << "\n"
       << "#nf : "<<bb.nf <<"\n"
       << "#mu : "<<bb.mu <<"\n"
       << "#L0r  L1r  L2r :" <<bb.L0r <<' '<<bb.L1r<<' '<<bb.L2r<<"\n"
       << "#L3r  L4r  L5r :" <<bb.L3r <<' '<<bb.L4r<<' '<<bb.L5r<<"\n"
       << "#L6r  L7r  L8r :" <<bb.L6r <<' '<<bb.L7r<<' '<<bb.L8r<<"\n"
       << "#L9r  L10r L11r:" <<bb.L9r<<' '<<bb.L10r<<' '<<bb.L11r<<"\n"
       << "#H1r  H2r      :" <<bb.H1r<<' '<<bb.H2r<<"\n";
  os.unsetf(std::ios_base::floatfield);
  os << std::setprecision(6);
  return os;
}

std::istream & operator>>(std::istream & is, Linf & Linfout){
  std::string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  std::string temps2;
  double l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,l11r,h1r,h2r,mu;
  int nf;
  std::string name;
  std::getline(is,temps2,':');// reads in characters up to ':'
  std::getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# Linf set ") std::cout << "trouble reading in Linf\n";
  std::getline(is,temps2,':');   is >> nf;  getline(is,temps);
  if(temps2 != "#nf ") std::cout<<"trouble reading in nf, number of flavours\n";
  std::getline(is,temps2,':');   is >> mu;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l0r >> l1r >> l2r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l3r >> l4r >> l5r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l6r >> l7r >> l8r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> l9r >> l10r >> l11r;  getline(is,temps);
  std::getline(is,temps2,':');   is >> h1r >> h2r;  getline(is,temps);
  Linfout = Linf(l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,l11r,h1r,h2r,mu,name);
  return is;
}


void Linf::setlinf(const int n, const double lin){
  switch(n){
  case 0: L0r = lin; break;
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
  case 11: L11r = lin; break;  
  case 12: H1r = lin; break;  
  case 13: H2r = lin; break;  
  default: std::cout << "Trying to set Linf number "<<n<<" in set "<<name<<'\n';
  }
}

void Linf::setlinf(const double lin,const int n){
  setlinf(n,lin);
}

void Linf::setnf(const int nfin){
  nf = nfin;
}

void Linf::setmu(const double muin){
  mu = muin;
}

void Linf::setname(const std::string inputname){
  name = inputname;
}

// summing two sets of Linf
Linf Linf::operator+(const Linf & bb) const{
  double l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l9r,l8r,l10r,l11r,h1r,h2r,muout;
  l0r = L0r+bb.L0r;
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
  l11r = L11r+bb.L11r;
  h1r = H1r+bb.H1r;
  h2r = H2r+bb.H2r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: summing Linf \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return Linf(l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,l11r,h1r,h2r,muout);
}
// difference of two sets of Linf
Linf Linf::operator-(const Linf & bb) const{
  double l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l9r,l8r,l10r,l11r,h1r,h2r,muout;
  l0r = L0r-bb.L0r;
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
  l11r = L11r-bb.L11r;
  h1r = H1r-bb.H1r;
  h2r = H2r-bb.H2r;
  muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
    std::cout << "WARNING: difference of Linf \""<< name << "\" and \"" << bb.name 
	   << "\" with different mu " << "\n";}
  return Linf(l0r,l1r,l2r,l3r,l4r,l5r,l6r,l7r,l8r,l9r,l10r,l11r,h1r,h2r,muout);
}

// multiplying an Linf with a number

Linf Linf::operator*(const double & aa) const{
  return Linf(aa*L0r,aa*L1r,aa*L2r,aa*L3r,aa*L4r,aa*L5r,aa*L6r,aa*L7r,aa*L8r,
	      aa*L9r,aa*L10r,aa*L11r,aa*H1r,aa*H2r,mu);
}

Linf operator*(const double & aa,const Linf & bb){
  return Linf(aa*bb.L0r,aa*bb.L1r,aa*bb.L2r,aa*bb.L3r,aa*bb.L4r,aa*bb.L5r,
	      aa*bb.L6r,aa*bb.L7r,aa*bb.L8r,aa*bb.L9r,aa*bb.L10r,aa*bb.L11r,
	      aa*bb.H1r,aa*bb.H2r,bb.mu);
}

// changing scale

void Linf::changescale(const double newmu){
  // note -2 of what is in the setLi.hf files
  double nq = double(nf);
  const double g0 = nq/48.;
  const double g1 = 1./16.;
  const double g2 = 1./8.;
  const double g3 = nq/24.;
  const double g4 = 1./8.;
  const double g5 = nq/8.;
  const double g6 = (nq*nq+2.)/(16.*nq*nq);
  const double g7 = 0.;
  const double g8 = (nq*nq-4.)/(16.*nq);
  const double g9 = nq/12.;
  const double g10= -nq/12.;
  const double g11=  0.;
  const double g12= -nq/24.; 
  const double g13= (nq*nq-4.)/(16.*nq);

  double logmm = log(mu/newmu);

  L0r =  L0r+pi16*g0*logmm;
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
  L11r=L11r+pi16*g11*logmm;
  H1r = H1r+pi16*g12*logmm;
  H2r = H2r+pi16*g13*logmm;
  mu = newmu;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & L11t,double & H1t,double & H2t,
	       double & mut, std::string nameout, int & nft) const{
  L0t = L0r;
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
  L11t= L11r;
  H1t = H1r;
  H2t = H2r;
  mut =mu;
  nameout = name;
  nft = nf;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & L11t,double & H1t,double & H2t,
	       double & mut, int & nft) const{
  L0t = L0r;
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
  L11t= L11r;
  H1t = H1r;
  H2t = H2r;
  mut =mu;
  nft = nf;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	       double & L4t,double & L5t,double & L6t,
	       double & L7t,double & L8t,double & L9t,
	       double & L10t,double & L11t,double & H1t,double & H2t,
	       int & nft) const{
  L0t = L0r;
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
  L11t= L11r;
  H1t = H1r;
  H2t = H2r;
  nft = nf;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	       double & L4t,double & L5t,double & L6t,
	       double & L7t,double & L8t,double & L9t,
	       double & L10t,double & L11t,int & nft) const{
  L0t = L0r;
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
  L11t= L11r;
  nft = nf;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	       double & L4t,double & L5t,double & L6t,
	       double & L7t,double & L8t,double & L9t,
	       double & L10t,double & L11t) const{
  L0t = L0r;
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
  L11t= L11r;
}

void Linf::out(double & L0t,double & L1t,double & L2t,double & L3t,
	       double & L4t,double & L5t,double & L6t,
	       double & L7t,double & L8t,double & L9t,
	       double & L10t, double & L11t,double & mut,
	       int & nft) const{
  L0t = L0r;
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
  L11t= L11r;
  mut =mu;
  nft = nf;
}

double Linf::out(const int n) const{
  switch(n){
  case 0:  return L0r ;
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
  case 11: return L11r;
  case 12: return  H1r;
  case 13: return  H2r;  
  default:
      std::cout << "Trying to output Linf number "<<n<<" in set "<<name<<'\n';
  }
  return 0.;
}

int Linf::getnf(void) const{
  return nf;
}

Linf Linfrandom(void){
  double Linfr[14];
  for(int i=0;i<=13;i++){
     Linfr[i] = pi16*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Linf(Linfr[0],Linfr[1],Linfr[2],Linfr[3],Linfr[4],Linfr[5],
	     Linfr[6],Linfr[7],Linfr[8],Linfr[9],Linfr[10],Linfr[11],
	     Linfr[12],Linfr[13],mu,"Linf_random");
}

