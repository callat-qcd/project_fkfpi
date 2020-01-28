// Ci.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// implementation of the Ci class and associated functions

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>  // for rand


#include "Ci.h"
const double pi16 = 1./(16.*M_PI*M_PI);
const double pi162 = pi16*pi16;

Ci::Ci(const double muin,const std::string Name){
    for(int i=0;i<=94;i++){
	Cr[i] = 0.;}
    mu = muin;
    name = Name;
}

Ci::Ci(const double Crin[95],const double muin,const std::string Name){
    for(int i=0;i<=94;i++){
	Cr[i] = Crin[i];}
    mu = muin;
    name = Name;
}

// Ci with resonances constructor: see further down
// Ci Chinese nonlocal quark model:    "
// Ci random
// Ci random large Nc (suppressed zero)
// Ci random large Nc2 (suppressed 1/3 smaller)
Ci::~Ci(void){}

std::ostream & operator<<(std::ostream & os, const Ci & bb){
  os   
    << std::setprecision(7) << std::scientific << std::showpos
    << "# Ci set : " << bb.name << "\n"
    << "# mu  : "<<bb.mu <<"\n"
    <<"# C1r  C2r  C3r  C4r :"
    <<bb.Cr[1]<<' '<<bb.Cr[2]<<' '<<bb.Cr[3]<<' '<<bb.Cr[4]<<' '<<'\n'
    <<"# C5r  C6r  C7r  C8r :"
    <<bb.Cr[5]<<' '<<bb.Cr[6]<<' '<<bb.Cr[7]<<' '<<bb.Cr[8]<<' '<<'\n'
    <<"# C9r  C10r C11r C12r:"
    <<bb.Cr[9]<<' '<<bb.Cr[10]<<' '<<bb.Cr[11]<<' '<<bb.Cr[12]<<' '<<'\n'
    <<"# C13r C14r C15r C16r:"
    <<bb.Cr[13]<<' '<<bb.Cr[14]<<' '<<bb.Cr[15]<<' '<<bb.Cr[16]<<' '<<'\n'
    <<"# C17r C18r C19r C20r:"
    <<bb.Cr[17]<<' '<<bb.Cr[18]<<' '<<bb.Cr[19]<<' '<<bb.Cr[20]<<' '<<'\n'
    <<"# C21r C22r C23r C24r:"
    <<bb.Cr[21]<<' '<<bb.Cr[22]<<' '<<bb.Cr[23]<<' '<<bb.Cr[24]<<' '<<'\n'
    <<"# C25r C26r C27r C28r:"
    <<bb.Cr[25]<<' '<<bb.Cr[26]<<' '<<bb.Cr[27]<<' '<<bb.Cr[28]<<' '<<'\n'
    <<"# C29r C30r C31r C32r:"
    <<bb.Cr[29]<<' '<<bb.Cr[30]<<' '<<bb.Cr[31]<<' '<<bb.Cr[32]<<' '<<'\n'
    <<"# C33r C34r C35r C36r:"
    <<bb.Cr[33]<<' '<<bb.Cr[34]<<' '<<bb.Cr[35]<<' '<<bb.Cr[36]<<' '<<'\n'
    <<"# C37r C38r C39r C40r:"
    <<bb.Cr[37]<<' '<<bb.Cr[38]<<' '<<bb.Cr[39]<<' '<<bb.Cr[40]<<' '<<'\n'
    <<"# C41r C42r C43r C44r:"
    <<bb.Cr[41]<<' '<<bb.Cr[42]<<' '<<bb.Cr[43]<<' '<<bb.Cr[44]<<' '<<'\n'
    <<"# C45r C46r C47r C48r:"
    <<bb.Cr[45]<<' '<<bb.Cr[46]<<' '<<bb.Cr[47]<<' '<<bb.Cr[48]<<' '<<'\n'
    <<"# C49r C50r C51r C52r:"
    <<bb.Cr[49]<<' '<<bb.Cr[50]<<' '<<bb.Cr[51]<<' '<<bb.Cr[52]<<' '<<'\n'
    <<"# C53r C54r C55r C56r:"
    <<bb.Cr[53]<<' '<<bb.Cr[54]<<' '<<bb.Cr[55]<<' '<<bb.Cr[56]<<' '<<'\n'
    <<"# C57r C58r C59r C60r:"
    <<bb.Cr[57]<<' '<<bb.Cr[58]<<' '<<bb.Cr[59]<<' '<<bb.Cr[60]<<' '<<'\n'
    <<"# C61r C62r C63r C64r:"
    <<bb.Cr[61]<<' '<<bb.Cr[62]<<' '<<bb.Cr[63]<<' '<<bb.Cr[64]<<' '<<'\n'
    <<"# C65r C66r C67r C68r:"
    <<bb.Cr[65]<<' '<<bb.Cr[66]<<' '<<bb.Cr[67]<<' '<<bb.Cr[68]<<' '<<'\n'
    <<"# C69r C70r C71r C72r:"
    <<bb.Cr[69]<<' '<<bb.Cr[70]<<' '<<bb.Cr[71]<<' '<<bb.Cr[72]<<' '<<'\n'
    <<"# C73r C74r C75r C76r:"
    <<bb.Cr[73]<<' '<<bb.Cr[74]<<' '<<bb.Cr[75]<<' '<<bb.Cr[76]<<' '<<'\n'
    <<"# C77r C78r C79r C80r:"
    <<bb.Cr[77]<<' '<<bb.Cr[78]<<' '<<bb.Cr[79]<<' '<<bb.Cr[80]<<' '<<'\n'
    <<"# C81r C82r C83r C84r:"
    <<bb.Cr[81]<<' '<<bb.Cr[82]<<' '<<bb.Cr[83]<<' '<<bb.Cr[84]<<' '<<'\n'
    <<"# C85r C86r C87r C88r:"
    <<bb.Cr[85]<<' '<<bb.Cr[86]<<' '<<bb.Cr[87]<<' '<<bb.Cr[88]<<' '<<'\n'
    <<"# C89r C90r C91r C92r:"
    <<bb.Cr[89]<<' '<<bb.Cr[90]<<' '<<bb.Cr[91]<<' '<<bb.Cr[92]<<' '<<'\n'
    <<"# C93r C94r          :"
    <<bb.Cr[93]<<' '<<bb.Cr[94]<<' '<<'\n';
  os.unsetf(std::ios_base::floatfield);
  os << std::setprecision(6);  
  return os;
}

std::istream & operator>>(std::istream & is, Ci & Ciout){
  std::string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  std::string temps2;
  double CC[95],mu;
  std::string name;
  getline(is,temps2,':');// reads in characters up to ':'
  getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# Ci set ") std::cout << "trouble reading in Ci\n";
  getline(is,temps2,':');   is >> mu;  getline(is,temps);
  for (int i = 0; i<92; i += 4){
    getline(is,temps2,':');   is >> CC[i+1] >> CC[i+2] >> CC[i+3] >> CC[i+4];
    getline(is,temps);
  }
  getline(is,temps2,':');   is >> CC[93] >> CC[94];
  getline(is,temps);
  Ciout = Ci(CC,mu,name);
  return is;
}

void Ci::setci(const int n, const double cin){
    Cr[n] = cin;
}

void Ci::setci(const double cin,const int n){
    Cr[n] = cin;
}

void Ci::setmu(const double muin){
  mu = muin;
}

void Ci::setname(const std::string inputname){
  name = inputname;
}

// summing two sets of Ci
Ci Ci::operator+(const Ci & bb) const{
  double Crout[95];
  for(int i=0;i<=94;i++){
      Crout[i] = Cr[i]+bb.Cr[i];}
  double muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
      std::cout << "WARNING: summing Ci "<< name << " and " << bb.name 
	   << " with different mu " << "\n";}
  return Ci(Crout,muout);
}
// difference
Ci Ci::operator-(const Ci & bb) const{
  double Crout[95];
  for(int i=0;i<=94;i++){
      Crout[i] = Cr[i]-bb.Cr[i];}
  double muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
      std::cout << "WARNING: difference of Ci "<< name << " and " << bb.name 
	   << " with different mu " << "\n";}
  return Ci(Crout,muout);
}

// multiplying a Ci with a number

Ci Ci::operator*(const double & aa) const{
    double Crin[95];
    for(int i=0;i<=94;i++){
	Crin[i] = aa*Cr[i];}
  return Ci(Crin,mu);
}

Ci operator*(const double & aa,const Ci & bb){
    double Crin[95];
    for(int i=0;i<=94;i++){
	Crin[i] = aa*bb.Cr[i];}
  return Ci(Crin,bb.mu);
}
// changing scale
void Ci::changescale(const double newmu) const{
  std::cout << "WARNING: attempt to change scale to : " << newmu<<'\n';
  std::cout << "WARNING:makes no sense to change scale in Ci without Li specified"
	 <<'\n';
}

void Ci::changescale(Li & li,const double newmu){
    changescale(newmu,li);
}

void Ci::changescale(const double newmu,Li & li){
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mut;
  std::string namet;
  li.out(L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r,mut,namet);
  if(fabs(1.-mu/mut)>=1e-6){
      std::cout << "ERROR: cannot change scale when Ci and Li different scale"
	   << '\n' << " using Li " << namet << " and Ci " << name
	   <<'\n';
      assert(0);
  }
// first the Ci, then the Li, see notes
  double logm = log(mu/newmu);
  double    logmm = -4.*log(mu/newmu);
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
  //const double g11=-1./8.; 
  //const double g12= 5./24.; 

// do this BEFORE resetting the lir
  double L1t =  pi16*((-2.)*logm*L1r -pi16*g1*logm*logm);
  double L2t =  pi16*((-2.)*logm*L2r -pi16*g2*logm*logm);
  double L3t =  pi16*((-2.)*logm*L3r -pi16*g3*logm*logm);
  double L4t =  pi16*((-2.)*logm*L4r -pi16*g4*logm*logm);
  double L5t =  pi16*((-2.)*logm*L5r -pi16*g5*logm*logm);
  double L6t =  pi16*((-2.)*logm*L6r -pi16*g6*logm*logm);
  double L7t =  pi16*((-2.)*logm*L7r -pi16*g7*logm*logm);
  double L8t =  pi16*((-2.)*logm*L8r -pi16*g8*logm*logm);
  double L9t =  pi16*((-2.)*logm*L9r -pi16*g9*logm*logm);
  double L10t=  pi16*((-2.)*logm*L10r-pi16*g10*logm*logm);
  li.changescale(newmu);

  double La1 = pi16*pi16*logmm;

  Cr[1] = Cr[1]
      + 7./6.*L1t + 35./12.*L2t + 9./8.*L3t - 1./2.*L9t + 71./
      1536.*La1;
  Cr[2] = Cr[2]
      + 3./2.*L1t + 1./2.*L2t + 5./12.*L3t + 1./256.*La1;
  Cr[3] = Cr[3]
      - 1./2.*L1t - 5./4.*L2t - 1./4.*L3t + 1./8.*L9t - 43./1536.
      *La1;
  Cr[4] = Cr[4]
      + 2./3.*L1t + 5./3.*L2t + 3./8.*L3t - 3./8.*L9t + 79./1536.
      *La1;
  Cr[5] = Cr[5]
      - 4./3.*L1t - 19./6.*L2t - 59./72.*L3t + L4t - 3./4.*L5t
      - 85./768.*La1;
  Cr[6] = Cr[6]
      - 5./3.*L1t - 8./9.*L2t - 65./108.*L3t - 2*L4t - 1./4.*L5t
      + 89./1152.*La1;
  Cr[7] = Cr[7]
      + 3*L1t + 3*L2t + 3./2.*L3t - 3./2.*L4t - 1./4.*L5t + 9./
      64.*La1;
  Cr[8] = Cr[8]
      - 1./3.*L1t - 13./6.*L2t - 23./72.*L3t + 1./2.*L4t - 3./4.
      *L5t + 1./24.*La1;
  Cr[9] = Cr[9]
      + L1t + 5./2.*L2t - 1./2.*L4t + 3./128.*La1;
  Cr[10] = Cr[10]
      - 4./3.*L1t - 13./6.*L2t + 1./72.*L3t + 3./4.*L5t - 41./
      384.*La1;
  Cr[11] = Cr[11]
      + 8./3.*L1t + 7./18.*L2t + 157./216.*L3t + 3./2.*L4t + 1./
      4.*L5t + 19./1152.*La1;
  Cr[12] = Cr[12]
      + 1./3.*L1t + 2./3.*L2t + 4./9.*L3t + 3./8.*L5t - 13./
      1536.*La1;
  Cr[13] = Cr[13]
      + 4./3.*L1t + 4./9.*L2t + 23./54.*L3t + 3./4.*L4t + 1./8.
      *L5t - 1./4608.*La1;
  Cr[14] = Cr[14]
      - 79./36.*L1t + 1./72.*L2t - 37./144.*L3t + 2./9.*L4t + 2
      *L6t - L7t - 3./4.*L8t - 479./9216.*La1;
  Cr[15] = Cr[15]
      + 22./9.*L1t + 11./18.*L2t + 53./72.*L3t + 19./36.*L4t -
      1./18.*L5t - 7./2.*L6t - 1./2.*L8t + 31./384.*La1;
  Cr[16] = Cr[16]
      + 13./9.*L1t + 13./36.*L2t + 61./144.*L3t - 17./36.*L4t
      + 1./8.*L5t - 1./2.*L6t - 1./4.*L8t + 43./1536.*La1;
  Cr[17] = Cr[17]
      - 35./36.*L1t + 23./72.*L2t + 1./36.*L3t + 1./9.*L4t - 1./
      24.*L5t + L6t + L7t + 55./9216.*La1;
  Cr[18] = Cr[18]
      + 19./18.*L1t - 1./9.*L2t + 1./72.*L3t - 1./9.*L4t - 1./4.
      *L5t - L6t + 1./2.*L8t + 235./20736.*La1;
  Cr[19] = Cr[19]
      + 22./81.*L1t + 4./81.*L2t + 5./81.*L3t + 1./6.*L4t - 1./
      8.*L5t + 2./27.*L6t - 16./27.*L7t - 1./12.*L8t + 1517./186624.
      *La1;
  Cr[20] = Cr[20]
      - 11./27.*L1t - 2./27.*L2t - 5./54.*L3t + 1./6.*L4t + 13./
      24.*L5t - 17./18.*L6t + 1./18.*L7t - 31./36.*L8t - 1517./
      124416.*La1;
  Cr[21] = Cr[21]
      + 11./81.*L1t + 2./81.*L2t + 5./162.*L3t + 7./18.*L4t - 1.
      /54.*L5t - 31./54.*L6t + 5./54.*L7t + 1517./373248.*La1;
  Cr[22] = Cr[22]
      - 8./3.*L1t - 7./6.*L2t - L3t - 1./3.*L4t + 1./8.*L5t - 5.
      /288.*La1;
  Cr[23] = Cr[23]
      + 3./2.*L1t + 3./4.*L2t + 1./2.*L3t - 1./4.*L5t + 15./512.
      *La1;
  Cr[24] = Cr[24]
      + 31./9.*L1t + 47./18.*L2t + 3./2.*L3t + 1./4.*L5t - 1./4.
      *L9t + 49./576.*La1;
  Cr[25] = Cr[25]
      - 5*L1t - 11./2.*L2t - 5./2.*L3t + 2./3.*L4t + 1./2.*L5t
      + 3./4.*L9t - 173./1152.*La1;
  Cr[26] = Cr[26]
      + 191./54.*L1t + 379./108.*L2t + 15./8.*L3t - 2./3.*L4t
      + 5./12.*L5t - L6t + 2*L7t - 3./4.*L8t - 1./3.*L9t + 1397./
      13824.*La1;
  Cr[27] = Cr[27]
      - 67./27.*L1t - 83./54.*L2t - 19./18.*L3t + 1./3.*L4t - 1.
      /2.*L5t - 7./2.*L7t - 1./2.*L8t + 1./12.*L9t - 67./1296.*La1;
  Cr[28] = Cr[28]
      + 29./27.*L1t + 14./27.*L2t + 3./8.*L3t - 1./2.*L7t - 1./
      4.*L8t - 1./24.*L9t + 1./13824.*La1;
  Cr[29] = Cr[29]
      + 37./54.*L1t - 49./108.*L2t + 1./6.*L3t - 1./3.*L4t + 2./
      3.*L5t + L6t + L7t + 5./24.*L9t - 59./3456.*La1;
  Cr[30] = Cr[30]
      + 49./27.*L1t + 37./27.*L2t + 8./9.*L3t + 17./6.*L4t + 7./
      36.*L5t - L7t + 1./2.*L8t - 1./12.*L9t + 91./6912.*La1;
  Cr[31] = Cr[31]
      + 31./27.*L1t + 22./27.*L2t + 17./27.*L3t + 1./2.*L4t + 3.
      /8.*L5t - 7./9.*L6t + 11./9.*L7t - 1./12.*L8t + 2359./124416.*
      La1;
  Cr[32] = Cr[32]
      + 25./27.*L1t + 10./27.*L2t + 1./3.*L3t + 1./2.*L4t + 1./
      8.*L5t - 1./9.*L6t - 11./18.*L7t - 13./36.*L8t - 851./124416.*
      La1;
  Cr[33] = Cr[33]
      - 28./27.*L1t - 16./27.*L2t - 13./27.*L3t - 1./2.*L4t - 1.
      /4.*L5t + 4./9.*L6t - 37./18.*L7t - 1./2.*L8t - 451./31104.*
      La1;
  Cr[34] = Cr[34]

      - 1./2.*L1t - 3./4.*L2t - 11./24.*L3t + 1./2.*L4t - 3./4.
      *L5t - 31./4608.*La1;
  Cr[35] = Cr[35]
      + 1./3.*L1t + 1./2.*L2t + 11./36.*L3t - 1./3.*L4t + 1./2.
      *L5t + 43./20736.*La1;
  Cr[36] = Cr[36]
      - 8./3.*L1t - 4./3.*L2t - 35./36.*L3t - 25./6.*L4t - 1./2.
      *L5t - 1./24.*La1;
  Cr[37] = Cr[37]
      + 173./41472.*La1;
  Cr[38] = Cr[38]
      - 5./12.*L5t + 49./4608.*La1;
  Cr[39] = Cr[39]
      - 4./3.*L4t - 11./36.*L5t;
  Cr[40] = Cr[40]
      - 1./3.*L1t - 11./2.*L2t - 31./24.*L3t + 1./2.*L9t - 379./
      2304.*La1;
  Cr[41] = Cr[41]
      - 3*L1t + 5./4.*L2t - 1./3.*L3t + 71./512.*La1;
  Cr[42] = Cr[42]
      - 13./3.*L2t - 11./16.*L3t + 1./4.*L9t - 233./2304.*La1;
  Cr[43] = Cr[43]
      + 9./2.*L2t + 27./256.*La1;
  Cr[44] = Cr[44]
      - 1./3.*L1t - 1./2.*L2t + 5./8.*L3t - L9t + 53./2304.*La1
      ;
  Cr[45] = Cr[45]
      + 3*L1t + 17./8.*L2t + 13./12.*L3t + 47./1024.*La1;
  Cr[46] = Cr[46]
      + 1./3.*L1t + 7./6.*L2t - 7./48.*L3t - 1./4.*L9t + 281./
      2304.*La1;
  Cr[47] = Cr[47]
      + 1./3.*L1t - 13./3.*L2t - 3./16.*L3t + 1./2.*L9t - 451./
      2304.*La1;
  Cr[48] = Cr[48]
      - 1./6.*L1t - 11./12.*L2t + 11./24.*L3t + 11./16.*L9t -
      65./384.*La1;
  Cr[49] = Cr[49]
      - 6*L1t - 7./2.*L2t - 13./6.*L3t + 5./256.*La1;
  Cr[50] = Cr[50]
      - 1./3.*L1t - 5./2.*L2t + 2./3.*L3t - 1./3.*L4t - 1./4.*
      L5t + 3./8.*L9t - 107./576.*La1;
  Cr[51] = Cr[51]
      - 2./3.*L1t - 13./3.*L2t - 17./12.*L3t - 1./3.*L4t - 1./4.
      *L5t + 7./4.*L9t - 43./288.*La1;
  Cr[52] = Cr[52]
      - 2./3.*L1t + 4*L2t + 5./12.*L3t + 1./3.*L4t + 1./4.*L5t
      + 1./8.*L9t + 197./2304.*La1;
  Cr[53] = Cr[53]
      - 1./3.*L1t - 1./12.*L3t + 1./6.*L4t + 1./8.*L5t - 7./16.
      *L9t - 7./8.*L10t - 7./1152.*La1;
  Cr[54] = Cr[54]
      - 1./4.*L2t - 1./12.*L3t - 7./512.*La1;
  Cr[55] = Cr[55]
      + 1./3.*L1t + 5./24.*L3t - 1./6.*L4t - 1./8.*L5t - 11./16.
      *L9t - 1./4.*L10t + 43./1152.*La1;
  Cr[56] = Cr[56]
      + 1./3.*L2t - 1./3.*L3t + 1./3.*L4t + 1./4.*L5t - 5./8.*
      L9t + 1./12.*La1;
  Cr[57] = Cr[57]
      + L1t - 1./6.*L2t + 5./12.*L3t - 2./3.*L4t - 1./2.*L5t -
      1./2.*L9t + 7./64.*La1;
  Cr[58] = Cr[58]
      + L2t + 1./3.*L3t - 11./128.*La1;
  Cr[59] = Cr[59]
      - 1./2.*L1t - 1./12.*L2t - 7./24.*L3t + 1./6.*L4t + 1./8.
      *L5t + 9./16.*L9t + 5./96.*La1;
  Cr[60] = Cr[60]
      - 3./64.*La1;
  Cr[61] = Cr[61]
      - 3./8.*L9t - 3./8.*L10t;
  Cr[62] = Cr[62]
      - 1./8.*L9t - 1./8.*L10t;
  Cr[63] = Cr[63]
      + 1./6.*L1t - 11./12.*L2t - 5./72.*L3t - 3./8.*L5t + 3./
      16.*L9t + 5./768.*La1;
  Cr[64] = Cr[64]
      - 16./3.*L1t - 16./9.*L2t - 157./108.*L3t - 3*L4t - 1./2.
      *L5t + 1./4.*L9t - 5./288.*La1;
  Cr[65] = Cr[65]
      + 1./3.*L1t - 11./6.*L2t - 8./9.*L3t - 3./4.*L5t + 3./8.*
      L9t + 5./384.*La1;
  Cr[66] = Cr[66]
      + 1./3.*L1t + 1./2.*L2t - 5./12.*L3t - 5./8.*L9t + 115./
      2304.*La1;
  Cr[67] = Cr[67];
  Cr[68] = Cr[68]
      - 6*L1t - L2t - 4./3.*L3t - 7./128.*La1;
  Cr[69] = Cr[69]
      - 1./3.*L1t - 1./2.*L2t - 5./6.*L3t - 1./2.*L9t + 29./
      2304.*La1;
  Cr[70] = Cr[70]
      - 4./3.*L1t - 17./6.*L2t - 5./4.*L3t + 1./6.*L4t + 1./8.*
      L5t + 5./16.*L9t + 7./8.*L10t - 389./4608.*La1;
  Cr[71] = Cr[71]
      - 3./2.*L1t - 1./2.*L2t - 5./12.*L3t - 7./256.*La1;
  Cr[72] = Cr[72]
      + 2./3.*L1t + 7./6.*L2t + 3./8.*L3t - 1./6.*L4t - 1./8.*
      L5t + 1./16.*L9t + 1./4.*L10t + 251./4608.*La1;
  Cr[73] = Cr[73]
      + 5./3.*L1t + 9./2.*L2t + 11./12.*L3t - 1./3.*L4t - 1./4.
      *L5t - 7./8.*L9t + 385./2304.*La1;
  Cr[74] = Cr[74]
      - 4./3.*L1t - 10./3.*L2t - 2*L3t - 1./2.*L9t - 13./256.*
      La1;
  Cr[75] = Cr[75];
  Cr[76] = Cr[76]
      - 5./6.*L1t - 9./4.*L2t - 13./12.*L3t + 1./6.*L4t + 1./8.
      *L5t - 1./16.*L9t - 217./4608.*La1;
  Cr[77] = Cr[77]
      - 6*L1t - L2t - 4./3.*L3t - 7./128.*La1;
  Cr[78] = Cr[78]
      + 1./6.*L1t - 1./12.*L2t + 1./8.*L3t - 1./6.*L4t - 1./8.*
      L5t - 9./16.*L9t + 49./2304.*La1;
  Cr[79] = Cr[79]
      - 1./6.*L1t + 1./12.*L2t - 1./8.*L3t + 1./6.*L4t + 1./8.*
      L5t - 3./16.*L9t - 11./1152.*La1;
  Cr[80] = Cr[80]
      - 1./3.*L1t - 2./3.*L2t - 4./9.*L3t - 3./8.*L5t + 3./8.*
      L10t - 17./1536.*La1;
  Cr[81] = Cr[81]
      - 4./3.*L1t - 4./9.*L2t - 23./54.*L3t - 3./4.*L4t - 1./8.
      *L5t + 1./8.*L10t - 101./4608.*La1;
  Cr[82] = Cr[82]
      - 1./12.*L1t + 1./24.*L2t - 1./16.*L3t + 1./12.*L4t + 1./
      16.*L5t + 5./32.*L9t + 5./288.*La1;
  Cr[83] = Cr[83]
      + 3./2.*L1t + 11./4.*L2t + 97./72.*L3t - 5./6.*L4t + 1./8.
      *L5t + 359./4608.*La1;
  Cr[84] = Cr[84]
      + 25./3.*L1t + 43./18.*L2t + 133./54.*L3t + 3*L4t + 1./4.
      *L9t + 215./2304.*La1;
  Cr[85] = Cr[85]
      + 1./3.*L1t + 7./6.*L2t + 31./36.*L3t + 3./4.*L5t + 3./8.
      *L9t - 1./384.*La1;
  Cr[86] = Cr[86]
      + 16./3.*L1t + 8./9.*L2t + 79./54.*L3t + 3*L4t + 1./2.*
      L5t + 1./4.*L9t + 5./144.*La1;
  Cr[87] = Cr[87]
      - 1./4.*L9t + 3./256.*La1;
  Cr[88] = Cr[88]
      - 1./3.*L1t + 1./6.*L2t - 1./4.*L3t + 1./3.*L4t + 1./4.*
      L5t + 1./8.*L9t - 31./1152.*La1;
  Cr[89] = Cr[89]
      + 1./3.*L1t - 1./6.*L2t + 1./4.*L3t - 7./8.*L9t + 1./64.*
      La1;
  Cr[90] = Cr[90]
      + 1./3.*L4t + 1./4.*L5t - 37./1152.*La1;
  Cr[91] = Cr[91]
      - 49./1152.*La1;
  Cr[92] = Cr[92]
      - 2*L9t;
  Cr[93] = Cr[93]
      + 1./2.*L9t;
  Cr[94] = Cr[94]
      - 88./27.*L1t - 16./27.*L2t - 20./27.*L3t - 2*L4t - 8./9.
      *L6t - 44./9.*L7t - 1517./15552.*La1;
  mu = newmu;
}

void Ci::out(double Cit[95],double & mut,std::string namet) const{
  for(int i =0;i<=94;i++){
      Cit[i]=Cr[i];}
  mut = mu;
  namet = name;
}

void Ci::out(double Cit[95],double & mut) const{
  for(int i =0;i<=94;i++){
      Cit[i]=Cr[i];}
  mut = mu;
}

void Ci::out(double Cit[95]) const{
  for(int i =0;i<=94;i++){
      Cit[i]=Cr[i];}
}

double Ci::out(const int n) const{
  return Cr[n];
}

// constructor with resonance input

Ci::Ci(const double fpi, const double muin
       ,const double MV, const double fV, const double gV
       ,const double fChi, const double aV
       ,const double MS, const double cd, const double cm, const double l3ss
       ,const double MEP, const double dmt
       ,const double MP, const double dm){
  for(int i =0;i<=94;i++){
    Cr[i] = 0.;
  }
  mu = muin;
  name = "Ciresonances"; 

  double fpi2 = fpi*fpi;

  // singlet pseudoscalar
  double mepfpi = fpi2/pow(MEP,4);
  double dmt2 = dmt*dmt;

  Cr[18] += mepfpi * (  - 1./2.*dmt2 );
  Cr[19] += mepfpi * (  - 1./9.*dmt2 );
  Cr[20] += mepfpi * (    1./6.*dmt2 );
  Cr[21] += mepfpi * (  - 1./18.*dmt2 );
  Cr[27] += mepfpi * (  - dmt2 );
  Cr[31] += mepfpi * (  - 1./3.*dmt2 );
  Cr[32] += mepfpi * (    1./6.*dmt2 );
  Cr[33] += mepfpi * (  - 1./6.*dmt2 );
  Cr[35] += mepfpi * (    dmt2 );
  Cr[37] += mepfpi * (  - 1./2.*dmt2 );
  Cr[94] += mepfpi * (    4./3.*dmt2 );

  // scalars
  double cd2 = cd*cd;
  double cm2 = cm*cm;
  double msfpi = fpi2/pow(MS,4); 
  double l3ssb = cm2*l3ss;
  //double rs = cd/cm;
  Cr[1 ] += msfpi * (  - 1./4.*cd2 );
  Cr[5 ] += msfpi * ( 1./2.*cd*cm +cd2*l3ss);
  Cr[8 ] += msfpi * ( 1./2.*cd*cm );
  Cr[10] += msfpi * (  - cd*cm );
  Cr[12] += msfpi * (  - 1./2.*cd*cm );
  Cr[14] += msfpi * ( 2.*cd*cm*l3ss);
  Cr[19] += msfpi * ( 1./27.*cd*cm +l3ssb);
  Cr[20] += msfpi * (  - 1./18.*cd*cm );
  Cr[21] += msfpi * ( 1./54.*cd*cm );
  Cr[22] += msfpi * ( 1./8.*cd2 );
  Cr[24] += msfpi * (  - 1./6.*cd2 );
  Cr[25] += msfpi * ( 1./4.*cd2 );
  Cr[26] += msfpi * (  - 1./4.*cm2 - 1./2.*cd*cm - 5./36.*cd2 );
  Cr[27] += msfpi * ( 1./3.*cd*cm + 1./18.*cd2 );
  Cr[28] += msfpi * (  - 1./36.*cd2 );
  Cr[29] += msfpi * (  - 1./4.*cm2 - 1./2.*cd*cm + 1./18.*cd2 );
  Cr[30] += msfpi * (  - 1./18.*cd2 );
  Cr[31] += msfpi * (  - 7./18.*cd*cm );
  Cr[32] += msfpi * (  - 1./18.*cd*cm );
  Cr[33] += msfpi * ( 2./9.*cd*cm );
  Cr[34] += msfpi * ( 1./2.*cm2 + 1./2.*cd*cm );
  Cr[35] += msfpi * (  - 1./3.*cd*cm );
  Cr[38] += msfpi * ( 1./2.*cm2 );
  Cr[40] += msfpi * ( 1./4.*cd2 );
  Cr[42] += msfpi * ( 1./4.*cd2 );
  Cr[44] += msfpi * (  - 1./2.*cd2);
  Cr[48] += msfpi * ( 1./4.*cd2 );
  Cr[51] += msfpi * ( 1./2.*cd2 );
  Cr[63] += msfpi * ( 1./2.*cd*cm );
  Cr[65] += msfpi * ( cd*cm );
  Cr[66] += msfpi * ( 1./4.*cd2 );
  Cr[69] += msfpi * ( 1./4.*cd2 );
  Cr[70] += msfpi * ( 1./4.*cd2 );
  Cr[74] += msfpi * ( 1./2.*cd2 );
  Cr[76] += msfpi * ( 1./4.*cd2 );
  Cr[80] += msfpi * ( 1./2.*cd*cm );
  Cr[83] += msfpi * (  - 1./2.*cd*cm - 1./8.*cd2);
  Cr[85] += msfpi * (  - cd*cm );
  Cr[94] += msfpi * (  - 4./9.*cd*cm );


  // pseudoscalar nonet
  double mpfpi = fpi2/pow(MP,4);
  double dm2 = dm*dm;
  Cr[14] += mpfpi* (-1./4.*dm2);
  Cr[17] += mpfpi* (-1./4.*dm2);
  Cr[26] += mpfpi* (-1./2.*dm2);
  Cr[29] += mpfpi* (-1./2.*dm2);
  Cr[31] += mpfpi* (-1./2.*dm2);
  Cr[33] += mpfpi* ( 1./6.*dm2);
  Cr[34] += mpfpi* ( 1./2.*dm2);
  Cr[38] += mpfpi* (-1./2.*dm2);
  Cr[91] += mpfpi* (    2.*dm2);

  // vectors
  double mvfpi = fpi2/pow(MV,2);
  double gV2 = gV*gV;
  double fV2 = fV*fV; 
  double sqrt2 = sqrt(2.);
  Cr[1 ] += mvfpi * ( 1./8.*gV2 );
  Cr[4 ] += mvfpi * ( 1./8.*gV2 );
  Cr[22] += mvfpi * ( 1./16.*gV2 +1./(2.*sqrt2)*gV*fChi);
  Cr[24] += mvfpi * ( 1./12.*gV2 );
  Cr[25] += mvfpi * (  - 3./8.*gV2    -1./sqrt2*gV*fChi);
  Cr[26] += mvfpi * ( 7./36.*gV2    +fChi*fChi+1./sqrt2*gV*fChi);
  Cr[27] += mvfpi * (  - 1./36.*gV2 );
  Cr[28] += mvfpi * ( 1./72.*gV2 );
  Cr[29] += mvfpi * (  - 11./72.*gV2    -fChi*fChi-1./sqrt2*gV*fChi);
  Cr[30] += mvfpi * ( 1./36.*gV2 );
  Cr[40] += mvfpi * (  - 1./8.*gV2 );
  Cr[42] += mvfpi * (  - 1./8.*gV2 );
  Cr[44] += mvfpi * ( 1./4.*gV2 );
  Cr[48] += mvfpi * (  - 1./8.*gV2 );
  Cr[50] += mvfpi * (  + 1./4.*gV*fV    +1./sqrt2*fV*fChi);
  Cr[51] += mvfpi * (  + 1./4.*gV*fV - 1./4.*gV2   +1./sqrt2*fV*fChi);
  Cr[52] += mvfpi * ( - 1./4.*gV*fV    -1./sqrt2*fV*fChi);
  Cr[53] += mvfpi * (  - 3./16.*fV2 - 1./8.*gV*fV  -1./(2.*sqrt2)*fV*fChi  );
  Cr[55] += mvfpi * (3./16.*fV2 + 1./8.*gV*fV  +1./(2.*sqrt2)*fV*fChi  );
  Cr[56] += mvfpi * (3./8.*fV2 - 1./4.*gV*fV  -1./sqrt2*fV*fChi  );
  Cr[57] += mvfpi * (1./8.*fV2 + 1./2.*gV*fV   +2./sqrt2*fV*fChi  );
  Cr[59] += mvfpi * (-1./4.*fV2 - 1./8.*gV*fV  -1./(2.*sqrt2)*fV*fChi  );
  Cr[66] += mvfpi * ( 1./8.*gV2+1./(2.*sqrt2)*gV*aV );
  Cr[69] += mvfpi * (  - 1./8.*gV2 -1./(2.*sqrt2)*gV*aV );
  Cr[70] += mvfpi * ( 1./8.*fV2 - 1./8.*gV*fV-1./8.*gV2-1./(2.*sqrt2)*fV*fChi );
  Cr[72] += mvfpi * (  - 1./8.*fV2 + 1./8.*gV*fV+1./(2.*sqrt2)*fV*fChi );
  Cr[73] += mvfpi * (  - 1./8.*fV2 + 1./4.*gV*fV+1./sqrt2*fV*fChi );
  Cr[74] += mvfpi * (  - 1./4.*gV2 -aV*aV-1./sqrt2*gV*aV );
  Cr[76] += mvfpi * ( 1./16.*fV2 - 1./8.*gV*fV +1./2.*aV*aV+1./(2.*sqrt2)*gV*aV-1./(2.*sqrt2)*fV*fChi );
  Cr[78] += mvfpi * ( 1./4.*fV2 + 1./8.*gV*fV +1./(2.*sqrt2)*fV*fChi );
  Cr[79] += mvfpi * ( 1./8.*fV2 - 1./8.*gV*fV -1./(2.*sqrt2)*fV*fChi );
  Cr[82] += mvfpi * (  - 1./16.*fV2 -1./16.*gV*fV-1./(4.*sqrt2)*fV*fChi );
  Cr[83] += mvfpi * ( 3./16.*gV2+1./(2.*sqrt2)*gV*aV+aV*fChi+1./(2.*sqrt2)*gV*fChi );
  Cr[87] += mvfpi * ( 1./8.*fV2 );
  Cr[88] += mvfpi * ( -1./4.*gV*fV -1./sqrt2*fV*fChi );
  Cr[89] += mvfpi * ( 1./2.*fV2 + 1./4.*gV*fV    +1./sqrt2*fV*aV );
  Cr[90] += mvfpi * ( -1./sqrt2*fV*fChi );
  Cr[92] += mvfpi * ( fV2 );
  Cr[93] += mvfpi * (  - 1./4.*fV2 );
}

Ci Cirandom(void){
 double CC[95];
 CC[0]=0.;
 for(int i=1;i<=94;++i){
     CC[i] = pi162*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Ci(CC,mu,"Ci_random");

}

Ci CirandomlargeNc(void){
  const int largeNc[54]={1,3,4,5,8,10,12,14,17,19,22,25,26,29,31,34,38,40,42,44,
			46,47,48,50,51,52,53,55,56,57,59,61,63,65,66,69,70,
			72,73,74,76,78,79,80,82,83,85,87,88,89,90,91,92,93};
 double CC[95];
 CC[0]=0.;
 for(int i=0;i<=94;++i){
     CC[i] = 0;
 }
 for(int i=0;i<54;++i){
     CC[largeNc[i]] = pi162*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Ci(CC,mu,"Ci_random largeNc");
}

Ci CirandomlargeNc2(void){
  const int largeNc[54]={1,3,4,5,8,10,12,14,17,19,22,25,26,29,31,34,38,40,42,44,
			46,47,48,50,51,52,53,55,56,57,59,61,63,65,66,69,70,
			72,73,74,76,78,79,80,82,83,85,87,88,89,90,91,92,93};
 double CC[95];
 CC[0]=0.;
 for(int i=0;i<=94;++i){
     CC[i] =  1./3.*pi162*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 for(int i=0;i<54;++i){
   CC[largeNc[i]] *= 3.;
 }
 const double mu = 0.77;
 return Ci(CC,mu,"Ci_random largeNc2");
}

// produces a nonlocal quark model from arXiv:0907.5229
//Ci Ci_JZLW(const physmass mass){
//  double CC[95];
//  CC[0]=0.;
//  CC[ 1]=  3.79	; CC[31]= -0.63	; CC[61]=  2.88	;
//  CC[ 2]=  0.00	; CC[32]=  0.18	; CC[62]=  0.00	;
//  CC[ 3]= -0.05	; CC[33]=  0.09	; CC[63]=  2.99	;
//  CC[ 4]=  3.10	; CC[34]=  1.59	; CC[64]=  0.00	;
//  CC[ 5]= -1.01	; CC[35]=  0.17	; CC[65]= -2.43	;
//  CC[ 6]=  0.00	; CC[36]=  0.00	; CC[66]=  1.71	;
//  CC[ 7]=  0.00	; CC[37]= -0.56	; CC[67]=  0.00	;
//  CC[ 8]=  2.31	; CC[38]=  0.41	; CC[68]=  0.00	;
//  CC[ 9]=  0.00	; CC[39]=  0.00	; CC[69]= -0.86	;
//  CC[10]= -1.05	; CC[40]= -6.35	; CC[70]=  1.73	;
//  CC[11]=  0.00	; CC[41]=  0.00	; CC[71]=  0.00	;
//  CC[12]= -0.34	; CC[42]=  0.60	; CC[72]= -3.30	;
//  CC[13]=  0.00	; CC[43]=  0.00	; CC[73]=  0.50	;
//  CC[14]= -0.83	; CC[44]=  6.32	; CC[74]= -5.07	;
//  CC[15]=  0.00	; CC[45]=  0.00	; CC[75]=  0.00	;
//  CC[16]=  0.00	; CC[46]= -0.60 ; CC[76]= -1.44	;
//  CC[17]=  0.01	; CC[47]=  0.08	; CC[77]=  0.00	;
//  CC[18]= -0.56	; CC[48]=  3.41	; CC[78]= 17.51	;
//  CC[19]= -0.48	; CC[49]=  0.00	; CC[79]= -0.56	;
//  CC[20]=  0.18	; CC[50]=  8.71	; CC[80]=  0.87	;
//  CC[21]= -0.06	; CC[51]=-11.49 ; CC[81]=  0.00	;
//  CC[22]=  0.27	; CC[52]= -5.04	; CC[82]= -7.13	;
//  CC[23]=  0.00	; CC[53]=-11.99 ; CC[83]=  0.07	;
//  CC[24]=  1.62	; CC[54]=  0.00	; CC[84]=  0.00	;
//  CC[25]= -5.98	; CC[55]= 16.79 ; CC[85]= -0.82	;
//  CC[26]=  3.35	; CC[56]= 19.34 ; CC[86]=  0.00	;
//  CC[27]= -1.54	; CC[57]=  7.92	; CC[87]=  7.57	;
//  CC[28]=  0.30	; CC[58]=  0.00	; CC[88]= -5.47	;
//  CC[29]= -3.08	; CC[59]=-22.49 ; CC[89]= 34.74	;
//  CC[30]=  0.60	; CC[60]=  0.00	; CC[90]=  2.44	;
//
//for(int i=91; i<=94; i++) CC[i]=0.;
//
//for(int i=0; i<=94; i++) CC[i] *= 1.e-3*pow(mass.fpi,2);
//
// return Ci(CC,mass.mu,"Ci_JZLW");
//}
