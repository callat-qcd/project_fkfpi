// Ki.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// implementation of the Ki class and associated functions

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>  // for rand


#include "Ki.h"
const double pi16 = 1./(16.*M_PI*M_PI);
const double pi162 = pi16*pi16;

Ki::Ki(const double muin, const std::string Name, const int nfin){
    for(int i=0;i<=115;i++){
	Kr[i] = 0.;}
    mu = muin;
    name = Name;
    nf = nfin;
}

Ki::Ki(const double Krin[116],const double muin,const std::string Name,
       const int nfin){
    for(int i=0;i<=115;i++){
	Kr[i] = Krin[i];}
    mu = muin;
    name = Name;
    nf = nfin;
}

// Ki random below
Ki::~Ki(void){}

std::ostream & operator<<(std::ostream & os, const Ki & bb){
  os   
    << std::setprecision(7) << std::scientific << std::showpos
    << "# Ki set : " << bb.name << "\n"
    << "# nf  : "<<bb.nf <<"\n"
    << "# mu  : "<<bb.mu <<"\n"
    <<"# K1r  K2r  K3r  K4r :"
    <<bb.Kr[1]<<' '<<bb.Kr[2]<<' '<<bb.Kr[3]<<' '<<bb.Kr[4]<<' '<<'\n'
    <<"# K5r  K6r  K7r  K8r :"
    <<bb.Kr[5]<<' '<<bb.Kr[6]<<' '<<bb.Kr[7]<<' '<<bb.Kr[8]<<' '<<'\n'
    <<"# K9r  K10r K11r K12r:"
    <<bb.Kr[9]<<' '<<bb.Kr[10]<<' '<<bb.Kr[11]<<' '<<bb.Kr[12]<<' '<<'\n'
    <<"# K13r K14r K15r K16r:"
    <<bb.Kr[13]<<' '<<bb.Kr[14]<<' '<<bb.Kr[15]<<' '<<bb.Kr[16]<<' '<<'\n'
    <<"# K17r K18r K19r K20r:"
    <<bb.Kr[17]<<' '<<bb.Kr[18]<<' '<<bb.Kr[19]<<' '<<bb.Kr[20]<<' '<<'\n'
    <<"# K21r K22r K23r K24r:"
    <<bb.Kr[21]<<' '<<bb.Kr[22]<<' '<<bb.Kr[23]<<' '<<bb.Kr[24]<<' '<<'\n'
    <<"# K25r K26r K27r K28r:"
    <<bb.Kr[25]<<' '<<bb.Kr[26]<<' '<<bb.Kr[27]<<' '<<bb.Kr[28]<<' '<<'\n'
    <<"# K29r K30r K31r K32r:"
    <<bb.Kr[29]<<' '<<bb.Kr[30]<<' '<<bb.Kr[31]<<' '<<bb.Kr[32]<<' '<<'\n'
    <<"# K33r K34r K35r K36r:"
    <<bb.Kr[33]<<' '<<bb.Kr[34]<<' '<<bb.Kr[35]<<' '<<bb.Kr[36]<<' '<<'\n'
    <<"# K37r K38r K39r K40r:"
    <<bb.Kr[37]<<' '<<bb.Kr[38]<<' '<<bb.Kr[39]<<' '<<bb.Kr[40]<<' '<<'\n'
    <<"# K41r K42r K43r K44r:"
    <<bb.Kr[41]<<' '<<bb.Kr[42]<<' '<<bb.Kr[43]<<' '<<bb.Kr[44]<<' '<<'\n'
    <<"# K45r K46r K47r K48r:"
    <<bb.Kr[45]<<' '<<bb.Kr[46]<<' '<<bb.Kr[47]<<' '<<bb.Kr[48]<<' '<<'\n'
    <<"# K49r K50r K51r K52r:"
    <<bb.Kr[49]<<' '<<bb.Kr[50]<<' '<<bb.Kr[51]<<' '<<bb.Kr[52]<<' '<<'\n'
    <<"# K53r K54r K55r K56r:"
    <<bb.Kr[53]<<' '<<bb.Kr[54]<<' '<<bb.Kr[55]<<' '<<bb.Kr[56]<<' '<<'\n'
    <<"# K57r K58r K59r K60r:"
    <<bb.Kr[57]<<' '<<bb.Kr[58]<<' '<<bb.Kr[59]<<' '<<bb.Kr[60]<<' '<<'\n'
    <<"# K61r K62r K63r K64r:"
    <<bb.Kr[61]<<' '<<bb.Kr[62]<<' '<<bb.Kr[63]<<' '<<bb.Kr[64]<<' '<<'\n'
    <<"# K65r K66r K67r K68r:"
    <<bb.Kr[65]<<' '<<bb.Kr[66]<<' '<<bb.Kr[67]<<' '<<bb.Kr[68]<<' '<<'\n'
    <<"# K69r K70r K71r K72r:"
    <<bb.Kr[69]<<' '<<bb.Kr[70]<<' '<<bb.Kr[71]<<' '<<bb.Kr[72]<<' '<<'\n'
    <<"# K73r K74r K75r K76r:"
    <<bb.Kr[73]<<' '<<bb.Kr[74]<<' '<<bb.Kr[75]<<' '<<bb.Kr[76]<<' '<<'\n'
    <<"# K77r K78r K79r K80r:"
    <<bb.Kr[77]<<' '<<bb.Kr[78]<<' '<<bb.Kr[79]<<' '<<bb.Kr[80]<<' '<<'\n'
    <<"# K81r K82r K83r K84r:"
    <<bb.Kr[81]<<' '<<bb.Kr[82]<<' '<<bb.Kr[83]<<' '<<bb.Kr[84]<<' '<<'\n'
    <<"# K85r K86r K87r K88r:"
    <<bb.Kr[85]<<' '<<bb.Kr[86]<<' '<<bb.Kr[87]<<' '<<bb.Kr[88]<<' '<<'\n'
    <<"# K89r K90r K91r K92r:"
    <<bb.Kr[89]<<' '<<bb.Kr[90]<<' '<<bb.Kr[91]<<' '<<bb.Kr[92]<<' '<<'\n'
    <<"# K89r K90r K91r K92r:"
    <<bb.Kr[89]<<' '<<bb.Kr[90]<<' '<<bb.Kr[91]<<' '<<bb.Kr[92]<<' '<<'\n'
    <<"# K93r K94r K95r K96r:"
    <<bb.Kr[93]<<' '<<bb.Kr[94]<<' '<<bb.Kr[95]<<' '<<bb.Kr[96]<<' '<<'\n'
    <<"# K97r K98r K99r K100r:"
    <<bb.Kr[97]<<' '<<bb.Kr[98]<<' '<<bb.Kr[99]<<' '<<bb.Kr[100]<<' '<<'\n'
    <<"# K101r K102r K103r K104r:"
    <<bb.Kr[101]<<' '<<bb.Kr[102]<<' '<<bb.Kr[103]<<' '<<bb.Kr[104]<<' '<<'\n'
    <<"# K105r K106r K107r K108r:"
    <<bb.Kr[105]<<' '<<bb.Kr[106]<<' '<<bb.Kr[107]<<' '<<bb.Kr[108]<<' '<<'\n'
    <<"# K109r K110r K111r K112r:"
    <<bb.Kr[109]<<' '<<bb.Kr[110]<<' '<<bb.Kr[111]<<' '<<bb.Kr[112]<<' '<<'\n'
    <<"# K113r K114r KK115r     :"
    <<bb.Kr[113]<<' '<<bb.Kr[114]<<' '<<bb.Kr[115]<<' '<<'\n';
  os.unsetf(std::ios_base::floatfield);
  os << std::setprecision(6);  
  return os;
}

std::istream & operator>>(std::istream & is, Ki & Kiout){
  std::string temps;// used to sneak to the next line and check correct input
  // strings should be the same as the ones in mass output
  // up to the :
  std::string temps2;
  double KK[116],mu;
  std::string name;
  int nf;
  getline(is,temps2,':');// reads in characters up to ':'
  getline(is,temps); // to remove end of line
  name = temps;
  if (temps2 != "# Ki set ") std::cout << "trouble reading in Ki\n";
  getline(is,temps2,':');   is >> nf;  getline(is,temps);
  getline(is,temps2,':');   is >> mu;  getline(is,temps);
  for (int i = 0; i<112; i += 4){
    getline(is,temps2,':');   is >> KK[i+1] >> KK[i+2] >> KK[i+3] >> KK[i+4];
    getline(is,temps);
  }
  getline(is,temps2,':');   is >> KK[113] >> KK[114] >> KK[115];
  getline(is,temps);
  Kiout = Ki(KK,mu,name,nf);
  return is;
}

void Ki::setki(const int n, const double kin){
    Kr[n] = kin;
}

void Ki::setki(const double kin,const int n){
    Kr[n] = kin;
}

void Ki::setmu(const double muin){
  mu = muin;
}

void Ki::setname(const std::string inputname){
  name = inputname;
}

// summing two sets of Ki
Ki Ki::operator+(const Ki & bb) const{
  double Krout[116];
  for(int i=0;i<=114;i++){
      Krout[i] = Kr[i]+bb.Kr[i];}
  double muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
      std::cout << "WARNING: summing Ki "<< name << " and " << bb.name 
	   << " with different mu " << "\n";}
  return Ki(Krout,muout);
}
// difference
Ki Ki::operator-(const Ki & bb) const{
  double Krout[116];
  for(int i=0;i<=114;i++){
      Krout[i] = Kr[i]-bb.Kr[i];}
  double muout = mu;
  if ( fabs(1.-mu/bb.mu)>1e-8){
      std::cout << "WARNING: difference of Ki "<< name << " and " << bb.name 
	   << " with different mu " << "\n";}
  return Ki(Krout,muout);
}

// multiplying a Ki with a number

Ki Ki::operator*(const double & aa) const{
    double Krin[116];
    for(int i=0;i<=115;i++){
	Krin[i] = aa*Kr[i];}
  return Ki(Krin,mu);
}

Ki operator*(const double & aa,const Ki & bb){
    double Krin[116];
    for(int i=0;i<=114;i++){
	Krin[i] = aa*bb.Kr[i];}
  return Ki(Krin,bb.mu);
}

// changing scale
void Ki::changescale(const double newmu) const{
  std::cout << "WARNING, attempt to change scale to :"<<newmu<<'\n';
  std::cout << "WARNING:makes no sense to change scale in Ki without Linf specified"
	    <<'\n';
}

void Ki::changescale(Linf & linf,const double newmu){
    changescale(newmu,linf);
}



void Ki::changescale(const double newmu,Linf & linf){
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,mut;
  std::string namet;
  int nft;
  linf.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,
	   mut,namet,nft);
  std::cout << "nft = " <<nft<<'\n';
  if(fabs(1.-mu/mut)>=1e-6){
      std::cout << "ERROR: cannot change scale when Ki and Linf different scale"
	   << '\n' << " using Linf " << namet << " and Ki " << name
	   <<'\n';
      assert(0);
  }
  if( nft != nf){
      std::cout << "ERROR: cannot change scale when Ki and Linf different number of flavours"
	   << '\n' << " using Linf " << namet << " and Ki " << name
	   <<'\n';
      assert(0);
  }
// first the Ki, then the Linf, see notes
  double logm = log(mu/newmu);
  double logmm = -4.*log(mu/newmu);
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
  //const double g11=  0.;

// do this BEFORE resetting the lir
  double L0t =  pi16*((-2.)*logm*L0r -pi16*g0*logm*logm);
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
  //double L11t=  pi16*((-2.)*logm*L11r-pi16*g11*logm*logm);
  linf.changescale(newmu);

  double La1 = pi16*pi16*logmm;
                    
  Kr[1] +=  + 7./6.*L1t + 35./12.*L2t + 1./12.*nq*L0t + 3./8.*nq*L3t
          + 5./192.*La1 - 5./13824.*La1*pow(nq,2);


  Kr[2] +=  + 5./12.*L0t + 3./8.*L3t - 1./16.*L9t + 1./2.*nq*L1t + 1./
         6.*nq*L2t + 11./6912.*La1*nq;

  Kr[3] +=  - 1./2.*L1t - 5./4.*L2t - 1./128.*La1 - 1./4608.*La1*pow(
         nq,2);

  Kr[4] +=  - 1./12.*L3t - 1./8.*L9t + 1./1728.*La1*nq;

  Kr[5] +=  + 2./3.*L1t + 5./3.*L2t - 1./12.*nq*L0t + 1./24.*nq*L3t + 
         7./384.*La1 + 5./13824.*La1*pow(nq,2);

  Kr[6] +=  + 3./2.*L0t + 1./4.*L3t - 3./8.*L9t + 23./2304.*La1*nq;

  Kr[7] +=  - 4./3.*pow(nq,-1)*L0t - 4./3.*pow(nq,-1)*L3t + 2./3.*L1t
          + 1./3.*L2t + 1./12.*nq*L0t + 1./24.*nq*L3t - 1./4.*nq*L5t - 
         13./6912.*La1*pow(nq,2);

  Kr[8] +=  - 5./3.*pow(nq,-2)*L0t - 5./3.*pow(nq,-2)*L3t + pow(nq,-1)
         *L1t + 1./3.*pow(nq,-1)*L2t - 1./3.*L0t - 5./12.*L3t - 1./4.*
         L5t - nq*L1t - 1./3.*nq*L2t - 1./2.*nq*L4t + 35./3456.*La1*nq;

  Kr[9] +=  + 1./2.*L0t + 5./4.*L3t - 1./4.*L5t + 1./2.*nq*L1t + 1./4.
         *nq*L2t - 1./4.*nq*L4t + 5./192.*La1*nq;

  Kr[10] +=  + 1./2.*L1t + 1./4.*L2t - 1./4.*L4t + 5./256.*La1;

  Kr[11] +=  + 14./3.*pow(nq,-1)*L0t + 14./3.*pow(nq,-1)*L3t + 2./3.*
         L1t + 1./3.*L2t - 5./12.*nq*L0t - 11./24.*nq*L3t - 1./4.*nq*
         L5t + 1./16.*La1 + 1./216.*La1*pow(nq,2);

  Kr[12] +=  + L0t - 1./2.*L3t - 5./384.*La1*nq;

  Kr[13] +=  - 10./3.*pow(nq,-1)*L0t - 10./3.*pow(nq,-1)*L3t - 4./3.*
         L1t - 2./3.*L2t + 1./12.*nq*L0t + 13./24.*nq*L3t + 1./4.*nq*
         L5t - 1./16.*La1 - 5./3456.*La1*pow(nq,2);

  Kr[14] +=  + 5./3.*pow(nq,-2)*L0t + 5./3.*pow(nq,-2)*L3t - pow(
         nq,-1)*L1t - 1./3.*pow(nq,-1)*L2t + 1./12.*L0t + 13./24.*L3t
          + 1./4.*L5t + nq*L1t + 1./3.*nq*L2t + 1./2.*nq*L4t + 5./1728.
         *La1*nq;

  Kr[15] +=  - 2*L0t - 1./2.*L3t - 1./2.*nq*L2t - 1./96.*La1*nq;

  Kr[16] +=  - 1./2.*L2t + 1./128.*La1;

  Kr[17] +=  - 5./3.*pow(nq,-1)*L0t - 5./3.*pow(nq,-1)*L3t + 1./3.*L1t
          + 2./3.*L2t + 1./6.*nq*L0t + 1./3.*nq*L3t + 1./8.*nq*L5t - 13.
         /13824.*La1*pow(nq,2);

  Kr[18] +=  + 5./6.*pow(nq,-2)*L0t + 5./6.*pow(nq,-2)*L3t - 1./2.*
         pow(nq,-1)*L1t - 1./6.*pow(nq,-1)*L2t + 1./6.*L0t + 1./3.*L3t
          + 1./8.*L5t + 1./2.*nq*L1t + 1./6.*nq*L2t + 1./4.*nq*L4t - 1./
         13824.*La1*nq;

  Kr[19] +=  - 7./2.*pow(nq,-1)*L0t - 17./4.*pow(nq,-1)*L3t + 1./4.*
         L1t + 5./8.*L2t - L7t + 1./8.*nq*L0t + 5./16.*nq*L3t - 1./4.*
         nq*L8t + 1./64.*La1*pow(nq,-2) - 3./256.*La1 - 1./3072.*La1*
         pow(nq,2);

  Kr[20] +=  + 3*pow(nq,-2)*L0t + 3*pow(nq,-2)*L3t - 1./2.*pow(nq,-2)*
         L5t + 1./4.*L0t + 5./8.*L3t - 1./2.*L8t + 1./4.*nq*L4t - 1./2.
         *nq*L6t + 1./128.*La1*pow(nq,-1) + 5./384.*La1*nq;

  Kr[21] +=  + 3./2.*pow(nq,-2)*L0t + 3./2.*pow(nq,-2)*L3t - 2*pow(
         nq,-1)*L1t - 1./2.*pow(nq,-1)*L2t + pow(nq,-1)*L4t + 1./8.*L0t
          + 5./16.*L3t + 1./8.*L5t - 1./4.*L8t + 1./2.*nq*L1t + 1./8.*
         nq*L2t - 1./4.*nq*L4t - 1./256.*La1*pow(nq,-1) + 5./768.*La1*
         nq;

  Kr[22] +=  - 3./2.*pow(nq,-3)*L0t - 3./2.*pow(nq,-3)*L3t + pow(
         nq,-2)*L1t + 1./4.*pow(nq,-2)*L2t - 1./2.*pow(nq,-2)*L4t + 1./
         2.*L1t + 1./8.*L2t - 1./2.*L6t + 5./512.*La1;

  Kr[23] +=  - pow(nq,-1)*L0t - 1./4.*pow(nq,-1)*L3t + pow(nq,-1)*L5t
          + 1./4.*L1t + 5./8.*L2t + L7t - 1./8.*nq*L5t + 1./64.*La1*
         pow(nq,-2) - 1./384.*La1 + 3./1024.*La1*pow(nq,2);

  Kr[24] +=  - 1./2.*pow(nq,-1)*L1t - 5./4.*pow(nq,-1)*L2t + 1./2.*L0t
          + 1./8.*L3t - 1./4.*L5t + 1./2.*L8t - 1./96.*La1*pow(nq,-3)
          - 3./128.*La1*pow(nq,-1);

  Kr[25] +=  - 3./2.*pow(nq,-1)*L5t + 2*pow(nq,-1)*L8t - L7t + 1./8.*
         nq*L5t - 1./4.*nq*L8t;

  Kr[26] +=  + 3./2.*pow(nq,-2)*L5t - pow(nq,-2)*L8t - pow(nq,-1)*L4t
          + 2*pow(nq,-1)*L6t + 2*pow(nq,-1)*L7t + 3./8.*L5t - 3./4.*L8t
          + 1./4.*nq*L4t - 1./2.*nq*L6t;

  Kr[27] +=  - 1./2.*pow(nq,-3)*L5t + 1./2.*pow(nq,-2)*L4t - pow(
         nq,-2)*L6t - pow(nq,-2)*L7t + 1./4.*L4t - 1./2.*L6t;


  Kr[28] +=  + 1./3.*L1t + 1./3.*L2t - 1./3.*L4t - 1./6.*nq*L0t - 1./6.
         *nq*L3t + 1./24.*nq*L5t + 17./2304.*La1 + 1./864.*La1*pow(
         nq,2);


  Kr[29] +=  + 1./2.*L0t + 1./4.*L3t - 1./4.*L5t + 1./256.*La1*nq;


  Kr[30] +=  + 4./3.*pow(nq,-1)*L1t + 10./3.*pow(nq,-1)*L2t + 1./3.*
         L0t + 5./6.*L3t + 1./4.*L5t + 7./192.*La1*pow(nq,-1) + 41./
         6912.*La1*nq;


  Kr[31] +=  - 2*L1t - 4*L2t + 2./3.*L4t - 1./2.*nq*L3t + 1./6.*nq*L5t
          - 59./1152.*La1 - 1./2304.*La1*pow(nq,2);
 

  Kr[32] +=  - 3./2.*L0t - 3./4.*L3t + 3./8.*L9t - nq*L1t - 1./2.*nq*
         L2t - 25./1152.*La1*nq;


  Kr[33] +=  - 2*pow(nq,-1)*L0t - 2*pow(nq,-1)*L3t + 1./2.*pow(nq,-1)*
         L5t + 11./6.*L1t + 43./12.*L2t - 2./3.*L4t - L6t + 1./4.*nq*
         L0t + 19./24.*nq*L3t + 1./12.*nq*L5t - 1./4.*nq*L8t + 25./576.
         *La1 + 1./1536.*La1*pow(nq,2);


  Kr[34] +=  + 4*pow(nq,-2)*L0t + 4*pow(nq,-2)*L3t - 7./3.*pow(nq,-1)*
         L1t - 29./6.*pow(nq,-1)*L2t + pow(nq,-1)*L4t - 2./3.*L0t - 3./
         2.*L3t - 1./2.*L5t - 1./2.*L8t - 1./2.*nq*L7t - 1./48.*La1*
         pow(nq,-3) - 1./48.*La1*pow(nq,-1) - 7./1728.*La1*nq;


  Kr[35] +=  + 5./12.*L0t + 3./8.*L3t - 1./4.*L8t - 1./16.*L9t + 1./2.
         *nq*L1t + 1./6.*nq*L2t + 1./64.*La1*pow(nq,-1) + 13./13824.*
         La1*nq;


  Kr[36] +=  + 2./3.*pow(nq,-2)*L1t + 5./3.*pow(nq,-2)*L2t - 1./3.*
         pow(nq,-1)*L0t + 1./16.*pow(nq,-1)*L9t - 1./2.*L1t - 1./6.*L2t
          - 1./2.*L7t + 7./384.*La1*pow(nq,-2) - 23./2304.*La1;


  Kr[37] +=  - 2*pow(nq,-1)*L0t - 2*pow(nq,-1)*L3t + 1./2.*pow(nq,-1)*
         L5t - 1./6.*L1t - 5./12.*L2t - 1./3.*L4t + L6t + 1./3.*nq*L0t
          + 1./3.*nq*L3t + 1./6.*nq*L5t + 1./16.*La1*pow(nq,-2) - 17./
         1152.*La1 - 1./1728.*La1*pow(nq,2);


  Kr[38] +=  + 2*pow(nq,-2)*L0t + 2*pow(nq,-2)*L3t - 1./2.*pow(nq,-2)*
         L5t - pow(nq,-1)*L1t - 1./2.*pow(nq,-1)*L2t - 1./2.*pow(nq,-1)
         *L4t + 1./2.*L0t + 2./3.*L3t + 1./4.*L5t + 1./2.*L8t - 1./8.*
         L9t + nq*L1t + 1./2.*nq*L2t + nq*L4t - 3./64.*La1*pow(nq,-1)
          + 103./6912.*La1*nq;


  Kr[39] +=  - 5./3.*pow(nq,-1)*L0t - 5./3.*pow(nq,-1)*L3t + 2*pow(
         nq,-1)*L8t + 1./3.*L1t + 2./3.*L2t - L6t + 1./6.*nq*L0t + 1./3.
         *nq*L3t + 1./8.*nq*L5t - 1./4.*nq*L8t - 1./8.*La1*pow(nq,-2)
          + 1./192.*La1 + 5./13824.*La1*pow(nq,2);


  Kr[40] +=  + 5./6.*pow(nq,-2)*L0t + 5./6.*pow(nq,-2)*L3t - pow(
         nq,-2)*L8t - 1./2.*pow(nq,-1)*L1t - 1./6.*pow(nq,-1)*L2t + 1./
         6.*L0t + 1./3.*L3t + 1./8.*L5t - 1./4.*L8t + 1./2.*nq*L1t + 1./
         6.*nq*L2t + 1./4.*nq*L4t + 1./16.*La1*pow(nq,-3) + 1./64.*La1*
         pow(nq,-1) - 5./6912.*La1*nq;


  Kr[41] +=  + 10./3.*pow(nq,-2)*L0t + 10./3.*pow(nq,-2)*L3t - 2./3.*
         pow(nq,-1)*L1t - 4./3.*pow(nq,-1)*L2t + 2*pow(nq,-1)*L6t + 2*
         pow(nq,-1)*L7t - 1./3.*L0t - 2./3.*L3t - 1./4.*L5t - 1./2.*L8t
          - 1./2.*nq*L7t + 1./48.*La1*pow(nq,-3) + 11./192.*La1*pow(
         nq,-1) - 23./6912.*La1*nq;


  Kr[42] +=  - 5./2.*pow(nq,-3)*L0t - 5./2.*pow(nq,-3)*L3t + 5./6.*
         pow(nq,-2)*L1t + 5./6.*pow(nq,-2)*L2t - pow(nq,-2)*L6t - pow(
         nq,-2)*L7t - 1./2.*L1t - 1./6.*L2t - 1./4.*L4t - 1./2.*L7t - 1.
         /48.*La1*pow(nq,-4) - 55./4608.*La1;


  Kr[43] +=  + 2*pow(nq,-1)*L0t + 2*pow(nq,-1)*L3t - 1./2.*L1t - 3./4.
         *L2t + 1./2.*L4t - 1./4.*nq*L0t - 3./8.*nq*L3t - 1./4.*nq*L5t
          - 1./32.*La1*pow(nq,-2) + 1./384.*La1 - 1./1536.*La1*pow(
         nq,2);


  Kr[44] +=  - 4*pow(nq,-2)*L0t - 4*pow(nq,-2)*L3t + pow(nq,-1)*L1t + 
         3./2.*pow(nq,-1)*L2t - pow(nq,-1)*L4t + 1./2.*L0t + 3./4.*L3t
          + 1./2.*L5t + 1./48.*La1*pow(nq,-3) + 1./64.*La1*pow(nq,-1)
          - 1./768.*La1*nq;


  Kr[45] +=  - 2*pow(nq,-2)*L0t - 2*pow(nq,-2)*L3t + pow(nq,-1)*L1t + 
         1./2.*pow(nq,-1)*L2t + pow(nq,-1)*L4t - 1./2.*L0t - 3./4.*L3t
          - 1./2.*L5t - nq*L1t - 1./2.*nq*L2t - 3./2.*nq*L4t + 1./64.*
         La1*pow(nq,-1) - 1./64.*La1*nq;


  Kr[46] +=  - 1./96.*La1*pow(nq,-3) - 1./64.*La1*pow(nq,-1) + 5./1536.
         *La1*nq;


  Kr[47] +=  + pow(nq,-1)*L5t - 1./4.*nq*L5t + 1./32.*La1*pow(nq,-2)
          - 1./96.*La1 + 1./512.*La1*pow(nq,2);


  Kr[48] +=  - 1./2.*pow(nq,-2)*L5t + 1./2.*pow(nq,-1)*L4t - 1./4.*L5t
          - 1./2.*nq*L4t;


  Kr[49] +=  - 1./3.*L1t - 5./2.*L2t + 1./12.*nq*L0t - 7./24.*nq*L3t
          - 37./1152.*La1 - 5./3456.*La1*pow(nq,2);


  Kr[50] +=  - 3./4.*L0t - 5./24.*L3t + 1./8.*L9t - nq*L1t - 1./4.*nq*
         L2t + 41./3456.*La1*nq;


  Kr[51] +=  + 1./8.*L2t + 5./512.*La1;


  Kr[52] +=  - 4./3.*L2t - 3./8.*nq*L0t - 25./48.*nq*L3t + 7./288.*La1
          + 1./512.*La1*pow(nq,2);


  Kr[53] +=  + 1./2.*L0t - 5./8.*L3t - 5./512.*La1*nq;



  Kr[54] +=  - 1./3.*L1t + 11./2.*L2t - 1./12.*nq*L0t + 5./8.*nq*L3t
          + 41./1152.*La1 + 5./3456.*La1*pow(nq,2);


  Kr[55] +=  - 1./3.*L0t + 1./6.*L3t - 1./3.*nq*L2t - 1./144.*La1*nq;


  Kr[56] +=  - 1./4.*L2t + 3./256.*La1;


  Kr[57] +=  + 5./8.*L0t + 31./48.*L3t - 1./8.*L9t + nq*L1t + 3./8.*nq
         *L2t + 53./6912.*La1*nq;


  Kr[58] +=  + 1./3.*L1t + 7./6.*L2t + 1./24.*nq*L0t + 1./48.*nq*L3t
          + 65./1152.*La1 + 5./13824.*La1*pow(nq,2);


  Kr[59] +=  + 1./2.*L0t - 3./8.*L3t + 65./4608.*La1*nq;


  Kr[60] +=  + 1./3.*L1t - 17./6.*L2t - 1./24.*nq*L0t - 1./48.*nq*L3t
          - 97./1152.*La1 - 23./13824.*La1*pow(nq,2);


  Kr[61] +=  + 1./2.*L0t + 5./8.*L3t - 23./4608.*La1*nq;


  Kr[62] +=  - 13./6.*L0t - 5./12.*L3t - 1./6.*nq*L2t - 5./576.*La1*nq
         ;


  Kr[63] +=  - L2t - 1./64.*La1;


  Kr[64] +=  - 1./6.*L1t - 23./12.*L2t - 5./12.*nq*L0t - 7./24.*nq*L3t
          + 1./16.*nq*L9t - 11./384.*La1 - 1./1728.*La1*pow(nq,2);


  Kr[65] +=  - 7./3.*L0t - 4./3.*L3t + 1./2.*L9t - 2*nq*L1t - 5./6.*nq
         *L2t - 7./384.*La1*nq;


  Kr[66] +=  - 1./3.*L1t - 5./2.*L2t - 1./3.*L4t - 1./3.*nq*L0t - 1./
         12.*nq*L5t - 1./24.*nq*L9t - 7./192.*La1 + 5./3456.*La1*pow(
         nq,2);


  Kr[67] +=  - 2./3.*L1t - 19./3.*L2t - 1./3.*L4t - 1./6.*nq*L0t - 7./
         12.*nq*L3t - 1./12.*nq*L5t + 1./12.*nq*L9t - 3./32.*La1 + 1./
         864.*La1*pow(nq,2);


  Kr[68] +=  - 2./3.*L1t + 3*L2t + 1./3.*L4t + 1./6.*nq*L0t - 1./12.*
         nq*L3t + 1./12.*nq*L5t + 1./24.*nq*L9t + 11./192.*La1 + 1./
         6912.*La1*pow(nq,2);


  Kr[69] +=  - L0t + L3t - 47./1152.*La1*nq;


  Kr[70] +=  + 4./3.*L0t + 1./2.*L3t + 1./4.*L9t + 1./3.*nq*L2t + 1./
         432.*La1*nq;


  Kr[71] +=  - 1./3.*L1t + 1./6.*L4t - 1./12.*nq*L3t + 1./24.*nq*L5t
          + 1./48.*nq*L9t - 1./8.*nq*L10t - 1./96.*La1 - 1./576.*La1*
         pow(nq,2);


  Kr[72] +=  - 1./12.*L0t - 1./24.*L3t - 1./8.*L9t - 1./8.*L10t - 1./
         12.*nq*L2t - 5./1728.*La1*nq;


  Kr[73] +=  + 1./3.*L1t - 1./6.*L4t - 1./12.*nq*L0t + 1./24.*nq*L3t
          - 1./24.*nq*L5t - 7./48.*nq*L9t + 1./96.*La1 + 13./6912.*La1*
         pow(nq,2);


  Kr[74] +=  + 1./12.*L3t - 1./4.*L9t - 1./4.*L10t + 23./6912.*La1*nq;


  Kr[75] +=  + 1./3.*L2t + 1./3.*L4t + 1./3.*nq*L0t + 1./12.*nq*L5t - 
         5./24.*nq*L9t - 1./144.*La1 - 1./1728.*La1*pow(nq,2);


  Kr[76] +=  + L1t - 1./6.*L2t - 2./3.*L4t - 1./6.*nq*L0t + 1./4.*nq*
         L3t - 1./6.*nq*L5t - 1./6.*nq*L9t - 7./576.*La1 + 5./1728.*La1
         *pow(nq,2);


  Kr[77] +=  + 1./3.*L0t + 1./6.*L3t + 1./3.*nq*L2t - 11./864.*La1*nq;


  Kr[78] +=  - 1./2.*L1t - 1./12.*L2t + 1./6.*L4t + 1./12.*nq*L0t - 1./
         24.*nq*L3t + 1./24.*nq*L5t + 3./16.*nq*L9t + 11./1152.*La1 - 1.
         /1728.*La1*pow(nq,2);


  Kr[79] +=  - 1./6.*L3t + 55./3456.*La1*nq;


  Kr[80] +=  - 1./6.*L3t + 1./3456.*La1*nq;


  Kr[81] +=  - 1./8.*nq*L9t - 1./8.*nq*L10t;


  Kr[82] +=  - 1./8.*L9t - 1./8.*L10t;


  Kr[83] +=  + 5./3.*pow(nq,-1)*L0t + 5./3.*pow(nq,-1)*L3t + 1./6.*L1t
          - 11./12.*L2t - 5./12.*nq*L0t - 5./24.*nq*L3t - 1./8.*nq*L5t
          + 1./16.*nq*L9t + 5./6912.*La1*pow(nq,2);


  Kr[84] +=  - 10./3.*pow(nq,-2)*L0t - 10./3.*pow(nq,-2)*L3t + 2*pow(
         nq,-1)*L1t + 2./3.*pow(nq,-1)*L2t - 7./6.*L0t - 13./12.*L3t - 
         1./2.*L5t + 1./4.*L9t - 2*nq*L1t - 2./3.*nq*L2t - nq*L4t - 5./
         864.*La1*nq;


  Kr[85] +=  + 10./3.*pow(nq,-1)*L0t + 10./3.*pow(nq,-1)*L3t + 1./3.*
         L1t - 11./6.*L2t - 1./3.*nq*L0t - 2./3.*nq*L3t - 1./4.*nq*L5t
          + 1./8.*nq*L9t + 5./3456.*La1*pow(nq,2);


  Kr[86] +=  + 1./3.*L1t + 1./2.*L2t - 1./6.*nq*L0t - 1./4.*nq*L3t - 1.
         /24.*nq*L9t + 7./1152.*La1 + 1./2304.*La1*pow(nq,2);


  Kr[87] +=  + L0t + 1./6.*L3t - 1./4.*L9t + 23./3456.*La1*nq;


  Kr[88] +=  + 2./3.*L0t - L3t - 1./2.*L9t - 2*nq*L1t - 1./3.*nq*L2t
          - 17./3456.*La1*nq;


  Kr[89] +=  - 1./3.*L1t - 1./2.*L2t - 1./3.*nq*L3t - 1./12.*nq*L9t - 
         7./1152.*La1 - 1./6912.*La1*pow(nq,2);


  Kr[90] +=  - 4./3.*L1t - 17./6.*L2t + 1./6.*L4t - 5./12.*nq*L3t + 1./
         24.*nq*L5t - 1./16.*nq*L9t + 1./8.*nq*L10t - 47./1152.*La1 - 
         13./4608.*La1*pow(nq,2);


  Kr[91] +=  - 5./12.*L0t - 3./8.*L3t + 1./16.*L9t + 1./8.*L10t - 1./2.
         *nq*L1t - 1./6.*nq*L2t - 31./3456.*La1*nq;


  Kr[92] +=  + 2./3.*L1t + 7./6.*L2t - 1./6.*L4t - 1./12.*nq*L0t + 1./
         24.*nq*L3t - 1./24.*nq*L5t + 1./16.*nq*L9t + 13./576.*La1 + 17.
         /13824.*La1*pow(nq,2);


  Kr[93] +=  + 1./12.*L3t + 1./8.*L9t + 1./4.*L10t + 1./3456.*La1*nq;


  Kr[94] +=  + 5./3.*L1t + 9./2.*L2t - 1./3.*L4t - 1./6.*nq*L0t + 1./
         12.*nq*L3t - 1./12.*nq*L5t + 1./24.*nq*L9t + 11./192.*La1 + 23.
         /6912.*La1*pow(nq,2);


  Kr[95] +=  - 4./3.*L1t - 10./3.*L2t - 2./3.*nq*L3t - 1./6.*nq*L9t - 
         7./192.*La1 - 11./6912.*La1*pow(nq,2);


  Kr[96] +=  + L0t + 1./6.*L3t - 1./4.*L9t + 23./3456.*La1*nq;


  Kr[97] +=  - 5./6.*L1t - 9./4.*L2t + 1./6.*L4t - 1./3.*nq*L3t + 1./
         24.*nq*L5t - 1./16.*nq*L9t - 11./384.*La1 - 13./13824.*La1*
         pow(nq,2);


  Kr[98] +=  + 1./2.*L0t + 1./12.*L3t - 1./8.*L9t + 23./6912.*La1*nq;


  Kr[99] +=  - 4./3.*L0t - 4./3.*L3t - 2*nq*L1t - 1./3.*nq*L2t - 7./
         384.*La1*nq;


  Kr[100] +=  + 1./6.*L1t - 1./12.*L2t - 1./6.*L4t - 1./12.*nq*L0t + 1.
         /24.*nq*L3t - 1./24.*nq*L5t - 3./16.*nq*L9t + 5./1152.*La1 + 
         13./6912.*La1*pow(nq,2);


  Kr[101] +=  - 1./6.*L1t + 1./12.*L2t + 1./6.*L4t + 1./12.*nq*L0t - 1.
         /24.*nq*L3t + 1./24.*nq*L5t - 1./16.*nq*L9t - 5./1152.*La1 - 1.
         /1728.*La1*pow(nq,2);


  Kr[102] +=  + 5./3.*pow(nq,-1)*L0t + 5./3.*pow(nq,-1)*L3t - 1./3.*
         L1t - 2./3.*L2t - 1./6.*nq*L0t - 1./3.*nq*L3t - 1./8.*nq*L5t
          + 1./8.*nq*L10t - 17./13824.*La1*pow(nq,2);


  Kr[103] +=  - 5./6.*pow(nq,-2)*L0t - 5./6.*pow(nq,-2)*L3t + 1./2.*
         pow(nq,-1)*L1t + 1./6.*pow(nq,-1)*L2t - 1./6.*L0t - 1./3.*L3t
          - 1./8.*L5t + 1./8.*L10t - 1./2.*nq*L1t - 1./6.*nq*L2t - 1./4.
         *nq*L4t - 101./13824.*La1*nq;


  Kr[104] +=  - 1./12.*L1t + 1./24.*L2t + 1./12.*L4t + 1./24.*nq*L0t
          - 1./48.*nq*L3t + 1./48.*nq*L5t + 5./96.*nq*L9t + 19./2304.*
         La1 + 7./6912.*La1*pow(nq,2);


  Kr[105] +=  - 4./3.*pow(nq,-1)*L0t - 4./3.*pow(nq,-1)*L3t + 3./2.*
         L1t + 11./4.*L2t - 5./6.*L4t + 1./12.*nq*L0t + 13./24.*nq*L3t
          + 1./24.*nq*L5t + 1./12.*nq*L9t + 31./768.*La1 + 1./512.*La1*
         pow(nq,2);


  Kr[106] +=  + 8./3.*pow(nq,-2)*L0t + 8./3.*pow(nq,-2)*L3t - 2*pow(
         nq,-1)*L1t - 1./3.*pow(nq,-1)*L2t + 5./6.*L0t + 25./12.*L3t + 
         3./8.*L9t + 3*nq*L1t + 5./6.*nq*L2t + nq*L4t + 1./36.*La1*nq;


  Kr[107] +=  - 8./3.*pow(nq,-1)*L0t - 8./3.*pow(nq,-1)*L3t + 1./3.*
         L1t + 7./6.*L2t + 1./6.*nq*L0t + 7./12.*nq*L3t + 1./4.*nq*L5t
          + 1./8.*nq*L9t - 1./3456.*La1*pow(nq,2);


  Kr[108] +=  + 8./3.*pow(nq,-2)*L0t + 8./3.*pow(nq,-2)*L3t - 2*pow(
         nq,-1)*L1t - 1./3.*pow(nq,-1)*L2t + 1./3.*L0t + 7./6.*L3t + 1./
         2.*L5t + 1./4.*L9t + 2*nq*L1t + 1./3.*nq*L2t + nq*L4t + 5./432.
         *La1*nq;


  Kr[109] +=  - 1./12.*nq*L9t + 1./768.*La1*pow(nq,2);


  Kr[110] +=  - 1./3.*L1t + 1./6.*L2t + 1./3.*L4t + 1./6.*nq*L0t - 1./
         12.*nq*L3t + 1./12.*nq*L5t + 1./24.*nq*L9t - 5./576.*La1 - 7./
         3456.*La1*pow(nq,2);


  Kr[111] +=  + 1./3.*L1t - 1./6.*L2t - 1./6.*nq*L0t + 1./12.*nq*L3t
          - 7./24.*nq*L9t + 1./576.*La1*pow(nq,2);


  Kr[112] +=  + 1./3.*L4t + 1./12.*nq*L5t - 5./576.*La1 - 1./384.*La1*
         pow(nq,2);


  Kr[113] +=  - 1./8.*La1*pow(nq,-2) + 1./24.*La1 - 1./128.*La1*pow(
         nq,2);


  Kr[114] +=  - 2./3.*nq*L9t;


  Kr[115] +=  + 1./6.*nq*L9t;

  mu = newmu;
}


void Ki::out(double Kit[116],double & mut,std::string namet) const{
  for(int i =0;i<=115;i++){
      Kit[i]=Kr[i];}
  mut = mu;
  namet = name;
}

void Ki::out(double Kit[116],double & mut) const{
  for(int i =0;i<=115;i++){
      Kit[i]=Kr[i];}
  mut = mu;
}

void Ki::out(double Kit[116]) const{
  for(int i =0;i<=115;i++){
      Kit[i]=Kr[i];}
}

double Ki::out(const int n) const{
  return Kr[n];
}

int Ki::getnf(void) const{
  return nf;
}

Ki Kirandom(void){
 double KK[116];
 KK[0]=0.;
 for(int i=1;i<=115;++i){
     KK[i] = pi162*2.*(double(rand())/(double(RAND_MAX) + 1.)-0.5);
 }
 const double mu = 0.77;
 return Ki(KK,mu,"Ki_random");
}

