// massdecayvevnf.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results for QCD like theories for masses, decay constants and
// the qbarq vacuume expectation value
// derived in
// J.~Bijnens and J.~Lu,
// Technicolor and other QCD-like theories at next-to-next-to-leading order,''
// JHEP {\bf 0911} (2009) 116
// [arXiv:0910.5424 [hep-ph]].
// %%CITATION = ARXIV:0910.5424;%%

// in terms of the lowest order mass and the chiral limit
// decay constant

#include "oneloopintegrals.h"
#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"
#include "massdecayvevnf.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////
/////////////// SUN ///////////////////////////////////////////////////////
// mass SUN
double mnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfSUNp4L(nf,mass,Liin)+mnfSUNp4R(nf,mass);
}

double mnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSUNp4L\n";
  }
  double m112 = pow(mass.getB0mq(1)*2.,2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double MM12p4x11i =
    (+ 16.*L6r*nf
     - 8.*L4r*nf
     + 16.*L8r
     - 8.*L5r)*m112;
  return MM12p4x11i/pow(mass.getf0(),2);
}

double mnfSUNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return -m11*Ab(m11,mu2)*nfm/pow(mass.getf0(),2);
}

double mnfSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  mnfSUNp6K(nf,mass,Kiin)+mnfSUNp6L(nf,mass,Liin)+mnfSUNp6R(nf,mass);
}

double mnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  double  MM12p6x11iK =  KK[17] * (  - 32*pow(m11,3) );

  MM12p6x11iK +=  + KK[18] * (  - 32*pow(m11,3)*nf );

  MM12p6x11iK +=  + KK[19] * (  - 16*pow(m11,3) );

  MM12p6x11iK +=  + KK[20] * (  - 16*pow(m11,3)*nf );

  MM12p6x11iK +=  + KK[21] * (  - 16*pow(m11,3)*nf );

  MM12p6x11iK +=  + KK[22] * (  - 16*pow(m11,3)*pow(nf,2) );

  MM12p6x11iK +=  + KK[23] * (  - 16*pow(m11,3) );

  MM12p6x11iK +=  + KK[25] * ( 48*pow(m11,3) );

  MM12p6x11iK +=  + KK[26] * ( 48*pow(m11,3)*nf );

  MM12p6x11iK +=  + KK[27] * ( 48*pow(m11,3)*pow(nf,2) );

  MM12p6x11iK +=  + KK[39] * ( 32*pow(m11,3) );

  MM12p6x11iK +=  + KK[40] * ( 32*pow(m11,3)*nf );
  return MM12p6x11iK/pow(mass.getf0(),4);
}



double mnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MM12p6x11iL =
       + pi16 * ( 16*pow(m11,3)*L8r*nfm + 16*pow(m11,3)*L6r - 8*pow(
         m11,3)*L5r*nfm - 8*pow(m11,3)*L4r - 4*pow(m11,3)*L3r*nfm + 2*
         pow(m11,3)*L3r*nf + 2*pow(m11,3)*L2r*pow(nf,2) + 4*pow(m11,3)*
         L1r - 4*pow(m11,3)*L0r*nfm + 4*pow(m11,3)*L0r*nf );

      MM12p6x11iL +=  + Ab(m11,mu2) * (  - 80*pow(m11,2)*L8r*nfm + 32*
         pow(m11,2)*L8r*nf + 40*pow(m11,2)*L5r*nfm - 16*pow(m11,2)*L5r*
         nf - 24*pow(m11,2)*L3r*nfm + 20*pow(m11,2)*L3r*nf - 24*pow(
         m11,2)*L0r*nfm + 8*pow(m11,2)*L0r*nf - 8*( - 5 + 2*pow(nf,2))*
         pow(m11,2)*L4r + 16*( - 1 + nf)*(1 + nf)*pow(m11,2)*L6r + 8*(
          - 1 + 2*pow(nf,2))*pow(m11,2)*L1r + 4*(4 + pow(nf,2))*pow(
         m11,2)*L2r );

      MM12p6x11iL +=  - 128*pow(m11,3)*L5r*L8r - 128*pow(m11,3)*L5r*L6r
         *nf + 64*pow(m11,3)*pow(L5r,2) - 128*pow(m11,3)*L4r*L8r*nf - 
         128*pow(m11,3)*L4r*L6r*pow(nf,2) + 128*pow(m11,3)*L4r*L5r*nf
          + 64*pow(m11,3)*pow(L4r,2)*pow(nf,2);

  return MM12p6x11iL/pow(mass.getf0(),4);
}

double mnfSUNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSUNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double MM12p6x11iR =
       + pow(pi16,2) * ( 3./4.*pow(m11,3)*pow(nfm,2) + 1./384.*( - 96
          + 169*pow(nf,2))*pow(m11,3) );

      MM12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 4*pow(m11,2)*pow(nfm,2) + 
         1./48.*( - 80 + 57*pow(nf,2))*pow(m11,2) );

      MM12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 9./2.*m11*pow(nfm,2) + 1.
         /8.*( - 4 + 3*pow(nf,2))*m11 );

  return  MM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// decay constant SUN
double fnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfSUNp4L(nf,mass,Liin)+fnfSUNp4R(nf,mass);
}

double fnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSUNp4L\n";
  }
  double m11 = 2.*B0mq[0];
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double    MF12p4x11iL =
       + 4*m11*L5r + 4*m11*L4r*nf;
  return MF12p4x11iL/pow(mass.getf0(),2);
}

double fnfSUNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return + Ab(m11,mu2) * ( 1./2.*nf )/pow(mass.getf0(),2);
}

double fnfSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  fnfSUNp6K(nf,mass,Kiin)+fnfSUNp6L(nf,mass,Liin)+fnfSUNp6R(nf,mass);
}

double fnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  double  MF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );

      MF12p6x11iK +=  + KK[20] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[21] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[22] * ( 8*pow(m11,2)*pow(nf,2) );
      MF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );
  return MF12p6x11iK/pow(mass.getf0(),4);
}



double fnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  MF12p6x11iL =
       + pi16 * (  - 8*pow(m11,2)*L8r*nf - 8*pow(m11,2)*L6r*pow(nf,2)
          + 4*pow(m11,2)*L5r*nf + 4*pow(m11,2)*L4r*pow(nf,2) + 2*pow(
         m11,2)*L3r*nfm - pow(m11,2)*L3r*nf - pow(m11,2)*L2r*pow(nf,2)
          - 2*pow(m11,2)*L1r + 2*pow(m11,2)*L0r*nfm - 2*pow(m11,2)*L0r*
         nf );

      MF12p6x11iL +=  + Ab(m11,mu2) * ( 8*m11*L8r*nf + 8*m11*L6r*pow(
         nf,2) - 4*m11*L5r*nfm - 2*m11*L5r*nf + 12*m11*L3r*nfm - 10*m11
         *L3r*nf + 12*m11*L0r*nfm - 4*m11*L0r*nf - 4*( - 1 + 2*pow(
         nf,2))*m11*L1r - 2*(2 + pow(nf,2))*m11*L4r - 2*(4 + pow(nf,2))
         *m11*L2r );

      MF12p6x11iL +=  - 8*pow(m11,2)*pow(L5r,2) - 16*pow(m11,2)*L4r*L5r
         *nf - 8*pow(m11,2)*pow(L4r,2)*pow(nf,2);

  return MF12p6x11iL/pow(mass.getf0(),4);
}

double fnfSUNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSUNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double MF12p6x11iR =
       + pow(pi16,2) * ( 7./8.*pow(m11,2)*pow(nfm,2) + 1./768.*( - 224
          + pow(nf,2))*pow(m11,2) );

      MF12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 1./2.*m11*pow(nfm,2) - 
         1./96.*( - 64 + 59*pow(nf,2))*m11 );

      MF12p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 1./16.*(8 + 3*pow(
         nf,2)) );

  return  MF12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// vacuum expectation value SUN
double qnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfSUNp4L(nf,mass,Liin)+qnfSUNp4R(nf,mass);
}

double qnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSUNp4L\n";
  }
  double m11 = 2.*B0mq[0];
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double   MQ11p4x11iL =
       + 4*m11*H2r + 8*m11*L8r + 16*m11*L6r*nf;
  return MQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfSUNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return  + Ab(m11,mu2) * (  - nfm + nf )/pow(mass.getf0(),2);
}

double qnfSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  qnfSUNp6K(nf,mass,Kiin)+qnfSUNp6L(nf,mass,Liin)+qnfSUNp6R(nf,mass);
}

double qnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);

  double   MQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );
      MQ11p6x11iK +=  + KK[26] * ( 48*pow(m11,2)*nf );
      MQ11p6x11iK +=  + KK[27] * ( 48*pow(m11,2)*pow(nf,2) );

  return MQ11p6x11iK/pow(mass.getf0(),4);
}



double qnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MQ11p6x11iL =
       + pi16 * ( 16*pow(m11,2)*L8r*nfm - 16*pow(m11,2)*L8r*nf - 8*pow(
         m11,2)*L5r*nfm + 8*pow(m11,2)*L5r*nf - 16*( - 1 + nf)*(1 + nf)
         *pow(m11,2)*L6r + 8*( - 1 + nf)*(1 + nf)*pow(m11,2)*L4r );

      MQ11p6x11iL +=  + Ab(m11,mu2) * (  - 48*m11*L8r*nfm + 48*m11*L8r*
         nf + 24*m11*L5r*nfm - 24*m11*L5r*nf + 48*( - 1 + nf)*(1 + nf)*
         m11*L6r - 24*( - 1 + nf)*(1 + nf)*m11*L4r );

  return MQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfSUNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSUNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   MQ11p6x11iR =
       + Ab(m11,mu2)*pi16 * ( m11 - m11*pow(nfm,2) );

      MQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 3./2. + 3./2.*pow(
         nfm,2) );

  return  MQ11p6x11iR/pow(mass.getf0(),4);
}

///////////////////////////////////////////////////////////////////////////
/////////////// SON ///////////////////////////////////////////////////////
// mass SON
double mnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfSONp4L(nf,mass,Liin)+mnfSONp4R(nf,mass);
}

double mnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSONp4L\n";
  }
  double m11 = 2.*B0mq[0];
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double MM12p4x11iL =
       + 16*pow(m11,2)*L8r + 16*pow(m11,2)*L6r*nf - 8*pow(m11,2)*L5r - 
         8*pow(m11,2)*L4r*nf;

  return MM12p4x11iL/pow(mass.getf0(),2);
}

double mnfSONp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return + Ab(m11,mu2) * ( 1./2.*m11 - m11*nfm )/pow(mass.getf0(),2);
}

double mnfSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return mnfSONp6K(nf,mass,Kiin)+mnfSONp6L(nf,mass,Liin)+mnfSONp6R(nf,mass);
}

double mnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  */
  return 0.;//MM12p6x11iK/pow(mass.getf0(),4);
}



double mnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  MM12p6x11iL =
       + pi16 * (  - 8*pow(m11,3)*L8r + 16*pow(m11,3)*L8r*nfm + 4*pow(
         m11,3)*L5r - 8*pow(m11,3)*L5r*nfm - 4*pow(m11,3)*L3r*nfm + 4*
         pow(m11,3)*L1r - 4*pow(m11,3)*L0r*nfm - 8*( - 2 + nf)*pow(
         m11,3)*L6r + 4*( - 2 + nf)*pow(m11,3)*L4r + (1 + nf)*pow(
         m11,3)*L2r*nf + 2*(1 + nf)*pow(m11,3)*L0r + (2 + nf)*pow(
         m11,3)*L3r );

      MM12p6x11iL +=  + Ab(m11,mu2) * (  - 80*pow(m11,2)*L8r*nfm + 40*
         pow(m11,2)*L5r*nfm - 24*pow(m11,2)*L3r*nfm - 24*pow(m11,2)*L0r
         *nfm - 4*( - 10 + 5*nf + 2*pow(nf,2))*pow(m11,2)*L4r + 8*( - 2
          + 3*nf + pow(nf,2))*pow(m11,2)*L6r + 8*( - 1 + nf + pow(nf,2)
         )*pow(m11,2)*L1r + 4*(3 + nf)*pow(m11,2)*L0r + 8*(5 + 2*nf)*
         pow(m11,2)*L8r - 4*(5 + 2*nf)*pow(m11,2)*L5r + 2*(6 + 5*nf)*
         pow(m11,2)*L3r + 2*(8 + nf + pow(nf,2))*pow(m11,2)*L2r );

      MM12p6x11iL +=  - 128*pow(m11,3)*L5r*L8r - 128*pow(m11,3)*L5r*L6r
         *nf + 64*pow(m11,3)*pow(L5r,2) - 128*pow(m11,3)*L4r*L8r*nf - 
         128*pow(m11,3)*L4r*L6r*pow(nf,2) + 128*pow(m11,3)*L4r*L5r*nf
          + 64*pow(m11,3)*pow(L4r,2)*pow(nf,2);

  return MM12p6x11iL/pow(mass.getf0(),4);
}

double mnfSONp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSONp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double  MM12p6x11iR =
       + pow(pi16,2) * (  - 3./8.*pow(m11,3)*nfm + 3./4.*pow(m11,3)*
         pow(nfm,2) + 1./1536.*( - 96 + 386*nf + 169*pow(nf,2))*pow(
         m11,3) );

      MM12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 3./2.*pow(m11,2)*nfm + 
         4*pow(m11,2)*pow(nfm,2) + 1./192.*( - 2 + 3*nf)*(64 + 19*nf)*
         pow(m11,2) );

      MM12p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 3*m11*nfm + 9./2.*m11
         *pow(nfm,2) + 3./32.*(4 + 2*nf + pow(nf,2))*m11 );

  return  MM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// decay constant SON

double fnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfSONp4L(nf,mass,Liin)+fnfSONp4R(nf,mass);
}

double fnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSONp4L\n";
  }
  double m11 = 2.*B0mq[0];
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  MF12p4x11iL =
       + 4*m11*L5r + 4*m11*L4r*nf;
  return MF12p4x11iL/pow(mass.getf0(),2);
}

double fnfSONp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return + Ab(m11,mu2) * ( 1./4.*nf )/pow(mass.getf0(),2);
}

double fnfSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  fnfSONp6K(nf,mass,Kiin)+fnfSONp6L(nf,mass,Liin)+fnfSONp6R(nf,mass);
}

double fnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  double  MF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );

      MF12p6x11iK +=  + KK[20] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[21] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[22] * ( 8*pow(m11,2)*pow(nf,2) );
      MF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );
  return MF12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double fnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MF12p6x11iL =
       + pi16 * (  - 4*pow(m11,2)*L8r*nf - 4*pow(m11,2)*L6r*pow(nf,2)
          + 2*pow(m11,2)*L5r*nf + 2*pow(m11,2)*L4r*pow(nf,2) + 2*pow(
         m11,2)*L3r*nfm - 2*pow(m11,2)*L1r + 2*pow(m11,2)*L0r*nfm - 1./
         2.*(1 + nf)*pow(m11,2)*L2r*nf - (1 + nf)*pow(m11,2)*L0r - 1./2.
         *(2 + nf)*pow(m11,2)*L3r );

      MF12p6x11iL +=  + Ab(m11,mu2) * ( 4*m11*L8r*nf + 4*m11*L6r*pow(
         nf,2) - 4*m11*L5r*nfm + 12*m11*L3r*nfm + 12*m11*L0r*nfm - ( - 
         2 + nf)*m11*L5r - 4*( - 1 + nf + pow(nf,2))*m11*L1r - 2*(3 + 
         nf)*m11*L0r - (4 - 2*nf + pow(nf,2))*m11*L4r - (6 + 5*nf)*m11*
         L3r - (8 + nf + pow(nf,2))*m11*L2r );

      MF12p6x11iL +=  - 8*pow(m11,2)*pow(L5r,2) - 16*pow(m11,2)*L4r*L5r
         *nf - 8*pow(m11,2)*pow(L4r,2)*pow(nf,2);

  return MF12p6x11iL/pow(mass.getf0(),4);
}

double fnfSONp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSONp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   MF12p6x11iR =
       + pow(pi16,2) * (  - 7./16.*pow(m11,2)*nfm + 7./8.*pow(m11,2)*
         pow(nfm,2) + 1./3072.*( - 224 + 114*nf + pow(nf,2))*pow(m11,2)
          );

      MF12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 1./4.*m11*nfm - 1./2.*m11*
         pow(nfm,2) - 1./384.*( - 112 + 174*nf + 59*pow(nf,2))*m11 );

      MF12p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 1./64.*(16 - 6*nf + 3
         *pow(nf,2)) );

  return  MF12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// vacuum expectation value SON

double qnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfSONp4L(nf,mass,Liin)+qnfSONp4R(nf,mass);
}

double qnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSONp4L\n";
  }
  double m11 = 2.*B0mq[0];
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double  MQ11p4x11iL =
       + 4*m11*H2r + 8*m11*L8r + 16*m11*L6r*nf;
 
 return MQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfSONp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return  + Ab(m11,mu2) * (  - nfm + 1./2.*(1 + nf) )/pow(mass.getf0(),2);
}

double qnfSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  qnfSONp6K(nf,mass,Kiin)+qnfSONp6L(nf,mass,Liin)+qnfSONp6R(nf,mass);
}

double qnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);

  double   MQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );
      MQ11p6x11iK +=  + KK[26] * ( 48*pow(m11,2)*nf );
      MQ11p6x11iK +=  + KK[27] * ( 48*pow(m11,2)*pow(nf,2) );
 
  return MQ11p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double qnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MQ11p6x11iL =
       + pi16 * ( 16*pow(m11,2)*L8r*nfm - 8*pow(m11,2)*L5r*nfm - 8*( - 
         1 + nf)*(2 + nf)*pow(m11,2)*L6r + 4*( - 1 + nf)*(2 + nf)*pow(
         m11,2)*L4r - 8*(1 + nf)*pow(m11,2)*L8r + 4*(1 + nf)*pow(m11,2)
         *L5r );

      MQ11p6x11iL +=  + Ab(m11,mu2) * (  - 48*m11*L8r*nfm + 24*m11*L5r*
         nfm + 24*( - 1 + nf)*(2 + nf)*m11*L6r - 12*( - 1 + nf)*(2 + nf
         )*m11*L4r + 24*(1 + nf)*m11*L8r - 12*(1 + nf)*m11*L5r );

  return MQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfSONp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSONp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double  MQ11p6x11iR =
       + Ab(m11,mu2)*pi16 * ( m11*nfm - m11*pow(nfm,2) - 1./4.*( - 1 + 
         nf)*m11 );

      MQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 3./2.*nfm + 3./2.*
         pow(nfm,2) + 3./8.*( - 1 + nf) );

  return  MQ11p6x11iR/pow(mass.getf0(),4);
}

///////////////////////////////////////////////////////////////////////////
/////////////// SPN ///////////////////////////////////////////////////////
// mass SPN
double mnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfSPNp4L(nf,mass,Liin)+mnfSPNp4R(nf,mass);
}

double mnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSPNp4L\n";
  }
  double m11 = 2.*B0mq[0];
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MM12p4x11iL =
       + 16*pow(m11,2)*L8r + 32*pow(m11,2)*L6r*nf - 8*pow(m11,2)*L5r - 
         16*pow(m11,2)*L4r*nf;

  return MM12p4x11iL/pow(mass.getf0(),2);
}

double mnfSPNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return + Ab(m11,mu2) * (  - 1./2.*m11 - 1./2.*m11*nfm )/pow(mass.getf0(),2);
}

double mnfSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return mnfSPNp6K(nf,mass,Kiin)+mnfSPNp6L(nf,mass,Liin)+mnfSPNp6R(nf,mass);
}

double mnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: mnfSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  */
  return 0.;//MM12p6x11iK/pow(mass.getf0(),4);
}



double mnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  MM12p6x11iL =
       + pi16 * ( 8*pow(m11,3)*L8r + 8*pow(m11,3)*L8r*nfm - 4*pow(
         m11,3)*L5r - 4*pow(m11,3)*L5r*nfm - 2*pow(m11,3)*L3r*nfm + 4*
         pow(m11,3)*L1r - 2*pow(m11,3)*L0r*nfm + 2*( - 1 + nf)*pow(
         m11,3)*L3r + 2*( - 1 + 2*nf)*pow(m11,3)*L2r*nf + 2*( - 1 + 2*
         nf)*pow(m11,3)*L0r + 16*(1 + nf)*pow(m11,3)*L6r - 8*(1 + nf)*
         pow(m11,3)*L4r );

      MM12p6x11iL +=  + Ab(m11,mu2) * (  - 40*pow(m11,2)*L8r*nfm + 20*
         pow(m11,2)*L5r*nfm - 12*pow(m11,2)*L3r*nfm - 12*pow(m11,2)*L0r
         *nfm - 8*( - 5 - 5*nf + 4*pow(nf,2))*pow(m11,2)*L4r + 8*( - 5
          + 4*nf)*pow(m11,2)*L8r - 4*( - 5 + 4*nf)*pow(m11,2)*L5r + 4*(
          - 3 + 2*nf)*pow(m11,2)*L0r + 4*( - 3 + 5*nf)*pow(m11,2)*L3r
          + 16*( - 1 - 3*nf + 2*pow(nf,2))*pow(m11,2)*L6r + 8*( - 1 - 2
         *nf + 4*pow(nf,2))*pow(m11,2)*L1r + 4*(4 - nf + 2*pow(nf,2))*
         pow(m11,2)*L2r );

      MM12p6x11iL +=  - 128*pow(m11,3)*L5r*L8r - 256*pow(m11,3)*L5r*L6r
         *nf + 64*pow(m11,3)*pow(L5r,2) - 256*pow(m11,3)*L4r*L8r*nf - 
         512*pow(m11,3)*L4r*L6r*pow(nf,2) + 256*pow(m11,3)*L4r*L5r*nf
          + 256*pow(m11,3)*pow(L4r,2)*pow(nf,2);

  return MM12p6x11iL/pow(mass.getf0(),4);
}

double mnfSPNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: mnfSPNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double    MM12p6x11iR =
       + pow(pi16,2) * ( 3./16.*pow(m11,3)*nfm + 3./16.*pow(m11,3)*pow(
         nfm,2) + 1./384.*( - 24 - 193*nf + 169*pow(nf,2))*pow(m11,3) )
         ;

      MM12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 3./4.*pow(m11,2)*nfm + 
         pow(m11,2)*pow(nfm,2) + 1./48.*(1 + 3*nf)*( - 32 + 19*nf)*pow(
         m11,2) );

      MM12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 3./2.*m11*nfm + 9./8.*
         m11*pow(nfm,2) + 3./8.*(1 - nf + pow(nf,2))*m11 );

  return  MM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// decay constant SPN

double fnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfSPNp4L(nf,mass,Liin)+fnfSPNp4R(nf,mass);
}

double fnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSPNp4L\n";
  }
  double m11 = 2.*B0mq[0];
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  MF12p4x11iL =
       + 4*m11*L5r + 8*m11*L4r*nf;
  return MF12p4x11iL/pow(mass.getf0(),2);
}

double fnfSPNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return + Ab(m11,mu2) * ( 1./2.*nf )/pow(mass.getf0(),2);
}

double fnfSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  fnfSPNp6K(nf,mass,Kiin)+fnfSPNp6L(nf,mass,Liin)+fnfSPNp6R(nf,mass);
}

double fnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: fnfSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);
  double  MF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );

      MF12p6x11iK +=  + KK[20] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[21] * ( 8*pow(m11,2)*nf );
      MF12p6x11iK +=  + KK[22] * ( 8*pow(m11,2)*pow(nf,2) );
      MF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );
  return MF12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double fnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MF12p6x11iL =
       + pi16 * (  - 8*pow(m11,2)*L8r*nf - 16*pow(m11,2)*L6r*pow(nf,2)
          + 4*pow(m11,2)*L5r*nf + 8*pow(m11,2)*L4r*pow(nf,2) + pow(
         m11,2)*L3r*nfm - 2*pow(m11,2)*L1r + pow(m11,2)*L0r*nfm - ( - 1
          + nf)*pow(m11,2)*L3r - ( - 1 + 2*nf)*pow(m11,2)*L2r*nf - ( - 
         1 + 2*nf)*pow(m11,2)*L0r );

      MF12p6x11iL +=  + Ab(m11,mu2) * ( 8*m11*L8r*nf + 16*m11*L6r*pow(
         nf,2) - 2*m11*L5r*nfm + 6*m11*L3r*nfm + 6*m11*L0r*nfm - 2*( - 
         3 + 2*nf)*m11*L0r - 2*( - 3 + 5*nf)*m11*L3r - 4*( - 1 - 2*nf
          + 4*pow(nf,2))*m11*L1r - 2*(1 + nf)*m11*L5r - 4*(1 + nf + 
         pow(nf,2))*m11*L4r - 2*(4 - nf + 2*pow(nf,2))*m11*L2r );

      MF12p6x11iL +=  - 8*pow(m11,2)*pow(L5r,2) - 32*pow(m11,2)*L4r*L5r
         *nf - 32*pow(m11,2)*pow(L4r,2)*pow(nf,2);

  return MF12p6x11iL/pow(mass.getf0(),4);
}

double fnfSPNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: fnfSPNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   MF12p6x11iR =
       + pow(pi16,2) * ( 7./32.*pow(m11,2)*nfm + 7./32.*pow(m11,2)*pow(
         nfm,2) + 1./768.*( - 56 - 57*nf + pow(nf,2))*pow(m11,2) );

      MF12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 1./8.*m11*nfm - 1./8.*
         m11*pow(nfm,2) - 1./96.*( - 28 - 87*nf + 59*pow(nf,2))*m11 );

      MF12p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 1./16.*(4 + 3*nf + 3*
         pow(nf,2)) );

  return  MF12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////
// vacuum expectation value SPN

double qnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfSPNp4L(nf,mass,Liin)+qnfSPNp4R(nf,mass);
}

double qnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){

  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSPNp4L\n";
  }
  double m11 = 2.*B0mq[0];
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double  MQ11p4x11iL =
       + 4*m11*H2r + 8*m11*L8r + 32*m11*L6r*nf;
 
 return MQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfSPNp4R(const int nf,const quarkmassnf mass){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  return  + Ab(m11,mu2) * (  - 1./2.*nfm + 1./2.*( - 1 + 2*nf) )/pow(mass.getf0(),2);
}

double qnfSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return  qnfSPNp6K(nf,mass,Kiin)+qnfSPNp6L(nf,mass,Liin)+qnfSPNp6R(nf,mass);
}

double qnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 1){
    cout << "ERROR: qnfSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double KK[116];
  Kiin.out(KK);

  double   MQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );
      MQ11p6x11iK +=  + KK[26] * ( 48*pow(m11,2)*nf );
      MQ11p6x11iK +=  + KK[27] * ( 48*pow(m11,2)*pow(nf,2) );
 
  return MQ11p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double qnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   MQ11p6x11iL =
       + pi16 * ( 8*pow(m11,2)*L8r*nfm - 4*pow(m11,2)*L5r*nfm - 16*( - 
         1 + nf)*(1 + 2*nf)*pow(m11,2)*L6r + 8*( - 1 + nf)*(1 + 2*nf)*
         pow(m11,2)*L4r - 8*( - 1 + 2*nf)*pow(m11,2)*L8r + 4*( - 1 + 2*
         nf)*pow(m11,2)*L5r );

      MQ11p6x11iL +=  + Ab(m11,mu2) * (  - 24*m11*L8r*nfm + 12*m11*L5r*
         nfm + 48*( - 1 + nf)*(1 + 2*nf)*m11*L6r - 24*( - 1 + nf)*(1 + 
         2*nf)*m11*L4r + 24*( - 1 + 2*nf)*m11*L8r - 12*( - 1 + 2*nf)*
         m11*L5r );

  return MQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfSPNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 1){
    cout << "ERROR: qnfSPNp6R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = B0mq[0]*2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double  MQ11p6x11iR =
       + Ab(m11,mu2)*pi16 * (  - 1./2.*m11*nfm - 1./4.*m11*pow(nfm,2)
          + 1./4.*(1 + 2*nf)*m11 );

      MQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 3./4.*nfm + 3./8.*pow(
         nfm,2) - 3./8.*(1 + 2*nf) );

  return  MQ11p6x11iR/pow(mass.getf0(),4);
}


