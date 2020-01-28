// massdecayvevnfPQ.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results for QCD like theories for masses, decay constants and
// the qbarq vacuum expectation value, partially quenched case
// derived in
// J.~Bijnens and T.~Rössler

// in terms of the lowest order mass and the chiral limit
// decay constant

#include "oneloopintegrals.h"
#include "quenchedsunsetintegrals.h"
#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"
#include "massdecayvevnfPQ.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////
/////////////SUN case//////////////////////////////////////////////////////
// mass SUN
double mnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfPQSUNp4L(nf,mass,Liin)+mnfPQSUNp4R(nf,mass);
}

double mnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSUNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double PQM12p4x11iL =
       + 16*m11*m44*L6r*nf - 8*m11*m44*L4r*nf + 16*pow(m11,2)*L8r - 8*
         pow(m11,2)*L5r;

  return PQM12p4x11iL/pow(mass.getf0(),2);
}

double mnfPQSUNp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double PQM12p4x11iR =
       + pi16 * (  - (m44 - m11)*m11*nfm );

      PQM12p4x11iR +=  + Ab(m11,mu2) * ( (m44 - 2*m11)*nfm );

  return PQM12p4x11iR/pow(mass.getf0(),2);
}

double mnfPQSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return mnfPQSUNp6K(nf,mass,Kiin)+mnfPQSUNp6L(nf,mass,Liin)
    +mnfPQSUNp6R(nf,mass);
}

double mnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfPQSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQM12p6x11iK = 0.;
      PQM12p6x11iK +=  + KK[17] * (  - 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[18] * (  - 32*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[19] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[20] * (  - 16*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[21] * (  - 16*m11*pow(m44,2)*nf );
      PQM12p6x11iK +=  + KK[22] * (  - 16*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[23] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[25] * ( 48*pow(m11,3) );
      PQM12p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m11*m44*nf );
      PQM12p6x11iK +=  + KK[27] * ( 48*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[39] * ( 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[40] * ( 32*pow(m11,2)*m44*nf );
 
  return PQM12p6x11iK/pow(mass.getf0(),4);
}



double mnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQM12p6x11iL =
       + pi16 * ( 4*pow(m11,3)*L1r - 16*(m44 - m11)*(m44 - m11)*m11*L7r
          + 1./2.*(m44 + m11)*(m44 + m11)*m11*L3r*nf + (m44 + m11)*(m44
          + m11)*m11*L0r*nf - 16*(3*m44 - 4*m11)*m11*m44*L6r + 8*(4*m44
          - 5*m11)*m11*m44*L4r - 4*(4*m44 - 3*m11)*pow(m11,2)*L3r*nfm
          - 4*(4*m44 - 3*m11)*pow(m11,2)*L0r*nfm + 8*(7*m44 - 8*m11)*
         pow(m11,2)*L5r*nfm + 2*( - pow(m44,2) + pow(m44,2)*pow(nf,2)
          + pow(m11,2))*m11*L2r - 16*(pow(m44,2) + 5*m11*m44 - 7*pow(
         m11,2))*m11*L8r*nfm );

      PQM12p6x11iL +=  + Ab(m11,mu2) * ( 20*pow(m11,2)*L2r + 8*pow(
         m11,2)*L1r - 24*(m44 - 2*m11)*m44*L4r + 16*(m44 - m11)*(m44 - 
         3*m11)*L7r + 32*(m44 - m11)*(m44 - m11)*L6r + 24*(2*m44 - 3*
         m11)*m11*L3r*nfm + 24*(2*m44 - 3*m11)*m11*L0r*nfm - 8*(8*m44
          - 13*m11)*m11*L5r*nfm + 16*(pow(m44,2) + 4*m11*m44 - 10*pow(
         m11,2))*L8r*nfm );

      PQM12p6x11iL +=  + Ab(m14,mu2) * ( 16*(m44 + m11)*m11*L8r*nf - 8*
         (m44 + m11)*m11*L5r*nf + 10*(m44 + m11)*m11*L3r*nf + 4*(m44 + 
         m11)*m11*L0r*nf );

      PQM12p6x11iL +=  + Ab(m44,mu2) * ( 16*( - 1 + nf)*(1 + nf)*m11*
         m44*L6r - 16*( - 1 + nf)*(1 + nf)*m11*m44*L4r + 4*( - 1 + nf)*
         (1 + nf)*m11*m44*L2r + 16*( - 1 + nf)*(1 + nf)*m11*m44*L1r );

      PQM12p6x11iL +=  - 128*m11*pow(m44,2)*L4r*L6r*pow(nf,2) + 64*m11*
         pow(m44,2)*pow(L4r,2)*pow(nf,2) - 128*pow(m11,2)*m44*L5r*L6r*
         nf - 128*pow(m11,2)*m44*L4r*L8r*nf + 128*pow(m11,2)*m44*L4r*
         L5r*nf - 128*pow(m11,3)*L5r*L8r + 64*pow(m11,3)*pow(L5r,2);
 
  return PQM12p6x11iL/pow(mass.getf0(),4);
}

double mnfPQSUNp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSUNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQM12p6x11iR =
       + pow(pi16,2) * (  - 1./48.*( - m44 + m44*nf + 3*m11)*(m44 + m44
         *nf - 3*m11)*m11*pow(pi,2) + 1./3.*(m44 - m11)*pow(m11,2)*pow(
         pi,2)*pow(nfm,2) - 1./4.*(m44*pow(nf,2) - 4*m11 + m11*pow(
         nf,2))*m11*m14 - 1./24.*(m44*pow(nf,2) - 4*m11 + m11*pow(nf,2)
         )*m11*m14*pow(pi,2) - 1./576.*( - 54*pow(m44,2) - 3*pow(m44,2)
         *pow(nf,2) + 396*m11*m44 - 26*m11*m44*pow(nf,2) - 702*pow(
         m11,2) + 4*pow(m11,2)*pow(nf,2))*m11 + 3./4.*(2*pow(m44,2) - 4
         *m11*m44 + 3*pow(m11,2))*m11*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 1./12.*(15*m44 - 17*m11)*
         m11 - (2*pow(m44,2) - 6*m11*m44 + pow(m11,2))*pow(nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 1./2.*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) - (2*m44 - 7*m11)*pow(nfm,2) - 1./8.*(6*m44
          - 25*m11) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 2*m44 );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( m44 - m44*pow(
         nfm,2) );

      PQM12p6x11iR +=  + Ab(m14,mu2)*pi16 * ( 1./3.*(6 + pow(nf,2))*m11
         *m44 );

      PQM12p6x11iR +=  + pow(Ab(m14,mu2),2) * (  - 1./4.*m11*pow(nf,2)
          - 1./4.*(m44*pow(nf,2) - 4*m11 + m11*pow(nf,2))*m11*pow(
         m14,-1) );

      PQM12p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - m11*m44 + m11*m44*pow(
         nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m44,mu2),2) * (  - 1./8.*( - 1 + nf)*(1
          + nf)*m11 );

      PQM12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * ( 1./4.*(m44 - 4*
         m11)*m11 );

      PQM12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * ( 1./4.*( - 1 + nf
         )*(1 + nf)*m11*m44 );

      PQM12p6x11iR +=  + hh(2,m11,m11,m11,m11,mu2) * (  - 4*(m44 - m11)
         *pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * ( 5./4.*(m44 - m11
         )*pow(m11,2) );

      PQM12p6x11iR +=  + hh(5,m11,m11,m11,m11,mu2) * ( 2*(m44 - m11)*(
         m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * (  - 2*(m44 - m11
         )*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(1,m11,m14,m14,m11,mu2) * ( 3./4.*pow(
         m11,2) );

      PQM12p6x11iR +=  + hh21(1,m44,m14,m14,m11,mu2) * ( 3./4.*( - 1 + 
         nf)*(1 + nf)*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(2,m11,m14,m14,m11,mu2) * ( 3./4.*(m44 - 
         m11)*pow(m11,2) );

  return  PQM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////////////
// decay SUN
double fnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfPQSUNp4L(nf,mass,Liin)+fnfPQSUNp4R(nf,mass);
}

double fnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSUNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double PQF12p4x11iL =
       + 4*m44*L4r*nf + 4*m11*L5r;

  return PQF12p4x11iL/pow(mass.getf0(),2);
}

double fnfPQSUNp4R(const int nf,const quarkmassnf mass){
  //const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQF12p4x11iR =
       + Ab(m14,mu2) * ( 1./2.*nf );

  return PQF12p4x11iR/pow(mass.getf0(),2);
}

double fnfPQSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return fnfPQSUNp6K(nf,mass,Kiin)+fnfPQSUNp6L(nf,mass,Liin)
    +fnfPQSUNp6R(nf,mass);
}

double fnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfPQSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double  PQF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );
      PQF12p6x11iK +=  + KK[20] * ( 8*m11*m44*nf );
      PQF12p6x11iK +=  + KK[21] * ( 8*pow(m44,2)*nf );
      PQF12p6x11iK +=  + KK[22] * ( 8*pow(m44,2)*pow(nf,2) );
      PQF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );

  return PQF12p6x11iK/pow(mass.getf0(),4);
}



double fnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQF12p6x11iL =
       + pi16 * (  - 2*pow(m11,2)*L1r - 4*(m44 - m11)*m11*L5r*nfm - 4*(
         m44 + m11)*m44*L6r*pow(nf,2) + 2*(m44 + m11)*m44*L4r*pow(nf,2)
          - 2*(m44 + m11)*(m44 + m11)*L8r*nf + (m44 + m11)*(m44 + m11)*
         L5r*nf - 1./4.*(m44 + m11)*(m44 + m11)*L3r*nf - 1./2.*(m44 + 
         m11)*(m44 + m11)*L0r*nf + 2*(4*m44 - 3*m11)*m11*L3r*nfm + 2*(4
         *m44 - 3*m11)*m11*L0r*nfm - ( - pow(m44,2) + pow(m44,2)*pow(
         nf,2) + pow(m11,2))*L2r );

      PQF12p6x11iL +=  + Ab(m11,mu2) * (  - 10*m11*L2r - 4*m11*L1r + 4*
         (m44 - 2*m11)*L5r*nfm - 12*(2*m44 - 3*m11)*L3r*nfm - 12*(2*m44
          - 3*m11)*L0r*nfm );

      PQF12p6x11iL +=  + Ab(m14,mu2) * (  - 2*m44*L4r*pow(nf,2) + 2*m11
         *L5r*nf + 4*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2) - 2*(m44
          + m11)*pow(m14,-1)*m44*L4r*pow(nf,2) - 5*(m44 + m11)*L3r*nf
          - 2*(m44 + m11)*L0r*nf + 2*(m44 + m11)*(m44 + m11)*pow(
         m14,-1)*L8r*nf - (m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r*nf );

      PQF12p6x11iL +=  + Ab(m44,mu2) * ( 4*( - 1 + nf)*(1 + nf)*m44*L4r
          - 2*( - 1 + nf)*(1 + nf)*m44*L2r - 8*( - 1 + nf)*(1 + nf)*m44
         *L1r );

      PQF12p6x11iL +=  - 8*pow(m44,2)*pow(L4r,2)*pow(nf,2) - 16*m11*m44
         *L4r*L5r*nf - 8*pow(m11,2)*pow(L5r,2);
 
  return PQF12p6x11iL/pow(mass.getf0(),4);
}

double fnfPQSUNp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSUNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQF12p6x11iR =
       + pow(pi16,2) * ( 1./8.*m14*m44*pow(nf,2) + 1./48.*m14*m44*pow(
         pi,2)*pow(nf,2) + 7./8.*pow(m11,2)*pow(nfm,2) + 1./96.*( - m44
          + m44*nf + m11)*(m44 + m44*nf - m11)*pow(pi,2) - 1./384.*(18*
         pow(m44,2) + pow(m44,2)*pow(nf,2) - 36*m11*m44 + 20*m11*m44*
         pow(nf,2) - 38*pow(m11,2)) );

      PQF12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 1./2.*m11*pow(nfm,2)
          + 1./24.*(3*m44 + 7*m11) );

      PQF12p6x11iR +=  + pow(Ab(m11,mu2),2) * (  - 1./16. + 1./8.*pow(
         m11,-1)*m44 );

      PQF12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1./4.*(m44 + 
         m11)*pow(m14,-1) );

      PQF12p6x11iR +=  + Ab(m14,mu2)*pi16 * (  - 1./24.*(4*m44 + m11)*
         pow(nf,2) );

      PQF12p6x11iR +=  + pow(Ab(m14,mu2),2) * ( 1./8.*pow(m14,-1)*m44*
         pow(nf,2) );

      PQF12p6x11iR +=  + pow(Ab(m44,mu2),2) * ( 1./16.*( - 1 + nf)*(1
          + nf) );

      PQF12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * (  - 1./8.*m44 );

      PQF12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * (  - 1./8.*( - 1
          + nf)*(1 + nf)*m44 );

      PQF12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * (  - 5./36.*(m44
          - m11)*m11 );

      PQF12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * ( 1./72.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hh1(3,m14,m11,m14,m11,mu2) * ( 1./36.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hhd(1,m11,m14,m14,m11,mu2) * ( 1./8.*(m44 - 4*
         m11)*m11 );

      PQF12p6x11iR +=  + hhd(1,m44,m14,m14,m11,mu2) * ( 1./8.*( - 1 + 
         nf)*(1 + nf)*m11*m44 );

      PQF12p6x11iR +=  + hhd(2,m11,m11,m11,m11,mu2) * (  - 2*(m44 - m11
         )*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11iR +=  + hhd(2,m11,m14,m14,m11,mu2) * ( 5./8.*(m44 - 
         m11)*pow(m11,2) );

      PQF12p6x11iR +=  + hhd(5,m11,m11,m11,m11,mu2) * ( (m44 - m11)*(
         m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11iR +=  + hh1d(2,m11,m14,m14,m11,mu2) * (  - (m44 - m11)
         *pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m11,m14,m14,m11,mu2) * ( 3./8.*pow(
         m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m44,m14,m14,m11,mu2) * ( 3./8.*( - 1
          + nf)*(1 + nf)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(2,m11,m14,m14,m11,mu2) * ( 3./8.*(m44 - 
         m11)*pow(m11,2) );

  return  PQF12p6x11iR/pow(mass.getf0(),4);
}

/////////////////////////////////////////////////////////////////////////////
// vev SUN
double qnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfPQSUNp4L(nf,mass,Liin)+qnfPQSUNp4R(nf,mass);
}

double qnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSUNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double   PQQ11p4x11iL =
       + 16*m44*L6r*nf + 4*m11*H2r + 8*m11*L8r;

  return PQQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfPQSUNp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQQ11p4x11iR =
       + pi16 * (  - (m44 - m11)*nfm );

      PQQ11p4x11iR +=  + Ab(m11,mu2) * ( pow(m11,-1)*m44*nfm - 2*nfm );

      PQQ11p4x11iR +=  + Ab(m14,mu2) * ( nf );

  return PQQ11p4x11iR/pow(mass.getf0(),2);
}

double qnfPQSUNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return qnfPQSUNp6K(nf,mass,Kiin)+qnfPQSUNp6L(nf,mass,Liin)
    +qnfPQSUNp6R(nf,mass);
}

double qnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfPQSUNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );

      PQQ11p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m44*nf );

      PQQ11p6x11iK +=  + KK[27] * ( 48*pow(m44,2)*pow(nf,2) );

  return PQQ11p6x11iK/pow(mass.getf0(),4);
}



double qnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSUNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQQ11p6x11iL =
       + pi16 * (  - 16*(m44 - m11)*(m44 - m11)*L7r - 4*(m44 + m11)*(
         m44 + m11)*L8r*nf + 2*(m44 + m11)*(m44 + m11)*L5r*nf + 8*(5*
         m44 - 6*m11)*m11*L5r*nfm - 8*(6*m44 + m44*pow(nf,2) - 8*m11 + 
         m11*pow(nf,2))*m44*L6r + 4*(6*m44 + m44*pow(nf,2) - 8*m11 + 
         m11*pow(nf,2))*m44*L4r - 16*(pow(m44,2) + 3*m11*m44 - 5*pow(
         m11,2))*L8r*nfm );

      PQQ11p6x11iL +=  + Ab(m11,mu2) * ( 16*pow(m11,-1)*pow(m44,2)*L8r*
         nfm + 16*pow(m11,-1)*pow(m44,2)*L7r + 32*pow(m11,-1)*pow(
         m44,2)*L6r - 16*pow(m11,-1)*pow(m44,2)*L4r - 64*m44*L6r + 32*
         m44*L4r + 32*(m44 - 3*m11)*L8r*nfm - 24*(2*m44 - 3*m11)*L5r*
         nfm - 16*(4*m44 - 3*m11)*L7r );

      PQQ11p6x11iL +=  + Ab(m14,mu2) * ( 16*m44*L6r*pow(nf,2) - 8*m44*
         L4r*pow(nf,2) + 8*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2) - 
         4*(m44 + m11)*pow(m14,-1)*m44*L4r*pow(nf,2) + 16*(m44 + m11)*
         L8r*nf - 8*(m44 + m11)*L5r*nf + 4*(m44 + m11)*(m44 + m11)*pow(
         m14,-1)*L8r*nf - 2*(m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r*nf
          );

      PQQ11p6x11iL +=  + Ab(m44,mu2) * ( 16*( - 1 + nf)*(1 + nf)*m44*
         L6r - 8*( - 1 + nf)*(1 + nf)*m44*L4r );
 
  return PQQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfPQSUNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSUNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double PQQ11p6x11iR =
       + pow(pi16,2) * ( 1./2.*(m44 - m11)*(3*m44 - 7*m11)*pow(nfm,2) )
         ;

      PQQ11p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 2*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 1./2.*(m44 - m11) + 10*(m44 - m11)*pow(
         nfm,2) );

      PQQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 3./2. + 1./2.*pow(
         m11,-2)*pow(m44,2)*pow(nfm,2) - 4*pow(m11,-1)*m44*pow(nfm,2)
          + 6*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 2 - pow(m11,-1)
         *m44 - 1./2.*(m44 + m11)*pow(m14,-1) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( pow(m11,-1)*m44 - 
         pow(m11,-1)*m44*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m14,mu2)*pi16 * ( (m44 + m11) );

      PQQ11p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - m44 + m44*pow(nfm,2) )
         ;

  return  PQQ11p6x11iR/pow(mass.getf0(),4);
}

/////////////SON case//////////////////////////////////////////////////////
// mass SON
double mnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfPQSONp4L(nf,mass,Liin)+mnfPQSONp4R(nf,mass);
}

double mnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSONp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double  PQM12p4x11iL =
       + 16*m11*m44*L6r*nf - 8*m11*m44*L4r*nf + 16*pow(m11,2)*L8r - 8*
         pow(m11,2)*L5r;

  return PQM12p4x11iL/pow(mass.getf0(),2);
}

double mnfPQSONp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQM12p4x11iR =
       + pi16 * (  - (m44 - m11)*m11*nfm );

      PQM12p4x11iR +=  + Ab(m11,mu2) * ( 1./2.*m11 + (m44 - 2*m11)*nfm
          );

  return PQM12p4x11iR/pow(mass.getf0(),2);
}

double mnfPQSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return mnfPQSONp6K(nf,mass,Kiin)+mnfPQSONp6L(nf,mass,Liin)
    +mnfPQSONp6R(nf,mass);
}

double mnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfPQSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQM12p6x11iK = 0.;
      PQM12p6x11iK +=  + KK[17] * (  - 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[18] * (  - 32*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[19] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[20] * (  - 16*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[21] * (  - 16*m11*pow(m44,2)*nf );
      PQM12p6x11iK +=  + KK[22] * (  - 16*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[23] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[25] * ( 48*pow(m11,3) );
      PQM12p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m11*m44*nf );
      PQM12p6x11iK +=  + KK[27] * ( 48*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[39] * ( 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[40] * ( 32*pow(m11,2)*m44*nf );
 
  return PQM12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double mnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQM12p6x11iL =
       + pi16 * (  - 8*pow(m11,3)*L8r + 4*pow(m11,3)*L5r + 4*pow(m11,3)
         *L1r - 16*(m44 - m11)*(m44 - m11)*m11*L7r - 4*(4*m44 - 3*m11)*
         pow(m11,2)*L3r*nfm - 4*(4*m44 - 3*m11)*pow(m11,2)*L0r*nfm - 8*
         (6*m44 - 8*m11 + m11*nf)*m11*m44*L6r + 8*(7*m44 - 8*m11)*pow(
         m11,2)*L5r*nfm + 4*(8*m44 - 10*m11 + m11*nf)*m11*m44*L4r + (
          - 2*pow(m44,2) + pow(m44,2)*nf + pow(m44,2)*pow(nf,2) + 2*
         pow(m11,2))*m11*L2r - 16*(pow(m44,2) + 5*m11*m44 - 7*pow(
         m11,2))*m11*L8r*nfm + 1./2.*(pow(m44,2)*nf + 2*m11*m44*nf + 4*
         pow(m11,2) + pow(m11,2)*nf)*m11*L0r + 1./4.*(pow(m44,2)*nf + 2
         *m11*m44*nf + 8*pow(m11,2) + pow(m11,2)*nf)*m11*L3r );

      PQM12p6x11iL +=  + Ab(m11,mu2) * ( 40*pow(m11,2)*L8r - 20*pow(
         m11,2)*L5r + 12*pow(m11,2)*L3r + 20*pow(m11,2)*L2r + 8*pow(
         m11,2)*L1r + 12*pow(m11,2)*L0r + 16*(m44 - m11)*(m44 - 3*m11)*
         L7r - 12*(2*m44 - 4*m11 + m11*nf)*m44*L4r + 24*(2*m44 - 3*m11)
         *m11*L3r*nfm + 24*(2*m44 - 3*m11)*m11*L0r*nfm - 8*(8*m44 - 13*
         m11)*m11*L5r*nfm + 16*(pow(m44,2) + 4*m11*m44 - 10*pow(m11,2))
         *L8r*nfm + 16*(2*pow(m44,2) - 4*m11*m44 + m11*m44*nf + 2*pow(
         m11,2))*L6r );

      PQM12p6x11iL +=  + Ab(m14,mu2) * ( 8*(m44 + m11)*m11*L8r*nf - 4*(
         m44 + m11)*m11*L5r*nf + 5*(m44 + m11)*m11*L3r*nf + 2*(m44 + 
         m11)*m11*L0r*nf );

      PQM12p6x11iL +=  + Ab(m44,mu2) * ( 8*( - 1 + nf)*(2 + nf)*m11*m44
         *L6r - 8*( - 1 + nf)*(2 + nf)*m11*m44*L4r + 2*( - 1 + nf)*(2
          + nf)*m11*m44*L2r + 8*( - 1 + nf)*(2 + nf)*m11*m44*L1r );

      PQM12p6x11iL +=  - 128*m11*pow(m44,2)*L4r*L6r*pow(nf,2) + 64*m11*
         pow(m44,2)*pow(L4r,2)*pow(nf,2) - 128*pow(m11,2)*m44*L5r*L6r*
         nf - 128*pow(m11,2)*m44*L4r*L8r*nf + 128*pow(m11,2)*m44*L4r*
         L5r*nf - 128*pow(m11,3)*L5r*L8r + 64*pow(m11,3)*pow(L5r,2);
 
  return PQM12p6x11iL/pow(mass.getf0(),4);
}

double mnfPQSONp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSONp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQM12p6x11iR =
       + pow(pi16,2) * (  - 3./8.*pow(m11,3)*nfm - 1./12.*(m44 - m11)*
         pow(m11,2)*pow(pi,2)*nfm + 1./3.*(m44 - m11)*pow(m11,2)*pow(
         pi,2)*pow(nfm,2) - 1./16.*(m44*nf + m44*pow(nf,2) - 8*m11 + 3*
         m11*nf + m11*pow(nf,2))*m11*m14 - 1./96.*(m44*nf + m44*pow(
         nf,2) - 8*m11 + 3*m11*nf + m11*pow(nf,2))*m11*m14*pow(pi,2) - 
         1./2304.*( - 108*pow(m44,2) - 6*pow(m44,2)*nf - 3*pow(m44,2)*
         pow(nf,2) + 792*m11*m44 - 76*m11*m44*nf - 26*m11*m44*pow(nf,2)
          - 1548*pow(m11,2) + 320*pow(m11,2)*nf + 4*pow(m11,2)*pow(
         nf,2))*m11 - 1./192.*( - 2*pow(m44,2) + pow(m44,2)*nf + pow(
         m44,2)*pow(nf,2) + 12*m11*m44 - 18*pow(m11,2) + 5*pow(m11,2)*
         nf)*m11*pow(pi,2) + 3./4.*(2*pow(m44,2) - 4*m11*m44 + 3*pow(
         m11,2))*m11*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - m11*m44*nfm + 1./24.*(
         15*m44 - 13*m11)*m11 - (2*pow(m44,2) - 6*m11*m44 + pow(m11,2))
         *pow(nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 1./2.*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 1./2.*(m44 - 8*m11)*nfm - (2*m44 - 7*m11)*
         pow(nfm,2) - 1./32.*(12*m44 - 70*m11 + 5*m11*nf) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1./4.*(4*m44 + 
         m11*nf) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( 1./2.*m44 + 1./2.*
         m44*nfm - m44*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m14,mu2)*pi16 * ( 1./48.*(48*m44 + 5*m44*nf
          + 4*m44*pow(nf,2) + 3*m11*nf)*m11 );

      PQM12p6x11iR +=  + pow(Ab(m14,mu2),2) * (  - 1./16.*m11*pow(nf,2)
          - 1./16.*(m44*nf + m44*pow(nf,2) - 8*m11 + 3*m11*nf + m11*
         pow(nf,2))*m11*pow(m14,-1) );

      PQM12p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - 1./2.*m11*m44 - 1./2.*
         m11*m44*nfm + m11*m44*pow(nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m44,mu2),2) * (  - 1./32.*( - 1 + nf)*(
         2 + nf)*m11 );

      PQM12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * ( 1./16.*(2*m44 - 
         8*m11 + 5*m11*nf)*m11 );

      PQM12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * ( 1./16.*( - 1 + 
         nf)*(2 + nf)*m11*m44 );

      PQM12p6x11iR +=  + hh(2,m11,m11,m11,m11,mu2) * ( (m44 - m11)*pow(
         m11,2)*nfm - 4*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * ( 5./8.*(m44 - m11
         )*pow(m11,2) );

      PQM12p6x11iR +=  + hh(5,m11,m11,m11,m11,mu2) * ( 2*(m44 - m11)*(
         m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh1(1,m11,m14,m14,m11,mu2) * (  - 1./2.*pow(
         m11,2)*nf );

      PQM12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * (  - (m44 - m11)*
         pow(m11,2) );

      PQM12p6x11iR +=  + hh21(1,m11,m14,m14,m11,mu2) * ( 3./16.*(2 + nf
         )*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(1,m44,m14,m14,m11,mu2) * ( 3./16.*( - 1
          + nf)*(2 + nf)*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(2,m11,m14,m14,m11,mu2) * ( 3./8.*(m44 - 
         m11)*pow(m11,2) );


  return  PQM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////////////
// decay SON
double fnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfPQSONp4L(nf,mass,Liin)+fnfPQSONp4R(nf,mass);
}

double fnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSONp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   PQF12p4x11iL =
       + 4*m44*L4r*nf + 4*m11*L5r;

  return PQF12p4x11iL/pow(mass.getf0(),2);
}

double fnfPQSONp4R(const int nf,const quarkmassnf mass){
  //const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQF12p4x11iR =
       + Ab(m14,mu2) * ( 1./4.*nf );

  return PQF12p4x11iR/pow(mass.getf0(),2);
}

double fnfPQSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return fnfPQSONp6K(nf,mass,Kiin)+fnfPQSONp6L(nf,mass,Liin)
    +fnfPQSONp6R(nf,mass);
}

double fnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfPQSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double  PQF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );
      PQF12p6x11iK +=  + KK[20] * ( 8*m11*m44*nf );
      PQF12p6x11iK +=  + KK[21] * ( 8*pow(m44,2)*nf );
      PQF12p6x11iK +=  + KK[22] * ( 8*pow(m44,2)*pow(nf,2) );
      PQF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );
      
  return PQF12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double fnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQF12p6x11iL =
       + pi16 * (  - 2*pow(m11,2)*L1r - 4*(m44 - m11)*m11*L5r*nfm - 2*(
         m44 + m11)*m44*L6r*pow(nf,2) + (m44 + m11)*m44*L4r*pow(nf,2)
          - (m44 + m11)*(m44 + m11)*L8r*nf + 1./2.*(m44 + m11)*(m44 + 
         m11)*L5r*nf + 2*(4*m44 - 3*m11)*m11*L3r*nfm + 2*(4*m44 - 3*m11
         )*m11*L0r*nfm - 1./2.*( - 2*pow(m44,2) + pow(m44,2)*nf + pow(
         m44,2)*pow(nf,2) + 2*pow(m11,2))*L2r - 1./4.*(pow(m44,2)*nf + 
         2*m11*m44*nf + 4*pow(m11,2) + pow(m11,2)*nf)*L0r - 1./8.*(pow(
         m44,2)*nf + 2*m11*m44*nf + 8*pow(m11,2) + pow(m11,2)*nf)*L3r )
         ;

      PQF12p6x11iL +=  + Ab(m11,mu2) * ( 2*m11*L5r - 6*m11*L3r - 10*m11
         *L2r - 4*m11*L1r - 6*m11*L0r + 4*(m44 - 2*m11)*L5r*nfm - 12*(2
         *m44 - 3*m11)*L3r*nfm - 12*(2*m44 - 3*m11)*L0r*nfm );

      PQF12p6x11iL +=  + Ab(m14,mu2) * (  - m44*L4r*pow(nf,2) + m11*L5r
         *nf + 2*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2) - (m44 + m11
         )*pow(m14,-1)*m44*L4r*pow(nf,2) - 5./2.*(m44 + m11)*L3r*nf - (
         m44 + m11)*L0r*nf + (m44 + m11)*(m44 + m11)*pow(m14,-1)*L8r*nf
          - 1./2.*(m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r*nf );

      PQF12p6x11iL +=  + Ab(m44,mu2) * ( 2*( - 1 + nf)*(2 + nf)*m44*L4r
          - ( - 1 + nf)*(2 + nf)*m44*L2r - 4*( - 1 + nf)*(2 + nf)*m44*
         L1r );

      PQF12p6x11iL +=  - 8*pow(m44,2)*pow(L4r,2)*pow(nf,2) - 16*m11*m44
         *L4r*L5r*nf - 8*pow(m11,2)*pow(L5r,2);
 
  return PQF12p6x11iL/pow(mass.getf0(),4);
}

double fnfPQSONp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSONp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQF12p6x11iR =
       + pow(pi16,2) * (  - 7./16.*pow(m11,2)*nfm + 7./8.*pow(m11,2)*
         pow(nfm,2) + 1./32.*(m44 + m44*nf + m11)*m14*nf + 1./192.*(m44
          + m44*nf + m11)*m14*pow(pi,2)*nf + 1./384.*( - 2*pow(m44,2)
          + pow(m44,2)*nf + pow(m44,2)*pow(nf,2) + 4*m11*m44 - 2*pow(
         m11,2) + pow(m11,2)*nf)*pow(pi,2) - 1./1536.*(36*pow(m44,2) + 
         2*pow(m44,2)*nf + pow(m44,2)*pow(nf,2) - 72*m11*m44 + 40*m11*
         m44*nf + 20*m11*m44*pow(nf,2) - 188*pow(m11,2)) );

      PQF12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 1./4.*m11*nfm - 1./2.*m11
         *pow(nfm,2) + 1./48.*(3*m44 + 5*m11) );

      PQF12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 1./16.*pow(m11,-1)*m44
          + 1./64.*( - 2 + nf) );

      PQF12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1./8.*(m44 + 
         m11)*pow(m14,-1) );

      PQF12p6x11iR +=  + Ab(m14,mu2)*pi16 * (  - 1./96.*(11*m44 + 4*m44
         *nf + 11*m11 + m11*nf)*nf );

      PQF12p6x11iR +=  + pow(Ab(m14,mu2),2) * ( 1./16.*nf + 1./32.*(3*
         m44 + m44*nf + 3*m11)*pow(m14,-1)*nf );

      PQF12p6x11iR +=  + pow(Ab(m44,mu2),2) * ( 1./64.*( - 1 + nf)*(2
          + nf) );

      PQF12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * (  - 1./32.*(2*m44
          + m11*nf) );

      PQF12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * (  - 1./32.*( - 1
          + nf)*(2 + nf)*m44 );

      PQF12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * (  - 5./72.*(m44
          - m11)*m11 );

      PQF12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * ( 1./144.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hh1(3,m14,m11,m14,m11,mu2) * ( 1./72.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hhd(1,m11,m14,m14,m11,mu2) * ( 1./32.*(2*m44
          - 8*m11 + 5*m11*nf)*m11 );

      PQF12p6x11iR +=  + hhd(1,m44,m14,m14,m11,mu2) * ( 1./32.*( - 1 + 
         nf)*(2 + nf)*m11*m44 );

      PQF12p6x11iR +=  + hhd(2,m11,m11,m11,m11,mu2) * ( 1./2.*(m44 - 
         m11)*pow(m11,2)*nfm - 2*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11iR +=  + hhd(2,m11,m14,m14,m11,mu2) * ( 5./16.*(m44 - 
         m11)*pow(m11,2) );

      PQF12p6x11iR +=  + hhd(5,m11,m11,m11,m11,mu2) * ( (m44 - m11)*(
         m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11iR +=  + hh1d(1,m11,m14,m14,m11,mu2) * (  - 1./4.*pow(
         m11,2)*nf );

      PQF12p6x11iR +=  + hh1d(2,m11,m14,m14,m11,mu2) * (  - 1./2.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m11,m14,m14,m11,mu2) * ( 3./32.*(2 + 
         nf)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m44,m14,m14,m11,mu2) * ( 3./32.*( - 1
          + nf)*(2 + nf)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(2,m11,m14,m14,m11,mu2) * ( 3./16.*(m44
          - m11)*pow(m11,2) );

  return  PQF12p6x11iR/pow(mass.getf0(),4);
}

/////////////////////////////////////////////////////////////////////////////
// vev SON
double qnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfPQSONp4L(nf,mass,Liin)+qnfPQSONp4R(nf,mass);
}

double qnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSONp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double  PQQ11p4x11iL =
       + 16*m44*L6r*nf + 4*m11*H2r + 8*m11*L8r;

  return PQQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfPQSONp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQQ11p4x11iR =
       + pi16 * (  - (m44 - m11)*nfm );

      PQQ11p4x11iR +=  + Ab(m11,mu2) * ( 1./2. + pow(m11,-1)*m44*nfm - 
         2*nfm );

      PQQ11p4x11iR +=  + Ab(m14,mu2) * ( 1./2.*nf );

  return PQQ11p4x11iR/pow(mass.getf0(),2);
}

double qnfPQSONp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return qnfPQSONp6K(nf,mass,Kiin)+qnfPQSONp6L(nf,mass,Liin)
    +qnfPQSONp6R(nf,mass);
}

double qnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfPQSONp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );

      PQQ11p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m44*nf );

      PQQ11p6x11iK +=  + KK[27] * ( 48*pow(m44,2)*pow(nf,2) );

  return PQQ11p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double qnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSONp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQQ11p6x11iL =
       + pi16 * (  - 16*(m44 - m11)*(m44 - m11)*L7r + 8*(5*m44 - 6*m11)
         *m11*L5r*nfm - 4*(12*m44 + m44*pow(nf,2) - 16*m11 + 2*m11*nf
          + m11*pow(nf,2))*m44*L6r + 2*(12*m44 + m44*pow(nf,2) - 16*m11
          + 2*m11*nf + m11*pow(nf,2))*m44*L4r - 16*(pow(m44,2) + 3*m11*
         m44 - 5*pow(m11,2))*L8r*nfm - 2*(pow(m44,2)*nf + 2*m11*m44*nf
          + 4*pow(m11,2) + pow(m11,2)*nf)*L8r + (pow(m44,2)*nf + 2*m11*
         m44*nf + 4*pow(m11,2) + pow(m11,2)*nf)*L5r );

      PQQ11p6x11iL +=  + Ab(m11,mu2) * ( 16*pow(m11,-1)*pow(m44,2)*L8r*
         nfm + 16*pow(m11,-1)*pow(m44,2)*L7r + 32*pow(m11,-1)*pow(
         m44,2)*L6r - 16*pow(m11,-1)*pow(m44,2)*L4r + 24*m11*L8r - 12*
         m11*L5r + 16*( - 4 + nf)*m44*L6r - 8*( - 4 + nf)*m44*L4r + 32*
         (m44 - 3*m11)*L8r*nfm - 24*(2*m44 - 3*m11)*L5r*nfm - 16*(4*m44
          - 3*m11)*L7r );

      PQQ11p6x11iL +=  + Ab(m14,mu2) * ( 8*m44*L6r*pow(nf,2) - 4*m44*
         L4r*pow(nf,2) + 4*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2) - 
         2*(m44 + m11)*pow(m14,-1)*m44*L4r*pow(nf,2) + 8*(m44 + m11)*
         L8r*nf - 4*(m44 + m11)*L5r*nf + 2*(m44 + m11)*(m44 + m11)*pow(
         m14,-1)*L8r*nf - (m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r*nf );

      PQQ11p6x11iL +=  + Ab(m44,mu2) * ( 8*( - 1 + nf)*(2 + nf)*m44*L6r
          - 4*( - 1 + nf)*(2 + nf)*m44*L4r );
 
  return PQQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfPQSONp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSONp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQQ11p6x11iR =
       + pow(pi16,2) * ( 1./2.*(m44 - m11)*m11*nfm + 1./2.*(m44 - m11)*
         (3*m44 - 7*m11)*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 2*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 1./4.*(m44 - 2*m11) + 10*(m44 - m11)*pow(
         nfm,2) - 1./2.*(4*m44 - 7*m11)*nfm );

      PQQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 9./8. + 1./2.*pow(
         m11,-2)*pow(m44,2)*pow(nfm,2) + pow(m11,-1)*m44*nfm - 4*pow(
         m11,-1)*m44*pow(nfm,2) - 3*nfm + 6*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1 - 1./2.*pow(
         m11,-1)*m44 - 1./4.*(m44 + m11)*pow(m14,-1) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( 1./2.*pow(m11,-1)*
         m44 + 1./2.*pow(m11,-1)*m44*nfm - pow(m11,-1)*m44*pow(nfm,2) )
         ;

      PQQ11p6x11iR +=  + Ab(m14,mu2)*pi16 * (  - 1./8.*( - 4 + nf)*(m44
          + m11) );

      PQQ11p6x11iR +=  + pow(Ab(m14,mu2),2) * ( 1./8.*nf + 1./8.*(m44
          + m11)*pow(m14,-1)*nf );

      PQQ11p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - 1./2.*m44 - 1./2.*m44*
         nfm + m44*pow(nfm,2) );

  return  PQQ11p6x11iR/pow(mass.getf0(),4);
}

/////////////SPN case//////////////////////////////////////////////////////
// mass SPN
double mnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  mnfPQSPNp4L(nf,mass,Liin)+mnfPQSPNp4R(nf,mass);
}

double mnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSPNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   PQM12p4x11iL =
       + 32*m11*m44*L6r*nf - 16*m11*m44*L4r*nf + 16*pow(m11,2)*L8r - 8*
         pow(m11,2)*L5r;

  return PQM12p4x11iL/pow(mass.getf0(),2);
}

double mnfPQSPNp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double   PQM12p4x11iR =
       + pi16 * (  - 1./2.*(m44 - m11)*m11*nfm );

      PQM12p4x11iR +=  + Ab(m11,mu2) * (  - 1./2.*m11 + 1./2.*(m44 - 2*
         m11)*nfm );

  return PQM12p4x11iR/pow(mass.getf0(),2);
}

double mnfPQSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return mnfPQSPNp6K(nf,mass,Kiin)+mnfPQSPNp6L(nf,mass,Liin)
    +mnfPQSPNp6R(nf,mass);
}

double mnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in mnfPQSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQM12p6x11iK = 0.;
      PQM12p6x11iK +=  + KK[17] * (  - 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[18] * (  - 32*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[19] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[20] * (  - 16*pow(m11,2)*m44*nf );
      PQM12p6x11iK +=  + KK[21] * (  - 16*m11*pow(m44,2)*nf );
      PQM12p6x11iK +=  + KK[22] * (  - 16*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[23] * (  - 16*pow(m11,3) );
      PQM12p6x11iK +=  + KK[25] * ( 48*pow(m11,3) );
      PQM12p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m11*m44*nf );
      PQM12p6x11iK +=  + KK[27] * ( 48*m11*pow(m44,2)*pow(nf,2) );
      PQM12p6x11iK +=  + KK[39] * ( 32*pow(m11,3) );
      PQM12p6x11iK +=  + KK[40] * ( 32*pow(m11,2)*m44*nf );
 
  return PQM12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double mnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in mnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQM12p6x11iL =
       + pi16 * ( 8*pow(m11,3)*L8r - 4*pow(m11,3)*L5r + 4*pow(m11,3)*
         L1r - 8*( - 4*m44 + 5*m11 + m11*nf)*m11*m44*L4r + 16*( - 3*m44
          + 4*m11 + m11*nf)*m11*m44*L6r - 16*(m44 - m11)*(m44 - m11)*
         m11*L7r - 2*(4*m44 - 3*m11)*pow(m11,2)*L3r*nfm - 2*(4*m44 - 3*
         m11)*pow(m11,2)*L0r*nfm + 4*(7*m44 - 8*m11)*pow(m11,2)*L5r*nfm
          + 2*( - pow(m44,2) - pow(m44,2)*nf + 2*pow(m44,2)*pow(nf,2)
          + pow(m11,2))*m11*L2r - 8*(pow(m44,2) + 5*m11*m44 - 7*pow(
         m11,2))*m11*L8r*nfm + 1./2.*(pow(m44,2)*nf + 2*m11*m44*nf - 4*
         pow(m11,2) + pow(m11,2)*nf)*m11*L3r + (pow(m44,2)*nf + 2*m11*
         m44*nf - 2*pow(m11,2) + pow(m11,2)*nf)*m11*L0r );

      PQM12p6x11iL +=  + Ab(m11,mu2) * (  - 40*pow(m11,2)*L8r + 20*pow(
         m11,2)*L5r - 12*pow(m11,2)*L3r + 20*pow(m11,2)*L2r + 8*pow(
         m11,2)*L1r - 12*pow(m11,2)*L0r + 24*( - m44 + 2*m11 + m11*nf)*
         m44*L4r + 16*(m44 - m11)*(m44 - 3*m11)*L7r + 12*(2*m44 - 3*m11
         )*m11*L3r*nfm + 12*(2*m44 - 3*m11)*m11*L0r*nfm - 4*(8*m44 - 13
         *m11)*m11*L5r*nfm - 32*( - pow(m44,2) + 2*m11*m44 + m11*m44*nf
          - pow(m11,2))*L6r + 8*(pow(m44,2) + 4*m11*m44 - 10*pow(m11,2)
         )*L8r*nfm );

      PQM12p6x11iL +=  + Ab(m14,mu2) * ( 16*(m44 + m11)*m11*L8r*nf - 8*
         (m44 + m11)*m11*L5r*nf + 10*(m44 + m11)*m11*L3r*nf + 4*(m44 + 
         m11)*m11*L0r*nf );

      PQM12p6x11iL +=  + Ab(m44,mu2) * ( 16*( - 1 + nf)*(1 + 2*nf)*m11*
         m44*L6r - 16*( - 1 + nf)*(1 + 2*nf)*m11*m44*L4r + 4*( - 1 + nf
         )*(1 + 2*nf)*m11*m44*L2r + 16*( - 1 + nf)*(1 + 2*nf)*m11*m44*
         L1r );

      PQM12p6x11iL +=  - 512*m11*pow(m44,2)*L4r*L6r*pow(nf,2) + 256*m11
         *pow(m44,2)*pow(L4r,2)*pow(nf,2) - 256*pow(m11,2)*m44*L5r*L6r*
         nf - 256*pow(m11,2)*m44*L4r*L8r*nf + 256*pow(m11,2)*m44*L4r*
         L5r*nf - 128*pow(m11,3)*L5r*L8r + 64*pow(m11,3)*pow(L5r,2);
 
  return PQM12p6x11iL/pow(mass.getf0(),4);
}

double mnfPQSPNp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: mnfPQSPNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQM12p6x11iR =
       + pow(pi16,2) * ( 3./16.*pow(m11,3)*nfm + 1./24.*(m44 - m11)*
         pow(m11,2)*pow(pi,2)*nfm + 1./12.*(m44 - m11)*pow(m11,2)*pow(
         pi,2)*pow(nfm,2) - 1./8.*( - m44*nf + 2*m44*pow(nf,2) - 4*m11
          - 3*m11*nf + 2*m11*pow(nf,2))*m11*m14 - 1./48.*( - m44*nf + 2
         *m44*pow(nf,2) - 4*m11 - 3*m11*nf + 2*m11*pow(nf,2))*m11*m14*
         pow(pi,2) - 1./576.*( - 27*pow(m44,2) + 3*pow(m44,2)*nf - 3*
         pow(m44,2)*pow(nf,2) + 198*m11*m44 + 38*m11*m44*nf - 26*m11*
         m44*pow(nf,2) - 387*pow(m11,2) - 160*pow(m11,2)*nf + 4*pow(
         m11,2)*pow(nf,2))*m11 - 1./96.*( - pow(m44,2) - pow(m44,2)*nf
          + 2*pow(m44,2)*pow(nf,2) + 6*m11*m44 - 9*pow(m11,2) - 5*pow(
         m11,2)*nf)*m11*pow(pi,2) + 3./16.*(2*pow(m44,2) - 4*m11*m44 + 
         3*pow(m11,2))*m11*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*pi16 * ( 1./2.*m11*m44*nfm + 1./24.
         *(15*m44 - 13*m11)*m11 - 1./4.*(2*pow(m44,2) - 6*m11*m44 + 
         pow(m11,2))*pow(nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 1./8.*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 1./16.*( - 6*m44 + 35*m11 + 5*m11*nf) - 1./
         4.*(m44 - 8*m11)*nfm - 1./4.*(2*m44 - 7*m11)*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * ( 1./2.*( - 2*m44 + 
         m11*nf) );

      PQM12p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( 1./2.*m44 - 1./4.*
         m44*nfm - 1./4.*m44*pow(nfm,2) );

      PQM12p6x11iR +=  + Ab(m14,mu2)*pi16 * ( 1./24.*(24*m44 - 5*m44*nf
          + 8*m44*pow(nf,2) - 3*m11*nf)*m11 );

      PQM12p6x11iR +=  + pow(Ab(m14,mu2),2) * (  - 1./4.*m11*pow(nf,2)
          - 1./8.*( - m44*nf + 2*m44*pow(nf,2) - 4*m11 - 3*m11*nf + 2*
         m11*pow(nf,2))*m11*pow(m14,-1) );

      PQM12p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - 1./2.*m11*m44 + 1./4.*
         m11*m44*nfm + 1./4.*m11*m44*pow(nfm,2) );

      PQM12p6x11iR +=  + pow(Ab(m44,mu2),2) * (  - 1./16.*( - 1 + nf)*(
         1 + 2*nf)*m11 );

      PQM12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * (  - 1./8.*( - m44
          + 4*m11 + 5*m11*nf)*m11 );

      PQM12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * ( 1./8.*( - 1 + nf
         )*(1 + 2*nf)*m11*m44 );

      PQM12p6x11iR +=  + hh(2,m11,m11,m11,m11,mu2) * (  - 1./2.*(m44 - 
         m11)*pow(m11,2)*nfm - (m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * ( 5./8.*(m44 - m11
         )*pow(m11,2) );

      PQM12p6x11iR +=  + hh(5,m11,m11,m11,m11,mu2) * ( 1./2.*(m44 - m11
         )*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQM12p6x11iR +=  + hh1(1,m11,m14,m14,m11,mu2) * ( pow(m11,2)*nf )
         ;

      PQM12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * (  - (m44 - m11)*
         pow(m11,2) );

      PQM12p6x11iR +=  + hh21(1,m11,m14,m14,m11,mu2) * (  - 3./8.*( - 1
          + nf)*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(1,m44,m14,m14,m11,mu2) * ( 3./8.*( - 1 + 
         nf)*(1 + 2*nf)*pow(m11,2) );

      PQM12p6x11iR +=  + hh21(2,m11,m14,m14,m11,mu2) * ( 3./8.*(m44 - 
         m11)*pow(m11,2) );

  return  PQM12p6x11iR/pow(mass.getf0(),4);
}
/////////////////////////////////////////////////////////////////////////////
// decay SPN
double fnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  fnfPQSPNp4L(nf,mass,Liin)+fnfPQSPNp4R(nf,mass);
}

double fnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSPNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double   PQF12p4x11iL =
       + 8*m44*L4r*nf + 4*m11*L5r;

  return PQF12p4x11iL/pow(mass.getf0(),2);
}

double fnfPQSPNp4R(const int nf,const quarkmassnf mass){
  //const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  //double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double  PQF12p4x11iR =
       + Ab(m14,mu2) * ( 1./2.*nf );

  return PQF12p4x11iR/pow(mass.getf0(),2);
}

double fnfPQSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return fnfPQSPNp6K(nf,mass,Kiin)+fnfPQSPNp6L(nf,mass,Liin)
    +fnfPQSPNp6R(nf,mass);
}

double fnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in fnfPQSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double  PQF12p6x11iK =
       + KK[19] * ( 8*pow(m11,2) );
      PQF12p6x11iK +=  + KK[20] * ( 8*m11*m44*nf );
      PQF12p6x11iK +=  + KK[21] * ( 8*pow(m44,2)*nf );
      PQF12p6x11iK +=  + KK[22] * ( 8*pow(m44,2)*pow(nf,2) );
      PQF12p6x11iK +=  + KK[23] * ( 8*pow(m11,2) );
      
  return PQF12p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double fnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in fnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQF12p6x11iL =
       + pi16 * (  - 2*pow(m11,2)*L1r - 2*(m44 - m11)*m11*L5r*nfm - 8*(
         m44 + m11)*m44*L6r*pow(nf,2) + 4*(m44 + m11)*m44*L4r*pow(nf,2)
          - 2*(m44 + m11)*(m44 + m11)*L8r*nf + (m44 + m11)*(m44 + m11)*
         L5r*nf + (4*m44 - 3*m11)*m11*L3r*nfm + (4*m44 - 3*m11)*m11*L0r
         *nfm - ( - pow(m44,2) - pow(m44,2)*nf + 2*pow(m44,2)*pow(nf,2)
          + pow(m11,2))*L2r - 1./4.*(pow(m44,2)*nf + 2*m11*m44*nf - 4*
         pow(m11,2) + pow(m11,2)*nf)*L3r - 1./2.*(pow(m44,2)*nf + 2*m11
         *m44*nf - 2*pow(m11,2) + pow(m11,2)*nf)*L0r );

      PQF12p6x11iL +=  + Ab(m11,mu2) * (  - 2*m11*L5r + 6*m11*L3r - 10*
         m11*L2r - 4*m11*L1r + 6*m11*L0r + 2*(m44 - 2*m11)*L5r*nfm - 6*
         (2*m44 - 3*m11)*L3r*nfm - 6*(2*m44 - 3*m11)*L0r*nfm );

      PQF12p6x11iL +=  + Ab(m14,mu2) * (  - 4*m44*L4r*pow(nf,2) + 2*m11
         *L5r*nf + 8*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2) - 4*(m44
          + m11)*pow(m14,-1)*m44*L4r*pow(nf,2) - 5*(m44 + m11)*L3r*nf
          - 2*(m44 + m11)*L0r*nf + 2*(m44 + m11)*(m44 + m11)*pow(
         m14,-1)*L8r*nf - (m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r*nf );

      PQF12p6x11iL +=  + Ab(m44,mu2) * ( 4*( - 1 + nf)*(1 + 2*nf)*m44*
         L4r - 2*( - 1 + nf)*(1 + 2*nf)*m44*L2r - 8*( - 1 + nf)*(1 + 2*
         nf)*m44*L1r );

      PQF12p6x11iL +=  - 32*pow(m44,2)*pow(L4r,2)*pow(nf,2) - 32*m11*
         m44*L4r*L5r*nf - 8*pow(m11,2)*pow(L5r,2);
 
  return PQF12p6x11iL/pow(mass.getf0(),4);
}

double fnfPQSPNp6R(const int nf,const quarkmassnf mass){
  const double pi = M_PI;
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: fnfPQSPNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQF12p6x11iR =
       + pow(pi16,2) * ( 7./32.*pow(m11,2)*nfm + 7./32.*pow(m11,2)*pow(
         nfm,2) + 1./16.*( - m44 + 2*m44*nf - m11)*m14*nf + 1./96.*( - 
         m44 + 2*m44*nf - m11)*m14*pow(pi,2)*nf + 1./192.*( - pow(
         m44,2) - pow(m44,2)*nf + 2*pow(m44,2)*pow(nf,2) + 2*m11*m44 - 
         pow(m11,2) - pow(m11,2)*nf)*pow(pi,2) - 1./384.*(9*pow(m44,2)
          - pow(m44,2)*nf + pow(m44,2)*pow(nf,2) - 18*m11*m44 - 20*m11*
         m44*nf + 20*m11*m44*pow(nf,2) - 47*pow(m11,2)) );

      PQF12p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 1./8.*m11*nfm - 1./8.*
         m11*pow(nfm,2) + 1./48.*(3*m44 + 5*m11) );

      PQF12p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 1./16.*pow(m11,-1)*m44
          - 1./32.*(1 + nf) );

      PQF12p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1./8.*(m44 + 
         m11)*pow(m14,-1) );

      PQF12p6x11iR +=  + Ab(m14,mu2)*pi16 * (  - 1./48.*( - 11*m44 + 8*
         m44*nf - 11*m11 + 2*m11*nf)*nf );

      PQF12p6x11iR +=  + pow(Ab(m14,mu2),2) * (  - 1./8.*nf + 1./16.*(
          - 3*m44 + 2*m44*nf - 3*m11)*pow(m14,-1)*nf );

      PQF12p6x11iR +=  + pow(Ab(m44,mu2),2) * ( 1./32.*( - 1 + nf)*(1
          + 2*nf) );

      PQF12p6x11iR +=  + hh(1,m11,m14,m14,m11,mu2) * ( 1./16.*( - m44
          + m11*nf) );

      PQF12p6x11iR +=  + hh(1,m44,m14,m14,m11,mu2) * (  - 1./16.*( - 1
          + nf)*(1 + 2*nf)*m44 );

      PQF12p6x11iR +=  + hh(2,m11,m14,m14,m11,mu2) * (  - 5./72.*(m44
          - m11)*m11 );

      PQF12p6x11iR +=  + hh1(2,m11,m14,m14,m11,mu2) * ( 1./144.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hh1(3,m14,m11,m14,m11,mu2) * ( 1./72.*(m44 - 
         m11)*m11 );

      PQF12p6x11iR +=  + hhd(1,m11,m14,m14,m11,mu2) * (  - 1./16.*( - 
         m44 + 4*m11 + 5*m11*nf)*m11 );

      PQF12p6x11iR +=  + hhd(1,m44,m14,m14,m11,mu2) * ( 1./16.*( - 1 + 
         nf)*(1 + 2*nf)*m11*m44 );

      PQF12p6x11iR +=  + hhd(2,m11,m11,m11,m11,mu2) * (  - 1./4.*(m44
          - m11)*pow(m11,2)*nfm - 1./2.*(m44 - m11)*pow(m11,2)*pow(
         nfm,2) );

      PQF12p6x11iR +=  + hhd(2,m11,m14,m14,m11,mu2) * ( 5./16.*(m44 - 
         m11)*pow(m11,2) );

      PQF12p6x11iR +=  + hhd(5,m11,m11,m11,m11,mu2) * ( 1./4.*(m44 - 
         m11)*(m44 - m11)*pow(m11,2)*pow(nfm,2) );

      PQF12p6x11iR +=  + hh1d(1,m11,m14,m14,m11,mu2) * ( 1./2.*pow(
         m11,2)*nf );

      PQF12p6x11iR +=  + hh1d(2,m11,m14,m14,m11,mu2) * (  - 1./2.*(m44
          - m11)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m11,m14,m14,m11,mu2) * (  - 3./16.*(
          - 1 + nf)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(1,m44,m14,m14,m11,mu2) * ( 3./16.*( - 1
          + nf)*(1 + 2*nf)*pow(m11,2) );

      PQF12p6x11iR +=  + hh21d(2,m11,m14,m14,m11,mu2) * ( 3./16.*(m44
          - m11)*pow(m11,2) );

  return  PQF12p6x11iR/pow(mass.getf0(),4);
}

/////////////////////////////////////////////////////////////////////////////
// vev SPN
double qnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin){
  return  qnfPQSPNp4L(nf,mass,Liin)+qnfPQSPNp4R(nf,mass);
}

double qnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin){
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp4L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSPNp4L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  int nnf;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r,nnf);
  double   PQQ11p4x11iL =
       + 32*m44*L6r*nf + 4*m11*H2r + 8*m11*L8r;

  return PQQ11p4x11iL/pow(mass.getf0(),2);
}

double qnfPQSPNp4R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp4R probably called with wrong number"
	 << "of quark masses\n";}
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);
  double  PQQ11p4x11iR =
       + pi16 * (  - 1./2.*(m44 - m11)*nfm );

      PQQ11p4x11iR +=  + Ab(m11,mu2) * (  - 1./2. + 1./2.*pow(m11,-1)*
         m44*nfm - nfm );

      PQQ11p4x11iR +=  + Ab(m14,mu2) * ( nf );

  return PQQ11p4x11iR/pow(mass.getf0(),2);
}

double qnfPQSPNp6(const int nf, const quarkmassnf mass, const Linf Liin,
		const Ki Kiin){
  return qnfPQSPNp6K(nf,mass,Kiin)+qnfPQSPNp6L(nf,mass,Liin)
    +qnfPQSPNp6R(nf,mass);
}

double qnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin){
  /*
  vector<double> B0mq = mass.getB0mq();
  if (B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp6K probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Kiin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Ki in qnfPQSPNp6K\n";
  }
  double m11 = 2.*B0mq[0];
  double m44 = 2.*B0mq[1];
  double KK[116];
  Kiin.out(KK);
  double   PQQ11p6x11iK =
       + KK[25] * ( 48*pow(m11,2) );

      PQQ11p6x11iK +=  + KK[26] * ( 16*(m44 + 2*m11)*m44*nf );

      PQQ11p6x11iK +=  + KK[27] * ( 48*pow(m44,2)*pow(nf,2) );

  return PQQ11p6x11iK/pow(mass.getf0(),4);
  */
  return 0.;
}



double qnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp6L probably called with wrong number"
	 << "of quark masses\n";}
  if ( nf != Liin.getnf()){
    cout << "ERROR: conflict in direct nf and nf in Li in qnfPQSPNp6L\n";
  }

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r;
  Liin.out(L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r);
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double  PQQ11p6x11iL =
       + pi16 * (  - 16*(m44 - m11)*(m44 - m11)*L7r - 16*(3*m44 + m44*
         pow(nf,2) - 4*m11 - m11*nf + m11*pow(nf,2))*m44*L6r + 8*(3*m44
          + m44*pow(nf,2) - 4*m11 - m11*nf + m11*pow(nf,2))*m44*L4r + 4
         *(5*m44 - 6*m11)*m11*L5r*nfm - 8*(pow(m44,2) + 3*m11*m44 - 5*
         pow(m11,2))*L8r*nfm - 4*(pow(m44,2)*nf + 2*m11*m44*nf - 2*pow(
         m11,2) + pow(m11,2)*nf)*L8r + 2*(pow(m44,2)*nf + 2*m11*m44*nf
          - 2*pow(m11,2) + pow(m11,2)*nf)*L5r );

      PQQ11p6x11iL +=  + Ab(m11,mu2) * ( 8*pow(m11,-1)*pow(m44,2)*L8r*
         nfm + 16*pow(m11,-1)*pow(m44,2)*L7r + 32*pow(m11,-1)*pow(
         m44,2)*L6r - 16*pow(m11,-1)*pow(m44,2)*L4r - 24*m11*L8r + 12*
         m11*L5r - 32*(2 + nf)*m44*L6r + 16*(2 + nf)*m44*L4r + 16*(m44
          - 3*m11)*L8r*nfm - 12*(2*m44 - 3*m11)*L5r*nfm - 16*(4*m44 - 3
         *m11)*L7r );

      PQQ11p6x11iL +=  + Ab(m14,mu2) * ( 32*m44*L6r*pow(nf,2) - 16*m44*
         L4r*pow(nf,2) + 16*(m44 + m11)*pow(m14,-1)*m44*L6r*pow(nf,2)
          - 8*(m44 + m11)*pow(m14,-1)*m44*L4r*pow(nf,2) + 16*(m44 + m11
         )*L8r*nf - 8*(m44 + m11)*L5r*nf + 4*(m44 + m11)*(m44 + m11)*
         pow(m14,-1)*L8r*nf - 2*(m44 + m11)*(m44 + m11)*pow(m14,-1)*L5r
         *nf );

      PQQ11p6x11iL +=  + Ab(m44,mu2) * ( 16*( - 1 + nf)*(1 + 2*nf)*m44*
         L6r - 8*( - 1 + nf)*(1 + 2*nf)*m44*L4r );
 
  return PQQ11p6x11iL/pow(mass.getf0(),4);
}

double qnfPQSPNp6R(const int nf,const quarkmassnf mass){
  const double pi16 = 1./pow(4.*M_PI,2);
  vector<double> B0mq = mass.getB0mq();
  if ( B0mq.size() != 2){
    cout << "ERROR: qnfPQSPNp6R probably called with wrong number"
	 << "of quark masses\n";}

  double m11 = B0mq[0]*2.;
  double m44 = B0mq[1]*2.;
  double m14 = (m11+m44)/2.;
  double nfm = 1./double(nf);
  double mu2 = pow(mass.getmu(),2);

  double   PQQ11p6x11iR =
       + pow(pi16,2) * (  - 1./4.*(m44 - m11)*m11*nfm + 1./8.*(m44 - 
         m11)*(3*m44 - 7*m11)*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*pi16 * (  - 1./2.*pow(m11,-1)*pow(
         m44,2)*pow(nfm,2) + 1./4.*(m44 - 2*m11) + 5./2.*(m44 - m11)*
         pow(nfm,2) + 1./4.*(4*m44 - 7*m11)*nfm );

      PQQ11p6x11iR +=  + pow(Ab(m11,mu2),2) * ( 9./8. + 1./8.*pow(
         m11,-2)*pow(m44,2)*pow(nfm,2) - 1./2.*pow(m11,-1)*m44*nfm - 
         pow(m11,-1)*m44*pow(nfm,2) + 3./2.*nfm + 3./2.*pow(nfm,2) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m14,mu2) * (  - 1 - 1./2.*pow(
         m11,-1)*m44 - 1./4.*(m44 + m11)*pow(m14,-1) );

      PQQ11p6x11iR +=  + Ab(m11,mu2)*Ab(m44,mu2) * ( 1./2.*pow(m11,-1)*
         m44 - 1./4.*pow(m11,-1)*m44*nfm - 1./4.*pow(m11,-1)*m44*pow(
         nfm,2) );

      PQQ11p6x11iR +=  + Ab(m14,mu2)*pi16 * ( 1./4.*(2 + nf)*(m44 + m11
         ) );

      PQQ11p6x11iR +=  + pow(Ab(m14,mu2),2) * (  - 1./4.*nf - 1./4.*(
         m44 + m11)*pow(m14,-1)*nf );

      PQQ11p6x11iR +=  + Ab(m44,mu2)*pi16 * (  - 1./2.*m44 + 1./4.*m44*
         nfm + 1./4.*m44*pow(nfm,2) );

  return  PQQ11p6x11iR/pow(mass.getf0(),4);
}

