// quenchedsunsetintegrals.h is part of the CHIRON ChPT
// at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains all the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]
// has the extra ones needed for the partially quenched case with the
// extra first integer index indicating the power of the propagators added
// Derived in:
//
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 73} (2006) 074509
//  [hep-lat/0602003].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 72} (2005) 074502
//  [hep-lat/0506004].
//  J.~Bijnens and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 71} (2005) 094502
//  [hep-lat/0501014].
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  Phys.\ Rev.\ D {\bf 70} (2004) 111503
//  [hep-lat/0406017].

#ifndef QUENCHEDSUNSETINTEGRALS_H
#define QUENCHEDSUNSETINTEGRALS_H


void setprecisionquenchedsunsetintegrals(const double eps);
double getprecisionquenchedsunsetintegrals(void);


double hh(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);
double hh1(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);
double hh21(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);
double hhd(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);
double hh1d(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);
double hh21d(const int iprop, const double m1sq, const double m2sq,
	  const double m3sq,const double qsq, const double mu2);



#endif  // QUENCHEDSUNSETINTEGRALS_H
