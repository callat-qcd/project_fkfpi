// sunsetintegrals.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// contains all the sunset functions
// Method is described in detail in 
// G. Amoros, J. Bijnens and P. Talavera,Nucl. Phys. B568 (2000) 319-363
// [hep-ph/9907264]

#ifndef SUNSETINTEGRALS_H
#define SUNSETINTEGRALS_H

#include<complex>

void setprecisionsunsetintegrals(const double eps);
double getprecisionsunsetintegrals(void);

typedef std::complex<double> dcomplex;

double hh(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
double hh1(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
double hh21(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
double hh31(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);

dcomplex zhh(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
dcomplex zhh1(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
dcomplex zhh21(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
dcomplex zhh31(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);

double hhd(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
double hh1d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
double hh21d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);

dcomplex zhhd(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
dcomplex zhh1d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);
dcomplex zhh21d(const double m1sq, const double m2sq,const double m3sq,
	  const double qsq, const double mu2);

#endif //  SUNSETINTEGRALS_H
