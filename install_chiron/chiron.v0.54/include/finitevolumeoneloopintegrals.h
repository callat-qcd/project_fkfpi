// finitevolumeoneloopintegrals.h is part of the CHIRON ChPT at two loops
// program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// uses the methods as derived in
//  J.~Bijnens, E.~Boström and T.~A.~Lähde,
//  %``Two-loop Sunset Integrals at Finite Volume,''
//  JHEP {\bf 1401} (2014) 019
//  [arXiv:1311.3531 [hep-lat]].
// and references therein

#ifndef FINITEVOLUMEONELOOPINTEGRALS_H
#define FINITEVOLUMEONELOOPINTEGRALS_H

// tadpoles
// theta functions
double AbVt(const double msq, const double xl);
double A22bVt(const double msq, const double xl);
double A23bVt(const double msq, const double xl);
double BbVt(const double msq, const double xl);
double B22bVt(const double msq, const double xl);
double B23bVt(const double msq, const double xl);
double CbVt(const double msq, const double xl);
double C22bVt(const double msq, const double xl);
double C23bVt(const double msq, const double xl);
double DbVt(const double msq, const double xl);
double D22bVt(const double msq, const double xl);
double D23bVt(const double msq, const double xl);

// Bessel functions
double AbVb(const double msq, const double xl);
double A22bVb(const double msq, const double xl);
double A23bVb(const double msq, const double xl);
double BbVb(const double msq, const double xl);
double B22bVb(const double msq, const double xl);
double B23bVb(const double msq, const double xl);
double CbVb(const double msq, const double xl);
double C22bVb(const double msq, const double xl);
double C23bVb(const double msq, const double xl);
double DbVb(const double msq, const double xl);
double D22bVb(const double msq, const double xl);
double D23bVb(const double msq, const double xl);

// tadpoles with different masses
// theta functions
double BbVt(const double m1sq,const double m2sq, const double xl);
//
//Bessel functions
double BbVb(const double m1sq,const double m2sq, const double xl);

////////////// accuracy setting //////////////////////////////////////
// theta functions
void setprecisionfinitevolumeoneloopt(const double Abacc=1e-10,
				      const double Bbacc=1e-9,
				      const bool printout=true);
// Bessel functions
void setprecisionfinitevolumeoneloopb(const int maxonesum=100,
				      const double Bbacc=1e-5,
				      const bool printout=true);

#endif // FINITEVOLUMEONELOOPINTEGRALS_H
