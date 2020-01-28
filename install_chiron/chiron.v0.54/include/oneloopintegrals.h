// oneloopintegrals.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of the one-loop integrals


#ifndef ONELOOPINTEGRALS_H
#define ONELOOPINTEGRALS_H

#include <complex>
typedef std::complex<double> dcomplex;

// tadpoles
double Ab(const double msq,const double xmu2);
double Bb(const double msq,const double xmu2);
double Cb(const double msq,const double xmu2);
double Ab(const int n, const double msq, const double mu2);
double Abeps(const double msq,const double xmu2);
double Bbeps(const double msq,const double xmu2);
double Cbeps(const double msq,const double xmu2);


// precision for those done with numerical integration over x
void setprecisiononeloopintegrals(const double eps=1e-10);
double getprecisiononeloopintegrals(void);

// bare bubble integral
dcomplex Bb(const double m1sq,const double m2sq,const double qsq,
	    const double mu2);
dcomplex Bb(const double msq,const double qsq,const double mu2);
dcomplex Bbnum(const double m1sq,const double m2sq,const double qsq,
	    const double mu2);


dcomplex B1b(const double m1sq,const double m2sq,const double qsq,
	    const double mu2);
// implemented for equal mass, remove analytically
dcomplex B1b(const double msq,const double qsq,const double mu2);
dcomplex B1bnum(const double m1sq,const double m2sq,const double qsq,
	    const double mu2);

dcomplex B21b(const double m1sq, const double m2sq, const double qsq,
	     const double xmu2);
dcomplex B21bnum(const double m1sq, const double m2sq, const double qsq,
	     const double xmu2);


dcomplex B22b(const double m1sq, const double m2sq, const double qsq,
	     const double xmu2);
dcomplex B22b(const double msq, const double qsq, const double xmu2);
dcomplex B22bnum(const double sq1,const double sq2,const double qsq,
	      const double xmu2);

dcomplex B31bnum(const double sq1,const double sq2,const double qsq,
	    const double xmu2);
dcomplex B32bnum(const double sq1,const double sq2,const double qsq,
	    const double xmu2);
//
//dcomplex c0int(const double sq1,const double sq2,const double sq3, 
//	       const double qsq,const double xmu2);
//dcomplex c0intana(const double sq1,const double sq2,const double sq3, 
//	       const double qsq,const double xmu2);
//
//
//double aaeps(const double msq,const double xmu2);
//double bbeps(const double msq,const double xmu2);
//double cceps(const double msq,const double xmu2);
//
//dcomplex b0eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//dcomplex b0epsana(const double sq1,const double sq2,
//	       const double qsq,const double xmu2);
//
//dcomplex b1eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//dcomplex b21eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//dcomplex b22eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//dcomplex b31eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//dcomplex b32eps(const double sq1,const double sq2,const double qsq,
//	    const double xmu2);
//
//dcomplex c0eps(const double sq1,const double sq2,const double sq3, 
//	       const double qsq,const double xmu2);
//
//// incomplete ones:
//dcomplex c0epsana(const double sq1,const double sq2,const double sq3, 
//	       const double qsq,const double xmu2);

#endif //  ONELOOPINTEGRALS_H

