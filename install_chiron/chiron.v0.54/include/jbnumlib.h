// jbnumlib.h is part of the numerical library jbdnumlib
// included with the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#ifndef JBNUMLIB_H
#define JBNUMLIB_H
#include <cmath>
#include <complex>

typedef std::complex<double> dcomplex; 

//**************************************************************************
//  integration routines
// eps is error but can be relative or absolute 
// real
double jbdgauss(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdgauss2(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdquad15(double (*f)(const double),const double a,const double b,
		const double eps);
double jbdquad21(double (*f)(const double),const double a,const double b,
		const double eps);

// real with singularity
double jbdcauch(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdcauch2(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdsing15(double (*f)(const double),const double a,const double b,
		const double s,const double eps);
double jbdsing21(double (*f)(const double),const double a,const double b,
		const double s,const double eps);

// complex
dcomplex jbwgauss(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwgauss2(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwquad15(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);
dcomplex jbwquad21(dcomplex (*f)(const dcomplex),const dcomplex a,
		  const dcomplex b,const double eps);

// complex with singularity
//dcomplex jbwcauch2(dcomplex (*f)(const dcomplex),const dcomplex a,
//		 const dcomplex b,const dcomplex s,const double eps);

// dimension 2 and 3 real integration
double jbdad2(double (*fcn)(double x[]),double a[],double b[],
	      const double releps, double &relerr,int &ifail);
double jbdad3(double (*fcn)(double x[]),double a[],double b[],
	      const double releps, double &relerr,int &ifail);

//**************************************************************************
// special functions
dcomplex jbdli2(const dcomplex x);
double jbdbesi0(const double x); 
double jbdbesi1(const double x); 
double jbdbesk0(const double x); 
double jbdbesk1(const double x); 
// implemented via recursion:
double jbdbesk2(const double x); 
double jbdbesk3(const double x); 
double jbdbesk4(const double x); 
double jbdtheta3(const double u,const double q);
double jbderivutheta3(const double u,const double q);
double jbderiv2utheta3(const double u,const double q);
double jbdtheta3(const double u,const double q);
double jbdtheta30(const double q);
double jbdtheta30m1(const double q);
double jbdtheta32(const double q);
double jbdtheta34(const double q);
double jbdtheta2d0(const double alpha, const double beta, const double gamma);
double jbdtheta2d0m1(const double alpha, const double beta, const double gamma);
double jbdtheta2d02(const double alpha, const double beta, const double gamma);

#endif // JBNUMLIB_H
