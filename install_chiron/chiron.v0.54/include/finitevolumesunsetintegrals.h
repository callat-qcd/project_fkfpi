// finitevolumeoneloopintegrals.h is part of the CHIRON ChPT at two loops
// program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.01
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// uses the methods derived in
//  J.~Bijnens, E.~Boström and T.~A.~Lähde,
//  Two-loop Sunset Integrals at Finite Volume,
//  JHEP {\bf 1401} (2014) 019
//  [arXiv:1311.3531 [hep-lat]].
 
#ifndef FINITEVOLUMESUNSETINTEGRALS_H
#define FINITEVOLUMESUNSETINTEGRALS_H

double hhVt(const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1Vt(const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21Vt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22Vt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27Vt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

double hhdVt(const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1dVt(const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21dVt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22dVt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27dVt(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

double hhVb(const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1Vb(const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21Vb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22Vb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27Vb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

double hhdVb(const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1dVb(const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21dVb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22dVb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27dVb(const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);


double hhVt(const int nprop,
	    const double m1sq,const double m2sq,const double m3sq,
	    const double qsq,const double xl,const double mu2);
double hh1Vt(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,
const double xl,const double mu2);
double hh21Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh22Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh27Vt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);

double hhdVt(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1dVt(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27dVt(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

double hhVb(const int nprop,
	    const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1Vb(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27Vb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

double hhdVb(const int nprop,
	     const double m1sq,const double m2sq,const double m3sq,
	     const double qsq,const double xl,const double mu2);
double hh1dVb(const int nprop,
	      const double m1sq,const double m2sq,const double m3sq,
	      const double qsq,const double xl,const double mu2);
double hh21dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh22dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);
double hh27dVb(const int nprop,
	       const double m1sq,const double m2sq,const double m3sq,
	       const double qsq,const double xl,const double mu2);

////////////// accuracy setting //////////////////////////////////////
void setprecisionfinitevolumesunsett(const double sunsetracc=1e-5,
				     const double sunsetrsacc=1e-4,
				     const bool out=true);

void setprecisionfinitevolumesunsetb(const int maxonesum=100,
				     const int maxtwosum=40,
				     const double sunsetracc=1e-5,
				     const double sunsetrsacc=1e-4,
				     const bool out=true);

#endif // FINITEVOLUMESUNSETINTEGRALS_H
