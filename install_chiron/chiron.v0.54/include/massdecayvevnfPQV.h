// massdecayvevnfPQV.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results for QCD like theories for masses, decay constants and
// the qbarq vacuume expectation value
// derived in
// J.~Bijnens and T.~Rössler

// in terms of the lowest order mass and the chiral limit
// decay constant

#ifndef MASSDECAYVEVNFPQV_H
#define  MASSDECAYVEVNFPQV_H

#include "inputsnf.h"
#include "Linf.h"

// theta functions
// masses
double mnfPQSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfPQSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfPQSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfPQSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfPQSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfPQSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfPQSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfPQSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfPQSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfPQSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfPQSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfPQSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSPNp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfPQSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSPNp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfPQSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfPQSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSPNp6RVt(const int nf,const quarkmassnf mass, const double L);


// Bessel functions
// masses
double mnfPQSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfPQSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfPQSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double qnfPQSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfPQSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfPQSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfPQSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfPQSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double qnfPQSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfPQSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfPQSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfPQSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfPQSPNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfPQSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfPQSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfPQSPNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double qnfPQSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfPQSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfPQSPNp6RVb(const int nf,const quarkmassnf mass, const double L);

#endif  //  MASSDECAYVEVNFPQV_H
