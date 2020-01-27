// massdecayvevnfV.h is part of the CHIRON ChPT at two loops program collection
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

#ifndef MASSDECAYVEVNFV_H
#define  MASSDECAYVEVNFV_H

#include "inputsnf.h"
#include "Linf.h"

// theta functions
// masses
double mnfSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSUNp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfSUNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSUNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSUNp6RVt(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSONp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfSONp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSONp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSONp6RVt(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double mnfSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSPNp6RVt(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double fnfSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSPNp6RVt(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSPNp4Vt(const int nf, const quarkmassnf mass, const double L);
double qnfSPNp6Vt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSPNp6LVt(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSPNp6RVt(const int nf,const quarkmassnf mass, const double L);


// Bessel functions
// masses
double mnfSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSUNp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfSUNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSUNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSUNp6RVb(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSONp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfSONp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSONp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSONp6RVb(const int nf,const quarkmassnf mass, const double L);
// masses
double mnfSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double mnfSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double mnfSPNp6RVb(const int nf,const quarkmassnf mass, const double L);
// decay
double fnfSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double fnfSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double fnfSPNp6RVb(const int nf,const quarkmassnf mass, const double L);
// vev
double qnfSPNp4Vb(const int nf, const quarkmassnf mass, const double L);
double qnfSPNp6Vb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSPNp6LVb(const int nf,const quarkmassnf mass, const Linf Liin,
		  const double L);
double qnfSPNp6RVb(const int nf,const quarkmassnf mass, const double L);

#endif  //  MASSDECAYVEVNFV_H
