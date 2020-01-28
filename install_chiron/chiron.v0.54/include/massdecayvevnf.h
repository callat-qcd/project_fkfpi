// massdecayvevnf.h is part of the CHIRON ChPT at two loops program collection
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

#ifndef MASSDECAYVEVNF_H
#define  MASSDECAYVEVNF_H

#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"

// SUN case
// masses
double mnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSUNp4R(const int nf,const quarkmassnf mass);
double mnfSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfSUNp6R(const int nf,const quarkmassnf mass);
// decay constant
double fnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSUNp4R(const int nf,const quarkmassnf mass);
double fnfSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfSUNp6R(const int nf,const quarkmassnf mass);
// qbar
double qnfSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSUNp4R(const int nf,const quarkmassnf mass);
double qnfSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfSUNp6R(const int nf,const quarkmassnf mass);

// SON case
// masses
double mnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSONp4R(const int nf,const quarkmassnf mass);
double mnfSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfSONp6R(const int nf,const quarkmassnf mass);
// decay constant
double fnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSONp4R(const int nf,const quarkmassnf mass);
double fnfSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfSONp6R(const int nf,const quarkmassnf mass);
// qbar
double qnfSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSONp4R(const int nf,const quarkmassnf mass);
double qnfSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfSONp6R(const int nf,const quarkmassnf mass);

// SPN case
// masses
double mnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSPNp4R(const int nf,const quarkmassnf mass);
double mnfSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfSPNp6R(const int nf,const quarkmassnf mass);
// decay constant
double fnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSPNp4R(const int nf,const quarkmassnf mass);
double fnfSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfSPNp6R(const int nf,const quarkmassnf mass);
// qbar
double qnfSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSPNp4R(const int nf,const quarkmassnf mass);
double qnfSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfSPNp6R(const int nf,const quarkmassnf mass);

#endif  //  MASSDECAYVEVNF_H
