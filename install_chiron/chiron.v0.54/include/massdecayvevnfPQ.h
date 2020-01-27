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

#ifndef MASSDECAYVEVNFPQ_H
#define  MASSDECAYVEVNFPQ_H

#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"

//////////// SUN //////////////////////////////////////////////////////
// masses

double mnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSUNp4R(const int nf,const quarkmassnf mass);
double mnfPQSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfPQSUNp6R(const int nf,const quarkmassnf mass);

// decay constant

double fnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSUNp4R(const int nf,const quarkmassnf mass);
double fnfPQSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfPQSUNp6R(const int nf,const quarkmassnf mass);

// vev

double qnfPQSUNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfPQSUNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSUNp4R(const int nf,const quarkmassnf mass);
double qnfPQSUNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfPQSUNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSUNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfPQSUNp6R(const int nf,const quarkmassnf mass);

//////////// SON //////////////////////////////////////////////////////
// masses

double mnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSONp4R(const int nf,const quarkmassnf mass);
double mnfPQSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfPQSONp6R(const int nf,const quarkmassnf mass);

// decay constant

double fnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSONp4R(const int nf,const quarkmassnf mass);
double fnfPQSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfPQSONp6R(const int nf,const quarkmassnf mass);

// vev

double qnfPQSONp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfPQSONp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSONp4R(const int nf,const quarkmassnf mass);
double qnfPQSONp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfPQSONp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSONp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfPQSONp6R(const int nf,const quarkmassnf mass);

//////////// SPN //////////////////////////////////////////////////////
// masses

double mnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double mnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSPNp4R(const int nf,const quarkmassnf mass);
double mnfPQSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double mnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double mnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double mnfPQSPNp6R(const int nf,const quarkmassnf mass);

// decay constant

double fnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double fnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSPNp4R(const int nf,const quarkmassnf mass);
double fnfPQSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double fnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double fnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double fnfPQSPNp6R(const int nf,const quarkmassnf mass);

// vev

double qnfPQSPNp4(const int nf, const quarkmassnf mass, const Linf Liin);
double qnfPQSPNp4L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSPNp4R(const int nf,const quarkmassnf mass);
double qnfPQSPNp6(const int nf,const quarkmassnf mass, const Linf Liin,
		const Ki Kiin);
double qnfPQSPNp6L(const int nf,const quarkmassnf mass, const Linf Liin);
double qnfPQSPNp6K(const int nf,const quarkmassnf mass, const Ki Kiin);
double qnfPQSPNp6R(const int nf,const quarkmassnf mass);

#endif  //  MASSDECAYVEVNFPQ_H
