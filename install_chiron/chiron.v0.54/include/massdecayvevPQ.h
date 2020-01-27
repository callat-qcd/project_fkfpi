// massdecayvevPQ.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results for partially quenched masses and decay constants
// derived in
//J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  ``Three-flavor partially quenched chiral perturbation theory at NNLO
// for meson masses and decay constants,''
//  Phys.\ Rev.\ D {\bf 73} (2006) 074509
//  [hep-lat/0602003].
// J.~Bijnens and T.~A.~Lahde,
//  ``Decay constants of pseudoscalar mesons to two loops in three-flavor
// partially quenched (chi)PT,''
//  Phys.\ Rev.\ D {\bf 71} (2005) 094502
// [hep-lat/0501014].
//  J.~Bijnens, N.~Danielsson and T.~A.~Lahde,
//  ``The Pseudoscalar meson mass to two loops in three-flavor partially
//  quenched chiPT,''
//  Phys.\ Rev.\ D {\bf 70} (2004) 111503
//  [hep-lat/0406017].
//
// contains the partially quenched  expressions for masses and decay
// constants for three sea-quark flavours 
// in terms of the lowest order masses and the chiral limit
// decay constant

#include "inputsnf.h"
#include "Linf.h"
#include "Ki.h"

// masses
// 1 valence mass 1 sea mass *****************************************

double mv1s1nf3p4(const quarkmassnf mass, const Linf Liin);
double mv1s1nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv1s1nf3p4R(const quarkmassnf mass);
double mv1s1nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv1s1nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv1s1nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv1s1nf3p6R(const quarkmassnf mass);

// 2 valence mass 1 sea mass *****************************************

double mv2s1nf3p4(const quarkmassnf mass, const Linf Liin);
double mv2s1nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv2s1nf3p4R(const quarkmassnf mass);
double mv2s1nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv2s1nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv2s1nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv2s1nf3p6R(const quarkmassnf mass);

// 1 valence mass 2 sea mass *****************************************

double mv1s2nf3p4(const quarkmassnf mass, const Linf Liin);
double mv1s2nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv1s2nf3p4R(const quarkmassnf mass);
double mv1s2nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv1s2nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv1s2nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv1s2nf3p6R(const quarkmassnf mass);

// 2 valence mass 2 sea mass *****************************************

double mv2s2nf3p4(const quarkmassnf mass, const Linf Liin);
double mv2s2nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv2s2nf3p4R(const quarkmassnf mass);
double mv2s2nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv2s2nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv2s2nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv2s2nf3p6R(const quarkmassnf mass);

// 1 valence mass 3 sea mass *****************************************

double mv1s3nf3p4(const quarkmassnf mass, const Linf Liin);
double mv1s3nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv1s3nf3p4R(const quarkmassnf mass);
double mv1s3nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv1s3nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv1s3nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv1s3nf3p6R(const quarkmassnf mass);

// 2 valence mass 3 sea mass *****************************************

double mv2s3nf3p4(const quarkmassnf mass, const Linf Liin);
double mv2s3nf3p4L(const quarkmassnf mass, const Linf Liin);
double mv2s3nf3p4R(const quarkmassnf mass);
double mv2s3nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double mv2s3nf3p6L(const quarkmassnf mass, const Linf Liin);
double mv2s3nf3p6K(const quarkmassnf mass, const Ki Kiin);
double mv2s3nf3p6R(const quarkmassnf mass);

//********************************************************************
// decay constants
// 1 valence mass 1 sea mass *****************************************

double fv1s1nf3p4(const quarkmassnf mass, const Linf Liin);
double fv1s1nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv1s1nf3p4R(const quarkmassnf mass);
double fv1s1nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv1s1nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv1s1nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv1s1nf3p6R(const quarkmassnf mass);

// 2 valence mass 1 sea mass *****************************************

double fv2s1nf3p4(const quarkmassnf mass, const Linf Liin);
double fv2s1nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv2s1nf3p4R(const quarkmassnf mass);
double fv2s1nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv2s1nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv2s1nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv2s1nf3p6R(const quarkmassnf mass);

// 1 valence mass 2 sea mass *****************************************

double fv1s2nf3p4(const quarkmassnf mass, const Linf Liin);
double fv1s2nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv1s2nf3p4R(const quarkmassnf mass);
double fv1s2nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv1s2nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv1s2nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv1s2nf3p6R(const quarkmassnf mass);

// 2 valence mass 2 sea mass *****************************************

double fv2s2nf3p4(const quarkmassnf mass, const Linf Liin);
double fv2s2nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv2s2nf3p4R(const quarkmassnf mass);
double fv2s2nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv2s2nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv2s2nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv2s2nf3p6R(const quarkmassnf mass);

// 1 valence mass 3 sea mass *****************************************

double fv1s3nf3p4(const quarkmassnf mass, const Linf Liin);
double fv1s3nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv1s3nf3p4R(const quarkmassnf mass);
double fv1s3nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv1s3nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv1s3nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv1s3nf3p6R(const quarkmassnf mass);

// 2 valence mass 3 sea mass *****************************************

double fv2s3nf3p4(const quarkmassnf mass, const Linf Liin);
double fv2s3nf3p4L(const quarkmassnf mass, const Linf Liin);
double fv2s3nf3p4R(const quarkmassnf mass);
double fv2s3nf3p6(const quarkmassnf mass, const Linf Liin, const Ki Kiin);
double fv2s3nf3p6L(const quarkmassnf mass, const Linf Liin);
double fv2s3nf3p6K(const quarkmassnf mass, const Ki Kiin);
double fv2s3nf3p6R(const quarkmassnf mass);
