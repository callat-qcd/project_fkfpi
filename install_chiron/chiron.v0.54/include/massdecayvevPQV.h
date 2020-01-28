// massdecayvevPQV.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results for the finite volume correction
// to the partially quenched masses and decay constants
// derived in
// J.~Bijnens and T.~Rössler,
//  %``Finite Volume for Three-Flavour Partially Quenched Chiral Perturbation Theory through NNLO in the Meson Sector,''
//  arXiv:1508.07238 [hep-lat].
//  %%CITATION = ARXIV:1508.07238;%%

//
// contains the finite volume corrections from the
// partially quenched  expressions for masses and decay
// constants for three sea-quark flavours 
// in terms of the lowest order masses and the chiral limit
// decay constant

#ifndef MASSDECAYVEVPQV_H
#define  MASSDECAYVEVPQV_H


#include "inputsnf.h"
#include "Linf.h"

// with theta functions

// masses
// 1 valence mass 1 sea mass *****************************************

double mv1s1nf3p4Vt(const quarkmassnf mass, const double L);
double mv1s1nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s1nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s1nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 1 sea mass *****************************************

double mv2s1nf3p4Vt(const quarkmassnf mass, const double L);
double mv2s1nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s1nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s1nf3p6RVt(const quarkmassnf mass, const double L);

// 1 valence mass 2 sea mass *****************************************

double mv1s2nf3p4Vt(const quarkmassnf mass, const double L);
double mv1s2nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s2nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s2nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 2 sea mass *****************************************

double mv2s2nf3p4Vt(const quarkmassnf mass, const double L);
double mv2s2nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s2nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s2nf3p6RVt(const quarkmassnf mass, const double L);

// 1 valence mass 3 sea mass *****************************************

double mv1s3nf3p4Vt(const quarkmassnf mass, const double L);
double mv1s3nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s3nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s3nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 3 sea mass *****************************************

double mv2s3nf3p4Vt(const quarkmassnf mass, const double L);
double mv2s3nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s3nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s3nf3p6RVt(const quarkmassnf mass, const double L);

//********************************************************************
// decay constants
// 1 valence mass 1 sea mass *****************************************

double fv1s1nf3p4Vt(const quarkmassnf mass, const double L);
double fv1s1nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s1nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s1nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 1 sea mass *****************************************

double fv2s1nf3p4Vt(const quarkmassnf mass, const double L);
double fv2s1nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s1nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s1nf3p6RVt(const quarkmassnf mass, const double L);

// 1 valence mass 2 sea mass *****************************************

double fv1s2nf3p4Vt(const quarkmassnf mass, const double L);
double fv1s2nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s2nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s2nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 2 sea mass *****************************************

double fv2s2nf3p4Vt(const quarkmassnf mass, const double L);
double fv2s2nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s2nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s2nf3p6RVt(const quarkmassnf mass, const double L);

// 1 valence mass 3 sea mass *****************************************

double fv1s3nf3p4Vt(const quarkmassnf mass, const double L);
double fv1s3nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s3nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s3nf3p6RVt(const quarkmassnf mass, const double L);

// 2 valence mass 3 sea mass *****************************************

double fv2s3nf3p4Vt(const quarkmassnf mass, const double L);
double fv2s3nf3p6Vt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s3nf3p6LVt(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s3nf3p6RVt(const quarkmassnf mass, const double L);


// with Bessel functions functions

// masses
// 1 valence mass 1 sea mass *****************************************

double mv1s1nf3p4Vb(const quarkmassnf mass, const double L);
double mv1s1nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s1nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s1nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 1 sea mass *****************************************

double mv2s1nf3p4Vb(const quarkmassnf mass, const double L);
double mv2s1nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s1nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s1nf3p6RVb(const quarkmassnf mass, const double L);

// 1 valence mass 2 sea mass *****************************************

double mv1s2nf3p4Vb(const quarkmassnf mass, const double L);
double mv1s2nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s2nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s2nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 2 sea mass *****************************************

double mv2s2nf3p4Vb(const quarkmassnf mass, const double L);
double mv2s2nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s2nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s2nf3p6RVb(const quarkmassnf mass, const double L);

// 1 valence mass 3 sea mass *****************************************

double mv1s3nf3p4Vb(const quarkmassnf mass, const double L);
double mv1s3nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s3nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv1s3nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 3 sea mass *****************************************

double mv2s3nf3p4Vb(const quarkmassnf mass, const double L);
double mv2s3nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s3nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double mv2s3nf3p6RVb(const quarkmassnf mass, const double L);

//********************************************************************
// decay constants
// 1 valence mass 1 sea mass *****************************************

double fv1s1nf3p4Vb(const quarkmassnf mass, const double L);
double fv1s1nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s1nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s1nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 1 sea mass *****************************************

double fv2s1nf3p4Vb(const quarkmassnf mass, const double L);
double fv2s1nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s1nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s1nf3p6RVb(const quarkmassnf mass, const double L);

// 1 valence mass 2 sea mass *****************************************

double fv1s2nf3p4Vb(const quarkmassnf mass, const double L);
double fv1s2nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s2nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s2nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 2 sea mass *****************************************

double fv2s2nf3p4Vb(const quarkmassnf mass, const double L);
double fv2s2nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s2nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s2nf3p6RVb(const quarkmassnf mass, const double L);

// 1 valence mass 3 sea mass *****************************************

double fv1s3nf3p4Vb(const quarkmassnf mass, const double L);
double fv1s3nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s3nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv1s3nf3p6RVb(const quarkmassnf mass, const double L);

// 2 valence mass 3 sea mass *****************************************

double fv2s3nf3p4Vb(const quarkmassnf mass, const double L);
double fv2s3nf3p6Vb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s3nf3p6LVb(const quarkmassnf mass, const Linf Liin, const double L);
double fv2s3nf3p6RVb(const quarkmassnf mass, const double L);

#endif  //  MASSDECAYVEVPQV_H
