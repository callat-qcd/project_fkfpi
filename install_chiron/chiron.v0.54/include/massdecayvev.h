// massdecayvev.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results derived in
//  G.~Amoros, J.~Bijnens and P.~Talavera,
//  Two point functions at two loops in three flavor chiral perturbation theory,
//  Nucl.\ Phys.\ B {\bf 568} (2000) 319 [hep-ph/9907264]
//
// G.~Amoros, J.~Bijnens and P.~Talavera,
// ``K(lepton 4) form-factors and pi pi scattering,''
// Nucl.\ Phys.\ B {\bf 585} (2000) 293
// [Erratum-ibid.\ B {\bf 598} (2001) 665] [hep-ph/0003258].


#ifndef MASSDECAYVEV_H
#define MASSDECAYVEV_H

#include <complex>
typedef std::complex<double> dcomplex;

#include "inputs.h"
#include "Li.h"
#include "Ci.h"

double mpi4(const physmass mass, const Li liin);
double mpi4L(const physmass mass, const Li liin);
double mpi4R(const physmass mass);
double mpi6(const physmass mass, const Li liin, const Ci Ciin);
double mpi6LC(const physmass mass, const Li liin, const Ci Ciin);
double mpi6L(const physmass mass, const Li Liin);
double mpi6C(const physmass mass, const Ci Ciin);
double mpi6R(const physmass mass);

double mk4(const physmass mass, const Li liin);
double mk4L(const physmass mass, const Li liin);
double mk4R(const physmass mass);
double mk6(const physmass mass, const Li liin, const Ci Ciin);
double mk6LC(const physmass mass, const Li liin, const Ci Ciin);
double mk6L(const physmass mass, const Li Liin);
double mk6C(const physmass mass, const Ci Ciin);
double mk6R(const physmass mass);

double meta4(const physmass mass, const Li liin);
double meta4L(const physmass mass, const Li liin);
double meta4R(const physmass mass);
double meta6(const physmass mass, const Li liin, const Ci Ciin);
double meta6LC(const physmass mass, const Li liin, const Ci Ciin);
double meta6L(const physmass mass, const Li Liin);
double meta6C(const physmass mass, const Ci Ciin);
double meta6R(const physmass mass);

double fpi4(const physmass mass, const Li liin);
double fpi4L(const physmass mass, const Li liin);
double fpi4R(const physmass mass);
double fpi6(const physmass mass, const Li liin, const Ci Ciin);
double fpi6LC(const physmass mass, const Li liin, const Ci Ciin);
double fpi6L(const physmass mass, const Li Liin);
double fpi6C(const physmass mass, const Ci Ciin);
double fpi6R(const physmass mass);

double fk4(const physmass mass, const Li liin);
double fk4L(const physmass mass, const Li liin);
double fk4R(const physmass mass);
double fk6(const physmass mass, const Li liin, const Ci Ciin);
double fk6LC(const physmass mass, const Li liin, const Ci Ciin);
double fk6L(const physmass mass, const Li Liin);
double fk6C(const physmass mass, const Ci Ciin);
double fk6R(const physmass mass);

double feta4(const physmass mass, const Li liin);
double feta4L(const physmass mass, const Li liin);
double feta4R(const physmass mass);
double feta6(const physmass mass, const Li liin, const Ci Ciin);
double feta6LC(const physmass mass, const Li liin, const Ci Ciin);
double feta6L(const physmass mass, const Li Liin);
double feta6C(const physmass mass, const Ci Ciin);
double feta6R(const physmass mass);

double qqup4(const physmass mass, const Li liin);
double qqup4L(const physmass mass, const Li liin);
double qqup4R(const physmass mass);
double qqup6(const physmass mass, const Li liin, const Ci Ciin);
double qqup6LC(const physmass mass, const Li liin, const Ci Ciin);
double qqup6L(const physmass mass, const Li Liin);
double qqup6C(const physmass mass, const Ci Ciin);
double qqup6R(const physmass mass);

double qqstrange4(const physmass mass, const Li liin);
double qqstrange4L(const physmass mass, const Li liin);
double qqstrange4R(const physmass mass);
double qqstrange6(const physmass mass, const Li liin, const Ci Ciin);
double qqstrange6LC(const physmass mass, const Li liin, const Ci Ciin);
double qqstrange6L(const physmass mass, const Li Liin);
double qqstrange6C(const physmass mass, const Ci Ciin);
double qqstrange6R(const physmass mass);

#endif // MASSDECAYVEV_H

