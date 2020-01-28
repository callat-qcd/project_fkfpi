// massdecayvevlo.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// results derived in
//  G.~Amoros, J.~Bijnens and P.~Talavera,
//  Two point functions at two loops in three flavor chiral perturbation theory,
//  Nucl.\ Phys.\ B {\bf 568} (2000) 319 [hep-ph/9907264]


#ifndef MASSDECAYVEVLO_H
#define MASSDECAYVEVLO_H

#include "inputs.h"
#include "Li.h"
#include "Ci.h"

double mpi4lo(const lomass mass, const Li liin);
double mpi4Llo(const lomass mass, const Li liin);
double mpi4Rlo(const lomass mass);
double mpi6lo(const lomass mass, const Li liin, const Ci Ciin);
double mpi6LClo(const lomass mass, const Li liin, const Ci Ciin);
double mpi6Llo(const lomass mass, const Li Liin);
double mpi6Clo(const lomass mass, const Ci Ciin);
double mpi6Rlo(const lomass mass);

double mk4lo(const lomass mass, const Li liin);
double mk4Llo(const lomass mass, const Li liin);
double mk4Rlo(const lomass mass);
double mk6lo(const lomass mass, const Li liin, const Ci Ciin);
double mk6LClo(const lomass mass, const Li liin, const Ci Ciin);
double mk6Llo(const lomass mass, const Li Liin);
double mk6Clo(const lomass mass, const Ci Ciin);
double mk6Rlo(const lomass mass);

double meta4lo(const lomass mass, const Li liin);
double meta4Llo(const lomass mass, const Li liin);
double meta4Rlo(const lomass mass);
double meta6lo(const lomass mass, const Li liin, const Ci Ciin);
double meta6LClo(const lomass mass, const Li liin, const Ci Ciin);
double meta6Llo(const lomass mass, const Li Liin);
double meta6Clo(const lomass mass, const Ci Ciin);
double meta6Rlo(const lomass mass);

double fpi4lo(const lomass mass, const Li liin);
double fpi4Llo(const lomass mass, const Li liin);
double fpi4Rlo(const lomass mass);
double fpi6lo(const lomass mass, const Li liin, const Ci Ciin);
double fpi6LClo(const lomass mass, const Li liin, const Ci Ciin);
double fpi6Llo(const lomass mass, const Li Liin);
double fpi6Clo(const lomass mass, const Ci Ciin);
double fpi6Rlo(const lomass mass);

double fk4lo(const lomass mass, const Li liin);
double fk4Llo(const lomass mass, const Li liin);
double fk4Rlo(const lomass mass);
double fk6lo(const lomass mass, const Li liin, const Ci Ciin);
double fk6LClo(const lomass mass, const Li liin, const Ci Ciin);
double fk6Llo(const lomass mass, const Li Liin);
double fk6Clo(const lomass mass, const Ci Ciin);
double fk6Rlo(const lomass mass);

double feta4lo(const lomass mass, const Li liin);
double feta4Llo(const lomass mass, const Li liin);
double feta4Rlo(const lomass mass);
double feta6lo(const lomass mass, const Li liin, const Ci Ciin);
double feta6LClo(const lomass mass, const Li liin, const Ci Ciin);
double feta6Llo(const lomass mass, const Li Liin);
double feta6Clo(const lomass mass, const Ci Ciin);
double feta6Rlo(const lomass mass);

double qqup4lo(const lomass mass, const Li liin);
double qqup4Llo(const lomass mass, const Li liin);
double qqup4Rlo(const lomass mass);
double qqup6lo(const lomass mass, const Li liin, const Ci Ciin);
double qqup6LClo(const lomass mass, const Li liin, const Ci Ciin);
double qqup6Llo(const lomass mass, const Li Liin);
double qqup6Clo(const lomass mass, const Ci Ciin);
double qqup6Rlo(const lomass mass);

double qqstrange4lo(const lomass mass, const Li liin);
double qqstrange4Llo(const lomass mass, const Li liin);
double qqstrange4Rlo(const lomass mass);
double qqstrange6lo(const lomass mass, const Li liin, const Ci Ciin);
double qqstrange6LClo(const lomass mass, const Li liin, const Ci Ciin);
double qqstrange6Llo(const lomass mass, const Li Liin);
double qqstrange6Clo(const lomass mass, const Ci Ciin);
double qqstrange6Rlo(const lomass mass);

#endif // MASSDECAYVEVLO_H

