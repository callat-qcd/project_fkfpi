// getfpimeta.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#ifndef GETFPIMETA_H
#define GETFPIMETA_H

physmass getfpimeta6(const double mpiin, const double mkin,
		     const physmass massin, const Li liin, const Ci ciin);
physmass getfpimeta4(const double mpiin, const double mkin,
		     const physmass massin, const Li liin);

#endif  // GETFPIMETA_H
