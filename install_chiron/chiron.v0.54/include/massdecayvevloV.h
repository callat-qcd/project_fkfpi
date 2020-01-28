// massdecayvevloV.h is part of the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.1
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

#ifndef MASSDECAYVEVLOV_H
#define  MASSDECAYVEVLOV_H


#include "inputs.h"

// theta functions
double mpi4loVt(const lomass mass, const double L);
double mpi6LloVt(const lomass mass, const Li li,const double L);
double mpi6RloVt(const lomass mass, const double L);
double mpi6loVt(const lomass mass, const Li li,const double L);
double mk4loVt(const lomass mass, const double L);
double mk6LloVt(const lomass mass, const Li li,const double L);
double mk6RloVt(const lomass mass, const double L);
double mk6loVt(const lomass mass, const Li li,const double L);
double meta4loVt(const lomass mass, const double L);
double meta6LloVt(const lomass mass, const Li li,const double L);
double meta6RloVt(const lomass mass, const double L);
double meta6loVt(const lomass mass, const Li li,const double L);
double fpi4loVt(const lomass mass, const double L);
double fpi6LloVt(const lomass mass, const Li li,const double L);
double fpi6RloVt(const lomass mass, const double L);
double fpi6loVt(const lomass mass, const Li li,const double L);
double fk4loVt(const lomass mass, const double L);
double fk6LloVt(const lomass mass, const Li li,const double L);
double fk6RloVt(const lomass mass, const double L);
double fk6loVt(const lomass mass, const Li li,const double L);
double feta4loVt(const lomass mass, const double L);
double feta6LloVt(const lomass mass, const Li li,const double L);
double feta6RloVt(const lomass mass, const double L);
double feta6loVt(const lomass mass, const Li li,const double L);
double qqup4loVt(const lomass mass, const double L);
double qqup6loVt(const lomass mass, const Li liin, const double L);
double qqup6LloVt(const lomass mass, const Li Liin, const double L);
double qqup6RloVt(const lomass mass, const double L);
double qqstrange4loVt(const lomass mass, const double L);
double qqstrange6loVt(const lomass mass, const Li liin, const double L);
double qqstrange6LloVt(const lomass mass, const Li Liin, const double L);
double qqstrange6RloVt(const lomass mass, const double L);

// Bessel functions
double mpi4loVb(const lomass mass, const double L);
double mpi6LloVb(const lomass mass, const Li li,const double L);
double mpi6RloVb(const lomass mass, const double L);
double mpi6loVb(const lomass mass, const Li li,const double L);
double mk4loVb(const lomass mass, const double L);
double mk6LloVb(const lomass mass, const Li li,const double L);
double mk6RloVb(const lomass mass, const double L);
double mk6loVb(const lomass mass, const Li li,const double L);
double meta4loVb(const lomass mass, const double L);
double meta6LloVb(const lomass mass, const Li li,const double L);
double meta6RloVb(const lomass mass, const double L);
double meta6loVb(const lomass mass, const Li li,const double L);
double fpi4loVb(const lomass mass, const double L);
double fpi6LloVb(const lomass mass, const Li li,const double L);
double fpi6RloVb(const lomass mass, const double L);
double fpi6loVb(const lomass mass, const Li li,const double L);
double fk4loVb(const lomass mass, const double L);
double fk6LloVb(const lomass mass, const Li li,const double L);
double fk6RloVb(const lomass mass, const double L);
double fk6loVb(const lomass mass, const Li li,const double L);
double feta4loVb(const lomass mass, const double L);
double feta6LloVb(const lomass mass, const Li li,const double L);
double feta6RloVb(const lomass mass, const double L);
double feta6loVb(const lomass mass, const Li li,const double L);
double qqup4loVb(const lomass mass, const double L);
double qqup6loVb(const lomass mass, const Li liin, const double L);
double qqup6LloVb(const lomass mass, const Li Liin, const double L);
double qqup6RloVb(const lomass mass, const double L);
double qqstrange4loVb(const lomass mass, const double L);
double qqstrange6loVb(const lomass mass, const Li liin, const double L);
double qqstrange6LloVb(const lomass mass, const Li Liin, const double L);
double qqstrange6RloVb(const lomass mass, const double L);


#endif //  MASSDECAYVEVLOV_H
