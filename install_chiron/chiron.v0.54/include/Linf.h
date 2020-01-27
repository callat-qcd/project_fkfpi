// Linf.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of the class Li as well
// as a few associated functions


#ifndef LINF_H
#define LINF_H

#include <string>
#include <iostream>

// L11 added for the two-flavour case to be possible
class Linf{
 private:
  int nf;// number of flavours
  double L0r,L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,L11r,H1r,H2r;
  double mu; 
  std::string name;
public:
  Linf(const double l0r = 0., const double l1r = 0.,const double l2r = 0.,
       const double l3r = 0., const double l4r = 0.,const double l5r = 0.,
       const double l6r = 0., const double l7r = 0.,const double l8r = 0.,
       const double l9r = 0., const double l10r = 0., const double l11r=0.,
       const double h1r = 0.,const double h2r = 0., const double mu = 0.77,
       const std::string Name = "nameless Linf", const int nfin=3);
  ~Linf(void);
  void setnf(const int nfin);
  void setlinf(const int n, const double lin);
  void setlinf(const double lin, const int n);
  void setmu(const double muin);
  void setname(const std::string inputname);
  // out and input
  friend std::ostream & operator<<(std::ostream & os,const Linf & bb);
  friend std::istream & operator>>(std::istream & is, Linf & Liout);
  Linf operator+(const Linf & bb) const; // defines the sum of two sets of Li
  Linf operator-(const Linf & bb) const; // defines the difference of two sets of Li
  Linf operator*(const double & aa) const; // multiplies Li with a number: aa*Li
  friend Linf operator*(const double & aa,const  Linf & bb);// aa*Li
  void changescale(const double newmu);
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t,double & H1t,double & H2t,
	   double & mut,std::string nameout, int & nft) const;
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t,double & H1t,double & H2t,
	   double & mut, int & nft) const;
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t,double & H1t,double & H2t,
	   int & nft) const;
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t,
	   double & mut, int & nft) const;
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t, int & nft) const;
  void out(double & L0t,double & L1t,double & L2t,double & L3t,
	   double & L4t,double & L5t,double & L6t,
	   double & L7t,double & L8t,double & L9t,
	   double & L10t,double & L11t) const;
  double out(const int n) const;
  int getnf(void) const;
};
Linf Linfrandom(void);// all +-pi16
#endif // LINF_H
