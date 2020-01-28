// Li.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of the class Li as well
// as a few associated functions


#ifndef LI_H
#define LI_H

#include <string>
#include <iostream>

class Li{
 private:
  double L1r,L2r,L3r,L4r,L5r,L6r,L7r,L8r,L9r,L10r,H1r,H2r;
  double mu;
  std::string name;
public:
  Li(const double l1r = 0.,const double l2r = 0.,const double l3r = 0.,
     const double l4r = 0.,const double l5r = 0.,const double l6r = 0.,
     const double l7r = 0.,const double l8r = 0.,const double l9r = 0.,
     const double l10r = 0.,const double h1r = 0.,const double h2r = 0.,
     const double mu = 0.77,
     const std::string Name = "nameless Li");
  ~Li(void);
  void setli(const int n, const double lin);
  void setli(const double lin, const int n);
  void setmu(const double muin);
  void setname(const std::string inputname);
  // out and input
  friend std::ostream & operator<<(std::ostream & os,const Li & bb);
  friend std::istream & operator>>(std::istream & is, Li & Liout);
  Li operator+(const Li & bb) const; // defines the sum of two sets of Li
  Li operator-(const Li & bb) const; // defines the difference of two sets of Li
  Li operator*(const double & aa) const; // multiplies Li with a number: aa*Li
  friend Li operator*(const double & aa,const  Li & bb);// aa*Li
  void changescale(const double newmu);
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & H1t,double & H2t,
	     double & mut,std::string nameout) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & H1t,double & H2t,
	     double & mut) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,double & H1t,double & H2t) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t,
	     double & mut) const;
  void out(double & L1t,double & L2t,double & L3t,
	     double & L4t,double & L5t,double & L6t,
	     double & L7t,double & L8t,double & L9t,
	     double & L10t) const;
  double out(const int n) const;
};
Li Lirandom(void);// all +-pi16
Li LirandomlargeNc(void);// only large Nc ones
Li LirandomlargeNc2(void);// single trace +-pi16, others pi16/3
#endif // LI_H
