// Ci.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of class Ci as well
// as associated functions


#ifndef CI_H
#define CI_H

#include <string>
#include <iostream>

#include "Li.h"

class Ci{
 private:
    double Cr[95];
    double mu;
    std::string name;
 public:
    Ci(const double muin = 0.77,const std::string Name="nameless Ci");
    Ci(const double Crin[95],const double muin = 0.77,
       const std::string Name="nameless Ci");
    ~Ci(void);
  void setci(const int n, const double cin);
  void setci(const double cin, const int n);
  void setmu(const double muin);
  void setname(const std::string inputname);
  friend std::ostream & operator<<(std::ostream & os,const Ci & bb);
  friend std::istream & operator>>(std::istream & is,Ci & bb); // input
  Ci operator+(const Ci & bb) const; // defines the sum of two sets of Ci
  Ci operator-(const Ci & bb) const; // defines the difference of two sets of Ci
  Ci operator*(const double & xx) const; // multiplies Ci with a number
  friend Ci operator*(const double & aa,const  Ci & bb);// aa*Ci
  void changescale(const double newmu) const;
  void changescale(const double newmu, Li & li);
  void changescale(Li & li,const double newmu);
  void out(double Cit[95],double & mut,std::string namet) const;
  void out(double Cit[95],double & mut) const;
  void out(double Cit[95]) const;
  double out(const int n) const;
  Ci(const double fpi, const double mu
     ,const double MV, const double fV, const double gV
     ,const double fChi, const double aV
     ,const double MS, const double cd, const double cm, const double l3ss
     ,const double MEP, const double dmt
     ,const double MP, const double dm
	    );
};
Ci Cirandom(void);// all +-pi162
Ci CirandomlargeNc(void);// only large Nc ones
Ci CirandomlargeNc2(void);// single trace +-pi162, others pi162/3
#endif // CI_H
