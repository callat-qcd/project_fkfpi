// Ki.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of class Ki as well
// as associated functions


#ifndef KI_H
#define KI_H

#include <string>
#include <iostream>

#include "Linf.h"

class Ki{
 private:
    double Kr[116];
    double mu;
    std::string name;
    int nf;
 public:
    Ki(const double muin = 0.77,const std::string Name="nameless Ki",
       const int nfin = 3);
    Ki(const double Krin[116],const double muin = 0.77,
       const std::string Name="nameless Ki", const int nfin = 3);
    ~Ki(void);
  void setki(const int n, const double kin);
  void setki(const double kin, const int n);
  void setmu(const double muin);
  void setname(const std::string inputname);
  friend std::ostream & operator<<(std::ostream & os,const Ki & bb);
  friend std::istream & operator>>(std::istream & is,Ki & bb); // input
  Ki operator+(const Ki & bb) const; // defines the sum of two sets of Ki
  Ki operator-(const Ki & bb) const; // defines the difference of two sets of Ki
  Ki operator*(const double & xx) const; // multiplies Ki with a number
  friend Ki operator*(const double & aa,const Ki & bb);// aa*Ki
  void changescale(const double newmu) const;
  void changescale(const double newmu, Linf & linf);
  void changescale(Linf & linf,const double newmu);
  void out(double Kit[116],double & mut,std::string namet) const;
  void out(double Kit[116],double & mut) const;
  void out(double Kit[116]) const;
  double out(const int n) const;
  int getnf(void) const;
  Ki(const double fpi, const double mu
     ,const double MV, const double fV, const double gV
     ,const double fChi, const double aV
     ,const double MS, const double cd, const double cm, const double l3ss
     ,const double MEP, const double dmt
     ,const double MP, const double dm
	    );
};
Ki Kirandom(void);// all +-pi162
#endif // KI_H
