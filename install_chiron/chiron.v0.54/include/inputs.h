// inputs.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of classes physmass, Li and Ci as well
// as associated functions


#ifndef INPUTS_H
#define INPUTS_H

#include <string>
#include <iostream>
using namespace std;

class physmass{
 private:
  double mpi,mk,meta,fpi,mu;
 public:
  physmass(const double mpiin =0.135, const double mkin=0.495,
	   const double metain =0.548,  const double fpiin=0.0922,
	   const double muin =0.77);
  void setmpi(const double mpiin=0.135);
  void setmk(const double mkin=0.495);
  void setmeta(const double metain=0.548);
  void setfpi(const double fpiin=0.0922);
  void setmu(const double muiin=0.77);
  void out(double &mpiout, double &mkout, double &metaout, double &fpiout,
      double &muout) const;
  double getmpi(void) const;
  double getmk(void) const;
  double getmeta(void) const;
  double getfpi(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const physmass & bb); // output
  friend istream & operator>>(istream & is,physmass & bb); // input
  friend bool operator==(const physmass mass1,const physmass mass2);
  ~physmass(void);
};

class quarkmass;
class lomass{
 private:
  double mp0,mk0,f0,mu;
 public:
  lomass(const double mp0 =0.135, const double mk0=0.495,
	   const double f0=0.090,   const double muin =0.77);
  lomass(const quarkmass mass);
  void setmp0(const double mp0in=0.135);
  void setmk0(const double mk0in=0.495);
  void setf0(const double f0in=0.09);
  void setmu(const double muin=0.77);
  void out(double &mp0out, double &mk0out, double &f0out, double &muout) const;
  double getmp0(void) const;
  double getmk0(void) const;
  double getf0(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const lomass & bb); // output
  friend istream & operator>>(istream & is,lomass & bb); // input
  friend bool operator==(const lomass mass1,const lomass mass2);
  ~lomass(void);
};


class quarkmass{
 private:
  double B0mhat,B0ms,f0,mu;
 public:
  quarkmass(const double B0mhat =0.01, const double B0ms=0.25,
	   const double f0=0.090,   const double muin =0.77);
  quarkmass(const lomass mass);
  void setB0mhat(const double B0mhatin=0.01);
  void setB0ms(const double B0msin=0.425);
  void setf0(const double f0in=0.09);
  void setmu(const double muin=0.77);
  void out(double &B0mhatout, double &B0msout, double &f0out, double &muout) const;
  double getB0mhat(void) const;
  double getB0ms(void) const;
  double getf0(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const quarkmass & bb); // output
  friend istream & operator>>(istream & is,quarkmass & bb); // input
  friend bool operator==(const quarkmass mass1,const quarkmass mass2);
  ~quarkmass(void);
};
#endif // INPUTS_H
