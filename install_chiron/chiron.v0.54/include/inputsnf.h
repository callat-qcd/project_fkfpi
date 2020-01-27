// inputsnf.h is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014-2015 Johan Bijnens, v1.1
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// Header file for the definitions of classes physmass, Li and Ci as well
// as associated functions


#ifndef INPUTSNF_H
#define INPUTSNF_H

#include<string>
#include<vector>
#include<iostream>
using namespace std;

class lomassnf;

class quarkmassnf{
 private:
  int nq;// holds the number of quark species
  vector<double> B0mq;// length should be nq if initialized correctly
  double f0,mu;
 public:
  quarkmassnf(const double f0=0.090, const double muin =0.77, const int nqin=3);
  quarkmassnf(const vector<double> B0mqin,const double f0in=0.090,
	      const double muin =0.77);
  quarkmassnf(const lomassnf mass);
  void setf0(const double f0in=0.09);
  void setmu(const double muin=0.77);
  void setB0mq(const vector<double> B0mqin);
  void setB0mq(const double B0mqin, const int i);
  void setB0mq(const int i, const double B0mqin=0.);
  void out(vector<double> &B0mq) const;
  void out(vector<double> &B0mq, double &f0out, double &muout) const;
  void out(vector<double> &B0mq, double &f0out, double &muout, int &nq) const;
  int getnq(void) const;
  vector<double> getB0mq(void) const;
  double getB0mq(const int i) const;
  double getf0(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const quarkmassnf & bb); // output
  friend istream & operator>>(istream & is,quarkmassnf & bb); // input
  friend bool operator==(const quarkmassnf mass1,const quarkmassnf mass2);
  ~quarkmassnf(void);
};

class lomassnf{
 private:
  int nmass;// holds the number of quark species
  vector<double> mass;// length should be nq if initialized correctly
  double f0,mu;
 public:
  lomassnf(const double f0=0.090, const double muin =0.77, const int nqin=3);
  lomassnf(const vector<double> mass,const double f0in=0.090,
	      const double muin =0.77);
  lomassnf(const quarkmassnf qmass);
  void setf0(const double f0in=0.09);
  void setmu(const double muin=0.77);
  void setmass(const vector<double> B0mqin);
  void setmass(const double B0mqin, const int i);
  void setmass(const int i, const double B0mqin=0.);
  void out(vector<double> &massout) const;
  void out(vector<double> &massout, double &f0out, double &muout) const;
  void out(vector<double> &massout, double &f0out, double &muout, int &nq) const;
  int getnmass(void) const;
  vector<double> getmass(void) const;
  double getmass(const int i) const;
  double getf0(void) const;
  double getmu(void) const;
  friend ostream & operator<<(ostream & os,const lomassnf & bb); // output
  friend istream & operator>>(istream & is,lomassnf & bb); // input
  friend bool operator==(const lomassnf mass1,const lomassnf mass2);
  ~lomassnf(void);
};

#endif // INPUTSNF_H
