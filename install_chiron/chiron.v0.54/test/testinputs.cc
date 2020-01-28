// testinputs.cc is part of the CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// test the masses etc in and output routines in inputs.cc

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include "inputs.h"

int main(void){
  // a number of different input values
  double  mpi = 0.15;
  double  mk = 0.5;
  double  meta = 0.55;
  double fpi = 0.95;
  double  mu = 0.9;
  double mp0 = 0.12;
  double mk0 = 0.41;
  double f0 = 0.88;
  double B0mhat = 0.011;
  double B0ms = 0.28;
  // creating a set with standard input values and puttin it out;
  physmass mass1;
  lomass mass1l;
  quarkmass mass1q;
  cout << "First mass case\n" << mass1 << mass1l << mass1q;
  // creating a set with the other inputs above and then changin one
  physmass mass2(mpi,mk,meta,fpi,mu);
  lomass mass2l(mp0,mk0,f0,mu);
  quarkmass mass2q(B0mhat,B0ms,f0,mu);
  cout << "Second mass case\n" << mass2 << mass2l << mass2q;
  physmass mass3 = mass2;
  lomass mass3l = mass2l;
  quarkmass mass3q = mass2q;
  // checking the conversion
  lomass mass3lp = mass2q;
  quarkmass mass3qp = mass2l;
  mass3.setmk(0.4);
  cout << "Third mass case\n" << mass3 << mass3l << mass3q
       << mass3lp << mass3qp;
  // writing them out to a file and then reading back in
  ofstream ouf("test.dat");
  ouf << mass2 << mass2l << mass2q;
  ouf.close();
  physmass mass4; lomass mass4l; quarkmass mass4q;
  ifstream inf("test.dat");
  inf >> mass4 >> mass4l >> mass4q;
  inf.close();
  cout << "written to and read from file\n" << mass4 << mass4l << mass4q;
  // readin out the data
  mass4.out(mpi,mk,meta,fpi,mu);
  cout <<"mpi,mk,meta,fpi,mu\n"
       <<mpi<<' '<<mk<<' '<<meta<<' '<<fpi<<' '<<mu<<endl;
  mass4l.out(mp0,mk0,f0,mu);
  cout << "mp0, mk0, f0, mu\n"
       << mp0<<' '<<mk0<<' '<<f0<<' '<<mu<<'\n';
  mass4q.out(B0mhat,B0ms,f0,mu);
  cout << "B0mhat, B0ms, f0, mu\n"
       << B0mhat<<' '<<B0ms<<' '<<f0<<' '<<mu<<'\n';
  // same but now read out with the separate fnuctions, resetting fpi, f0
  mass4.setfpi(0.12);
  mass4l.setf0(0.095);
  mass4l.setmu(0.81);
  mass4q.setf0(0.084);
  mpi = mass4.getmpi();
  mk  = mass4.getmk();
  meta= mass4.getmeta();
  fpi = mass4.getfpi();
  mu  = mass4.getmu();
  cout <<"mpi,mk,meta,fpi,mu\n"
       <<mpi<<' '<<mk<<' '<<meta<<' '<<fpi<<' '<<mu<<endl;
  mp0 = mass4l.getmp0();
  mk0 = mass4l.getmk0();
  f0 = mass4l.getf0();
  mu = mass4l.getmu();
  cout << "mp0, mk0, f0, mu\n"
       << mp0<<' '<<mk0<<' '<<f0<<' '<<mu<<'\n';
  B0mhat = mass4q.getB0mhat();
  B0ms = mass4q.getB0ms();
  f0 = mass4q.getf0();
  mu = mass4q.getmu();
  cout << "B0mhat, B0ms, f0, mu\n"
       << B0mhat<<' '<<B0ms<<' '<<f0<<' '<<mu<<'\n';
  // checking equality of two masses (note only to about 7 digits
  // to avoid calculated masses to become unequal
  if( mass1 == mass4){
    cout << "mass4 equal to mass1\n";}
  else{
    cout << "mass4 not equal to mass1\n";}
  if( mass1l == mass4l){
    cout << "mass4l equal to mass1l\n";}
  else{
    cout << "mass4l not equal to mass1l\n";}
  if( mass1q == mass4q){
    cout << "mass4q equal to mass1q\n";}
  else{
    cout << "mass4q not equal to mass1q\n";}
  physmass mass5; lomass mass5l; quarkmass mass5q;
  if( mass1 == mass5){
    cout << "mass5 equal to mass1\n";}
  else{
    cout << "mass5 not equal to mass1\n";}
  if( mass1l == mass5l){
    cout << "mass5l equal to mass1l\n";}
  else{
    cout << "mass5l not equal to mass1l\n";}
  if( mass1q == mass5q){
    cout << "mass5q equal to mass1q\n";}
  else{
    cout << "mass5q not equal to mass1q\n";}

  return 0;
}
