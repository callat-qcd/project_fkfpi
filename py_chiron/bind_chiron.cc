// This is modified from a file from J. Bijnens
// This computes the FF(x) function defined in Eqs(8-17) of 1711.11328
// by Ananthanarayan, Bijnens, Friot, and Ghosh

#include<iostream>
#include<iomanip>
#include<cmath>
#include "inputs.h"
#include "Li.h"
#include "Ci.h"
#include "massdecayvev.h"
#include "oneloopintegrals.h"
#include "sunsetintegrals.h"
#include "quenchedsunsetintegrals.h"
#include "jbnumlib.h"

#include <pybind11/pybind11.h>

const double pi16 = 1./(16.*M_PI*M_PI);
const double pi162 = 1./pow(16.*M_PI*M_PI,2);
const double pi = M_PI;
const double pi2 = M_PI*M_PI;

// changing the hhb function to hh
double hhb(int prop,double mp2, double mk2, double me2, double qsq, double mu2){
    double abmp2 = Ab(mp2,mu2);
    double abmk2 = Ab(mk2,mu2);
    double abme2 = Ab(me2,mu2);
    double abqsq = Ab(qsq,mu2);
    switch(prop){
        case 1:
        return hh(1,mp2,mk2,me2,qsq,mu2)-(abmp2*abmp2/mp2+abme2*abme2/me2
            +abmk2*abmk2/mk2
            +pi16*(abmp2+abmk2+abme2-0.5*abqsq));
        case 2:
        return hh(2,mp2,mk2,me2,qsq,mu2)-(abmp2*abmp2/(mp2*mp2)
            -pi16*(abmp2)/mp2);
        case 3:
        return hh(3,mp2,mk2,me2,qsq,mu2)-(abmk2*abmk2/(mk2*mk2)
            -pi16*(abmk2)/mk2);
        }
        std::cout << "hhb called with wrong prop\n";
    }

// this is the approximate F_F in some draft of the paper
double FFKPapprox(const double x){
    //     double x = mp2/mk2;
    const double a1 = 0.8740, a2 = -2.172,  a3 = 0.8294, a4 = -0.4583,
                 a5 = 3.716,  a6 = -0.1113, a7 = 0.8776, a8 = -1.635, a9 = 1.4697,
                 a10 = -0.1406, a11 = -1.343, a12 = 0.2731, a13 = -0.2109;
    return ( a1
           +(a2  + a3 *log(x) +a4  *log(x)*log(x) )*x
           +(a5  + a6 *log(x) +a7  *log(x)*log(x) )*pow(x,2)
           +(a8  + a9 *log(x) +a10 *log(x)*log(x) )*pow(x,3)
           +(a11 + a12*log(x) +a13 *log(x)*log(x) )*pow(x,4));
       }

double FFKP(double rho){
    // the return value is independent of \mu^2 up to numerics
    // We will pass in mpi**2 / mK**2
    double mu2 = 0.77*0.77;
    double mp2 = rho;
    double mk2 = 1.0;
    double me2 = 4./3.*mk2 -mp2/3.;
    double FKP6test=0.;
    double lampi  = log(mp2/mu2);
    double lamk   = log(mk2/mu2);
    double lameta = log(me2/mu2);
    FKP6test = + pi162 * (
        27./8.*pow(mp2,-1)*pow(mk2,3)
        - 1./48.*pow(mp2,-1)*pow(mk2,3)*pi2
        -755./6912.*pow(mk2,2)
        +113./648.*pow(mk2,2)*pi2
        -271./216.*mp2*mk2
        -715./2592.*mp2*mk2*pi2
        -22757./6912. *pow(mp2,2)
        +335./2592.   *pow(mp2,2)*pi2
        -1./8.  *pow(mp2,3)*pow(mk2,-1)
        -1./144.*pow(mp2,3)*pow(mk2,-1)*pi2 );

    FKP6test +=  +pow(lameta,2)*pi162 * (
        - 1./8. *pow(mp2,-1)*pow(mk2,3)
        + 5./64.*pow(mp2,3) *pow(mk2,-1) );

    FKP6test +=  +lamk*pi162 * ( 9./16.*pow(mk2,2) );

    FKP6test +=  + lamk*lameta*pi162 * ( 1./4.*pow(mp2,-1)*pow(mk2,3) );

    FKP6test +=  + pow(lamk,2)*pi162 * (
        -1./8.*pow(mp2,-1)*pow(mk2,3) - 3./32.*pow(mk2,2) );

    FKP6test +=  + lampi*pi162 * (  - 9./16.*pow(mk2,2) );

    FKP6test +=  + lampi*lameta*pi162 * (  - 5./32.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + lampi*lamk*pi162 * ( 3./16.*pow(mk2,2) );

    FKP6test +=  + pow(lampi,2)*pi162 * (
        -3./32.*pow(mk2,2) + 5./64.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + hhb(1,mp2,mk2,mk2,mp2,mu2) * (
        -3./16.*pow(mp2,-1)*pow(mk2,2) - 17./12.*mk2 - 1./16.*mp2 );

    FKP6test +=  + hhb(1,mp2,me2,me2,mp2,mu2) * (  - 1./36.*mp2 );

    FKP6test +=  + hhb(1,mk2,mp2,mp2,mk2,mu2) * (
        3./64.*mk2 + 17./16. *mp2 + 9./64.*pow(mp2,2)*pow(mk2,-1) );

    FKP6test +=  + hhb(1,mk2,mp2,me2,mk2,mu2) * (
        5./96.*mk2 + 7./6.*mp2 - 7./32.*pow(mp2,2)*pow(mk2,-1) );

    FKP6test +=  + hhb(1,mk2,me2,me2,mk2,mu2) * (
        197./192.*mk2 - 83./144.*mp2 + 5./64.*pow(mp2,2)*pow(mk2,-1) );

    FKP6test +=  + hhb(1,me2,mk2,mk2,mp2,mu2) * (
        -15./16.*pow(mp2,-1)*pow(mk2,2) + 13./36.*mk2 - 13./144.*mp2 );

    FKP6test +=  + hhb(2,mp2,mk2,mk2,mp2,mu2) * ( 2*mp2*mk2 );

    FKP6test +=  + hhb(2,mp2,mk2,me2,mk2,mu2) * (
        1./4.*mp2*mk2 - 15./8.*pow(mp2,2) + 3./8.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + hhb(2,mp2,me2,me2,mp2,mu2) * ( 1./36.*pow(mp2,2) );

    FKP6test +=  + hhb(2,mk2,mp2,mp2,mk2,mu2) * (  - 3./2.*mp2*mk2 );

    FKP6test +=  + hhb(2,mk2,mp2,me2,mk2,mu2) * (  - 1./2.*pow(mk2,2));

    FKP6test +=  + hhb(2,mk2,me2,me2,mk2,mu2) * (  - 1./9.*pow(mk2,2));

    FKP6test +=  + hhb(2,me2,mk2,mp2,mk2,mu2) * (
        -11./18.*pow(mk2,2) - 41./72.*mp2*mk2 - 11./72.*pow(mp2,2)
        + 1./12.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + hhb(2,me2,mk2,mk2,mp2,mu2) * (
        pow(mp2,-1)*pow(mk2,3) - 91./108.*pow(mk2,2)
        +5./27.*mp2*mk2 - 1./108.*pow(mp2,2) );

    FKP6test +=  + hhb(3,mp2,mk2,mk2,mp2,mu2) * (
        3./4.*pow(mp2,-1)*pow(mk2,3) + 13./6.*pow(mk2,2) + 1./12.*mp2*mk2 );

    FKP6test +=  + hhb(3,mp2,me2,me2,mp2,mu2) * (
        2./27.*mp2*mk2 - 1./54.*pow(mp2,2) );

    FKP6test +=  + hhb(3,mk2,mp2,mp2,mk2,mu2) * (
        -1./16.*mp2*mk2 -13./8.*pow(mp2,2) - 9./16.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + hhb(3,mk2,me2,me2,mk2,mu2) * (
        -445./108.*pow(mk2,2) + 173./48.*mp2*mk2
        -229./216.*pow(mp2,2) + 5./48.*pow(mp2,3)*pow(mk2,-1) );

    FKP6test +=  + hhb(3,me2,mk2,mk2,mp2,mu2) * (
        2*pow(mp2,-1)*pow(mk2,3) - 1./2.*pow(mk2,2) + 1./6.*mp2*mk2 );

    return FKP6test/pi162/(mk2*mk2);
    }


PYBIND11_MODULE(chiron, m) {
    m.doc() = "Chiron python module"; // optional module docstring
    m.def("FF", &FFKP, "Exact F_F for fK/FPi");
    m.def("FFApprox", &FFKPapprox, "Approximate F_F for FK/FPi");
    }
