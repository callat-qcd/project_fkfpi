// jbdadmul.cc  is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2014 Johan Bijnens, v1.0
// except for the CERNLIB parts which have their own copyright (usually CERN)
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.

// naively converted to C++ from the CERNLIB Fortran code JB 2008
// provides wrapper routines for jbdadmul, with work space allocated
// automatically as needed

// a,b contain the boundaries, eps wanted error,relerr achieved error,
// ifail indicates something was wrong, see dadmul in cernlib

/*
F
    (type according to t) Name of a user-supplied FUNCTION subprogram,
      declared EXTERNAL in the calling program. 
N
    (INTEGER) Number n of dimensions . 
A,B
  One-dimensional arrays of length n ,
a,b contain the lower and upper limits of integration, respectively.

MINPTS
    (INTEGER) Minimum number of function evaluations requested. 
 Must not exceed MAXPTS. 
MAXPTS
    (INTEGER) Maximum number ( > 2**N + 2N(N + 1) + 1) of function evaluations
                to be allowed. 
EPS
    (type according to t) Specified relative accuracy. 
WK
    (type according to t) One-dimensional array of length IWK,
         used as working space. 
IWK
    (INTEGER) Length (  (2N + 3) * (1 + MAXPTS/(2**N + 2N(N + 1) + 1))/2)
      of WK. 
RESULT
    (type according to t) Contains, on exit, an approximate value of
    the integral tex2html_wrap_inline211 . 
RELERR
    (type according to t) Contains, on exit, an estimation of the
     relative accuray of RESULT. 
NFNEVL
    (INTEGER) Contains, on exit, the number of function evaluations performed.

IFAIL
    (INTEGER) On exit:

    0
        Normal exit. At least MINPTS and at most MAXPTS calls to the
               function F were performed. 
    1
        MAXPTS is too small for the specified accuracy EPS. RESULT and
        RELERR contain the values obtainable for the specified value of MAXPTS. 
    2
        IWK is too small for the specified number MAXPTS of function
        evaluations.
        RESULT and RELERR contain the values obtainable for the specified
        value of IRK. 
    3
        n < 2 or n> 15 MINPTS > MAXPTS, or MAXPTS < 2**N + 2N(N + 1) + 1.
        RESULT and RELERR are set equal to zero. 
*/

#include <cmath>
#include <iostream>

void jbdadmul(double (*f)(double x[]),const int n,double a[],double b[],
	    const int minpts,const int maxpts,const double eps,double wk[],
	    const int iwk,double &result,double &relerr,
	    int &nfnevl,int &ifail);

double jbdad2(double (*fcn)(double x[]),double a[],double b[],
	      const double releps, double &relerr,int &ifail){
  const int ndim =2;
  const int minpts=100;
  const int maxpts=1000000;
  const int iwk = ((2*ndim+3)*(1+maxpts/(pow(2,ndim)+2*ndim*(ndim+1)+1))/2);
  int nfnevl;
  double wk[iwk+1];
  double result;
  jbdadmul(fcn,ndim,a,b,minpts,maxpts,releps,wk,iwk,result,relerr,nfnevl,ifail);
  if(ifail!=0){
    std::cout << "#accuracy not reached in jbdad2 in "<<nfnevl 
	      <<" calls, relerr is "<<relerr<<std::endl;
  }
  return result;
}

// due to the very large array, only use
double jbdad3(double (*fcn)(double x[]),double a[],double b[],
	      const double releps, double &relerr,int &ifail){
  const int ndim = 3;
  const int minpts=300;
  const int maxpts=2000000;
  const int iwk = ((2*ndim+3)*(1+maxpts/(pow(2,ndim)+2*ndim*(ndim+1)+1))/2);
  int nfnevl;
  double wk[iwk+1];
  double result;
  jbdadmul(fcn,ndim,a,b,minpts,maxpts,releps,wk,iwk,result,relerr,nfnevl,ifail);
  if(ifail!=0){
    std::cout << "#accuracy not reached in jbdad3 in "<<nfnevl 
	      <<" calls, relerr is "<<relerr<<std::endl;
  }
  return result;
}

// due to the very large array, only use
//double jbdad4(double (*fcn)(double x[]),double a[],double b[],
//	      const double releps, double &relerr,int &ifail){
//  const int ndim = 4;
//  const int minpts=500;
//  const int maxpts=5000000;
//  const int iwk = ((2*ndim+3)*(1+maxpts/(pow(2,ndim)+2*ndim*(ndim+1)+1))/2);
//  int nfnevl;
//  double wk[iwk+1];
//  double result;
//  jbdadmul(fcn,ndim,a,b,minpts,maxpts,releps,wk,iwk,result,relerr,nfnevl,ifail);
//  if(ifail!=0){
//    std::cout << "#accuracy not reached in jbdad4 in "<<nfnevl 
//	      <<" calls, relerr is "<<relerr<<std::endl;
//  }
//  return result;
//}
//
//#define WORK 2000000
//static double wk5[WORK];// too big to be on the stack
//double jbdad5(double (*fcn)(double x[]),double a[],double b[],
//	      const double releps, double &relerr,int &ifail){
//  const int ndim = 5;
//  const int minpts=1000;
//  const int maxpts=20000000;// if too large, remember to allow for more memory
//  const int iwk = ((2*ndim+3)*(1+maxpts/(pow(2,ndim)+2*ndim*(ndim+1)+1))/2);
//  int nfnevl;
//  if(iwk > (WORK-1)){
//    std::cout << "toolittle workspace in jbdqad5 !"<<std::endl;
//    return 0.;
//  }
//  double result;
//  jbdadmul(fcn,ndim,a,b,minpts,maxpts,releps,wk5,iwk,result,relerr,nfnevl,ifail);
//  if(ifail!=0){
//    std::cout << "#accuracy not reached in jbdad5 in "<<nfnevl 
//	      <<" calls, relerr is "<<relerr<<std::endl;
//  }
//  return result;
//}

void jbdadmul(double (*f)(double x[]),const int n,double aa[],double bb[],
	    const int minpts,const int maxpts,const double eps,double wk[],
	    const int iwk,double &result,double &relerr,
	    int &nfnevl,int &ifail){
 
  bool ldv;
  double ctr[16],wth[16],wthl[16],z[16],a[16],b[16];
  //    dimension w(2:15,5),wp(2:15,3)

  const double r1 = 1;
  const double hf = r1/2.;
  const double xl2 =  0.358568582800318073;
  const double xl4 =  0.948683298050513796;
  const double xl5 =  0.688247201611685289;
 
  const double w2 =  980.*r1/6561.;
  const double w4 = 200.*r1/19683.;
  const double wp2 =  245.*r1/486.;
  const double wp4 = 25.*r1/729.;
 
  double w[16][6] = {
    {0.,0.                       ,0.,0.                       ,0.,0.                      },
    {0.,0.                       ,0.,0.                       ,0.,0.                      },
    {0.,-0.193872885230909911e+00,0., 0.518213686937966768e-01,0.,0.871183254585174982e-01},
    {0.,-0.555606360818980835e+00,0., 0.314992633236803330e-01,0.,0.435591627292587508e-01},
    {0.,-0.876695625666819078e+00,0., 0.111771579535639891e-01,0.,0.217795813646293754e-01},
    {0.,-0.115714067977442459e+01,0.,-0.914494741655235473e-02,0.,0.108897906823146873e-01},
    {0.,-0.139694152314179743e+01,0.,-0.294670527866686986e-01,0.,0.544489534115734364e-02},
    {0.,-0.159609815576893754e+01,0.,-0.497891581567850424e-01,0.,0.272244767057867193e-02},
    {0.,-0.175461057765584494e+01,0.,-0.701112635269013768e-01,0.,0.136122383528933596e-02},
    {0.,-0.187247878880251983e+01,0.,-0.904333688970177241e-01,0.,0.680611917644667955e-03},
    {0.,-0.194970278920896201e+01,0.,-0.110755474267134071e+00,0.,0.340305958822333977e-03},
    {0.,-0.198628257887517146e+01,0.,-0.131077579637250419e+00,0.,0.170152979411166995e-03},
    {0.,-0.198221815780114818e+01,0.,-0.151399685007366752e+00,0.,0.850764897055834977e-04},
    {0.,-0.193750952598689219e+01,0.,-0.171721790377483099e+00,0.,0.425382448527917472e-04},
    {0.,-0.185215668343240347e+01,0.,-0.192043895747599447e+00,0.,0.212691224263958736e-04},
    {0.,-0.172615963013768225e+01,0.,-0.212366001117715794e+00,0.,0.106345612131979372e-04}
  };

  double wp[16][4] = {
    {0.,0.                       ,0.,0.                       },
    {0.,0.                       ,0.,0.                       },
    {0.,-0.133196159122085045e+01,0., 0.445816186556927292e-01},
    {0.,-0.229218106995884763e+01,0.,-0.240054869684499309e-01},
    {0.,-0.311522633744855959e+01,0.,-0.925925925925925875e-01},
    {0.,-0.380109739368998611e+01,0.,-0.161179698216735251e+00},
    {0.,-0.434979423868312742e+01,0.,-0.229766803840877915e+00},
    {0.,-0.476131687242798352e+01,0.,-0.298353909465020564e+00},
    {0.,-0.503566529492455417e+01,0.,-0.366941015089163228e+00},
    {0.,-0.517283950617283939e+01,0.,-0.435528120713305891e+00},
    {0.,-0.517283950617283939e+01,0.,-0.504115226337448555e+00},
    {0.,-0.503566529492455417e+01,0.,-0.572702331961591218e+00},
    {0.,-0.476131687242798352e+01,0.,-0.641289437585733882e+00},
    {0.,-0.434979423868312742e+01,0.,-0.709876543209876532e+00},
    {0.,-0.380109739368998611e+01,0.,-0.778463648834019195e+00},
    {0.,-0.311522633744855959e+01,0.,-0.847050754458161859e+00}
  };

  double rgnvol,sum1,sum2,sum3,sum4,sum5,difmax,dif,f2,f3,rgncmp,rgnval,rgnerr;
  int j,j1,k,l,m,idvaxn,isbtmp,isbtpp,idvax0;

  result=0.;
  double abserr=0.;
  ifail=3;
  if((n < 2) || (n > 15)) return;
  if(minpts > maxpts) return;
  int ifncls=0;
  ldv=false;
  double twondm=pow(2,n);
  int irgnst=2*n+3;
  int irlcls=pow(2,n)+2*n*(n+1)+1;
  int isbrgn=irgnst;
  int isbrgs=irgnst;
  if(maxpts < irlcls) return;

  for(j=1;j<=n;j++){
    a[j]=aa[j-1];
    b[j]=bb[j-1];
  }

  for(j=1;j<=n;j++){
    ctr[j]=(b[j]+a[j])*hf;
    wth[j]=(b[j]-a[j])*hf;
  }
label20:
  rgnvol=twondm;
  for(j=1;j<=n;j++){
    rgnvol=rgnvol*wth[j];
    z[j]=ctr[j];
  }
  sum1=f(&z[1]);
  difmax=0.;
  sum2=0.;
  sum3=0.;
  for(j=1;j<=n;j++){
    z[j]=ctr[j]-xl2*wth[j];
    f2=f(&z[1]);
    z[j]=ctr[j]+xl2*wth[j];
    f2=f2+f(&z[1]);
    wthl[j]=xl4*wth[j];
    z[j]=ctr[j]-wthl[j];
    f3=f(&z[1]);
    z[j]=ctr[j]+wthl[j];
    f3=f3+f(&z[1]);
    sum2=sum2+f2;
    sum3=sum3+f3;
    dif = fabs(7.*f2-f3-12.*sum1);
    difmax = fmax(dif,difmax);
    if(difmax == dif) idvaxn=j;
    z[j]=ctr[j];
  } 
  sum4=0.;
  for(j=2;j<=n;j++){
    j1=j-1;
    for(k=j;k<=n;k++){
      for(l=1;l<=2;l++){
	wthl[j1]=-wthl[j1];
	z[j1]=ctr[j1]+wthl[j1];
	for (m=1;m<=2;m++){
	  wthl[k]=-wthl[k];
	  z[k]=ctr[k]+wthl[k];
	  sum4=sum4+f(&z[1]);
	}
      }
      z[k]=ctr[k];
    }
    z[j1]=ctr[j1];
  }
  sum5=0.;
  for(j=1;j<=n;j++){
    wthl[j]=-xl5*wth[j];
    z[j]=ctr[j]+wthl[j];
  }
label90:
  sum5=sum5+f(&z[1]);
  for(j=1;j<=n;j++){
    wthl[j]=-wthl[j];
    z[j]=ctr[j]+wthl[j];
    if((wthl[j]/(b[j]-a[j])) > 0.) goto label90;
  }
  rgncmp=rgnvol*(wp[n][1]*sum1+wp2*sum2+wp[n][3]*sum3+wp4*sum4);
  rgnval=w[n][1]*sum1+w2*sum2+w[n][3]*sum3+w4*sum4+w[n][5]*sum5;
  rgnval=rgnvol*rgnval;
  rgnerr=fabs(rgnval-rgncmp);
  result=result+rgnval;
  abserr=abserr+rgnerr;
  ifncls=ifncls+irlcls;
 
  if(ldv) {
label110:
    isbtmp=2*isbrgn;
    if(isbtmp > isbrgs) goto label160;
    if(isbtmp < isbrgs) {
      isbtpp=isbtmp+irgnst;
      if(wk[isbtmp] < wk[isbtpp]) isbtmp=isbtpp;
    }
    if(rgnerr >= wk[isbtmp]) goto label160;
    for(k=0;k<=irgnst-1;k++) wk[isbrgn-k]=wk[isbtmp-k];
    isbrgn=isbtmp;
    goto label110;
  }
label140:
  isbtmp=(isbrgn/(2*irgnst))*irgnst;
  if((isbtmp >= irgnst) && (rgnerr > wk[isbtmp])) {
    for(k=0;k <= irgnst-1;k++) wk[isbrgn-k]=wk[isbtmp-k];
    isbrgn=isbtmp;
    goto label140;
  }
label160:
  wk[isbrgn]=rgnerr;
  wk[isbrgn-1]=rgnval;
  wk[isbrgn-2]=idvaxn;

  for(j=1;j<=n;j++){
    isbtmp=isbrgn-2*j-2;
    wk[isbtmp+1]=ctr[j];
    wk[isbtmp]=wth[j];
  }

  if(ldv){
    ldv=false;
    ctr[idvax0]=ctr[idvax0]+2*wth[idvax0];
    isbrgs=isbrgs+irgnst;
    isbrgn=isbrgs;
    goto label20;
  }
  relerr=abserr/(1e-315+fabs(result));//1e-315 added for problem if function=0
  if((isbrgs+irgnst) > iwk) ifail=2;
  if( (ifncls+2*irlcls) > maxpts) ifail=1;
  if((relerr <= eps) && (ifncls >= minpts) ) ifail=0;
  if(ifail == 3) {
    ldv=true;
    isbrgn=irgnst;
    abserr=abserr-wk[isbrgn];
    result=result-wk[isbrgn-1];
    idvax0=wk[isbrgn-2];
    for(j=1;j<=n;j++){
      isbtmp=isbrgn-2*j-2;
      ctr[j]=wk[isbtmp+1];
      wth[j]=wk[isbtmp];
    }
    wth[idvax0]=hf*wth[idvax0];
    ctr[idvax0]=ctr[idvax0]-wth[idvax0];
    goto label20;
  }
  nfnevl=ifncls;
  return;
}
