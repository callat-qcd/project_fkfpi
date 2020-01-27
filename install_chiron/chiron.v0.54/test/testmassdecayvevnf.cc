// testmassdecayvevnf.cc is part of the 
// numerical library jbnumlib included with the 
// CHIRON ChPT at two loops program collection
// Copyright (C) 2015 Johan Bijnens, v1.0
// CHIRON is licenced under the GNU GPL version 2 or later,
// see COPYING for details.
// Please respect the Guidelines, see GUIDELINES for details.



#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "inputsnf.h"
#include "Li.h"
#include "Ci.h"
#include "Linf.h"
#include "Ki.h"
#include "massdecayvevlo.h"
#include "massdecayvevPQ.h"
#include "massdecayvevnf.h"
#include "massdecayvevnfPQ.h"
#include "massdecayvevloV.h"
#include "massdecayvevPQV.h"
#include "massdecayvevnfV.h"
#include "massdecayvevnfPQV.h"
#include "finitevolumeoneloopintegrals.h"
#include "finitevolumesunsetintegrals.h"
using namespace std;
int main(void){

  double f0 = 0.0877;
  double mu = 0.77;
  Li LiBE14=Li(    0.000528552266181,0.000814177210776,-0.003071901580820,
		   0.000300000000000,0.001014976780543,0.000144108338197,
		   -0.000342935304325,0.000473156165745,0.005930000000000,
		   0.000000000000000,0.000000000000000,0.000000000000000,
		   mu,"Linf BE14");
  Linf LinfBE14=Linf(0.0,
		   0.000528552266181,0.000814177210776,-0.003071901580820,
		   0.000300000000000,0.001014976780543,0.000144108338197,
		   -0.000342935304325,0.000473156165745,0.005930000000000,
		   0.000000000000000,0.000000000000000,0.000000000000000,
		   0.,mu,"Linf BE14",3);
  Ci Ci0;
  Ki Ki0;

  setprecisionfinitevolumeoneloopt(1e-9,1e-8);
  setprecisionfinitevolumesunsett(1e-5,1e-4);
  setprecisionfinitevolumeoneloopb(50,1e-7);
  setprecisionfinitevolumesunsetb(50,30,1e-5,1e-4);
  cout << LinfBE14;

  double m11 = pow(0.140,2);
  double m44 = pow(0.180,2);

  double mq1 = m11/2;
  double mq4 = m44/2.;
  vector<double> Bmq1 = {mq1};
  quarkmassnf qm1(Bmq1,f0,mu),qm2({mq1,mq4},f0,mu);
  cout << qm1 << qm2;
  quarkmass qm1r(mq1,mq1,f0,mu);
  quarkmassnf qm3({mq1,0.99*mq1},f0,mu);
  quarkmassnf qm4({mq1,1.01*mq1},f0,mu);
  double L = 2./sqrt(m11);
  int nf = 3;

  // 3 flavour case here versus usual ChPT with equal masses (older results)
  cout << "Comparing with three flavour ChPT equal mass if possible\n";

  cout << "\nmass SUN infinite volume\n";
  cout << mnfSUNp4L(nf,qm1,LinfBE14)<<' '<<mpi4Llo(qm1r,LiBE14)<<'\n';
  cout << mnfSUNp4R(nf,qm1)<<' '<<mpi4Rlo(qm1r)<<'\n';
  cout << mnfSUNp6L(nf,qm1,LinfBE14)<<' '<<mpi6Llo(qm1r,LiBE14)<<'\n';
  cout << mnfSUNp6R(nf,qm1)<<' '<<mpi6Rlo(qm1r)<<'\n';
  cout << "mass PQ SUN infinite volume\n";
  cout << mnfPQSUNp4L(nf,qm2,LinfBE14)<<' '<<mv1s1nf3p4L(qm2,LinfBE14)<<'\n';
  cout << mnfPQSUNp4R(nf,qm2)<<' '<<mv1s1nf3p4R(qm2)<<'\n';
  cout << mnfPQSUNp6L(nf,qm2,LinfBE14)<<' '<<mv1s1nf3p6L(qm2,LinfBE14)<<'\n';
  cout << mnfPQSUNp6R(nf,qm2)<<' '<<mv1s1nf3p6R(qm2)<<'\n';
  cout << "mass SUN finite volume, theta\n";
  cout << mnfSUNp4Vt(nf,qm1,L)<<' '<<mpi4loVt(qm1r,L)<<'\n';
  cout << mnfSUNp6LVt(nf,qm1,LinfBE14,L)<<' '<<mpi6LloVt(qm1r,LiBE14,L)<<'\n';
  cout << mnfSUNp6RVt(nf,qm1,L)<<' '<<mpi6RloVt(qm1r,L)<<'\n';
  cout << "mass SUN finite volume, Bessel\n";
  cout << mnfSUNp4Vb(nf,qm1,L)<<' '<<mpi4loVb(qm1r,L)<<'\n';
  cout << mnfSUNp6LVb(nf,qm1,LinfBE14,L)<<' '<<mpi6LloVb(qm1r,LiBE14,L)<<'\n';
  cout << mnfSUNp6RVb(nf,qm1,L)<<' '<<mpi6RloVb(qm1r,L)<<'\n';
  cout << "mass PQ SUN finite volume theta\n";
  cout << mnfPQSUNp4Vt(nf,qm2,L)<<' '<<mv1s1nf3p4Vt(qm2,L)<<'\n';
  cout << mnfPQSUNp6LVt(nf,qm2,LinfBE14,L)<<' '<<mv1s1nf3p6LVt(qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSUNp6RVt(nf,qm2,L)<<' '<<mv1s1nf3p6RVt(qm2,L)<<'\n';
  cout << "mass PQ SUN finite volume Bessel\n";
  cout << mnfPQSUNp4Vb(nf,qm2,L)<<' '<<mv1s1nf3p4Vb(qm2,L)<<'\n';
  cout << mnfPQSUNp6LVb(nf,qm2,LinfBE14,L)<<' '<<mv1s1nf3p6LVb(qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSUNp6RVb(nf,qm2,L)<<' '<<mv1s1nf3p6RVb(qm2,L)<<'\n';

  cout << "\ndecay SUN infinite volume\n";
  cout << fnfSUNp4L(nf,qm1,LinfBE14)<<' '<<fpi4Llo(qm1r,LiBE14)<<'\n';
  cout << fnfSUNp4R(nf,qm1)<<' '<<fpi4Rlo(qm1r)<<'\n';
  cout << fnfSUNp6L(nf,qm1,LinfBE14)<<' '<<fpi6Llo(qm1r,LiBE14)<<'\n';
  cout << fnfSUNp6R(nf,qm1)<<' '<<fpi6Rlo(qm1r)<<'\n';
  cout << "decay PQ SUN infinite volume\n";
  cout << fnfPQSUNp4L(nf,qm2,LinfBE14)<<' '<<fv1s1nf3p4L(qm2,LinfBE14)<<'\n';
  cout << fnfPQSUNp4R(nf,qm2)<<' '<<fv1s1nf3p4R(qm2)<<'\n';
  cout << fnfPQSUNp6L(nf,qm2,LinfBE14)<<' '<<fv1s1nf3p6L(qm2,LinfBE14)<<'\n';
  cout << fnfPQSUNp6R(nf,qm2)<<' '<<fv1s1nf3p6R(qm2)<<'\n';
  cout << "decay SUN finite volume theta\n";
  cout << fnfSUNp4Vt(nf,qm1,L)<<' '<<fpi4loVt(qm1r,L)<<'\n';
  cout << fnfSUNp6LVt(nf,qm1,LinfBE14,L)<<' '<<fpi6LloVt(qm1r,LiBE14,L)<<'\n';
  cout << fnfSUNp6RVt(nf,qm1,L)<<' '<<fpi6RloVt(qm1r,L)<<'\n';
  cout << "decay SUN finite volume Bessel\n";
  cout << fnfSUNp4Vb(nf,qm1,L)<<' '<<fpi4loVb(qm1r,L)<<'\n';
  cout << fnfSUNp6LVb(nf,qm1,LinfBE14,L)<<' '<<fpi6LloVb(qm1r,LiBE14,L)<<'\n';
  cout << fnfSUNp6RVb(nf,qm1,L)<<' '<<fpi6RloVb(qm1r,L)<<'\n';
  cout << "decay PQ SUN finite volume theta\n";
  cout << fnfPQSUNp4Vt(nf,qm2,L)<<' '<<fv1s1nf3p4Vt(qm2,L)<<'\n';
  cout << fnfPQSUNp6LVt(nf,qm2,LinfBE14,L)<<' '<<fv1s1nf3p6LVt(qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSUNp6RVt(nf,qm2,L)<<' '<<fv1s1nf3p6RVt(qm2,L)<<'\n';
  cout << "decay PQ SUN finite volume Bessel\n";
  cout << fnfPQSUNp4Vb(nf,qm2,L)<<' '<<fv1s1nf3p4Vb(qm2,L)<<'\n';
  cout << fnfPQSUNp6LVb(nf,qm2,LinfBE14,L)<<' '<<fv1s1nf3p6LVb(qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSUNp6RVb(nf,qm2,L)<<' '<<fv1s1nf3p6RVb(qm2,L)<<'\n';

  cout << "\nvev SUN infinite volume\n";
  cout << qnfSUNp4L(nf,qm1,LinfBE14)<<' '<<qqup4Llo(qm1r,LiBE14)<<'\n';
  cout << qnfSUNp4R(nf,qm1)<<' '<<qqup4Rlo(qm1r)<<'\n';
  cout << qnfSUNp6L(nf,qm1,LinfBE14)<<' '<<qqup6Llo(qm1r,LiBE14)<<'\n';
  cout << qnfSUNp6R(nf,qm1)<<' '<<qqup6Rlo(qm1r)<<'\n';
  cout << "vev PQ SUN infinite volume\n";
  cout << qnfPQSUNp4L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSUNp4R(nf,qm2)<<'\n';
  cout << qnfPQSUNp6L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSUNp6R(nf,qm2)<<'\n';
  cout << "vev SUN finite volume theta\n";
  cout << qnfSUNp4Vt(nf,qm1,L)<<' '<<qqup4loVt(qm1r,L)
       <<' '<<qqstrange4loVt(qm1r,L)<<'\n';
  cout << qnfSUNp6LVt(nf,qm1,LinfBE14,L)<<' '<<qqup6LloVt(qm1r,LiBE14,L)
       <<' '<<qqstrange6LloVt(qm1r,LiBE14,L)<<'\n';
  cout << qnfSUNp6RVt(nf,qm1,L)<<' '<<qqup6RloVt(qm1r,L)
       <<' '<<qqstrange6RloVt(qm1r,L)<<'\n';
  cout << "vev SUN finite volume Bessel\n";
  cout << qnfSUNp4Vb(nf,qm1,L)<<' '<<qqup4loVb(qm1r,L)
       <<' '<<qqstrange4loVb(qm1r,L)<<'\n';
  cout << qnfSUNp6LVb(nf,qm1,LinfBE14,L)<<' '<<qqup6LloVb(qm1r,LiBE14,L)
       <<' '<<qqstrange6LloVb(qm1r,LiBE14,L)<<'\n';
  cout << qnfSUNp6RVb(nf,qm1,L)<<' '<<qqup6RloVb(qm1r,L)
       <<' '<<qqstrange6RloVb(qm1r,L)<<'\n';
  cout << "vev PQ SUN finite volume theta\n";
  cout << qnfPQSUNp4Vt(nf,qm2,L)<<'\n';
  cout << qnfPQSUNp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSUNp6RVt(nf,qm2,L)<<'\n';
  cout << "vev PQ SUN finite volume Bessel\n";
  cout << qnfPQSUNp4Vb(nf,qm2,L)<<'\n';
  cout << qnfPQSUNp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSUNp6RVb(nf,qm2,L)<<'\n';

  cout << "\nmass SON infinite volume\n";
  cout << mnfSONp4L(nf,qm1,LinfBE14)<<'\n';
  cout << mnfSONp4R(nf,qm1)<<'\n';
  cout << mnfSONp6L(nf,qm1,LinfBE14)<<'\n';
  cout << mnfSONp6R(nf,qm1)<<'\n';
  cout << "mass PQ SON infinite volume\n";
  cout << mnfPQSONp4L(nf,qm2,LinfBE14)<<'\n';
  cout << mnfPQSONp4R(nf,qm2)<<'\n';
  cout << mnfPQSONp6L(nf,qm2,LinfBE14)<<'\n';
  cout << mnfPQSONp6R(nf,qm2)<<'\n';
  cout << "mass SON finite volume theta\n";
  cout << mnfSONp4Vt(nf,qm1,L)<<'\n';
  cout << mnfSONp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << mnfSONp6RVt(nf,qm1,L)<<'\n';
  cout << "mass SON finite volume Bessel\n";
  cout << mnfSONp4Vb(nf,qm1,L)<<'\n';
  cout << mnfSONp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << mnfSONp6RVb(nf,qm1,L)<<'\n';
  cout << "mass PQ SON finite volume theta\n";
  cout << mnfPQSONp4Vt(nf,qm2,L)<<'\n';
  cout << mnfPQSONp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSONp6RVt(nf,qm2,L)<<'\n';
  cout << "mass PQ SON finite volume Bessel\n";
  cout << mnfPQSONp4Vb(nf,qm2,L)<<'\n';
  cout << mnfPQSONp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSONp6RVb(nf,qm2,L)<<'\n';

  cout << "\ndecay SON infinite volume\n";
  cout << fnfSONp4L(nf,qm1,LinfBE14)<<'\n';
  cout << fnfSONp4R(nf,qm1)<<'\n';
  cout << fnfSONp6L(nf,qm1,LinfBE14)<<'\n';
  cout << fnfSONp6R(nf,qm1)<<'\n';
  cout << "decay PQ SON infinite volume\n";
  cout << fnfPQSONp4L(nf,qm2,LinfBE14)<<'\n';
  cout << fnfPQSONp4R(nf,qm2)<<'\n';
  cout << fnfPQSONp6L(nf,qm2,LinfBE14)<<'\n';
  cout << fnfPQSONp6R(nf,qm2)<<'\n';
  cout << "decay SON finite volume theta\n";
  cout << fnfSONp4Vt(nf,qm1,L)<<'\n';
  cout << fnfSONp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << fnfSONp6RVt(nf,qm1,L)<<'\n';
  cout << "decay SON finite volume Bessel\n";
  cout << fnfSONp4Vb(nf,qm1,L)<<'\n';
  cout << fnfSONp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << fnfSONp6RVb(nf,qm1,L)<<'\n';
  cout << "decay PQ SON finite volume theta\n";
  cout << fnfPQSONp4Vt(nf,qm2,L)<<'\n';
  cout << fnfPQSONp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSONp6RVt(nf,qm2,L)<<'\n';
   cout << "decay PQ SON finite volume Bessel\n";
  cout << fnfPQSONp4Vb(nf,qm2,L)<<'\n';
  cout << fnfPQSONp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSONp6RVb(nf,qm2,L)<<'\n';
 
  cout << "\nvev SON infinite volume\n";
  cout << qnfSONp4L(nf,qm1,LinfBE14)<<'\n';
  cout << qnfSONp4R(nf,qm1)<<'\n';
  cout << qnfSONp6L(nf,qm1,LinfBE14)<<'\n';
  cout << qnfSONp6R(nf,qm1)<<'\n';
  cout << "vev PQ SON infinite volume\n";
  cout << qnfPQSONp4L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSONp4R(nf,qm2)<<'\n';
  cout << qnfPQSONp6L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSONp6R(nf,qm2)<<'\n';
  cout << "vev SON finite volume theta\n";
  cout << qnfSONp4Vt(nf,qm1,L)<<'\n';
  cout << qnfSONp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << qnfSONp6RVt(nf,qm1,L)<<'\n';
  cout << "vev SON finite volume Bessel\n";
  cout << qnfSONp4Vb(nf,qm1,L)<<'\n';
  cout << qnfSONp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << qnfSONp6RVb(nf,qm1,L)<<'\n';
  cout << "vev PQ SON finite volume theta\n";
  cout << qnfPQSONp4Vt(nf,qm2,L)<<'\n';
  cout << qnfPQSONp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSONp6RVt(nf,qm2,L)<<'\n';
  cout << "vev PQ SON finite volume Bessel\n";
  cout << qnfPQSONp4Vb(nf,qm2,L)<<'\n';
  cout << qnfPQSONp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSONp6RVb(nf,qm2,L)<<'\n';

  cout << "\nmass SPN infinite volume\n";
  cout << mnfSPNp4L(nf,qm1,LinfBE14)<<'\n';
  cout << mnfSPNp4R(nf,qm1)<<'\n';
  cout << mnfSPNp6L(nf,qm1,LinfBE14)<<'\n';
  cout << mnfSPNp6R(nf,qm1)<<'\n';
  cout << "mass PQ SPN infinite volume\n";
  cout << mnfPQSPNp4L(nf,qm2,LinfBE14)<<'\n';
  cout << mnfPQSPNp4R(nf,qm2)<<'\n';
  cout << mnfPQSPNp6L(nf,qm2,LinfBE14)<<'\n';
  cout << mnfPQSPNp6R(nf,qm2)<<'\n';
  cout << "mass SPN finite volume theta\n";
  cout << mnfSPNp4Vt(nf,qm1,L)<<'\n';
  cout << mnfSPNp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << mnfSPNp6RVt(nf,qm1,L)<<'\n';
  cout << "mass SPN finite volume Bessel\n";
  cout << mnfSPNp4Vb(nf,qm1,L)<<'\n';
  cout << mnfSPNp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << mnfSPNp6RVb(nf,qm1,L)<<'\n';
  cout << "mass PQ SPN finite volume theta\n";
  cout << mnfPQSPNp4Vt(nf,qm2,L)<<'\n';
  cout << mnfPQSPNp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSPNp6RVt(nf,qm2,L)<<'\n';
  cout << "mass PQ SPN finite volume Bessel\n";
  cout << mnfPQSPNp4Vb(nf,qm2,L)<<'\n';
  cout << mnfPQSPNp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << mnfPQSPNp6RVb(nf,qm2,L)<<'\n';

  cout << "\ndecay SPN infinite volume\n";
  cout << fnfSPNp4L(nf,qm1,LinfBE14)<<'\n';
  cout << fnfSPNp4R(nf,qm1)<<'\n';
  cout << fnfSPNp6L(nf,qm1,LinfBE14)<<'\n';
  cout << fnfSPNp6R(nf,qm1)<<'\n';
  cout << "decay PQ SPN infinite volume\n";
  cout << fnfPQSPNp4L(nf,qm2,LinfBE14)<<'\n';
  cout << fnfPQSPNp4R(nf,qm2)<<'\n';
  cout << fnfPQSPNp6L(nf,qm2,LinfBE14)<<'\n';
  cout << fnfPQSPNp6R(nf,qm2)<<'\n';
  cout << "decay SPN finite volume theta\n";
  cout << fnfSPNp4Vt(nf,qm1,L)<<'\n';
  cout << fnfSPNp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << fnfSPNp6RVt(nf,qm1,L)<<'\n';
  cout << "decay SPN finite volume Bessel\n";
  cout << fnfSPNp4Vb(nf,qm1,L)<<'\n';
  cout << fnfSPNp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << fnfSPNp6RVb(nf,qm1,L)<<'\n';
  cout << "decay PQ SPN finite volume theta\n";
  cout << fnfPQSPNp4Vt(nf,qm2,L)<<'\n';
  cout << fnfPQSPNp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSPNp6RVt(nf,qm2,L)<<'\n';
  cout << "decay PQ SPN finite volume Bessel\n";
  cout << fnfPQSPNp4Vb(nf,qm2,L)<<'\n';
  cout << fnfPQSPNp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << fnfPQSPNp6RVb(nf,qm2,L)<<'\n';

  cout << "\nvev SPN infinite volume\n";
  cout << qnfSPNp4L(nf,qm1,LinfBE14)<<'\n';
  cout << qnfSPNp4R(nf,qm1)<<'\n';
  cout << qnfSPNp6L(nf,qm1,LinfBE14)<<'\n';
  cout << qnfSPNp6R(nf,qm1)<<'\n';
  cout << "vev PQ SPN infinite volume\n";
  cout << qnfPQSPNp4L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSPNp4R(nf,qm2)<<'\n';
  cout << qnfPQSPNp6L(nf,qm2,LinfBE14)<<'\n';
  cout << qnfPQSPNp6R(nf,qm2)<<'\n';
  cout << "vev SPN finite volume theta\n";
  cout << qnfSPNp4Vt(nf,qm1,L)<<'\n';
  cout << qnfSPNp6LVt(nf,qm1,LinfBE14,L)<<'\n';
  cout << qnfSPNp6RVt(nf,qm1,L)<<'\n';
  cout << "vev SPN finite volume Bessel\n";
  cout << qnfSPNp4Vb(nf,qm1,L)<<'\n';
  cout << qnfSPNp6LVb(nf,qm1,LinfBE14,L)<<'\n';
  cout << qnfSPNp6RVb(nf,qm1,L)<<'\n';
  cout << "vev PQ SPN finite volume theta\n";
  cout << qnfPQSPNp4Vt(nf,qm2,L)<<'\n';
  cout << qnfPQSPNp6LVt(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSPNp6RVt(nf,qm2,L)<<'\n';
  cout << "vev PQ SPN finite volume Bessel\n";
  cout << qnfPQSPNp4Vb(nf,qm2,L)<<'\n';
  cout << qnfPQSPNp6LVb(nf,qm2,LinfBE14,L)<<'\n';
  cout << qnfPQSPNp6RVb(nf,qm2,L)<<'\n';


  return 0;
}
