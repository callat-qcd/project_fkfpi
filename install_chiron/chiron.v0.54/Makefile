# this Makefile produces the numerical library jbnumlib
# and the library libchiron.a
# it also allows you to produce a number of testing/example programs
# it is part of the CHIRON ChPT at two loops program collection
#
# Copyright (C) 2014-2015 Johan Bijnens, v1.02
# CHIRON is licenced under the GNU GPL version 2 or later,
# see COPYING for details.
# Please respect the Guidelines, see GUIDELINES for details.

# in both testing and library
CC = g++ -O3  -I./include
# added for the compilation of the libraries
CFLAGS = -fPIC -Wall -Wextra -Wconversion
#added for the testing files
CFLAGSTEST = -std=c++0x -Wall -Wextra -Wconversion

OBJECTSJBNUMLIB = jbdgauss.o jbdcauch.o jbwgauss.o \
  jbdgauss2.o jbdquad15.o jbdquad21.o \
  jbdcauch2.o jbdsing15.o jbdsing21.o \
  jbwgauss2.o jbwquad15.o jbwquad21.o \
  jbdadmul.o \
  jbdli2.o jbdbesik.o jbdtheta30.o jbdtheta30m1.o jbdtheta32.o jbdtheta34.o \
  jbdtheta2d0.o jbdtheta2d0m1.o jbdtheta2d02.o \
  jbdtheta3.o jbderivutheta3.o jbderiv2utheta3.o

OBJECTSCHIRON = inputs.o Li.o Ci.o inputsnf.o Linf.o Ki.o \
       massdecayvev.o getfpimeta.o \
       massdecayvevlo.o massdecayvevPQ.o \
       massdecayvevnf.o massdecayvevnfPQ.o \
       oneloopintegrals.o sunsetintegrals.o finitevolumeoneloopintegrals.o \
       finitevolumesunsetintegrals.o quenchedsunsetintegrals.o

OBJECTSCHIRONTHETA =  massdecayvevVt.o massdecayvevloVt.o massdecayvevPQVt.o \
       massdecayvevnfVt.o massdecayvevnfPQVt.o
OBJECTSCHIRONBESSEL =  massdecayvevVb.o massdecayvevloVb.o massdecayvevPQVb.o \
       massdecayvevnfVb.o massdecayvevnfPQVb.o
OBJECTSCHIRON2 = $(OBJECTSCHIRONTHETA) $(OBJECTSCHIRONBESSEL)

INCLUDECHIRON =  include/inputs.h include/Li.h include/Ci.h\
       include/inputsnf.h include/Linf.h include/Ki.h \
       include/massdecayvev.h \
       include/massdecayvevlo.h \
       include/massdecayvevPQ.h\
       include/massdecayvevnf.h \
       include/massdecayvevnfPQ.h \
       include/oneloopintegrals.h \
       include/sunsetintegrals.h\
       include/finitevolumeoneloopintegrals.h\
       include/finitevolumesunsetintegrals.h\
       include/massdecayvevV.h \
       include/massdecayvevloV.h \
       include/massdecayvevPQV.h\
       include/quenchedsunsetintegrals.h

TESTCHIRON = testintegralsreal testintegralsrealsingular \
       testintegralscomplex \
       testinputs testinputsnf testLi testCi testmassdecayvev \
       testgetfpimeta testoneloopintegrals testsunsetintegrals \
       testfinitevolumeoneloopintegrals testfinitevolumesunsetintegrals \
       testmassdecayvevV testquenchedsunsetintegrals testLinf testKi \
       testmassdecayvevlo testmassdecayvevPQ testmassdecayvevPQV \
       testmassdecayvevloV testmassdecayvevnf

all: libchiron.a libjbnumlib.a

libchiron.a: $(OBJECTSCHIRON) $(OBJECTSCHIRON2)
	ar r libchiron.a $(OBJECTSCHIRON) $(OBJECTSCHIRON2) \
        ; cp libchiron.a lib

libjbnumlib.a: $(OBJECTSJBNUMLIB) include/jbnumlib.h
	ar r libjbnumlib.a $(OBJECTSJBNUMLIB); cp libjbnumlib.a lib

massdecayvevVt.o: src/massdecayvevV.cc include/massdecayvevV.h
	$(CC) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevVt.o src/massdecayvevV.cc
massdecayvevloVt.o: src/massdecayvevloV.cc include/massdecayvevloV.h
	$(CC) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevloVt.o src/massdecayvevloV.cc
massdecayvevPQVt.o: src/massdecayvevPQV.cc include/massdecayvevPQV.h
	$(CC) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevPQVt.o src/massdecayvevPQV.cc

massdecayvevnfVt.o: src/massdecayvevnfV.cc include/massdecayvevnfV.h
	$(CC) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevnfVt.o src/massdecayvevnfV.cc

massdecayvevnfPQVt.o: src/massdecayvevnfPQV.cc include/massdecayvevnfPQV.h
	$(CC) -c $(CFLAGS) -D CHIRONTHETA -o massdecayvevnfPQVt.o src/massdecayvevnfPQV.cc

massdecayvevVb.o: src/massdecayvevV.cc include/massdecayvevV.h
	$(CC) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevVb.o src/massdecayvevV.cc

massdecayvevloVb.o: src/massdecayvevloV.cc include/massdecayvevloV.h
	$(CC) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevloVb.o src/massdecayvevloV.cc

massdecayvevPQVb.o: src/massdecayvevPQV.cc include/massdecayvevPQV.h
	$(CC) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevPQVb.o src/massdecayvevPQV.cc

massdecayvevnfVb.o: src/massdecayvevnfV.cc include/massdecayvevnfV.h
	$(CC) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevnfVb.o src/massdecayvevnfV.cc

massdecayvevnfPQVb.o: src/massdecayvevnfPQV.cc include/massdecayvevnfPQV.h
	$(CC) -c $(CFLAGS) -D CHIRONBESSEL -o massdecayvevnfPQVb.o src/massdecayvevnfPQV.cc

$(OBJECTSJBNUMLIB): %.o: src/%.cc
	$(CC) -c $(CFLAGS) $< -o $@

$(OBJECTSCHIRON): %.o: src/%.cc include/%.h
	$(CC) -c $(CFLAGS) $< -o $@

# making testing programs, output is always a.out

$(TESTCHIRON): %: test/%.cc  libjbnumlib.a libchiron.a
	$(CC) $(CFLAGSTEST)  -o a.out $<  -lchiron -ljbnumlib -L./lib


# some installing and cleaning up bits
.PHONY: clean installjbnumlib installchiron doc

installchiron: libjbnumlib.a libchiron.a
	cp libjbnumlib.a libchiron.a ~/lib ; cp include/jbnumlib.h $(INCLUDECHIRON) ~/include

installjbnumlib: libjbnumlib.a
	cp libjbnumlib.a ~/lib ; cp include/jbnumlib.h ~/include

clean:
	rm *.o *.a a.out lib/*.a \
          temp.dat test.dat *.log *.out *.toc *.aux \
          *.idx *.ilg *.ind *.pdf

doc:
	pdflatex doc/manual.tex ;pdflatex doc/manual.tex ; makeindex manual ;\
          pdflatex doc/manual.tex ; makeindex manual ; \
          pdflatex doc/manual.tex ; mv manual.pdf doc
