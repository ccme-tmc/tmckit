#include ../make.inc 
F90=ifort 
FFLAGS=-O2

OBJS= cdft_init.o dlatgen.o dm_init.o ioutils.o lapack.o latdist.o lsqfit.o randspin.o 
BINDIR = ../

all: cdft_init dm_init 


cdft_init: cdft_init.f90
	$(F90) $(FFLAGS) $(LDFLAGS) -o $(BINDIR)cdft_init.x cdft_init.f90

dm_init: dm_init.f90
	$(F90) $(FFLAGS) $(LDFLAGS) -o $(BINDIR)dm_init.x dm_init.f90

clean:
	rm -f *.o 

