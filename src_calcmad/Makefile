FC=ifort 
DOPT=-g  -traceback
NOPT=-O3 -mp1 -w -prec_div -pc80 -nopad -ip -fpe0 
#Common
OPT=-r8 

#Switch to debug
#FFLAGS=$(OPT) $(DOPT)
FFLAGS=$(OPT) $(NOPT)

TARGET = ../calcmad

calcmad: calcmadModule.f90 structMod.f90 EwaldModule.f90 DirectMadModule.f90 calcmad.f90 
	$(FC) -o $(TARGET) $(FFLAGS) $^ 

clean:
	rm -f *.o *.mod $(TARGET) 

