FC =     ifort
FOPT     = -O3 -FR 
LDFLAGS  =      

EXEC   = ../cif2struct

OBJS =	module.o spacegroup.o getlattype.o getsgnum.o getsgname.o  getzz.o \
        structgen.o scan_cif.o ciftbx.o hash_funcs.o test_sgname.o scan_in.o

$(EXEC):	$(OBJS) 
	$(FC) -o $(EXEC) $(FFLAGS) $(LDFLAGS) $(OBJS)

clean:	  
	-rm -f $(OBJS)  *.mod
	@echo clean

.f.o:
	$(FC) -c $(FOPT) $<
.f90.o:
	$(FC) $(FOPT) -c $<
