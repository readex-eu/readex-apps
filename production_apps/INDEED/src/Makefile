include make.cfg

SUBDIRS	= Utils \
	  Dmsys Dafeld Driver Element \
	  Indini Input Machine License \
	  Output Rezoning Solver Trimming

OBJECT= time_measure.o main.o module.o

LNKLIB= $(LIBNAM) $(DMSLIB) $(LICLIB) $(SERLIB)

OBJEXE= ../bin/Indeed_$(EXEMOD)$(BSIZE).exe

Indeed: $(OBJECT)
	$(MAKE) all
	$(LD) $(LDFLAGS) -o $(OBJEXE) $(OBJECT) $(LNKLIB) $(SOLVLIB) $(MPILIB) $(SCOREPLIB)

linker: $(OBJECT)
	$(LD) $(LDFLAGS) -o $(OBJEXE) $(OBJECT) $(LNKLIB) $(SOLVLIB) $(MPILIB) $(SCOREPLIB)

clobber::
	$(MAKE) clean
	rm -f *.o *.mod *.optrpt
#	rm -f $(OBJEXE)

include make.sub

.f.o:
	$(FC) $(DEFINE) $(MPATH) $(FFLAGS) -c $< -o $*.o
