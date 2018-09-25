SHELL = /bin/sh
.SUFFIXES: .cc .o

LULESH_EXEC = lulesh2.0

#SERCXX = scorep --online-access --user --nocompiler --noopenmp icc -DUSE_MPI=0 
#MPICXX = mpicc -DUSE_MPI=1 -qopenmp
#MPICXX = scorep --online-access --user --nocompiler --thread=none  mpicc -DUSE_MPI=1 -qopenmp
#MPICXX = scorep --online-access --user --nocompiler mpicc -DUSE_MPI=1 -qopenmp
#MPICXX = scorep --online-access --nocompiler --thread=none --user  mpicc -DUSE_MPI=1 -qopenmp
#MPICXX = scorep --online-access --user --nocompiler --thread=none  mpic++  -I$(MPI_INC) -L$(MPI_LIB) -lmpi  -lmpi_cxx -DUSE_MPI=1 -fopenmp
#MPICXX = scorep --online-access  mpic++ -DUSE_MPI=1 -fopenmp

#MPICXX = ${SCOREP_PREP} mpic++ ${READEX_CXX_FLAGS} -I$(MPI_INC) -L$(MPI_LIB) -lmpi  -lmpi_cxx -DUSE_MPI=1 -fopenmp ${DSWITCH} 
MPICXX = ${SCOREP_PREP} ${READEX_CXX} ${READEX_CXX_FLAGS} ${FILTER_INTEL} ${MERIC_INC_PATH} -DUSE_MPI=1 -fopenmp $(MERIC_LIB_PATH) $(MERIC_LIBS) ${DSWITCH}
CXX = $(MPICXX)

SOURCES2.0 = \
         lulesh.cc \
         lulesh-comm.cc \
         lulesh-viz.cc \
         lulesh-util.cc \
         lulesh-init.cc
OBJECTS2.0 = $(SOURCES2.0:.cc=.o)

#Default build suggestions with OpenMP for g++ CXXFLAGS = -g -O3 -fopenmp -I. $(DEVEL_INC) -Wall LDFLAGS = -g -O3 -fopenmp

.cc.o: lulesh.h
	@echo "Building $<"
	$(CXX) -c $(CXXFLAGS) -o $@  $<

all: $(LULESH_EXEC)

lulesh2.0: $(OBJECTS2.0)
	@echo "Linking"
	$(CXX) $(OBJECTS2.0) $(LDFLAGS) -lm -o $@

clean:
	/bin/rm -f *.o *~ $(OBJECTS) $(LULESH_EXEC)
	/bin/rm -rf *.dSYM

tar: clean
	cd .. ; tar cvf lulesh-2.0.tar LULESH-2.0 ; mv lulesh-2.0.tar LULESH-2.0

run: lulesh2.0
	rm job.out advice* *.psc; sbatch  runLuleshpsc.scp
