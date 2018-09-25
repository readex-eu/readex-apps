
CXX = ${READEX_CXX}
CXX_SCOREP		= scorep --online-access --user --mpp=mpi --thread=omp --noopenmp $(FILTER_GNU) $(CXX)  
CXX_SCOREP_MANUAL	= scorep --online-access --user --nocompiler --mpp=mpi --thread=omp --noopenmp $(CXX)  
SCOREP_SWITCH		= -DUSE_SCOREP
SCOREP_SWITCH_MANUAL	= -DUSE_SCOREP_MANUAL
 
CXX_FLAGS   		= ${READEX_CXX_FLAGS} -O3 -g -std=c++11 ${READEX_OMP_FLAG}

INC_FOLDER_MKL 		= -I$(MKLROOT)/include/
LIB_FOLDER_MKL 		= -L$(MKLROOT)/lib/intel64/
LIBS        		= ${READEX_OMP_LIB} -lmpi
LIBS_MKL        	= -lmkl_intel_lp64 -lmkl_core ${READEX_MKL_THREAD_LIB}

EXEC = blasbench

all:
	$(CXX) $(CXX_FLAGS) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)

ptf:
	$(CXX_SCOREP) $(CXX_FLAGS) $(FILTER_INTEL) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)

ptf_manual:
	$(CXX_SCOREP_MANUAL) $(CXX_FLAGS) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH_MANUAL) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)

autofilter:
	scorep $(CXX) $(CXX_FLAGS) $(FILTER_INTEL) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)

scorepmanual: 
	$(CXX_SCOREP) $(CXX_FLAGS) $(FILTER_INTEL) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)
	
meric:
	$(CXX) $(CXX_FLAGS) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) -DUSE_MERIC src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(MERIC_LIB_PATH) $(LIBS) $(LIBS_MKL) $(MERIC_LIBS)

dyndetect:
	$(CXX_SCOREP) $(CXX_FLAGS) $(FILTER_INTEL) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)
	
dyndetect_manual:
	$(CXX_SCOREP_MANUAL) $(CXX_FLAGS) $(INC_FOLDER_MKL) $(MERIC_INC_PATH) $(SCOREP_SWITCH_MANUAL) src/main.cpp src/**/*.cpp -o $(EXEC) $(LIB_FOLDER_MKL) $(LIBS) $(LIBS_MKL)

clean:
	rm -f $(EXEC)
	rm -f *.o
