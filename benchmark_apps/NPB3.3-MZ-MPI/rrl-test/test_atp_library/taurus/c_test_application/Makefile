COMPILER = scorep --online-access --user --nocompiler --thread=omp --openmp --mpp=none  g++ -std=c++11 -pthread

all:
	$(COMPILER) tuning_test.c -o tuning_function -I/home/mian/install/readex-rrl/include -L/home/mian/install/readex-rrl/lib -lscorep_substrate_rrl

clean:
	rm -rf tuning_function
