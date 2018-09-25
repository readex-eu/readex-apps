test run on 1 node:
export MERIC_NUM_THREADS=12
srun -n 2 ./blasbench --matrixFile bmwcra-1-mtx --mpi_ata 100,256 --dgemv 100,1000 --dgemm 100,1000 --tan 100 --io --io_threads 12 --sparse0 100 --sparse1 100 --sparse2 100 --sparse3 100 --mpi_ring 100,256 --runtime 0.1

matrix download https://www.cise.ufl.edu/research/sparse/matrices/GHS_psdef/bmwcra_1.html in mtx format

arguments:

--matrixFile [matrix]    path + name of the input sparse matrix
--io                     Parallel OMP I/O, requires --matrixFile
--io_threads [nThreads]  Number of I/O threads
--sparse0 [reps]         SpMV multiplication in CSR format, requires --matrixFile and --io
--sparse1 [reps]         SpMV multiplication in IJV format, requires --matrixFile and --io
--sparse2 [reps]         SpMM multiplication, requires --matrixFile and --io
--sparse3 [reps]         SpMM addition, requires --matrixFile and --io
--dgemm [reps],[size]    C = A*B, matrices' size is size*size*sizeof(double)
--dgemv [reps],[size]    matrix-vector product M*Vi=Vo, matrix size size*size*sizeof(double), vector size*sizeof(double)
--mpi_ata [reps],[size]  MPI_Alltoall, a single message size is size*size*sizeof(MPI_DOUBLE)
--mpi_ring [reps],[size] MPI_Isend+MPI_Irecv in a ring, all sending and receiving at the same time
                         a single message size is size*size*MPIranks
--runtime [multiplier]   Increses or decreases the overall runtime (e.q reps), 0.0 < multiplier < inf
