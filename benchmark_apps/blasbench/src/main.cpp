
#include <cmath>

#include <omp.h>
#include "solver/SparseSolverMKL.h"
#include <readex.h>
#include <mpi.h>
#include <getopt.h>

bool IO_flag = false;
bool SPARSE0_flag = false;
int  SPARSE0_reps = 1;//2500
bool SPARSE1_flag = false;
int  SPARSE1_reps = 1;//400
bool SPARSE2_flag = false;
int  SPARSE2_reps = 1;//10
bool SPARSE3_flag = false;
int  SPARSE3_reps = 1;//50
bool DENSE0_flag = false;
int DENSE0_reps = 1;//400
bool DENSE1_flag = false;
int DENSE1_reps = 1;//62500
bool MPI_ata_flag = false;
int MPI_ata_reps = 1;//250
bool MPI_ring_flag = false;
int MPI_ring_reps = 1;//500
bool TAN_flag = false;
int TAN_reps = 1;//300000
bool PHASE_flag = false;
int PHASE_reps = 1;
double RUNTIME = 1.0;
int IO_THREADS = 24;
std::string matrixName;
int matrixSizeDGEMM = 1000;
int matrixSizeDGEMV = 1000;
int mpiATASize = 1000;
int mpiRingSize = 1000;

READEX_REGION_DEFINE(Main);
READEX_REGION_DEFINE(Parallel_IO);
READEX_REGION_DEFINE(SPARSE0);
READEX_REGION_DEFINE(SPARSE1);
READEX_REGION_DEFINE(SPARSE2);
READEX_REGION_DEFINE(SPARSE3);
READEX_REGION_DEFINE(MPI_ata);
READEX_REGION_DEFINE(MPI_ring);
READEX_REGION_DEFINE(DGEMM_);
READEX_REGION_DEFINE(DGEMV_);
READEX_REGION_DEFINE(TAN);


void ParseCommandLine(int argc, char** argv)
{
    char c;
    int longIndex;
    std::string delimeter = ",";

    const struct option longOpts[] = {
        { "io", no_argument, NULL, 0 },
        { "sparse0", required_argument, NULL, 0 },
        { "sparse1", required_argument, NULL, 0 },
        { "sparse2", required_argument, NULL, 0 },
        { "sparse3", required_argument, NULL, 0 },
        { "dgemm", required_argument, NULL, 0 },
        { "dgemv", required_argument, NULL, 0 },
        { "mpi_ata", required_argument, NULL, 0 },
        { "mpi_ring", required_argument, NULL, 0 },
        { "tan", required_argument, NULL, 0 },
        { "phase", required_argument, NULL, 0 },
        { "runtime", required_argument, NULL, 0 },
        { "io_threads", required_argument, NULL, 0 },
        { "matrixFile", required_argument, NULL, 0 },
        { NULL, no_argument, NULL, 0 }
    };

    // Short parameters //
    while ((c = getopt_long(argc, argv, "", longOpts, &longIndex)) != -1) {
        switch (c) {
        // long option without a short arg
        case 0: {
            if (strcmp("io", longOpts[longIndex].name) == 0) {
                IO_flag = true;
            }
            else if (strcmp("sparse0", longOpts[longIndex].name) == 0) {
                SPARSE0_flag = true;
                SPARSE0_reps = atol(optarg);
            }
            else if (strcmp("sparse1", longOpts[longIndex].name) == 0) {
                SPARSE1_flag = true;
                SPARSE1_reps = atol(optarg);
            }
            else if (strcmp("sparse2", longOpts[longIndex].name) == 0) {
                SPARSE2_flag = true;
                SPARSE2_reps = atol(optarg);
            }
            else if (strcmp("sparse3", longOpts[longIndex].name) == 0) {
                SPARSE3_flag = true;
                SPARSE3_reps = atol(optarg);
            }
            else if (strcmp("dgemm", longOpts[longIndex].name) == 0) {
                DENSE0_flag = true;
                std::string s(optarg);
                DENSE0_reps = atol(s.substr(0, s.find(delimeter)).c_str());
                matrixSizeDGEMM = atol(s.substr(s.find(delimeter)+1).c_str());
            }
            else if (strcmp("dgemv", longOpts[longIndex].name) == 0) {
                DENSE1_flag = true;
                std::string s(optarg);
                DENSE1_reps = atol(s.substr(0, s.find(delimeter)).c_str());
                matrixSizeDGEMV = atol(s.substr(s.find(delimeter)+1).c_str());
            }
            else if (strcmp("mpi_ata", longOpts[longIndex].name) == 0) {
                MPI_ata_flag = true;
                std::string s(optarg);
                MPI_ata_reps = atol(s.substr(0, s.find(delimeter)).c_str());
                mpiATASize = atol(s.substr(s.find(delimeter)+1).c_str());
            }
            else if (strcmp("mpi_ring", longOpts[longIndex].name) == 0) {
                MPI_ring_flag = true;
                std::string s(optarg);
                MPI_ring_reps = atol(s.substr(0, s.find(delimeter)).c_str());
                mpiRingSize = atol(s.substr(s.find(delimeter)+1).c_str());
            }
            else if (strcmp("tan", longOpts[longIndex].name) == 0) {
                TAN_flag = true;
                TAN_reps = atol(optarg);
            }
            else if (strcmp("phase", longOpts[longIndex].name) == 0) {
                PHASE_flag = true;
                PHASE_reps = atol(optarg);
            }
            else if (strcmp("runtime", longOpts[longIndex].name) == 0) {
                RUNTIME = atof(optarg);
            }
            else if (strcmp("io_threads", longOpts[longIndex].name) == 0) {
                IO_THREADS = atol(optarg);
            }
            else if (strcmp("matrixFile", longOpts[longIndex].name) == 0) {
                matrixName = optarg;
            }
            break;
        }
        default: {
            fprintf(stderr, "Wrong argument!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        }
    }

    if (RUNTIME <= 0.0 || matrixSizeDGEMV <= 0 || matrixSizeDGEMM <= 0 || IO_THREADS <= 0)
    {
      std::cerr << "Invalid input parameter!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if ((SPARSE0_flag || SPARSE1_flag || SPARSE2_flag || SPARSE3_flag) && (matrixName.empty() || !IO_flag))
    {
      std::cerr << "Sparse routines need a name of the input sparse matrix!" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (IO_flag && matrixName.empty())
    {
      std::cerr << "Need an input sparse matrix file to measure I/O" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

} // end of ParseCommandLine
//------------------------------------------------------------------------------

std::vector<std::streampos> getDistribution(std::ifstream& file, size_t parts)
{
    std::streampos start = file.tellg();
    file.seekg(0, file.end);
    std::streampos end = file.tellg();

    if (start > end) {
        std::cerr << "Distribution of interval <" << start << "," << end << "> is not possible.";
    }
    size_t size = end - start;
    std::vector<std::streampos> distribution(parts + 1, 0);
    size_t chunkSize = std::ceil(size / (double)parts);

    for (size_t t = 1; t < parts; t++) {
        distribution[t] = t * chunkSize + start;
        if (distribution[t] > end) {
            distribution[t] = end;
        }
    }
    distribution[0] = start;
    distribution[parts] = end;

    return distribution;
}

/*
void ReadMatrix(espreso::SparseMatrix matrix, char*)
{
}

int areWeDone(double start, int rank)
{
    int done = 0;

    if (!rank && (MPI_Wtime() - start > RUNTIME)) {
        done = 1;
    }

    MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);

    return done;
}
*/


void Parallel_IO_readex(espreso::SparseMatrix & matrix)
{        
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    READEX_REGION_START(Parallel_IO, "Parallel_IO", SCOREP_USER_REGION_TYPE_COMMON);
      
    matrix.type = 'G';

    std::ifstream file(matrixName.c_str());
    if (!file.is_open()) {
        std::cerr << "Cannot open the file '" << matrixName << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    else
    {
        while (file.good()) {
            std::string line;
            getline(file, line);
            if (line.size() && line[0] != '%') {
                std::stringstream ss(line);
                ss >> matrix.rows >> matrix.cols >> matrix.nnz;
                break;
            }
        }

        std::vector<std::streampos> distribution = getDistribution(file, IO_THREADS);
        std::vector<espreso::SparseMatrix> matrices(distribution.size() - 1);

        double start = omp_get_wtime();

        #pragma omp parallel for
        for (size_t t = 0; t < distribution.size() - 1; t++) {
            std::ifstream chunk(matrixName.c_str());
            chunk.seekg(distribution[t]);
            std::string line;
            int i, j;
            double v = 1;

            if (t) { // we cannot know if we have the full line -> skip to the next
                getline(chunk, line);
            }

            // guess of the chunk size in order to avoid reallocation
            matrices[t].I_row_indices.reserve(2 * matrix.nnz / distribution.size());
            matrices[t].J_col_indices.reserve(2 * matrix.nnz / distribution.size());
            matrices[t].V_values.reserve(2 * matrix.nnz / distribution.size());
            while (chunk.tellg() < distribution[t + 1]) {
                getline(chunk, line);
                std::stringstream values(line);
                values >> i >> j >> v;
                matrices[t].I_row_indices.push_back(i);
                matrices[t].J_col_indices.push_back(j);
                matrices[t].V_values.push_back(v);
            }

            if (t + 2 < distribution.size() && chunk.tellg() == distribution[t + 1]) { // read the skipped line
                getline(chunk, line);
                std::stringstream values(line);
                values >> i >> j >> v;
                matrices[t].I_row_indices.push_back(i);
                matrices[t].J_col_indices.push_back(j);
                matrices[t].V_values.push_back(v);
            }
        }

        double mid = omp_get_wtime();

        for (size_t t = 0; t < distribution.size() - 1; t++) {
            matrix.I_row_indices.insert(matrix.I_row_indices.end(), matrices[t].I_row_indices.begin(), matrices[t].I_row_indices.end());
            matrix.J_col_indices.insert(matrix.J_col_indices.end(), matrices[t].J_col_indices.begin(), matrices[t].J_col_indices.end());
            matrix.V_values.insert(matrix.V_values.end(), matrices[t].V_values.begin(), matrices[t].V_values.end());
        }

        if (matrix.nnz != matrix.V_values.size()) {
            std::cerr << "Broken parallelization of reading matrix. FIX IT :)\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        }

        double end = omp_get_wtime();

        #pragma omp parallel
        {
            #pragma omp master
            {
                if (!rank)
                    std::cout << "Number of running OMP threads (from OMP runtime)  : " << omp_get_num_threads() << "\n";
            }
        }

        if (!rank)
        {
            std::cout << "Parallel I/O - Number of OMP threads - chunks     : " << distribution.size() - 1 << "\n";
            std::cout << "Parallel I/O - loading matrix from file           : " << mid - start << "\n";
            std::cout << "Parallel I/O - Merging partial arrays into final  : " << end - mid << "\n";
        }       
    }

    READEX_REGION_STOP(Parallel_IO);
}

long long int TAN_wrapper(  int repsLimit ,double& min )
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    long long int reps = 0;
    
    READEX_REGION_START(TAN, "TAN", SCOREP_USER_REGION_TYPE_COMMON);
    
    #pragma omp parallel for reduction(+:reps) reduction(min:min)
    for (int i = 0; i < repsLimit; ++i)
    {
       #pragma omp simd reduction(+:reps) reduction(min:min)
       for (int j = 0; j < 15000; ++j)
       {
         double result = std::tan(0.333*reps*i*(rank+1));
         min = result < min ? result : min;
         reps++;
       }
    }

    READEX_REGION_STOP(TAN);
    return reps;
}
 
void TAN_readex(int repsLimit)
{   
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    double start = MPI_Wtime();

    double min  = 0.0;
   
    long long int reps = TAN_wrapper(repsLimit, min); 

    if (!rank)
        printf("TAN time: %f (reps %lld, min %f)\n", MPI_Wtime() - start, reps, min);
}

int SpMVmultCSR_readex_wrapper(espreso::SparseMatrix & matrix,
                               std::vector<double> & in_vector_db, 
                               std::vector<double> & out_vector_db,
                               int repsLimit)
{
    int reps = 0;
    READEX_REGION_START(SPARSE0, "SPARSE0", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        matrix.MatVec(in_vector_db, out_vector_db, 'N');
        ++reps;
    } while (reps < repsLimit); //!areWeDone(start, rank));

    READEX_REGION_STOP(SPARSE0);
    return reps;
}                               


void SpMVmultCSR_readex(int repsLimit, 
                        espreso::SparseMatrix & matrix,  
                        std::vector<double> & in_vector_db, 
                        std::vector<double> & out_vector_db)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    matrix.MatVec(in_vector_db, out_vector_db, 'N');
    matrix.MatVec(in_vector_db, out_vector_db, 'N');

    double start = MPI_Wtime();

    int reps = SpMVmultCSR_readex_wrapper(matrix, in_vector_db, out_vector_db, repsLimit);

    double end = MPI_Wtime();

    if (!rank)
        std::cout << "Sparse Matrix Vector Mult - CSR                   : " << (end - start) / reps << " (reps " << reps << ")\n";
}

int SpMVmultIJV_readex_wrapper(espreso::SparseMatrix & matrix,  
                                std::vector<double> & in_vector_db, 
                                std::vector<double> & out_vector_db,
                                int repsLimit)
{
    int reps = 0;
    READEX_REGION_START(SPARSE1, "SPARSE1", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        matrix.MatVecCOO(in_vector_db, out_vector_db, 'N');
        ++reps;
    } while (reps < repsLimit);//!areWeDone(start, rank));

    READEX_REGION_STOP(SPARSE1);
    return reps;
}                                


void SpMVmultIJV_readex(int repsLimit, 
                        espreso::SparseMatrix & matrix,  
                        std::vector<double> & in_vector_db, 
                        std::vector<double> & out_vector_db)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    matrix.MatVecCOO(in_vector_db, out_vector_db, 'N');
    matrix.MatVecCOO(in_vector_db, out_vector_db, 'N');

    double start = MPI_Wtime();
    
    int reps = SpMVmultIJV_readex_wrapper(matrix, in_vector_db, out_vector_db, repsLimit);

    double end = MPI_Wtime();

    if (!rank)
        std::cout << "Sparse Matrix Vector Mult - IJV/COO               : " << (end - start) / reps << " (reps " << reps << ")\n";
}

int SpMMmult_wrapper(int repsLimit,
                      espreso::SparseMatrix & matrix, 
                      espreso::SparseMatrix & matrix2, 
                      espreso::SparseMatrix & matrix3,  
                      std::vector<double> & in_vector_db, 
                      std::vector<double> & out_vector_db)
{
    int reps = 0;
    READEX_REGION_START(SPARSE2, "SPARSE2", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        matrix3.MatMat(matrix, 'N', matrix2);
        ++reps;
    } while (reps < repsLimit);//!areWeDone(start, rank));

    READEX_REGION_STOP(SPARSE2);

    return reps;
}

void SpMMmult(int repsLimit,
              espreso::SparseMatrix & matrix, 
              espreso::SparseMatrix & matrix2, 
              espreso::SparseMatrix & matrix3,  
              std::vector<double> & in_vector_db, 
              std::vector<double> & out_vector_db)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    matrix3.MatMat(matrix, 'N', matrix2);
    matrix3.MatMat(matrix, 'N', matrix2);

    double start = MPI_Wtime();

    int reps = SpMMmult_wrapper(RUNTIME*SPARSE2_reps, matrix, matrix2, matrix3, in_vector_db, out_vector_db);

    double end = MPI_Wtime();

    if (!rank)
        std::cout << "Sparse Matrix Matrix Mult - CSRxCSR               : " << (end - start) / reps << " (reps " << reps << ")\n";
}

int SpMMadd_wrapper(int repsLimit,
                     espreso::SparseMatrix & matrix, 
                     espreso::SparseMatrix & matrix2, 
                     espreso::SparseMatrix & matrix3,  
                     std::vector<double> & in_vector_db, 
                     std::vector<double> & out_vector_db)
{ 
    int reps = 0;

    READEX_REGION_START(SPARSE3, "SPARSE3", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        matrix3.MatAdd(matrix, matrix2, 'N', 1.0);
        ++reps;
    } while (reps < repsLimit);//!areWeDone(start, rank));

    READEX_REGION_STOP(SPARSE3);

    return reps;
}

void SpMMadd(int repsLimit,
             espreso::SparseMatrix & matrix, 
             espreso::SparseMatrix & matrix2, 
             espreso::SparseMatrix & matrix3,  
             std::vector<double> & in_vector_db, 
             std::vector<double> & out_vector_db)
{ 
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    matrix3.MatAdd(matrix, matrix2, 'N', 1.0);
    matrix3.MatAdd(matrix, matrix2, 'N', 1.0);


    double start = MPI_Wtime();

    int reps = SpMMadd_wrapper(RUNTIME*SPARSE3_reps, matrix, matrix2, matrix3, in_vector_db, out_vector_db);

    double end = MPI_Wtime();

    if (!rank)
        std::cout << "Sparse Matrix Matrix Add - CSR+CSR                : " << (end - start) / reps << " (reps " << reps << ")\n";
}

int MPI_A2A_readex_wrapper(int repsLimit, double *A, double *B)
{ 
    int reps = 0;
    READEX_REGION_START(MPI_ata, "MPI_ata", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        MPI_Alltoall(A, mpiATASize * mpiATASize, MPI_DOUBLE, B, mpiATASize * mpiATASize, MPI_DOUBLE, MPI_COMM_WORLD);
        reps++;
    } while (reps < repsLimit); //!areWeDone(start, rank));

    READEX_REGION_STOP(MPI_ata);
    return reps;

}

int MPI_A2A_readex(int repsLimit)
{   
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double* A = (double*)mkl_malloc(mpiATASize * mpiATASize * size * sizeof(double), 64);
    double* B = (double*)mkl_malloc(mpiATASize * mpiATASize * size * sizeof(double), 64);
    for (int i = 0; i < mpiATASize * mpiATASize * size; i++) {
        A[i] = (double)i; //(double)(i+1);
    }

    double start = MPI_Wtime();

    int reps = MPI_A2A_readex_wrapper(repsLimit, A, B);

    if (!rank)
        printf("MPI time: %f (reps %d)\n", MPI_Wtime() - start, reps);
}

int MPI_ring_readex_wrapper(int repsLimit, double *A, double *B)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int reps = 0;

    MPI_Request send_request, recv_request;
    MPI_Status status;

    READEX_REGION_START(MPI_ring, "MPI_ring", SCOREP_USER_REGION_TYPE_COMMON);

    do {
      MPI_Irecv(B, mpiRingSize * mpiRingSize * size, MPI_DOUBLE, (rank-1)<0 ? size-1 : rank-1, 0, MPI_COMM_WORLD, &recv_request);  
      MPI_Isend(A, mpiRingSize * mpiRingSize * size, MPI_DOUBLE, (rank+1)%size, 0, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&recv_request, &status);
      reps++;
    } while (reps < repsLimit); //!areWeDone(start, rank));

    READEX_REGION_STOP(MPI_ring);

    return reps;
}

void MPI_ring_readex(int repsLimit)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double* A = (double*)mkl_malloc(mpiRingSize * mpiRingSize * size * sizeof(double), 64);
    double* B = (double*)mkl_malloc(mpiRingSize * mpiRingSize * size * sizeof(double), 64);
    for (int i = 0; i < mpiRingSize * mpiRingSize * size; i++) {
        A[i] = (double)i; //(double)(i+1);
    }

    double start = MPI_Wtime();

    int reps = MPI_ring_readex_wrapper(repsLimit, A, B);

    if (!rank)
        printf("MPI time: %f (reps %d)\n", MPI_Wtime() - start, reps);
}

int DGEMM_wrapper(double* A, double* B, double* C, int repLimit)
{
    int reps = 0;
    
    READEX_REGION_START(DGEMM_, "DGEMM", SCOREP_USER_REGION_TYPE_COMMON);
    
    do {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            matrixSizeDGEMM, matrixSizeDGEMM, matrixSizeDGEMM, 1.0, A, matrixSizeDGEMM, B, matrixSizeDGEMM, 0.0, C, matrixSizeDGEMM);
        reps++;
    } while (reps < repLimit); //!areWeDone(start, rank));
    
    READEX_REGION_STOP(DGEMM_);
    
    return reps;
}

void DGEMM_readex(int repLimit)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
   
    double *A, *B, *C;
    int m, n, k, i, j, vMatN;
    double alpha, beta;
    time_t regionTime;

    m = matrixSizeDGEMM, k = matrixSizeDGEMM, n = matrixSizeDGEMM;

    alpha = 1.0;
    beta = 0.0;

    A = (double*)mkl_malloc(m * k * sizeof(double), 64);
    B = (double*)mkl_malloc(k * n * sizeof(double), 64);
    C = (double*)mkl_malloc(m * n * sizeof(double), 64);

    if (A == NULL || B == NULL || C == NULL) {
        printf("\n ERROR: Can't allocate memory. Aborting... \n\n");
        mkl_free(A);
        mkl_free(B);
        mkl_free(C);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    for (i = 0; i < (m * k); i++) {
        A[i] = (double)i; //(double)(i+1);
    }

    for (i = 0; i < (k * n); i++) {
        B[i] = (double)i; //(double)(-i-1);
    }

    for (i = 0; i < (m * n); i++) {
        C[i] = 0.0;
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        m, n, k, alpha, A, k, B, n, beta, C, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        m, n, k, alpha, A, k, B, n, beta, C, n);


    double start = MPI_Wtime();
    
    
    int reps = DGEMM_wrapper(A,B,C,repLimit);

    if (!rank)
        printf("DGEMM TIME: %f, (reps %d)\n", MPI_Wtime() - start, reps);

    mkl_free(A);
    mkl_free(B);
    mkl_free(C);
}

int DGEMV_wrapper(int vMatN, double alpha, double* Vmat, double* v, int incX, double beta, double* vOut, int incY, int repLimit)
{
    int reps = 0;
    READEX_REGION_START(DGEMV_, "DGEMV", SCOREP_USER_REGION_TYPE_COMMON);

    do {
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
            vMatN, vMatN, alpha, Vmat, vMatN, v, incX, beta, vOut, incY);
        reps++;
    } while (reps < repLimit); //!areWeDone(start, rank));

    READEX_REGION_STOP(DGEMV_);
    return reps;
}

int DGEMV_wrapper(int vMatN, double alpha, double* Vmat1, double* Vmat2, double* v, int incX, double beta, double* vOut, int incY, int repLimit)
{
	int reps = 0;
	READEX_REGION_START(DGEMV_, "DGEMV", SCOREP_USER_REGION_TYPE_COMMON);

	do {
		if (reps%2)
			cblas_dgemv(CblasRowMajor, CblasNoTrans,
				vMatN, vMatN, alpha, Vmat1, vMatN, v, incX, beta, vOut, incY);
		else
			cblas_dgemv(CblasRowMajor, CblasNoTrans,
				vMatN, vMatN, alpha, Vmat2, vMatN, v, incX, beta, vOut, incY);
	reps++;
	} while (reps < repLimit); //!areWeDone(start, rank));

	READEX_REGION_STOP(DGEMV_);
	return reps;
}


void DGEMV_readex(int repLimit)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
       
    double *Vmat1, *Vmat2, *v, *vOut;
    int m, n, k, i, j, vMatN;
    double alpha, beta;
    time_t regionTime;

    m = matrixSizeDGEMV, k = matrixSizeDGEMV, n = matrixSizeDGEMV, vMatN = matrixSizeDGEMV;

    alpha = 1.0;
    beta = 0.0;

    Vmat1 = (double*)mkl_malloc(vMatN * vMatN * sizeof(double), 64);
    Vmat2 = (double*)mkl_malloc(vMatN * vMatN * sizeof(double), 64);

    v = (double*)mkl_malloc(vMatN * sizeof(double), 64);
    vOut = (double*)mkl_malloc(vMatN * sizeof(double), 64);

    if (v == NULL || vOut == NULL) {
        printf("\n ERROR: Can't allocate memory. Aborting... \n\n");
        mkl_free(Vmat1);
        mkl_free(Vmat2);
        mkl_free(v);
        mkl_free(vOut);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    for (i = 0; i < (vMatN * vMatN); i++) {
        Vmat1[i] = 1.0;
        Vmat2[i] = i/2.0;
    }
    for (i = 0; i < vMatN; i++) {
        v[i] = (double)i;
        vOut[i] = 0.0;
    }

    cblas_dgemv(CblasRowMajor, CblasNoTrans,
        vMatN, vMatN, alpha, Vmat1, vMatN, v, 1, beta, vOut, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
        vMatN, vMatN, alpha, Vmat2, vMatN, v, 1, beta, vOut, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
        vMatN, vMatN, alpha, Vmat1, vMatN, v, 1, beta, vOut, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
        vMatN, vMatN, alpha, Vmat2, vMatN, v, 1, beta, vOut, 1);

    
    double start = MPI_Wtime();

    int reps = DGEMV_wrapper(vMatN, alpha, Vmat1, Vmat2, v, 1, beta, vOut, 1, repLimit);

    if (!rank)
        printf("DGEMV TIME: %f (reps %d)\n", MPI_Wtime() - start, reps);

    mkl_free(Vmat1);
    mkl_free(Vmat2);
    mkl_free(v);
    mkl_free(vOut);
}


/* ./blasbench MATRIX.mtx threadsUsedForIO whatToTest */
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    READEX_INIT();

    ParseCommandLine(argc, argv);

    for (int i = 0; i < PHASE_reps; ++i) {

    if( !rank )
      std::cout << "Starting phase " << i << std::endl;

    READEX_PHASE_START(Main, "Main", SCOREP_USER_REGION_TYPE_COMMON);

    espreso::SparseMatrix matrix;
    espreso::SparseMatrix matrix3;
    if (IO_flag)
    {
      //READEX_REGION_START(Parallel_IO, "Parallel_IO", SCOREP_USER_REGION_TYPE_COMMON);
      Parallel_IO_readex(matrix);
      //READEX_REGION_STOP(Parallel_IO);
    }
    espreso::SparseMatrix matrix2(matrix);
    matrix2.ConvertToCSR(0);

    double start = omp_get_wtime();
    matrix.ConvertToCSR(0); // Convert to CSR and keep IJV data as well
    double end = omp_get_wtime();
    if (!rank && IO_flag)
        std::cout << "Sparse Matrix - convert IJV to CSR format         : " << end - start << "\n";
    

    std::vector<double> in_vector_db(matrix.rows, 0.0);
    std::vector<float> in_vector_fl(matrix.rows, 0.0);

    std::vector<double> out_vector_db(matrix.rows);
    std::vector<float> out_vector_fl(matrix.rows);

    if (SPARSE0_flag)
    {
       //READEX_REGION_START(SPARSE0, "SPARSE0", SCOREP_USER_REGION_TYPE_COMMON);
         SpMVmultCSR_readex(RUNTIME*SPARSE0_reps, matrix, in_vector_db, out_vector_db);
         //READEX_REGION_STOP(SPARSE0);
    }
    if (SPARSE1_flag)
    {
         //READEX_REGION_START(SPARSE1, "SPARSE1", SCOREP_USER_REGION_TYPE_COMMON);
         SpMVmultIJV_readex(RUNTIME*SPARSE1_reps, matrix, in_vector_db, out_vector_db); 
         //READEX_REGION_STOP(SPARSE1);   
    }
    if (SPARSE2_flag)
    {
        //READEX_REGION_START(SPARSE2, "SPARSE2", SCOREP_USER_REGION_TYPE_COMMON);
        SpMMmult(RUNTIME*SPARSE2_reps, matrix, matrix2, matrix3, in_vector_db, out_vector_db);
        //READEX_REGION_STOP(SPARSE2);
    }
    if (SPARSE3_flag)
    {
        //READEX_REGION_START(SPARSE3, "SPARSE3", SCOREP_USER_REGION_TYPE_COMMON);
        SpMMadd(RUNTIME*SPARSE3_reps, matrix, matrix2, matrix3, in_vector_db, out_vector_db);
        //READEX_REGION_STOP(SPARSE3);
    }
    if (DENSE0_flag) 
    {
        //READEX_REGION_START(DGEMM_, "DGEMM", SCOREP_USER_REGION_TYPE_COMMON);
        DGEMM_readex(RUNTIME*DENSE0_reps);
        //READEX_REGION_STOP(DGEMM_);
    }
    if (DENSE1_flag)
    {
        //READEX_REGION_START(DGEMV_, "DGEMV", SCOREP_USER_REGION_TYPE_COMMON);
        DGEMV_readex(RUNTIME*DENSE1_reps);
        //READEX_REGION_STOP(DGEMV_);
    }
    if (TAN_flag)
    {
        //READEX_REGION_START(TAN, "TAN", SCOREP_USER_REGION_TYPE_COMMON);
        TAN_readex(RUNTIME*TAN_reps);
        //READEX_REGION_STOP(TAN);
    }    

    if (MPI_ata_flag) 
    {
        //READEX_REGION_START(MPI_ata, "MPI_ata", SCOREP_USER_REGION_TYPE_COMMON);
        MPI_A2A_readex(RUNTIME*MPI_ata_reps);
        //READEX_REGION_STOP(MPI_ata);
    }    
    if (MPI_ring_flag) 
    {
        //READEX_REGION_START(MPI_ring, "MPI_ring", SCOREP_USER_REGION_TYPE_COMMON);
        MPI_ring_readex(RUNTIME*MPI_ring_reps);
        //READEX_REGION_STOP(MPI_ring);
    }

    READEX_PHASE_STOP(Main);

    if( !rank )
      std::cout << "Finished phase " << i << std::endl;
    
    }

    READEX_CLOSE();
    MPI_Finalize();

    return 0;
}

