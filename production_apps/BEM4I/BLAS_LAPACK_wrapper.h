/*!
 * @file    BLAS_LAPACK_wrapper.h
 * @author  Michal Merta 
 * @date    July 6, 2013
 * @brief   The file contains the interface to several FORTRAN BLAS and LAPACK routines
 * 
 */

#ifndef BLAS_LAPACK_WRAPPER_H
#define BLAS_LAPACK_WRAPPER_H

#include <complex>
#include "Settings.h"

#define CHAR_MACRO(char_var) &char_var

#ifdef BLAS_INT

/// INT FOR INDEXING ///

//*** DOUBLE PRECISION IMPLEMENTATION WITH INT FOR INDEXING ***//

// LEVEL 1 BLAS

// scale a vector (matrix) by a scalar value
extern "C" {
void dscal_(
    int * N,
    double * ALPHA,
    double * X,
    int * INCX
    );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void daxpy_(
    int * N,
    double * A,
    double * X,
    int * INCX,
    double * Y,
    int * INCY
    );
}

// vector 2-norm
extern "C" {
double dnrm2_(
    int const * N,
    double * X,
    int * INCX
    );
}

// dot product
extern "C" {
double ddot_(
    int const * N,
    double * X,
    int * INCX,
    double * Y,
    int * INCY
    );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void dgemv_(
    char * TRANS,
    int * M,
    int * N,
    double * ALPHA,
    double * A,
    int * LDA,
    double * X,
    int * INCX,
    double * BETA,
    double * Y,
    int * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void dgemm_(
    char * TRANSA,
    char * TRANSB,
    int * M,
    int * N,
    int * K,
    double * ALPHA,
    double * A,
    int * LDA,
    double * B,
    int * LDB,
    double * BETA,
    double * C,
    int * LDC
    );
}

// LAPACK

// matrix norm
extern "C" {
double dlange_(
    char * norm,
    int const * M,
    int const * N,
    double const * A,
    int const * LDA,
    double * WORK
    );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void dgetrf_(
    int * M,
    int * N,
    double * A,
    int * LDA,
    int * IPIV,
    int * INFO
    );
}

extern "C" { // LAPACK LU solver
void dgetrs_(
    char * TRANS,
    int * N,
    int * NRHS,
    double * A,
    int * LDA,
    int * IPIV,
    double * B,
    int * LDB,
    int * INFO
    );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void dpotrf_(
    char * UPLO,
    int * M,
    double * A,
    int * LDA,
    int * INFO
    );
}

// LAPACK Choleski solver
extern "C" {
void dpotrs_(
    char * UPLO,
    int * N,
    int * NRHS,
    double * A,
    int * LDA,
    double * B,
    int * LDB,
    int * INFO
    );
}

// reduces real double symmetric matrix to tridiagonal form
extern "C" {
void dsytrd_( char *UPLO, int *N, double *A, int *LDA, double *D, double *E, double *TAU, double *WORK, int *LWORK, int *INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void dsyevx_( char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA,
    double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W,
    double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *IFAIL,
    int *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void dhseqr_( char *job, char *compz, int *n, int *ilo, int *ihi, double *h,
    int *ldh, double *wr, double *wi, double *z, int *ldz, double *work,
    int *lwork, int *info );
}

// computes Schur decomposition of a general matrix 
extern "C" {
void dgees_( char *jobvs, char *sort, void *select, int *n, double *h,
    int *ldh, int *sdim, double *wr, double *wi, double *z, int *ldz, double *work,
    int *lwork, void* bwork, int *info );
}


//*** SINGLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a scalar value
extern "C" {
void sscal_( int *N, float *ALPHA, float *X, int *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void saxpy_( int *N, float *A, float *X, int *INCX, float *Y, int *INCY );
}

// vector 2-norm
extern "C" {
float snrm2_( int const *N, float *X, int *INCX );
}

// dot product
extern "C" {
float sdot_( int const *N, float *X, int *INCX, float *Y, int *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void sgemv_(
    char * TRANS,
    int * M,
    int * N,
    float * ALPHA,
    float * A,
    int * LDA,
    float * X,
    int * INCX,
    float * BETA,
    float * Y,
    int * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void sgemm_( char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA, float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC );
}

// LAPACK

// matrix norm
extern "C" {
float slange_( char *norm, int const *M, int const *N, float const *A, int const *LDA, float *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void sgetrf_( int *M, int *N, float *A, int *LDA, int *IPIV, int *INFO );
}

extern "C" { // LAPACK LU solver
void sgetrs_( char *TRANS, int *N, int *NRHS, float *A, int *LDA, int *IPIV, float *B, int *LDB, int *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void spotrf_( char * UPLO, int * M, float * A, int * LDA, int * INFO );
}

// LAPACK Choleski solver
extern "C" {
void spotrs_( char * UPLO, int * N, int * NRHS, float * A, int * LDA,
    float * B, int * LDB, int * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void ssyevx_( char *JOBZ, char *RANGE, char *UPLO, int *N, float *A, int *LDA, float *VL, float *VU, int *IL, int *IU,
    float *ABSTOL, int *M, float *W, float *Z, int *LDZ, float *WORK, int *LWORK, int *IWORK,
    int *IFAIL, int *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void shseqr_( char *job, char *compz, int *n, int *ilo, int *ihi, float *h,
    int *ldh, float *wr, float *wi, float *z, int *ldz, float *work,
    int *lwork, int *info );
}

//*** COMPLEX SINGLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a complex value
extern "C" {
void cscal_( int *N, std::complex<float> *ALPHA, std::complex<float> *X, int *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void caxpy_( int *N, std::complex<float> *A, std::complex<float> *X, int *INCX, std::complex<float> *Y, int *INCY );
}

// vector 2-norm (cnrm2 returns only float but we need this because of templates)
extern "C" {
float scnrm2_( int const *N, std::complex<float> *X, int *INCX );
}

// conjugated dot product (cdotc returns only float but we need this because of templates)
extern "C" {
std::complex<float> cdotc_( int const *N, std::complex<float> *X, int *INCX, std::complex<float> *Y, int *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void cgemv_(
    char * TRANS,
    int * M,
    int * N,
    std::complex<float> * ALPHA,
    std::complex<float> * A,
    int * LDA,
    std::complex<float> * X,
    int * INCX,
    std::complex<float> * BETA,
    std::complex<float> * Y,
    int * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void cgemm_( char *TRANSA, char *TRANSB, int *M, int *N, int *K, std::complex<float> *ALPHA, std::complex<float> *A, int *LDA, std::complex<float> *B, int *LDB, std::complex<float> *BETA, std::complex<float> *C, int *LDC );
}

// LAPACK

// matrix norm
extern "C" {
float clange_( char *norm, int const *M, int const *N, std::complex<float> const *A, int const *LDA, float *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void cgetrf_( int *M, int *N, std::complex<float> *A, int *LDA, int *IPIV, int *INFO );
}

// LAPACK LU solver
extern "C" {
void cgetrs_( char *TRANS, int *N, int *NRHS, std::complex<float> *A, int *LDA, int *IPIV, std::complex<float> *B, int *LDB, int *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void cpotrf_( char * UPLO, int * M, std::complex< float > * A, int * LDA,
    int * INFO );
}

// LAPACK Choleski solver
extern "C" {
void cpotrs_( char * UPLO, int * N, int * NRHS, std::complex< float > * A,
    int * LDA, std::complex< float > * B, int * LDB, int * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void cheevx_( char *JOBZ, char *RANGE, char *UPLO, int *N, std::complex<float> *A, int *LDA, float *VL, float *VU, int *IL, int *IU,
    float *ABSTOL, int *M, float *W, std::complex<float> *Z, int *LDZ, std::complex<float> *WORK, int *LWORK, float *RWORK, int *IWORK,
    int *IFAIL, int *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void chseqr_( char *job, char *compz, int *n, int *ilo, int *ihi,
    std::complex<float> *h, int *ldh, std::complex<float> *w,
    std::complex<float> *z, int *ldz,
    std::complex<float> *work, int *lwork, int *info );
}
//*** COMPLEX DOUBLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a complex value
extern "C" {
void zscal_( int *N, std::complex<double> *ALPHA, std::complex<double> *X, int *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void zaxpy_( int *N, std::complex<double> *A, std::complex<double> *X, int *INCX, std::complex<double> *Y, int *INCY );
}

// vector 2-norm (cnrm2 returns only float but we need this because of templates)
extern "C" {
double dznrm2_( int const *N, std::complex<double> *X, int *INCX );
}

// conjugated dot product (cdotc returns only float but we need this because of templates)
extern "C" {
std::complex<double> zdotc_( int const *N, std::complex<double> *X, int *INCX, std::complex<double> *Y, int *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void zgemv_(
    char * TRANS,
    int * M,
    int * N,
    std::complex<double> * ALPHA,
    std::complex<double> * A,
    int * LDA,
    std::complex<double> * X,
    int * INCX,
    std::complex<double> * BETA,
    std::complex<double> * Y,
    int * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void zgemm_( char *TRANSA, char *TRANSB, int *M, int *N, int *K,
    std::complex<double> *ALPHA, std::complex<double> *A, int *LDA,
    std::complex<double> *B, int *LDB, std::complex<double> *BETA,
    std::complex<double> *C, int *LDC );
}

// LAPACK

// matrix norm
extern "C" {
double zlange_( char *norm, int const *M, int const *N, std::complex<double> const *A, int const *LDA, double *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void zgetrf_( int *M, int *N, std::complex<double> *A, int *LDA, int *IPIV, int *INFO );
}

// LAPACK LU solver
extern "C" {
void zgetrs_( char *TRANS, int *N, int *NRHS, std::complex<double> *A, int *LDA, int *IPIV, std::complex<double> *B, int *LDB, int *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void zpotrf_( char * UPLO, int * M, std::complex< double > * A, int * LDA, int * INFO );
}

// LAPACK Choleski solver
extern "C" {
void zpotrs_( char * UPLO, int * N, int * NRHS, std::complex< double > * A,
    int * LDA, std::complex< double > * B, int * LDB, int * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void zheevx_( char *JOBZ, char *RANGE, char *UPLO, int *N, std::complex<double> *A, int *LDA, double *VL, double *VU, int *IL, int *IU,
    double *ABSTOL, int *M, double *W, std::complex<double> *Z, int *LDZ, std::complex<double> *WORK, int *LWORK, double *RWORK, int *IWORK,
    int *IFAIL, int *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void zhseqr_( char *job, char *compz, int *n, int *ilo, int *ihi,
    std::complex<double> *h, int *ldh, std::complex<double> *w,
    std::complex<double> *z, int *ldz,
    std::complex<double> *work, int *lwork, int *info );
}
#endif

#ifdef BLAS_LONG

// LONG FOR INDEXINT ///

//*** DOUBLE PRECISION IMPLEMENTATION ***//

// LEVEL 1 BLAS

// scale a vector (matrix) by a scalar value
extern "C" {
void dscal_( long *N, double *ALPHA, double *X, long *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void daxpy_( long *N, double *A, double *X, long *INCX, double *Y, long *INCY );
}

// vector 2-norm
extern "C" {
double dnrm2_( long const *N, double *X, long *INCX );
}

// dot product
extern "C" {
double ddot_( long const *N, double *X, long *INCX, double *Y, long *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void dgemv_(
    char * TRANS,
    long * M,
    long * N,
    double * ALPHA,
    double * A,
    long * LDA,
    double * X,
    long * INCX,
    double * BETA,
    double * Y,
    long * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void dgemm_( char *TRANSA, char *TRANSB, long *M, long *N, long *K, double *ALPHA, double *A, long *LDA, double *B, long *LDB, double *BETA, double *C, long *LDC );
}

// LAPACK

// matrix norm
extern "C" {
double dlange_( char *norm, long const *M, long const *N, double const *A, long const *LDA, double *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void dgetrf_( long *M, long *N, double *A, long *LDA, long *IPIV, long *INFO );
}

extern "C" { // LAPACK LU solver
void dgetrs_( char *TRANS, long *N, long *NRHS, double *A, long *LDA, long *IPIV, double *B, long *LDB, long *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void dpotrf_( char * UPLO, long * M, double * A, long * LDA, long * INFO );
}

// LAPACK Choleski solver
extern "C" {
void dpotrs_( char * UPLO, long * N, long * NRHS, double * A, long * LDA,
    double * B, long * LDB, long * INFO );
}

// reduces real double symmetric matrix to tridiagonal form
extern "C" {
void dsytrd_( char *UPLO, long *N, double *A, long *LDA, double *D, double *E, double *TAU, double *WORK, long *LWORK, long *INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void dsyevx_( char *JOBZ, char *RANGE, char *UPLO, long *N, double *A, long *LDA, double *VL, double *VU, long *IL, long *IU,
    double *ABSTOL, long *M, double *W, double *Z, long *LDZ, double *WORK, long *LWORK, long *IWORK,
    long *IFAIL, long *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void dhseqr_( char *job, char *compz, long *n, long *ilo, long *ihi, double *h,
    long *ldh, double *wr, double *wi, double *z, long *ldz, double *work,
    long *lwork, long *info );
}

//*** SINGLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a scalar value
extern "C" {
void sscal_( long *N, float *ALPHA, float *X, long *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void saxpy_( long *N, float *A, float *X, long *INCX, float *Y, long *INCY );
}

// vector 2-norm
extern "C" {
float snrm2_( long const *N, float *X, long *INCX );
}

// dot product
extern "C" {
float sdot_( long const *N, float *X, long *INCX, float *Y, long *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void sgemv_(
    char * TRANS,
    long * M,
    long * N,
    float * ALPHA,
    float * A,
    long * LDA,
    float * X,
    long * INCX,
    float * BETA,
    float * Y,
    long * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void sgemm_( char *TRANSA, char *TRANSB, long *M, long *N, long *K, float *ALPHA, float *A, long *LDA, float *B, long *LDB, float *BETA, float *C, long *LDC );
}

// LAPACK

// matrix norm
extern "C" {
float slange_( char *norm, long const *M, long const *N, float const *A, long const *LDA, float *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void sgetrf_( long *M, long *N, float *A, long *LDA, long *IPIV, long *INFO );
}

extern "C" { // LAPACK LU solver
void sgetrs_( char *TRANS, long *N, long *NRHS, float *A, long *LDA, long *IPIV, float *B, long *LDB, long *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void spotrf_( char * UPLO, long * M, float * A, long * LDA, long * INFO );
}

// LAPACK Choleski solver
extern "C" {
void spotrs_( char * UPLO, long * N, long * NRHS, float * A, long * LDA,
    float * B, long * LDB, long * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void ssyevx_( char *JOBZ, char *RANGE, char *UPLO, long *N, float *A, long *LDA, float *VL, float *VU, long *IL, long *IU,
    float *ABSTOL, long *M, float *W, float *Z, long *LDZ, float *WORK, long *LWORK, long *IWORK,
    long *IFAIL, long *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void shseqr_( char *job, char *compz, long *n, long *ilo, long *ihi, float *h,
    long *ldh, float *wr, float *wi, float *z, long *ldz, float *work,
    long *lwork, long *info );
}

//*** COMPLEX SINGLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a complex value
extern "C" {
void cscal_( long *N, std::complex<float> *ALPHA, std::complex<float> *X, long *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void caxpy_( long *N, std::complex<float> *A, std::complex<float> *X, long *INCX, std::complex<float> *Y, long *INCY );
}

// vector 2-norm (cnrm2 returns only float but we need this because of templates)
extern "C" {
float scnrm2_( long const *N, std::complex<float> *X, long *INCX );
}

// conjugated dot product (cdotc returns only float but we need this because of templates)
extern "C" {
std::complex<float> cdotc_( long const *N, std::complex<float> *X, long *INCX, std::complex<float> *Y, long *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void cgemv_(
    char * TRANS,
    long * M,
    long * N,
    std::complex<float> * ALPHA,
    std::complex<float> * A,
    long * LDA,
    std::complex<float> * X,
    long * INCX,
    std::complex<float> * BETA,
    std::complex<float> * Y,
    long * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void cgemm_( char *TRANSA, char *TRANSB, long *M, long *N, long *K, std::complex<float> *ALPHA, std::complex<float> *A, long *LDA, std::complex<float> *B, long *LDB, std::complex<float> *BETA, std::complex<float> *C, long *LDC );
}

// LAPACK

// matrix norm
extern "C" {
float clange_( char *norm, long const *M, long const *N, std::complex<float> const *A, long const *LDA, float *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void cgetrf_( long *M, long *N, std::complex<float> *A, long *LDA, long *IPIV, long *INFO );
}

// LAPACK LU solver
extern "C" {
void cgetrs_( char *TRANS, long *N, long *NRHS, std::complex<float> *A, long *LDA, long *IPIV, std::complex<float> *B, long *LDB, long *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void cpotrf_( char * UPLO, long * M, std::complex< float > * A, long * LDA,
    long * INFO );
}

// LAPACK Choleski solver
extern "C" {
void cpotrs_( char * UPLO, long * N, long * NRHS, std::complex< float > * A,
    long * LDA, std::complex< float > * B, long * LDB, long * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void cheevx_( char *JOBZ, char *RANGE, char *UPLO, long *N, std::complex<float> *A, long *LDA, float *VL, float *VU, long *IL, long *IU,
    float *ABSTOL, long *M, float *W, std::complex<float> *Z, long *LDZ, std::complex<float> *WORK, long *LWORK, float *RWORK, long *IWORK,
    long *IFAIL, long *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void chseqr_( char *job, char *compz, long *n, long *ilo, long *ihi,
    std::complex<float> *h, long *ldh, std::complex<float> *w,
    std::complex<float> *z, long *ldz,
    std::complex<float> *work, long *lwork, long *info );
}

//*** COMPLEX DOUBLE PRECISION IMPLEMENTATION ***//


// LEVEL 1 BLAS

// scale a vector (matrix) by a complex value
extern "C" {
void zscal_( long *N, std::complex<double> *ALPHA, std::complex<double> *X, long *INCX );
}

// vector (matrix) addition y = y + alpha*x
extern "C" {
void zaxpy_( long *N, std::complex<double> *A, std::complex<double> *X, long *INCX, std::complex<double> *Y, long *INCY );
}

// vector 2-norm (cnrm2 returns only float but we need this because of templates)
extern "C" {
double dznrm2_( long const *N, std::complex<double> *X, long *INCX );
}

// conjugated dot product (cdotc returns only float but we need this because of templates)
extern "C" {
std::complex<double> zdotc_( long const *N, std::complex<double> *X, long *INCX, std::complex<double> *Y, long *INCY );
}

// LEVEL 2 BLAS

// BLAS matrix-vector multiplication
extern "C" {
#if N_MIC > 0
__attribute__( ( target( mic ) ) )
#endif
void zgemv_(
    char * TRANS,
    long * M,
    long * N,
    std::complex<double> * ALPHA,
    std::complex<double> * A,
    long * LDA,
    std::complex<double> * X,
    long * INCX,
    std::complex<double> * BETA,
    std::complex<double> * Y,
    long * INCY
    );
}


// LEVEL 3 BLAS

// matrix-matrix multiplication
extern "C" {
void zgemm_( char *TRANSA, char *TRANSB, long *M, long *N, long *K, std::complex<double> *ALPHA, std::complex<double> *A, long *LDA, std::complex<double> *B, long *LDB, std::complex<double> *BETA, std::complex<double> *C, long *LDC );
}

// LAPACK

// matrix norm
extern "C" {
double zlange_( char *norm, long const *M, long const *N, std::complex<double> const *A, long const *LDA, double *WORK );
}

// LAPACK LU factorization of a general matrix
extern "C" {
void zgetrf_( long *M, long *N, std::complex<double> *A, long *LDA, long *IPIV, long *INFO );
}

// LAPACK LU solver
extern "C" {
void zgetrs_( char *TRANS, long *N, long *NRHS, std::complex<double> *A, long *LDA, long *IPIV, std::complex<double> *B, long *LDB, long *INFO );
}

// LAPACK Choleski factorization of a general matrix
extern "C" {
void zpotrf_( char * UPLO, long * M, std::complex< double > * A, long * LDA,
    long * INFO );
}

// LAPACK Choleski solver
extern "C" {
void zpotrs_( char * UPLO, long * N, long * NRHS, std::complex< double > * A,
    long * LDA, std::complex< double > * B, long * LDB, long * INFO );
}

// computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
extern "C" {
void zheevx_( char *JOBZ, char *RANGE, char *UPLO, long *N, std::complex<double> *A, long *LDA, double *VL, double *VU, long *IL, long *IU,
    double *ABSTOL, long *M, double *W, std::complex<double> *Z, long *LDZ, std::complex<double> *WORK, long *LWORK, double *RWORK, long *IWORK,
    long *IFAIL, long *INFO );
}

// computes Schur decomposition of a Hessenberg matrix H
extern "C" {
void zhseqr_( char *job, char *compz, long *n, long *ilo, long *ihi,
    std::complex<double> *h, long *ldh, std::complex<double> *w,
    std::complex<double> *z, long *ldz,
    std::complex<double> *work, long *lwork, long *info );
}
#endif

/*
int dgetrf(const int &nMatRows,
        const int &nMatCols,
        double *dMatrix,
        const int &nMatrixLeadingDim,
        int *nPivotVector) {
  // Return status
  int nReturnCode;

  // Call method
  dgetrf_(
          nMatRows,
          nMatCols,
          dMatrix,
          nMatrixLeadingDim,
          nPivotVector,
          nReturnCode);

  // Return result;
  return nReturnCode;
};

void dgetrs(char TRANS, int n, int nrhs, double* A, int lda, int* IPIV, double* B, int ldb, int* info) {
  dgetrs_(CHAR_MACRO(TRANS), &n, &nrhs, A, &lda, IPIV, B, &ldb, info);
}


 */
#endif /* BLAS_LAPACK_WRAPPER_H */

