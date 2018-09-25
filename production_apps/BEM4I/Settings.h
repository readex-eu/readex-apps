/*!
 * @file    settings.h
 * @author  Michal Merta
 * @date    July 23, 2013
 * @brief   The file contains compilation-time settings
 *
 */
#ifndef SETTINGS_H
#define	SETTINGS_H

//#define OPTIMIZATION

// loads fenv.h and checks for floating-point exceptions
//#define FENV

/*
 * Intel Xeon Phi settings
 */
#define N_MIC ( 0 )
#define MIC_CHUNK ( 50.0e6 )
#define N_MIC_THREADS ( 240 )

// Host power relatively to the Intel Xeon Phi coprocessor
#define HOST_TO_MIC_PERF ( 1.0 )

#define DATA_ALIGN ( 64 )
#define DATA_WIDTH ( 4 )

/*!
 *  define the size of BLAS and LAPACK indexing
 *  this is necessary to load the right interface to lapack and blas
 * BLAS_INT     4 bytes
 * BLAS_LONG    8 bytes
 */
//#define BLAS_INT
#define BLAS_LONG

/*!
 * is PARDISO present at the system?
 */
#undef HAS_PARDISO

/*!
 *  whether to use experimental basis functions for time domain wave equation
 */
#undef EXPERIMENTAL_WAVE

/*!
 * is IPOPT present at the system?
 */
//#define IPOPT

//#define VERBOSE
//#define MORE_VERBOSE
//#define SILENT

#ifdef MORE_VERBOSE
#define VERBOSE
#endif

/*!
 *  define what technology to use to accelerate FMM matrix-vector multiplication
 * FMM_SCALAR     no acceleration
 * FMM_OMP        use OpenMP to accelerate matrix-vector multiplication
 */
//#define FMM_SCALAR
#define FMM_OMP

#endif	/* SETTINGS_H */
