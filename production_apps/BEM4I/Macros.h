/*!
 * @file    Macros.h
 * @author  Jan Zapletal
 * @date    July 16, 2013
 *
 */

#ifndef MACROS_H
#define	MACROS_H

#include "mpi.h"

#define EPS (0.000000000001)
//#define EPS (0.000001)
//#define EPS (eps)

#define PI_FACT ( (SCVT) 0.07957747154594766788 )

#define DOT3( a, b ) ( a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] )

#define NORM3( a ) (std::sqrt(DOT3(a,a)))

#define DIST3( a, b ) ( std::sqrt( ( a[ 0 ] - b[ 0 ] ) * ( a[ 0 ] - b[ 0 ] ) + \
( a[ 1 ] - b[ 1 ] ) * ( a[ 1 ] - b[ 1 ] ) + \
( a[ 2 ] - b[ 2 ] ) * ( a[ 2 ] - b[ 2 ] ) ) )

#define DIST3SQ( a, b ) ( ( a[ 0 ] - b[ 0 ] ) * ( a[ 0 ] - b[ 0 ] ) + \
( a[ 1 ] - b[ 1 ] ) * ( a[ 1 ] - b[ 1 ] ) + \
( a[ 2 ] - b[ 2 ] ) * ( a[ 2 ] - b[ 2 ] ) )

#define BUFFER_SIZE (40)

namespace bem4i {

template<class LO, class SC>
struct GetType {

};

template<>
struct GetType<int, std::complex<double> > {

  double eps = 1e-12;
  typedef double SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_INT;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_DOUBLE_COMPLEX;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_DOUBLE;
  }
};

template<>
struct GetType<int, std::complex<float> > {

  float eps = 1e-7;
  typedef float SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_INT;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_COMPLEX;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_FLOAT;
  }
};

template<>
struct GetType<int, double> {

  double eps = 1e-12;
  typedef double SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_INT;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_DOUBLE;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_DOUBLE;
  }
};

template<>
struct GetType<int, float> {

  double eps = 1e-7;
  typedef float SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_INT;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_FLOAT;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_FLOAT;
  }
};

template<>
struct GetType<long, std::complex<double> > {

  double eps = 1e-12;
  typedef double SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_LONG;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_DOUBLE_COMPLEX;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_DOUBLE;
  }
};

template<>
struct GetType<long, std::complex<float> > {

  float eps = 1e-7;
  typedef float SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_LONG;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_COMPLEX;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_FLOAT;
  }
};

template<>
struct GetType<long, double> {

  double eps = 1e-12;
  typedef double SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_LONG;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_DOUBLE;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_DOUBLE;
  }
};

template<>
struct GetType<long, float> {

  double eps = 1e-7;
  typedef float SCVT;

  static MPI_Datatype MPI_LO( ) {
    return MPI_LONG;
  }

  static MPI_Datatype MPI_SC( ) {
    return MPI_FLOAT;
  }

  static MPI_Datatype MPI_SCVT( ) {
    return MPI_FLOAT;
  }
};

}

#endif	/* MACROS_H */
