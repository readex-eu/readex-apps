#include "../../Settings.h"

#include <iostream>
#include <mpi.h>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../SparseMatrix.h"
#include "../../MPIBlockMatrix.h"
#include "../../Vector.h"

using namespace std;
using namespace bem4i;

void testMPIMatrix( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  testMPIMatrix( );

  return 0;
}

void testMPIMatrix( ) {
  MPI_Init( nullptr, nullptr );
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  typedef double SC;
  typedef int LO;

  // test block matrices
  SparseMatrix<LO, SC> A11( 3, 3 );
  SparseMatrix<LO, SC> A12( 3, 3 );
  SparseMatrix<LO, SC> A21( 3, 3 );
  FullMatrix<LO, SC> A22( 3, 3 );
  A11.set( 1, 1, 1.0 );
  A12.set( 0, 0, 1.0 );
  A21.set( 2, 1, 1.0 );
  A22.set( 1, 2, 1.0 );
  LO nRows[2] = { 3, 3 };
  LO nCols[2] = { 3, 3 };
  int ranks[4] = { 0, 1, 2, 3 };

  MPIBlockMatrix<LO, SC> matrix( 2, 2, &nRows[0], &nCols[0], ranks, MPI_COMM_WORLD );

  if ( rank == 0 ) {

    matrix.setBlock( 0, 0, &A11 );
  }
  matrix.setBlock( 0, 1, &A12 );
  matrix.setBlock( 1, 0, &A21 );
  matrix.setBlock( 1, 1, &A22 );
  Vector<LO, SC> x2( 6 );
  Vector<LO, SC> y2( 6 );
  x2.setAll( 1.0 );
  y2.setAll( 1.0 );
  matrix.apply( x2, y2, true, 2.0, 2.0 );
  y2.print( std::cout );

  MPI_Finalize( );
}