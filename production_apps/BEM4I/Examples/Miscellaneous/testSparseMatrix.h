#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../SparseMatrix.h"
#include "../../BlockMatrix.h"

using namespace std;
using namespace bem4i;

void testSparseMatrix( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  testSparseMatrix( );

  return 0;
}

void testSparseMatrix( ) {

  typedef double SC;
  typedef int LO;

  // test constructor using triplets (i,j,v)
  std::vector<LO> rowInd;
  std::vector<LO> colInd;
  std::vector<SC> values;
  rowInd.push_back( 1 );
  colInd.push_back( 3 );
  values.push_back( 3.1415 );
  rowInd.push_back( 3 );
  colInd.push_back( 4 );
  values.push_back( 2.71 );
  SparseMatrix<LO, SC> mat( 5, 5, rowInd, colInd, values );
  mat.print( std::cout );

  // test matrix-vector multiplication
  Vector<LO, SC> x( 5 );
  Vector<LO, SC> y( 5 );
  x.setAll( 2.0 );
  y.setAll( 1.0 );
  mat.apply( x, y, true, 2.0, 3.0 );
  y.print( std::cout );

  // test construction using setters
  SparseMatrix<LO, SC> mat2( 10, 10, 10 );
  mat2.set( 2, 2, 2.0 );
  mat2.set( 4, 2, 3.0 );
  mat2.makeCompressed( );
  mat2.add( 2, 2, 3.0 );

  // test scaling
  mat2.scale( 3.0 );
  SparseMatrix<LO, SC> mat3( 10, 10 );
  mat3.set( 1, 4, 2.5 );
  mat3.set( 2, 2, 3.14 );

  // test matrix-matrix addition
  mat2.add( mat3, 2.0 );
  mat2.multiply( mat3, mat2, true, true );
  mat2.print( std::cout );

  // test frobenius norm
  std::cout << "Fro-norm: " << mat2.normFro( ) << std::endl;

  // test copy constructor
  SparseMatrix<LO, SC> mat3_copy( mat3 );
  mat3_copy.print( std::cout );

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

  BlockMatrix<LO, SC> matrix( 2, 2, &nRows[0], &nCols[0] );
  matrix.setBlock( 0, 0, &A11 );
  matrix.setBlock( 0, 1, &A12 );
  matrix.setBlock( 1, 0, nullptr );
  matrix.setBlock( 1, 1, &A22 );
  Vector<LO, SC> x2( 6 );
  Vector<LO, SC> y2( 6 );
  x2.setAll( 1.0 );
  y2.setAll( 1.0 );
  matrix.apply( x2, y2, false, 2.0, 2.0 );
  y2.print( std::cout );
}