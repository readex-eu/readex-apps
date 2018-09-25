/*!
 * @file    BEBilinearFormLame1Layer.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * 
 */

#ifdef BEBILINEARFORMLAME2LAYER_H

namespace bem4i {

template<class LO, class SC>
BEBilinearFormLame2Layer<LO, SC>::BEBilinearFormLame2Layer( ) {
}

template<class LO, class SC>
BEBilinearFormLame2Layer<LO, SC>::BEBilinearFormLame2Layer(
    const BEBilinearFormLame2Layer& orig
    ) {
}

template<class LO, class SC>
BEBilinearFormLame2Layer<LO, SC>::~BEBilinearFormLame2Layer( ) {
}

template<class LO, class SC>
BEBilinearFormLame2Layer<LO, SC>::BEBilinearFormLame2Layer(
    BESpace<LO, SC>* space,
    int* quadratureOrder,
    quadratureType quadrature,
    int* quadratureOrderDisjointElems
    ) {
  this->space = space;
  this->quadrature = quadrature;

  if ( quadratureOrder ) {
    this->quadratureOrder = quadratureOrder;
  } else {
    switch ( quadrature ) {
      case SauterSchwab:
        this->quadratureOrder = defaultQuadraturesSauterSchwab;
        break;
      case Steinbach:
        this->quadratureOrder = defaultQuadraturesSteinbach;
    }
  }

  // default values = steel
  this->nu = 0.29;
  this->E = 2.0e11;

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
}

template<class LO, class SC>
void BEBilinearFormLame2Layer<LO, SC>::assemble(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & Vlaplace,
    FullMatrix<LO, SC> & Vlame,
    FullMatrix<LO, SC> & Klaplace,
    const std::vector< SparseMatrix<LO, SC>* > & T
    ) const {

  if ( this->space->getAnsatzFunctionType( ) == p1 &&
      this->space->getTestFunctionType( ) == p0 ) {

    this->assembleP0P1( matrix, Vlaplace, Vlame, Klaplace, T );
    return;
  }
}

template<class LO, class SC>
void BEBilinearFormLame2Layer<LO, SC>::assembleP0P1(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & Vlaplace,
    FullMatrix<LO, SC> & Vlame,
    FullMatrix<LO, SC> & Klaplace,
    const std::vector< SparseMatrix<LO, SC> * > & T
    ) const {


  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nNodes3 = 3 * nNodes;
  LO nElems = this->space->getMesh( )->getNElements( );
  LO nElems3 = 3 * nElems;

  matrix.resize( nElems3, nNodes3 );

  matrix.multiply( Vlame, *( T[3] ), false, false,
      this->E / ( (SCVT) 1.0 + this->nu ), (SCVT) 0.0 );

  // add the Laplace matrix K to the diagonal of the result
  SC val = 0.0;

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = Klaplace.get( i, j );
      matrix.add( i, j, val );
      matrix.add( i + nElems, j + nNodes, val );
      matrix.add( i + 2 * nElems, j + 2 * nNodes, val );
    }
  }

  FullMatrix<LO, SC> * VT = new FullMatrix<LO, SC>( nElems, nNodes );
  VT->multiply( Vlaplace, *( T[ 0 ] ) );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = -VT->get( i, j );
      matrix.add( i, j + nNodes, val );
      matrix.add( i + nElems, j, -val );
    }
  }

  VT->multiply( Vlaplace, *( T[ 1 ] ) );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = -VT->get( i, j );
      matrix.add( i, j + 2 * nNodes, val );
      matrix.add( i + 2 * nElems, j, -val );
    }
  }

  VT->multiply( Vlaplace, *( T[ 2 ] ) );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = -VT->get( i, j );
      matrix.add( i + nElems, j + 2 * nNodes, val );
      matrix.add( i + 2 * nElems, j + nNodes, -val );
    }
  }

  delete VT;
}

}

#endif
