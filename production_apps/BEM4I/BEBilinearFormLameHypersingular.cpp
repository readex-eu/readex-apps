/*!
 * @file    BEBilinearFormLameHypersingular.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    April 29, 2015
 * 
 */

#ifdef BEBILINEARFORMLAMEHYPERSINGULAR_H

namespace bem4i {

template <class LO, class SC>
BEBilinearFormLameHypersingular<LO, SC>::BEBilinearFormLameHypersingular( ) {
}

template<class LO, class SC>
BEBilinearFormLameHypersingular<LO, SC>::BEBilinearFormLameHypersingular(
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
  this->mu = this->E / ( 2 * ( 1 + this->nu ) );

  this->quadratureOrderDisjointElems = quadratureOrderDisjointElems;
}

template<class LO, class SC>
void BEBilinearFormLameHypersingular<LO, SC>::assemble(
    FullMatrix<LO, SC>& matrix,
    FullMatrix<LO, SC> & Vlaplace,
    FullMatrix<LO, SC> & Vlame,
    const std::vector< SparseMatrix<LO, SC>* > & T
    ) const {

  if ( this->space->getTestFunctionType( ) == p1 &&
      this->space->getAnsatzFunctionType( ) == p1 ) {

    this->assembleP1P1( matrix, Vlaplace, Vlame, T );
    return;
  }

}
/*
template<class LO, class SC>
void BEBilinearFormLameHypersingular<LO, SC>::assembleP1P1(
    FullMatrix<LO, SC>& matrix,
    FullMatrix<LO, SC> & Vlaplace,
    FullMatrix<LO, SC> & Vlame,
    const std::vector< SparseMatrix<LO, SC>* > & T
    ) const {

  LO dim = 3;

  LO nRowsScalar = this->space->getOuterDOFs( );
  LO nColsScalar = this->space->getInnerDOFs( );

  matrix.resize( dim * nRowsScalar, dim * nColsScalar );

  FullMatrix<LO, SC> *VT =
      new FullMatrix<LO, SC>( Vlame.getNRows( ), T[3]->getNCols( ) );

  VT->multiply( Vlame, *( T[3] ) );
  matrix.multiply( *( T[3] ), *VT, true, false, -4.0 * this->mu, 0.0 );

  FullMatrix<LO, SC> *TVT = new FullMatrix<LO, SC>( nRowsScalar, nColsScalar );

  VT->resize( Vlaplace.getNRows( ), T[2]->getNCols( ) );
  // VT23
  VT->multiply( Vlaplace, *( T[2] ) );

  // T23VT23
  TVT->multiply( *( T[2] ), *VT, true );

  SC val = 0.0;
  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i, j, val );
      matrix.add( i + nRowsScalar, j + nColsScalar, 4.0 * val );
      matrix.add( i + 2 * nRowsScalar, j + 2 * nColsScalar, 4.0 * val );
    }
  }

  // VT13
  VT->multiply( Vlaplace, *( T[1] ) );

  // T13VT13
  TVT->multiply( *( T[1] ), *VT, true );

  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i, j, 4.0 * val );
      matrix.add( i + nRowsScalar, j + nColsScalar, val );
      matrix.add( i + 2 * nRowsScalar, j + 2 * nColsScalar, 4.0 * val );
    }
  }

  // T23VT13
  TVT->multiply( *( T[2] ), *VT, true );
  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i + nRowsScalar, j, 2.0 * val );
      matrix.add( j + nRowsScalar, i, val );
      matrix.add( i, j + nColsScalar, val );
      matrix.add( j, i + nColsScalar, 2.0 * val );
    }
  }

  // VT12
  VT->multiply( Vlaplace, *( T[0] ) );

  // T12VT12
  TVT->multiply( *( T[0] ), *VT, true );

  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i, j, 4.0 * val );
      matrix.add( i + nRowsScalar, j + nColsScalar, 4.0 * val );
      matrix.add( i + 2 * nRowsScalar, j + 2 * nColsScalar, val );
    }
  }

  // T23VT12
  TVT->multiply( *( T[2] ), *VT, true );

  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i + 2 * nRowsScalar, j, -2.0 * val );
      matrix.add( j + 2 * nRowsScalar, i, -val );
      matrix.add( i, j + 2 * nColsScalar, -val );
      matrix.add( j, i + 2 * nColsScalar, -2.0 * val );
    }
  }

  // T13VT12
  TVT->multiply( *( T[1] ), *VT, true );

  for ( LO j = 0; j < nColsScalar; j++ ) {
    for ( LO i = 0; i < nRowsScalar; i++ ) {
      val = TVT->get( i, j );
      matrix.add( i + 2 * nRowsScalar, j + nColsScalar, 2.0 * val );
      matrix.add( j + 2 * nRowsScalar, i + nColsScalar, val );
      matrix.add( i + nRowsScalar, j + 2 * nColsScalar, val );
      matrix.add( j + nRowsScalar, i + 2 * nColsScalar, 2.0 * val );
    }
  }

  matrix.scale( this->mu );

  delete VT;
  delete TVT;

}
 */
///*

template<class LO, class SC>
void BEBilinearFormLameHypersingular<LO, SC>::assembleP1P1(
    FullMatrix<LO, SC> & matrix,
    FullMatrix<LO, SC> & Vlaplace,
    FullMatrix<LO, SC> & Vlame,
    const std::vector< SparseMatrix<LO, SC> * > & T
    ) const {

  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nNodes3 = 3 * nNodes;
  LO nElems = this->space->getMesh( )->getNElements( );
  LO nElems3 = 3 * nElems;

  matrix.resize( nNodes3, nNodes3 );

  FullMatrix<LO, SC> * VT = new FullMatrix<LO, SC>( nElems3, nNodes3 );

  VT->multiply( Vlame, *( T[ 3 ] ) );
  matrix.multiply( *( T[ 3 ] ), *VT, true, false, ( SCVT ) - 4.0 * this->mu,
      (SCVT) 0.0 );

  FullMatrix<LO, SC> * TVT = new FullMatrix<LO, SC>( nNodes, nNodes );

  VT->resize( nElems, nNodes );
  // VT23
  VT->multiply( Vlaplace, *( T[ 2 ] ) );

  // T23VT23
  TVT->multiply( *( T[ 2 ] ), *VT, true );

  SC val = 0.0;

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.add( i, j, val );
      matrix.add( i + nNodes, j + nNodes, (SCVT) 4.0 * val );
      matrix.add( i + 2 * nNodes, j + 2 * nNodes, (SCVT) 4.0 * val );
    }
  }

  // VT13
  VT->multiply( Vlaplace, *( T[ 1 ] ) );

  // T13VT13
  TVT->multiply( *( T[ 1 ] ), *VT, true );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.add( i, j, (SCVT) 4.0 * val );
      matrix.add( i + nNodes, j + nNodes, val );
      matrix.add( i + 2 * nNodes, j + 2 * nNodes, (SCVT) 4.0 * val );
    }
  }

  // T23VT13
  TVT->multiply( *( T[ 2 ] ), *VT, true );

#pragma omp parallel for schedule( dynamic, 32 ) private( val ) 
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.addAtomic( i + nNodes, j, (SCVT) 2.0 * val );
      matrix.addAtomic( j + nNodes, i, val );
      matrix.addAtomic( i, j + nNodes, val );
      matrix.addAtomic( j, i + nNodes, (SCVT) 2.0 * val );
    }
  }

  // VT12
  VT->multiply( Vlaplace, *( T[0] ) );

  // T12VT12
  TVT->multiply( *( T[0] ), *VT, true );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )  
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.add( i, j, (SCVT) 4.0 * val );
      matrix.add( i + nNodes, j + nNodes, (SCVT) 4.0 * val );
      matrix.add( i + 2 * nNodes, j + 2 * nNodes, val );
    }
  }

  // T23VT12
  TVT->multiply( *( T[2] ), *VT, true );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.addAtomic( i + 2 * nNodes, j, ( SCVT ) - 2.0 * val );
      matrix.addAtomic( j + 2 * nNodes, i, -val );
      matrix.addAtomic( i, j + 2 * nNodes, -val );
      matrix.addAtomic( j, i + 2 * nNodes, ( SCVT ) - 2.0 * val );
    }
  }

  // T13VT12
  TVT->multiply( *( T[1] ), *VT, true );

#pragma omp parallel for schedule( dynamic, 32 ) private( val )
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = TVT->get( i, j );
      matrix.addAtomic( i + 2 * nNodes, j + nNodes, (SCVT) 2.0 * val );
      matrix.addAtomic( j + 2 * nNodes, i + nNodes, val );
      matrix.addAtomic( i + nNodes, j + 2 * nNodes, val );
      matrix.addAtomic( j + nNodes, i + 2 * nNodes, (SCVT) 2.0 * val );
    }
  }

  matrix.scale( this->mu );

  delete VT;
  delete TVT;

}
//*/

template <class LO, class SC>
BEBilinearFormLameHypersingular<LO, SC>::BEBilinearFormLameHypersingular(
    const BEBilinearFormLameHypersingular& orig
    ) {
}

template <class LO, class SC>
BEBilinearFormLameHypersingular<LO, SC>::~BEBilinearFormLameHypersingular( ) {
}

}

#endif
