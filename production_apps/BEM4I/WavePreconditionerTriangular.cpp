/*!
 * @file    WavePreconditionerTriangular.cpp
 * @author  Michal Merta 
 * @date    March 14, 2014
 * 
 */

#ifdef WAVEPRECONDITIONERTRIANGULAR_H

namespace bem4i {

template<class LO, class SC>
WavePreconditionerTriangular<LO, SC>::WavePreconditionerTriangular( ) {
}

template<class LO, class SC>
WavePreconditionerTriangular<LO, SC>::WavePreconditionerTriangular( const WavePreconditionerTriangular& orig ) {
}

template<class LO, class SC>
WavePreconditionerTriangular<LO, SC>::WavePreconditionerTriangular(
    BlockMatrix<LO, SC>* A,
    int maxLevel
    ) {

  this->sysMatrix = A;
  preconditioner = new BlockMatrix<LO, SC>( *A );
  LO nRows = preconditioner->getNBlockRows( );
  for ( LO i = 0; i < A->getNBlockRows( ) - 1; i++ ) {
    for ( LO j = 0; j < A->getNBlockCols( ); j++ ) {
      if ( i != j ) {
        //  preconditioner->setBlock( i, j, NULL );
      }
    }

  }
  //}
}

template<class LO, class SC>
WavePreconditionerTriangular<LO, SC>::~WavePreconditionerTriangular( ) {
  delete this->preconditioner;
}

template<class LO, class SC>
void WavePreconditionerTriangular<LO, SC>::apply(
    const Vector<LO, SC>& x,
    Vector<LO, SC>& y,
    bool transA, SC alpha,
    SC beta
    ) {

  Vector<LO, SC> rhs( x );

  LO elemCount = 0;
  LO nBlockRows = preconditioner->getNBlockRows( );
  LO nLocalRows = preconditioner->getBlock( 0, 0 )->getNRows( );
  Vector<LO, SC> localUpdate( nLocalRows );

  // solve diagonal systems and update the right hand side
  for ( LO i = 0; i < nBlockRows; i++ ) {
    SparseMatrix<LO, SC>* diagMatrix = ( SparseMatrix<LO, SC>* ) preconditioner->getBlock( i, i );
    Vector<LO, SC> localX( diagMatrix->getNCols( ) );
    Vector<LO, SC> localRHS( diagMatrix->getNRows( ) );
    for ( LO j = 0; j < nLocalRows; j++ ) {
      localRHS.set( j, rhs.get( i * nLocalRows + j ) );
    }
    diagMatrix->LUSolve( localRHS, localX );
    for ( LO j = 0; j < nLocalRows; j++ ) {
      y.set( i * nLocalRows + j, localX.get( j ) );
    }

    // update rhs
//
//    for ( LO j = i + 1; j < nBlockRows; j++ ) {
//      if ( preconditioner->getBlock( j, i ) != NULL ) {
//        preconditioner->getBlock( j, i )->apply( localX, localUpdate );
//        for ( LO k = 0; k < nLocalRows; k++ ) {
//          rhs.add( k + j*nLocalRows, -localUpdate.get( k ) );
//        }
//      }
//    }
  }
}


//LeftIdentityPreconditioner<LO, SC>* M = new LeftIdentityPreconditioner<LO, SC>;
//preconditioner->GMRESSolve( x, y, 1e-14, 10, 999, M );
}



#endif
