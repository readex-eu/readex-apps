/*!
 * @file    WavePreconditioner.cpp
 * @author  Michal Merta 
 * @date    March 14, 2014
 * 
 */

#ifdef WAVEPRECONDITIONER_H

namespace bem4i {

template<class LO, class SC>
WavePreconditioner<LO, SC>::WavePreconditioner( ) {
  this->sysMatrix = nullptr;
  this->preconditioner = nullptr;
  this->sysMatrixMPI = nullptr;
  this->preconditionerMPI = nullptr;
}

template<class LO, class SC>
WavePreconditioner<LO, SC>::WavePreconditioner(
    const WavePreconditioner& orig
    ) {
}

template<class LO, class SC>
WavePreconditioner<LO, SC>::WavePreconditioner(
    BlockMatrix<LO, SC> * A,
    int maxLevel,
    int * maxIters
    ) {

  this->sysMatrix = A;
  preconditioner = new BlockMatrix<LO, SC>( *A );

  this->maxLevel = maxLevel;
  this->level = 1;

  if ( maxIters == nullptr ) {
    int size = maxLevel;
    if ( size < 1 ) {
      size = 1;
    }
    maxIters = new int[size];
    for ( int i = 0; i < size; i++ ) {
      maxIters[i] = defaultMaxIters;
    }
    deleteMaxIters = true;
  } else {
    deleteMaxIters = false;
  }
  this->maxIters = maxIters;

  // set MPI matrices to null
  this->sysMatrixMPI = nullptr;
  this->preconditionerMPI = nullptr;
}

template<class LO, class SC>
WavePreconditioner<LO, SC>::WavePreconditioner(
    MPIBlockMatrix<LO, SC> * A,
    int maxLevel,
    int * maxIters
    ) {

  this->sysMatrixMPI = A;
  preconditionerMPI = new MPIBlockMatrix<LO, SC>( *A );

  this->maxLevel = maxLevel;
  this->level = 1;

  if ( maxIters == nullptr ) {
    int size = maxLevel;
    if ( size < 1 ) {
      size = 1;
    }
    maxIters = new int[size];
    for ( int i = 0; i < size; i++ ) {
      maxIters[i] = defaultMaxIters;
    }
    deleteMaxIters = true;
  } else {
    deleteMaxIters = false;
  }
  this->maxIters = maxIters;

  // set not MPI matrices to null
  this->sysMatrix = nullptr;
  this->preconditioner = nullptr;
}

template<class LO, class SC>
WavePreconditioner<LO, SC>::~WavePreconditioner( ) {
  delete this->preconditioner;
  if ( deleteMaxIters == true ) {
    delete [] maxIters;
  }
}

template<class LO, class SC>
void WavePreconditioner<LO, SC>::apply(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  // call appropriate apply method depending on type of matrix
  if ( preconditioner ) {
    this->applyPreconditioner( *this->preconditioner, x, y );
  } else if ( preconditionerMPI ) {
    this->applyPreconditioner( *this->preconditionerMPI, x, y );
  }

}

template<class LO, class SC>
void WavePreconditioner<LO, SC>::applyPreconditioner(
    BlockMatrix<LO, SC> & matrix,
    const Vector<LO, SC>& x,
    Vector<LO, SC>& y
    ) {

  int nBlockRowsM1 = 0;
  int nBlockRowsM2 = 0;
  int nBlockColsM1 = 0;
  int nBlockColsM2 = 0;
  LO * numsOfRowsM1;
  LO * numsOfRowsM2;
  LO * numsOfColsM1;
  LO * numsOfColsM2;

  nBlockRowsM1 = (int) ceil( preconditioner->getNBlockRows( ) / 2.0 );
  nBlockRowsM2 = preconditioner->getNBlockRows( ) - nBlockRowsM1;
  nBlockColsM1 = (int) ceil( preconditioner->getNBlockCols( ) / 2.0 );
  nBlockColsM2 = preconditioner->getNBlockCols( ) - nBlockColsM1;

  numsOfRowsM1 = new LO[nBlockRowsM1];
  numsOfRowsM2 = new LO[nBlockRowsM2];
  numsOfColsM1 = new LO[nBlockColsM1];
  numsOfColsM2 = new LO[nBlockColsM2];
  for ( int i = 0; i < nBlockRowsM1; i++ ) {
    numsOfRowsM1[i] = preconditioner->getNRowsOfBlock( i );
    numsOfColsM1[i] = preconditioner->getNColsOfBlock( i );
  }

  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    numsOfRowsM2[i] = preconditioner->getNRowsOfBlock( nBlockRowsM1 + i );
    numsOfColsM2[i] = preconditioner->getNColsOfBlock( nBlockColsM1 + i );
  }

  BlockMatrix<LO, SC> * M11 = new BlockMatrix<LO, SC>( nBlockRowsM1,
      nBlockColsM1, numsOfRowsM1, numsOfColsM1 );
  BlockMatrix<LO, SC> * M21 = new BlockMatrix<LO, SC>( nBlockRowsM2,
      nBlockColsM1, numsOfRowsM2, numsOfColsM1 );
  BlockMatrix<LO, SC> * M22 = new BlockMatrix<LO, SC>( nBlockRowsM2,
      nBlockColsM2, numsOfRowsM2, numsOfColsM2 );

  for ( int i = 0; i < nBlockRowsM1; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      M11->setBlock( i, j, preconditioner->getBlock( i, j ) );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      M21->setBlock( i, j, preconditioner->getBlock( nBlockRowsM1 + i, j ) );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM2; j++ ) {
      M22->setBlock( i, j,
          preconditioner->getBlock( nBlockRowsM1 + i, nBlockColsM1 + j ) );
    }
  }

  Vector<LO, SC> * x1 = new Vector<LO, SC>( M11->getNRows( ) );
  Vector<LO, SC> * x2 = new Vector<LO, SC>( M21->getNRows( ) );
  Vector<LO, SC> * y1 = new Vector<LO, SC>( M11->getNRows( ) );
  Vector<LO, SC> * y2 = new Vector<LO, SC>( M21->getNRows( ) );
  y1->setAll( 0.0 );
  y2->setAll( 0.0 );

  for ( LO i = 0; i < M11->getNRows( ); i++ ) {
    x1->set( i, x.get( i ) );
  }

  for ( LO i = 0; i < M21->getNRows( ); i++ ) {
    x2->set( i, x.get( M11->getNRows( ) + i ) );
  }

  //x1->print();x2->print()
  LeftPreconditioner<LO, SC> * precM11; //( M11 );
  LeftPreconditioner<LO, SC> * precM22;

  if ( maxLevel != -1 && level >= maxLevel ) {
    // if we are deep enough do not use preconditioner any more
    precM11 = new LeftIdentityPreconditioner<LO, SC>; //WavePreconditionerTriangular<LO, SC>(M11);//
    precM22 = new LeftIdentityPreconditioner<LO, SC>; //new WavePreconditionerTriangular<LO, SC>(M22);//
  } else {
    //      precM11 = new LeftIdentityPreconditioner<LO, SC>; //WavePreconditionerTriangular<LO, SC>(M11);//
    //      precM22 = new LeftIdentityPreconditioner<LO, SC>;
    precM11 = new WavePreconditioner<LO, SC>( M11, this->maxLevel,
        this->maxIters );
    precM22 = new WavePreconditioner<LO, SC>( M22, this->maxLevel,
        this->maxIters );
    ( ( WavePreconditioner<LO, SC>* ) precM11 )->setCurrentLevel(
        this->level + 1 );
    ( ( WavePreconditioner<LO, SC>* ) precM22 )->setCurrentLevel(
        this->level + 1 );
  }

  if ( M11->getNBlockRows( ) > 1 ) {
    //WavePreconditioner<LO, SC> precM11( M11 );
    M11->FGMRESSolve( *x1, *y1, 1e-9, maxIters[level - 1], 200, precM11 );
    //M11->DGMRESSolve( *x1, *y1, 1e-9, 1000, 200, 200, 2, precM11 );
  } else {
    SparseMatrix<LO, SC>* diagMatrix =
        ( SparseMatrix<LO, SC>* ) M11->getBlock( 0, 0 );
    diagMatrix->LUSolve( *x1, *y1 );
  }

  M21->apply( *y1, *x2, false, -1.0, 1.0 );

  if ( M22->getNBlockRows( ) > 1 ) {
    //WavePreconditioner<LO, SC> precM22( M22 );
    M22->FGMRESSolve( *x2, *y2, 1e-9, maxIters[level - 1], 200, precM22 );
    //M22->DGMRESSolve( *x2, *y2, 1e-9, 1000, 200,200, 2, precM22 );
  } else {
    SparseMatrix<LO, SC>* diagMatrix =
        ( SparseMatrix<LO, SC>* ) M22->getBlock( 0, 0 );
    diagMatrix->LUSolve( *x2, *y2 );
  }


  for ( LO i = 0; i < y1->getLength( ); i++ ) {
    y.set( i, y1->get( i ) );
  }
  for ( LO i = 0; i < y2->getLength( ); i++ ) {
    y.set( M11->getNRows( ) + i, y2->get( i ) );
  }

  delete M11;
  delete M22;
  delete M21;
  delete precM11;
  delete precM22;
  delete x1;
  delete y1;
  delete x2;
  delete y2;
  delete [] numsOfRowsM1;
  delete [] numsOfRowsM2;
  delete [] numsOfColsM1;
  delete [] numsOfColsM2;

}

template<class LO, class SC>
void WavePreconditioner<LO, SC>::applyPreconditioner(
    MPIBlockMatrix<LO, SC> & matrix,
    const Vector<LO, SC> & x,
    Vector<LO, SC>& y
    ) {

  int nBlockRowsM1 = 0;
  int nBlockRowsM2 = 0;
  int nBlockColsM1 = 0;
  int nBlockColsM2 = 0;
  LO * numsOfRowsM1;
  LO * numsOfRowsM2;
  LO * numsOfColsM1;
  LO * numsOfColsM2;
  int * ranksM11, * ranksM21, * ranksM22;

  nBlockRowsM1 = (int) ceil( preconditionerMPI->getNBlockRows( ) / 2.0 );
  nBlockRowsM2 = preconditionerMPI->getNBlockRows( ) - nBlockRowsM1;
  nBlockColsM1 = (int) ceil( preconditionerMPI->getNBlockCols( ) / 2.0 );
  nBlockColsM2 = preconditionerMPI->getNBlockCols( ) - nBlockColsM1;

  numsOfRowsM1 = new LO[nBlockRowsM1];
  numsOfRowsM2 = new LO[nBlockRowsM2];
  numsOfColsM1 = new LO[nBlockColsM1];
  numsOfColsM2 = new LO[nBlockColsM2];
  ranksM11 = new int[nBlockRowsM1 * nBlockColsM1];
  ranksM21 = new int[nBlockRowsM2 * nBlockColsM1];
  ranksM22 = new int[nBlockRowsM2 * nBlockColsM2];

  for ( int i = 0; i < nBlockRowsM1; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      ranksM11[j * nBlockRowsM1 + i] = matrix.getOwner( i, j );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      ranksM21[j * nBlockRowsM2 + i] = matrix.getOwner(
          nBlockRowsM1 + i, j );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM2; j++ ) {
      ranksM22[j * nBlockRowsM2 + i] = matrix.getOwner(
          nBlockRowsM1 + i, nBlockColsM1 + j );
    }
  }

  for ( int i = 0; i < nBlockRowsM1; i++ ) {
    numsOfRowsM1[i] = preconditionerMPI->getNRowsOfBlock( i );
    numsOfColsM1[i] = preconditionerMPI->getNColsOfBlock( i );
  }

  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    numsOfRowsM2[i] = preconditionerMPI->getNRowsOfBlock( nBlockRowsM1 + i );
    numsOfColsM2[i] = preconditionerMPI->getNColsOfBlock( nBlockColsM1 + i );
  }

  MPIBlockMatrix<LO, SC> * M11 = new MPIBlockMatrix<LO, SC>( nBlockRowsM1,
      nBlockColsM1, numsOfRowsM1, numsOfColsM1, ranksM11,
      matrix.getCommunicator( ) );
  MPIBlockMatrix<LO, SC> * M21 = new MPIBlockMatrix<LO, SC>( nBlockRowsM2,
      nBlockColsM1, numsOfRowsM2, numsOfColsM1, ranksM21,
      matrix.getCommunicator( ) );
  MPIBlockMatrix<LO, SC> * M22 = new MPIBlockMatrix<LO, SC>( nBlockRowsM2,
      nBlockColsM2, numsOfRowsM2, numsOfColsM2, ranksM22,
      matrix.getCommunicator( ) );

  for ( int i = 0; i < nBlockRowsM1; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      M11->setBlock( i, j, preconditionerMPI->getBlock( i, j ) );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM1; j++ ) {
      M21->setBlock( i, j, preconditionerMPI->getBlock( nBlockRowsM1 + i, j ) );
    }
  }
  for ( int i = 0; i < nBlockRowsM2; i++ ) {
    for ( int j = 0; j < nBlockColsM2; j++ ) {
      M22->setBlock( i, j, 
          preconditionerMPI->getBlock( nBlockRowsM1 + i, nBlockColsM1 + j ) );
    }
  }

  Vector<LO, SC> * x1 = new Vector<LO, SC>( M11->getNRows( ) );
  Vector<LO, SC> * x2 = new Vector<LO, SC>( M21->getNRows( ) );
  Vector<LO, SC> * y1 = new Vector<LO, SC>( M11->getNRows( ) );
  Vector<LO, SC> * y2 = new Vector<LO, SC>( M21->getNRows( ) );
  y1->setAll( 0.0 );
  y2->setAll( 0.0 );

  for ( LO i = 0; i < M11->getNRows( ); i++ ) {
    x1->set( i, x.get( i ) );
  }

  for ( LO i = 0; i < M21->getNRows( ); i++ ) {
    x2->set( i, x.get( M11->getNRows( ) + i ) );
  }

  //x1->print();x2->print()
  LeftPreconditioner<LO, SC> * precM11; //( M11 );
  LeftPreconditioner<LO, SC> * precM22;

  if ( maxLevel != -1 && level >= maxLevel ) {
    // if we are deep enough do not use preconditionerMPI any more
    precM11 = new LeftIdentityPreconditioner<LO, SC>; //WavePreconditionerTriangular<LO, SC>(M11);//
    precM22 = new LeftIdentityPreconditioner<LO, SC>; //new WavePreconditionerTriangular<LO, SC>(M22);//
  } else {
    precM11 = new WavePreconditioner<LO, SC>( M11, this->maxLevel );
    precM22 = new WavePreconditioner<LO, SC>( M22, this->maxLevel );
    ( ( WavePreconditioner<LO, SC>* ) precM11 )->setCurrentLevel(
        this->level + 1 );
    ( ( WavePreconditioner<LO, SC>* ) precM22 )->setCurrentLevel(
        this->level + 1 );
  }

  MPI_Barrier( preconditionerMPI->getCommunicator( ) );
  if ( M11->getNBlockRows( ) > 1 ) {
    //WavePreconditioner<LO, SC> precM11( M11 );
    M11->FGMRESSolve( *x1, *y1, 1e-5, 10, 200, precM11 );
  } else {
    SparseMatrix<LO, SC>* diagMatrix = 
        ( SparseMatrix<LO, SC>* ) M11->getBlock( 0, 0 );
    diagMatrix->LUSolve( *x1, *y1 );
  }

  MPI_Barrier( preconditionerMPI->getCommunicator( ) );
  M21->apply( *y1, *x2, false, -1.0, 1.0 );

  MPI_Barrier( preconditionerMPI->getCommunicator( ) );
  if ( M22->getNBlockRows( ) > 1 ) {
    //WavePreconditioner<LO, SC> precM22( M22 );
    M22->FGMRESSolve( *x2, *y2, 1e-5, 10, 200, precM22 );
  } else {
    SparseMatrix<LO, SC>* diagMatrix = 
        ( SparseMatrix<LO, SC>* ) M22->getBlock( 0, 0 );
    diagMatrix->LUSolve( *x2, *y2 );
  }
  MPI_Barrier( preconditionerMPI->getCommunicator( ) );

  for ( LO i = 0; i < y1->getLength( ); i++ ) {
    y.set( i, y1->get( i ) );
  }
  for ( LO i = 0; i < y2->getLength( ); i++ ) {
    y.set( M11->getNRows( ) + i, y2->get( i ) );
  }

  MPI_Barrier( preconditionerMPI->getCommunicator( ) );
  delete M11;
  delete M22;
  delete M21;
  delete precM11;
  delete precM22;
  delete x1;
  delete y1;
  delete x2;
  delete y2;
  delete [] numsOfRowsM1;
  delete [] numsOfRowsM2;
  delete [] numsOfColsM1;
  delete [] numsOfColsM2;
  delete [] ranksM11;
  delete [] ranksM21;
  delete [] ranksM22;
}

//template<class LO, class SC>
//void WavePreconditioner<LO, SC>::applyGSPreconditioner(
//    BlockMatrix<LO, SC> & matrix,
//    const Vector<LO, SC>& x,
//    Vector<LO, SC>& y ) {
//
//  // x(k+1)= -(L+D)^(-1)*U*x(k) +  (L+D)^(-1)*b
//
//  LO nIters = 10;
//  //y.print();
//  y.setAll( 0.0 );
//  for ( LO i = 0; i < nIters; i++ ) {
//    Vector<LO, SC> Uxk( x.getLength( ) );
//
//    for ( LO j = 0; j < matrix.getNBlockRows( ) - 1; j++ ) {
//
//      SparseMatrix<LO, SC> *M = ( SparseMatrix<LO, SC>* )
//          matrix.getBlock( j, j + 1 );
//
//      LO n = M->getNRows( );
//      Vector<LO, SC> yLocal( n );
//      Vector<LO, SC> UxkLocal( n );
//      for ( LO k = 0; k < n; k++ ) {
//        yLocal.set( k, y.get( ( j + 1 ) * n + k ) );
//        //std::cout << y.get( (j) * n + k ) << std::endl;
//      }
//      M->apply( yLocal, UxkLocal );
//      for ( LO k = 0; k < n; k++ ) {
//        Uxk.set( j * n + k, -UxkLocal.get( k ) );
//      }
//    }
//    SC omega = 0.5;
//    Uxk.scale( omega );
//    Uxk.add( x, omega );
//    Vector<LO, SC> yy( y );
//
//    applyLInverse( matrix, Uxk, yy );
//    y.scale( 1.0 - omega );
//    y.add( yy, omega );
//
//    // y.print();
//  }
//
//}
//
//template<class LO, class SC>
//void WavePreconditioner<LO, SC>::applyLInverse(
//    BlockMatrix<LO, SC> & matrix,
//    const Vector<LO, SC>& x,
//    Vector<LO, SC>& y ) {
//
//  // apply only the diagonal and lower triangle of the matrix
//  Vector<LO, SC> xLocal( x );
//  LO currentRow = 0;
//
//  for ( LO i = 0; i < matrix.getNBlockRows( ); i++ ) {
//    SparseMatrix<LO, SC> *Mii = ( SparseMatrix<LO, SC>* )
//        matrix.getBlock( i, i );
//    LO ni = Mii->getNRows( );
//
//    Vector<LO, SC> xi( ni );
//    Vector<LO, SC> yi( ni );
//
//    yi.setAll( 0.0 );
//
//    for ( LO j = 0; j < ni; j++ ) {
//      xi.set( j, xLocal.get( ni * i + j ) );
//    }
//    //LeftIdentityPreconditioner<LO, SC>* prc = new LeftIdentityPreconditioner<LO, SC>;
//    //Mii->GMRESSolve( xi, yi, 1e-9, 2000, 200 , prc);
//
//    //if (i==19) {
//    //  Mii->print();
//    //  exit(0);
//    //}
//
//    Mii->LUSolve( xi, yi );
//
//    for ( LO j = 0; j < ni; j++ ) {
//      //      if (isnan(yi.get( j ))) {
//      //        std::cout << i << " " <<j <<std::endl;
//      //      }
//      y.set( ni * i + j, yi.get( j ) );
//    }
//
//    for ( LO j = i + 1; j < matrix.getNBlockRows( ); j++ ) {
//      SparseMatrix<LO, SC> *Mji = ( SparseMatrix<LO, SC>* )
//          matrix.getBlock( j, i );
//      if ( Mji != nullptr ) {
//        Vector<LO, SC> MjiXi( ni );
//        Mji->apply( yi, MjiXi );
//        for ( LO k = 0; k < ni; k++ ) {
//          xLocal.set( j * ni + k, xLocal.get( j * ni + k ) - 0.5 * MjiXi.get( k ) );
//        }
//      }
//    }
//
//    currentRow += ni;
//  }
//
//}


}

#endif
