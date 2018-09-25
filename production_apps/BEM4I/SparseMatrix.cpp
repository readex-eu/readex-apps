/*!
 * @file    SparseMatrix.cpp
 * @author  Michal Merta 
 * @date    November 22, 2013
 * 
 */

#ifdef SPARSEMATRIX_H

namespace bem4i {

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::SparseMatrix( ) {
  this->nRows = 0;
  this->nCols = 0;
  this->matrix = new Eigen::SparseMatrix<SC, storage, LO>( );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::SparseMatrix(
    const SparseMatrix<LO, SC, storage>& orig
    ) {
  this->nRows = orig.getNRows( );
  this->nCols = orig.getNCols( );
  this->matrix = new Eigen::SparseMatrix<SC, storage, LO>( *( orig.getEigenSparseMatrix( ) ) );
  this->factorized = orig.factorized;
  this->solver = orig.solver;
  this->QRFactorized = orig.QRFactorized;
  this->QRSolver = orig.QRSolver;
}

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::SparseMatrix( LO nRows, LO nCols ) {
  this->nRows = nRows;
  this->nCols = nCols;
  this->matrix = new Eigen::SparseMatrix<SC, storage, LO>( nRows, nCols );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::SparseMatrix( LO nRows, LO nCols, LO reserveSize ) {
  this->nRows = nRows;
  this->nCols = nCols;
  this->matrix = new Eigen::SparseMatrix<SC, storage, LO>( nRows, nCols );
  this->matrix->reserve( reserveSize );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::SparseMatrix(
    LO nRows,
    LO nCols,
    std::vector<LO> &rowInd,
    std::vector<LO> &colInd,
    std::vector<SC> &values ) {

  this->nRows = nRows;
  this->nCols = nCols;
  this->matrix = new Eigen::SparseMatrix<SC, storage, LO>( nRows, nCols );

  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( rowInd.size( ) );
  for ( LO i = 0; i < rowInd.size( ); i++ ) {
    tripletList.push_back( Eigen::Triplet<SC, LO>( rowInd[i], colInd[i], values[i] ) );
  }
  this->matrix->setFromTriplets( tripletList.begin( ), tripletList.end( ) );

  this->matrix->makeCompressed( );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::setFromTriplets(
    LO nRows,
    LO nCols,
    std::vector<LO> &rowInd,
    std::vector<LO> &colInd,
    std::vector<SC> &values
    ) {

  this->nRows = nRows;
  this->nCols = nCols;
  this->resize( nRows, nCols );

  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( rowInd.size( ) );
  for ( LO i = 0; i < rowInd.size( ); i++ ) {
    tripletList.push_back( Eigen::Triplet<SC, LO>( rowInd[i], colInd[i], values[i] ) );
  }
  this->matrix->setFromTriplets( tripletList.begin( ), tripletList.end( ) );

  this->matrix->makeCompressed( );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::setFromTriplets(
    LO nRows,
    LO nCols,
    LO nnz,
    LO *rowInd,
    LO *colInd,
    SC *values
    ) {

  this->nRows = nRows;
  this->nCols = nCols;
  this->resize( nRows, nCols );

  // create Eigen triplets and insert values into matrix (duplicated indices are summed!)
  std::vector<Eigen::Triplet<SC, LO> > tripletList;
  tripletList.reserve( nnz );
  for ( LO i = 0; i < nnz; i++ ) {
    tripletList.push_back( Eigen::Triplet<SC, LO>(
        rowInd[i], colInd[i], values[i] ) );
  }
  this->matrix->setFromTriplets( tripletList.begin( ), tripletList.end( ) );

  this->matrix->makeCompressed( );
  this->factorized = false;
  this->solver = nullptr;
  this->QRFactorized = false;
  this->QRSolver = nullptr;
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::apply(
    const Vector<LO, SC>& x,
    Vector<LO, SC>& y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  // create a map to convert from raw array in x vector to Eigen::Matrix
  typedef Eigen::Map<Eigen::Matrix<SC, Eigen::Dynamic, 1> > Map;
  Map x2map( x.getData( ), x.getLength( ) );
  Map y2map( y.getData( ), y.getLength( ) );

  if ( !transA ) {
    y2map = beta * y2map + alpha * ( *this->matrix ) * x2map;
  } else {
    y2map = beta * y2map + alpha * ( *this->matrix ).transpose( ) * x2map;
  }

}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::multiply( SparseMatrix<LO, SC> &A, SparseMatrix<LO, SC> &B, bool transA,
    bool transB, SC alpha, SC beta ) {

  typedef Eigen::SparseMatrix<SC, storage, LO> EigMat;

  EigMat *tmpA = A.getEigenSparseMatrix( );
  EigMat *tmpB = B.getEigenSparseMatrix( );

  // matrices have to have same storage, therefore we need to evaluate transposition to temp. matrix object
  if ( !transA && !transB ) {
    *( this->matrix ) = beta * ( *this->matrix ) + alpha * ( *tmpA )*( *tmpB );
  } else if ( transA && !transB ) {
    *( this->matrix ) = beta * ( *this->matrix ) + alpha * EigMat( ( *tmpA ).transpose( ) )*( *tmpB );
  } else if ( !transA && transB ) {
    *( this->matrix ) = beta * ( *this->matrix ) + alpha * ( *tmpA ) * EigMat( ( *tmpB ).transpose( ) );
  } else if ( transA && transB ) {
    *( this->matrix ) = beta * ( *this->matrix ) + alpha * EigMat( ( *tmpA ).transpose( ) ) * EigMat( ( *tmpB ).transpose( ) );
  }
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::addToPositions( std::vector<LO> const & rows, std::vector<LO> const & cols, FullMatrix<LO, SC> const & mat ) {
  LO iMax = mat.getNRows( );
  LO jMax = mat.getNCols( );

  for ( LO i = 0; i < iMax; i++ ) {
    for ( LO j = 0; j < jMax; j++ ) {
      this->add( rows[i], cols[j], mat.get( i, j ) );
    }
  }
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::LUSolve( const Vector<LO, SC>& rhs, Vector<LO, SC>& x ) {
  // create a map to convert from raw array in x vector to Eigen::Matrix
  typedef Eigen::Map<Eigen::Matrix<SC, Eigen::Dynamic, 1> > Map;
  Map rhs2map( rhs.getData( ), rhs.getLength( ) );
  Map x2map( x.getData( ), x.getLength( ) );

  if ( !factorized ) {
    this->solver = new Eigen::SparseLU<Eigen::SparseMatrix<SC, storage, LO>, Eigen::COLAMDOrdering<LO> >;
    solver->analyzePattern( *this->matrix );
    solver->factorize( *this->matrix );
    factorized = true;
  }
  x2map = solver->solve( rhs2map );

}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::QRSolve( const Vector<LO, SC>& rhs, Vector<LO, SC>& x ) {
  // create a map to convert from raw array in x vector to Eigen::Matrix
  typedef Eigen::Map<Eigen::Matrix<SC, Eigen::Dynamic, 1> > Map;
  Map rhs2map( rhs.getData( ), rhs.getLength( ) );
  Map x2map( x.getData( ), x.getLength( ) );

  if ( !QRFactorized ) {
    this->QRSolver = new Eigen::SparseQR<Eigen::SparseMatrix<SC, storage, LO>, Eigen::COLAMDOrdering<LO> >;
    this->QRSolver->analyzePattern( *this->matrix );
    this->QRSolver->factorize( *this->matrix );
    QRFactorized = true;
  }
  x2map = this->QRSolver->solve( rhs2map );
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::BiCGStabSolve( Vector<LO, SC>& rhs, Vector<LO, SC>& x ) {
  // create a map to convert from raw array in x vector to Eigen::Matrix
  typedef Eigen::Map<Eigen::Matrix<SC, Eigen::Dynamic, 1> > Map;
  Map rhs2map( rhs.getData( ), rhs.getLength( ) );
  Map x2map( x.getData( ), x.getLength( ) );

  Eigen::BiCGSTAB<Eigen::SparseMatrix<SC, storage> > solver( *this->matrix );
  solver.setTolerance( 1e-12 );

  x2map = solver.solve( rhs2map );
  std::cout << "#iterations: " << solver.iterations( ) << std::endl;
  std::cout << "estimated error: " << solver.error( ) << std::endl;

}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::saveTripletsBin( const std::string &fileName ) const {

  std::ofstream file( fileName.c_str( ), std::ios_base::out | std::ios_base::binary );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  LO xyn[3] = { this->matrix->rows( ), this->matrix->cols( ), this->matrix->nonZeros( ) };

  file.write( (char*) xyn, 3 * sizeof (LO ) );

  for ( LO k = 0; k < this->matrix->outerSize( ); ++k ) {
    typename Eigen::SparseMatrix<SC, storage, LO>::InnerIterator it( *this->matrix, k );
    for (; it; ++it ) {
      LO rc[2] = { it.row( ), it.col( ) };
      file.write( (char*) rc, 2 * sizeof (LO ) );
      SC v = it.value( );
      file.write( (char*) &v, sizeof (SC ) );
    }
  }
  file.close( );
}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::loadTripletsBin( const std::string &fileName ) const {

  std::ifstream file( fileName.c_str( ), std::ios_base::in | std::ios_base::binary );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  LO xyn[3];
  LO rc[2];
  SC v;
  file.read( (char*) xyn, sizeof (LO )*3 );
  this->resize( xyn[0], xyn[1] );
  std::vector < Eigen::Triplet<SC, LO> > trips;
  trips.reserve( xyn[2] );

  for ( int k = 0; k < xyn[2]; ++k ) {
    file.read( (char*) rc, 2 * sizeof (LO ) );
    file.read( (char*) &v, sizeof (SC ) );
    trips.push_back( Eigen::Triplet<SC, LO>( rc[0], rc[1], v ) );
  }

  file.close( );

  this->matrix->setFromTriplets( trips.begin( ), trips.end( ) );
  this->matrix->makeCompressed( );

}

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::loadTriplets(
    const std::string &fileName
    ) {

  std::ifstream file( fileName.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file >> this->nRows;
  file >> this->nCols;
  LO nnz;
  file >> nnz;

  std::vector<LO> rowInd;
  rowInd.resize( nnz );
  std::vector<LO> colInd;
  colInd.resize( nnz );
  std::vector<SC> val;
  val.resize( nnz );

  for ( int k = 0; k < nnz; ++k ) {
    file >> rowInd[ k ];
  }
  for ( int k = 0; k < nnz; ++k ) {
    file >> colInd[ k ];
  }
  for ( int k = 0; k < nnz; ++k ) {
    file >> val[ k ];
  }

  file.close( );

  this->setFromTriplets( this->nRows, this->nCols, rowInd, colInd, val );
}

template<class LO, class SC, int storage>
SparseMatrix<LO, SC, storage>::~SparseMatrix( ) {
  delete this->matrix;
  if ( factorized ) {
    delete this->solver;
  }
  if ( QRFactorized ) {
    delete this->QRSolver;
  }
}

#ifdef HAS_PARDISO

template<class LO, class SC, int storage>
void SparseMatrix<LO, SC, storage>::PARDISOSolve( Vector<LO, SC>& rhs, Vector<LO, SC>& x ) {
  // create a map to convert from raw array in x vector to Eigen::Matrix
  typedef Eigen::Map<Eigen::Matrix<SC, Eigen::Dynamic, 1> > Map;
  Map rhs2map( rhs.getData( ), rhs.getLength( ) );
  Map x2map( x.getData( ), x.getLength( ) );

  Eigen::PardisoLU<Eigen::SparseMatrix<SC, storage, LO> > solver;
  solver.analyzePattern( *this->matrix );
  solver.factorize( *this->matrix );
  x2map = solver.solve( rhs2map );
}
#endif

}

#endif
