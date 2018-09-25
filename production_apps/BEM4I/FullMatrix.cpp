/*!
 * @file    FullMatrix.cpp
 * @author  Michal Merta 
 * @date    July 5, 2013
 * 
 */

#ifdef FULLMATRIX_H


namespace bem4i {

  template<class LO, class SC>
  FullMatrix<LO, SC>::FullMatrix( ) {
    this->nRows = 0;
    this->nCols = 0;
    this->data = nullptr;
    this->work = nullptr;
    deleteData = true;
    //factorized = false;
    //reuseFact = false;
  }

  template<class LO, class SC>
  FullMatrix<LO, SC>::FullMatrix(
    const FullMatrix<LO, SC> & orig
    ) {

    this->nRows = orig.nRows;
    this->nCols = orig.nCols;

    // copy matrix data
    this->data = new SC[this->nRows * this->nCols];
    if ( orig.data ) {
      memcpy( this->data, orig.data, this->nRows * this->nCols * sizeof (SC ) );
    } else {
      this->data = nullptr;
    }
    if ( orig.work ) {
      this->work = new SC[ this->nRows ];
    } else {
      this->work = nullptr;
    }
    deleteData = true;
    //factorized = false;
    //reuseFact = false;
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::copy(
    FullMatrix<LO, SC> & copy
    ) const {

    copy.nRows = this->nRows;
    copy.nCols = this->nCols;
    if ( copy.data ) {
      delete [] copy.data;
    }

    if ( this->data ) {
      copy.data = new SC[ this->nRows * this->nCols ];
      memcpy( copy.data, this->data, this->nRows * this->nCols * sizeof ( SC ) );
    } else {
      copy.data = nullptr;
    }

    if ( this->work ) {
      if ( copy.work ) delete [] copy.work;
      copy.work = new SC[ this->nRows ];
    } else {
      copy.work = nullptr;
    }

    copy.deleteData = true;
    //copy.factorized = this->factorized;
    //copy.reuseFact = this->reuseFact;
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::copyToComplex(
    FullMatrix<LO, std::complex<SCVT> > & copy
    ) const {

    copy.nRows = this->nRows;
    copy.nCols = this->nCols;

    if ( copy.data ) {
      delete [] copy.data;
    }

    if ( this->data ) {
      copy.data = new std::complex<SCVT>[ this->nRows * this->nCols ];
      SCVT re, im;
      for ( LO j = 0; j < this->nCols; ++j ) {
        for ( LO i = 0; i < this->nRows; ++i ) {
          re = std::real( this->data[ i + this->nRows * j ] );
          im = std::imag( this->data[ i + this->nRows * j ] );
          copy.data[ i + this->nRows * j ] = std::complex< SCVT >( re, im );
        }
      }
    } else {
      copy.data = nullptr;
    }

    if ( this->work ) {
      if ( copy.work ) delete [] copy.work;
      copy.work = new std::complex< SCVT >[ this->nRows ];
    } else {
      copy.work = nullptr;
    }

    copy.deleteData = true;
    //copy.factorized = this->factorized;
    //copy.reuseFact = this->reuseFact;

  }

  template<class LO, class SC>
  FullMatrix<LO, SC>::FullMatrix(
    LO nRows,
    LO nCols,
    bool zeroOut,
    bool allocWork
    ) {

    this->nRows = nRows;
    this->nCols = nCols;
    this->data = new SC[nRows * nCols];
    deleteData = true;

    if ( zeroOut ) {
      this->setAll( 0.0 );
    }

    // allocating workspace for some LAPACK methods
    if ( allocWork ) {
      work = new SC[nRows];
    } else {
      work = nullptr;
    }
    //factorized = false;
    //reuseFact = false;
  }

  template<class LO, class SC>
  FullMatrix<LO, SC>::FullMatrix(
    LO nRows,
    LO nCols,
    SC * data,
    bool allocWork
    ) {

    this->nRows = nRows;
    this->nCols = nCols;

    this->data = data;
    deleteData = false;

    // allocating workspace for some LAPACK methods
    if ( allocWork ) {
      work = new SC[nRows];
    } else {
      work = nullptr;
    }
    //factorized = false;
    //reuseFact = false;

  }

  template<class LO, class SC>
  FullMatrix<LO, SC>::~FullMatrix( ) {

    if ( data && deleteData ) {
      delete [] data;
    }
    if ( work ) {
      delete [] work;
    }
    //  if ( factorized && reuseFact ) {
    //    delete [] ipiv;
    //  }
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::resize(
    LO nRows,
    LO nCols,
    bool zeroOut
    ) {

    this->nRows = nRows;
    this->nCols = nCols;

    if ( data ) {
      delete [] data;
    }

    data = new SC[nRows * nCols];

    if ( zeroOut ) {
      this->setAll( 0.0 );
    }

    if ( work ) {
      delete [] work;
      work = new SC[nRows];
    }

    // current factorization no longer useful
    //  if ( factorized && reuseFact ) {
    //    delete [] ipiv;
    //    factorized = false;
    //  }
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::addToPositions(
    std::vector<LO> const & rows,
    std::vector<LO> const & cols,
    FullMatrix<LO, SC> const & mat
    ) {

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );
    SC val;

    for ( LO j = 0; j < jMax; j++ ) {
      for ( LO i = 0; i < iMax; i++ ) {
        val = this->get( rows[i], cols[j] );
        this->set( rows[i], cols[j], val + mat.get( i, j ) );
      }
    }
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::addToPositionsAtomic(
    std::vector<LO> const & rows,
    std::vector<LO> const & cols,
    FullMatrix<LO, SC> const & mat
    ) {

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );

    SC * matData = mat.getData( );

    for ( LO j = 0; j < jMax; ++j ) {
      for ( LO i = 0; i < iMax; ++i ) {
#pragma omp atomic update
        this->data[ rows[ i ] + this->nRows * cols[ j ] ] +=
          matData[ i + j * iMax ];
      }
    }
  }

  template<>
  void FullMatrix<int, std::complex< double > >::addToPositionsAtomic(
    std::vector<int> const & rows,
    std::vector<int> const & cols,
    FullMatrix<int, std::complex< double > > const & mat
    ) {

    typedef int LO;
    typedef std::complex< double > SC;
    typedef typename SC::value_type SCVT;

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );

    SC * matData = mat.getData( );

    for ( LO j = 0; j < jMax; ++j ) {
      for ( LO i = 0; i < iMax; ++i ) {
#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) ] +=
          std::real( matData[ i + j * iMax ] );

#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) + 1 ] +=
          std::imag( matData[ i + j * iMax ] );
      }
    }
  }

  template<>
  void FullMatrix<int, std::complex< float > >::addToPositionsAtomic(
    std::vector<int> const & rows,
    std::vector<int> const & cols,
    FullMatrix<int, std::complex< float > > const & mat
    ) {

    typedef int LO;
    typedef std::complex< float > SC;
    typedef typename SC::value_type SCVT;

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );

    SC * matData = mat.getData( );

    for ( LO j = 0; j < jMax; ++j ) {
      for ( LO i = 0; i < iMax; ++i ) {
#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) ] +=
          std::real( matData[ i + j * iMax ] );

#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) + 1 ] +=
          std::imag( matData[ i + j * iMax ] );
      }
    }
  }

  template<>
  void FullMatrix<long, std::complex< double > >::addToPositionsAtomic(
    std::vector<long> const & rows,
    std::vector<long> const & cols,
    FullMatrix<long, std::complex< double > > const & mat
    ) {

    typedef long LO;
    typedef std::complex< double > SC;
    typedef typename SC::value_type SCVT;

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );

    SC * matData = mat.getData( );

    for ( LO j = 0; j < jMax; ++j ) {
      for ( LO i = 0; i < iMax; ++i ) {
#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) ] +=
          std::real( matData[ i + j * iMax ] );

#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) + 1 ] +=
          std::imag( matData[ i + j * iMax ] );
      }
    }
  }

  template<>
  void FullMatrix<long, std::complex< float > >::addToPositionsAtomic(
    std::vector<long> const & rows,
    std::vector<long> const & cols,
    FullMatrix<long, std::complex< float > > const & mat
    ) {

    typedef long LO;
    typedef std::complex< float > SC;
    typedef typename SC::value_type SCVT;

    LO iMax = mat.getNRows( );
    LO jMax = mat.getNCols( );

    SC * matData = mat.getData( );

    for ( LO j = 0; j < jMax; ++j ) {
      for ( LO i = 0; i < iMax; ++i ) {
#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) ] +=
          std::real( matData[ i + j * iMax ] );

#pragma omp atomic update
        reinterpret_cast < SCVT * > ( this->data )
          [ 2 * ( rows[ i ] + this->nRows * cols[ j ] ) + 1 ] +=
          std::imag( matData[ i + j * iMax ] );
      }
    }
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::addToPositions(
    std::vector<LO> const & rows,
    std::vector<LO> const & cols,
    std::vector<SC> const & values
    ) {

    SC val;

    auto it = rows.begin( );
    auto it2 = cols.begin( );
    auto it3 = values.begin( );

    for (; it != rows.end( ); ++it, ++it2, ++it3 ) {
      val = this->get( *it, *it2 );
      this->set( *it, *it2, val + *it3 );

    }
  }

  // GENERIC IMPLEMENTATION //

  template<class LO, class SC>
  void FullMatrix<LO, SC>::getCol(
    LO idx,
    SC *outCol ) const {

    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (SC ) );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::getRow(
    LO idx,
    SC *outCol
    ) const {

    for ( LO i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::scale( SC alpha ) {
    LO dataSize = this->nRows * this->nCols;
    LO incx = 1;
    dscal_( &dataSize, &alpha, data, &incx );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::conjugate( ) {

    SC val;

#pragma omp parallel for collapse( 2 ) private( val )
    for ( LO j = 0; j < this->nCols; ++j ) {
      for ( LO i = 0; i < this->nRows; ++i ) {
        val = ( SC ) std::conj( this->get( i, j ) );
        this->set( i, j, val );
      }
    }
  }

  template<class LO, class SC>
  SC FullMatrix<LO, SC>::norm1( ) const {
    char norm = '1';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<class LO, class SC>
  SC FullMatrix<LO, SC>::normFro( ) {
    char norm = 'F';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<class LO, class SC>
  SC FullMatrix<LO, SC>::normI( ) {
    char norm = 'I';
    if ( work == nullptr ) {
      work = new SC[this->nRows];
    }
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::add( FullMatrix<LO, SC> &A, SC alpha ) {
    LO size = this->nRows * this->nCols;
    LO inc = 1;
    daxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::add(
    SparseMatrix<LO, SC, Eigen::ColMajor> &A,
    SC alpha
    ) {

    Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> * Ae = A.getEigenSparseMatrix( );

    for ( LO j = 0; j < Ae->outerSize( ); ++j ) {
      typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
      it( *Ae, j );
      for (; it; ++it ) {
        this->add( it.row( ), it.col( ), alpha * it.value( ) );
      }
    }

  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::multiply(
    FullMatrix<LO, SC> &A,
    FullMatrix<LO, SC> &B,
    bool transA,
    bool transB,
    SC alpha,
    SC beta
    ) {

    char transAChar, transBChar;
    LO lda, ldb;
    LO nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    dgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::multiply(
    FullMatrix<LO, SC> &A,
    SparseMatrix<LO, SC, Eigen::ColMajor> &B,
    bool transA,
    bool transB,
    SC alpha,
    SC beta
    ) {

    // todo: implement all versions
    if ( transA != false || transB != false ) {
      std::cout << "Not implemented!" << std::endl;
      return;
    }

    if ( beta != 1.0 ) {
      this->scale( beta );
    }

    Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> * Be = B.getEigenSparseMatrix( );

#pragma omp parallel for shared( A, B, transA, transB, alpha, beta, Be ) schedule(dynamic,1)
    for ( LO i = 0; i < this->nCols; ++i ) {
      typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
      it( *Be, i );
      for (; it; ++it ) {
        for ( LO j = 0; j < this->nRows; ++j ) {
          this->add( j, i, alpha * it.value( ) * A.get( j, it.row( ) ) );
        } // over rows in A's ith column
      } // over rows in B's jth column    
    } // over columns of this

  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::multiply(
    SparseMatrix<LO, SC, Eigen::ColMajor> &A,
    FullMatrix<LO, SC> &B,
    bool transA,
    bool transB,
    SC alpha,
    SC beta
    ) {

    if ( beta != 1.0 ) {
      this->scale( beta );
    }

    Eigen::SparseMatrix<SC, Eigen::ColMajor, LO> * Ae = A.getEigenSparseMatrix( );

    if ( transA == true && transB == false ) {

#pragma omp parallel for schedule( dynamic, 32 )\
shared( A, B, transA, transB, alpha, beta, Ae )
      for ( LO i = 0; i < this->nCols; ++i ) {
        for ( LO j = 0; j < this->nRows; ++j ) {
          typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
          it( *Ae, j );
          for (; it; ++it ) {
            this->add( j, i, alpha * it.value( ) * B.get( it.row( ), i ) );
          } // over rows in A's jth column
        } // over rows of this
      } // over columns of this

    } else if ( transA == false && transB == false ) {

#pragma omp parallel for schedule( dynamic, 32 )\
shared( A, B, transA, transB, alpha, beta, Ae )    
      for ( LO i = 0; i < this->nCols; ++i ) {
        for ( LO j = 0; j < A.getNCols( ); ++j ) {
          typename Eigen::SparseMatrix<SC, Eigen::ColMajor, LO>::InnerIterator
          it( *Ae, j );
          for (; it; ++it ) {
            this->add( it.row( ), i, alpha * it.value( ) * B.get( j, i ) );
          } // over rows in A's jth column
        } // over columns of A
      } // over columns of this 

    } else {
      // todo: implement all versions
      std::cout << "Not implemented!" << std::endl;
      return;
    }

  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::multiply(
    FullMatrix<LO, SC> &A,
    FullMatrix<LO, SC> &B,
    LO ARows,
    LO ACols,
    LO BRows,
    LO BCols,
    bool transA,
    bool transB,
    SC alpha,
    SC beta
    ) {

    char transAChar, transBChar;
    LO lda, ldb;
    LO nACols, nARows, nBCols, nBRows;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    dgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;

    dgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::applyMIC(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta,
    int device
    ) {
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::applySubmatrix(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    LO ARows,
    LO ACols,
    bool transA,
    SC alpha,
    SC beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;

    dgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::LUSolve(
    Vector<LO, SC> & x,
    LO nRhs
    ) {

    LO ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    LO info = 0;
    char trans = 'N';

    // factorize a matrix
    LO *ipiv = new LO[ipivLength];
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );
    delete [] ipiv;
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::LUSolve(
    FullMatrix<LO, SC> & x,
    LO nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    LO ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    LO info = 0;
    char trans = 'N';

    // factorize a matrix
    LO *ipiv = new LO[ipivLength];
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::CholeskiSolve(
    Vector<LO, SC> & x,
    LO nRhs
    ) {

    LO info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::CholeskiSolve(
    FullMatrix<LO, SC> & x,
    LO nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    LO info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::backward(
    Vector<LO, SC> &x,
    LO nRhs,
    LO n
    ) {

    LO ipivLength;
    LO nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    LO *ipiv = new LO[ipivLength];
    for ( LO i = 0; i < nr; i++ ) {
      ipiv[ i ] = i + 1;
    }
    LO info;
    char trans = 'N';

    // solve a system
    dgetrs_( &trans, &nr, &nRhs, this->data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::schurOfHessenberg(
    FullMatrix<LO, SC> &Z,
    std::complex<SCVT> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    int ilo = 1; //
    int ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( LO i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    int info = 0;
    int workSize = 3 * this->getNRows( );
    SCVT *localWork = new SCVT[workSize];
    SCVT *wr = new SCVT[this->getNRows( )];
    SCVT *wi = new SCVT[this->getNRows( )];

    dhseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, wr, wi, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    for ( LO i = 0; i < this->getNRows( ); i++ ) {
      eigs[i] = std::complex<SCVT>( wr[i], wi[i] );
    }

    delete [] localWork;
    delete [] wr;
    delete [] wi;

    //   char jobvs = 'V';
    //  char sort = 'N';
    //
    //  int sdim = 0;
    // 
    //
    //  // if Z has not a good size, resize it
    //  if ( ( Z.getNRows( ) != this->getNRows( ) )
    //      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
    //    Z.resize( this->getNRows( ), this->getNCols( ) );
    //  }
    //  Z.setAll( 0.0 );
    //  for ( LO i = 0; i < this->getNRows( ); i++ ) {
    //    Z.set( i, i, 1.0 );
    //  }
    //
    //  int info = 0;
    //  int workSize = 3 * this->getNRows( );
    //  SCVT *localWork = new SCVT[workSize];
    //  SCVT *wr = new SCVT[this->getNRows( )];
    //  SCVT *wi = new SCVT[this->getNRows( )];
    //  
    //  dgees_( &jobvs, &sort, nullptr, &this->nRows, this->data,
    //      &this->nRows, &sdim, wr, wi, Z.data, &Z.nRows, localWork,
    //      &workSize, nullptr, &info );
    //
    //  for ( LO i = 0; i < this->getNRows( ); i++ ) {
    //    eigs[i] = std::complex<SCVT>( wr[i], wi[i] );
    //  }
    //
    //  delete [] localWork;
    //  delete [] wr;
    //  delete [] wi;

  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( LO i = 0; i < this->nRows; i++ ) {
      for ( LO j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }

  template<class LO, class SC>
  void FullMatrix<LO, SC>::apply(
    Vector<LO, SC> const &x,
    std::vector<LO>& indicesX,
    Vector<LO, SC> &y,
    std::vector<LO> &indicesY,
    bool transA,
    SC alpha,
    SC beta
    ) {

    Vector<LO, SC> localX( indicesX.size( ), false );
    Vector<LO, SC> localY( indicesY.size( ), false );

    for ( LO i = 0; i < indicesX.size( ); i++ ) {
      localX.set( i, x.get( indicesX[i] ) );
    }

    for ( LO i = 0; i < indicesY.size( ); i++ ) {
      localY.set( i, y.get( indicesY[i] ) );
    }

    this->apply( localX, localY, transA, alpha, beta );

    for ( LO i = 0; i < indicesY.size( ); i++ ) {
      y.set( indicesY[i], localY.get( i ) );
    }

    //  for ( LO i = 0; i < this->nCols; i++ ) {
    //    for ( LO j = 0; j < this->nRows; j++ ) {
    //      val = y.get( indicesY[j] );
    //      y.set( indicesY[j], val + this->get( j, i ) * x.get( indicesX[i] ) );
    //    }
    //  }
  }

#ifdef BLAS_INT

  template<>
  void FullMatrix< int, double >::applyMIC(
    Vector< int, double > const & x,
    Vector< int, double > & y,
    bool transA,
    double alpha,
    double beta,
    int device
    ) {

    typedef int LO;
    typedef double SC;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SC * data = this->data;
    SC * xdata = x.data;
    SC * ydata = y.data;

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : transAChar, nRows, nCols, alpha, one, beta ) 
// defaultmap( to : scalar )
#endif
    {
      dgemv_( &transAChar, &nRows, &nCols, &alpha, data, &nRows, xdata, &one,
        &beta, ydata, &one );
    }
  }

  template<>
  void FullMatrix<int, double>::getCol( int idx, double *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (double ) );
  }

  template<>
  void FullMatrix<int, double>::getRow( int idx, double *outCol ) const {
    for ( int i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<int, double>::scale( double alpha ) {
    int dataSize = this->nRows * this->nCols;
    int incx = 1;
    dscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  double FullMatrix<int, double>::norm1( ) const {
    char norm = '1';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  double FullMatrix<int, double>::normFro( ) {
    char norm = 'F';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  double FullMatrix<int, double>::normI( ) {
    char norm = 'I';
    if ( work == nullptr ) {
      work = new double[this->nRows];
    }
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  void FullMatrix<int, double>::add( FullMatrix<int, double> &A, double alpha ) {
    int size = this->nRows * this->nCols;
    int inc = 1;
    daxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<int, double>::multiply( FullMatrix<int, double> &A, FullMatrix<int, double> &B, bool transA,
    bool transB, double alpha, double beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nCols;
    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    dgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, double>::multiply(
    FullMatrix<int, double> &A,
    FullMatrix<int, double> &B,
    int ARows,
    int ACols,
    int BRows,
    int BCols,
    bool transA,
    bool transB,
    double alpha,
    double beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nACols, nARows, nBCols;
    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';

      nBCols = BRows;
    } else {
      transBChar = 'N';
      nBCols = BCols;
    }

    dgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, double>::apply(
    Vector<int, double> const &x,
    Vector<int, double> &y,
    bool transA,
    double alpha,
    double beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    dgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows, x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, double>::applySubmatrix(
    Vector<int, double> const &x,
    Vector<int, double> &y,
    int ARows,
    int ACols,
    bool transA,
    double alpha,
    double beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    dgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  //template<>
  //void FullMatrix<int, double>::solve( Vector<int, double> &x, int nRhs ) {
  //
  //  char trans = 'N';
  //  int info = 0;
  //
  //  // factorize a matrix
  //  if ( !factorized || !reuseFact ) {
  //    int ipivLength;
  //    if ( this->nRows < this->nCols ) {
  //      ipivLength = this->nRows;
  //    } else {
  //      ipivLength = this->nCols;
  //    }
  //    ipiv = new int[ipivLength];
  //    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );
  //    factorized = true;
  //  } else {
  //
  //  }
  //
  //  // solve a system
  //  dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );
  //
  //  if ( !reuseFact ) {
  //    delete [] ipiv;
  //  }
  //}

  template<>
  void FullMatrix<int, double>::LUSolve(
    Vector<int, double> & x,
    int nRhs
    ) {

    char trans = 'N';
    int info = 0;

    // factorize a matrix
    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int * ipiv = new int[ ipivLength ];
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, double>::LUSolve(
    FullMatrix<int, double> & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int info = 0;
    char trans = 'N';

    // factorize a matrix
    int * ipiv = new int[ipivLength];
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, double>::CholeskiSolve(
    Vector<int, double> & x,
    int nRhs
    ) {

    int info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, double>::CholeskiSolve(
    FullMatrix<int, double> & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );

  }

  template<>
  void FullMatrix<int, double>::backward( Vector<int, double> &x, int nRhs, int n ) {

    int ipivLength;
    int nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    for ( int i = 0; i < nr; i++ ) {
      ipiv[i] = i + 1;
    }
    int info;
    char trans = 'N';

    // solve a system
    dgetrs_( &trans, &nr, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, double>::schurOfHessenberg(
    FullMatrix<int, double> &Z,
    std::complex<double> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    int ilo = 1; //
    int ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( int i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    int info = 0;
    int workSize = 3 * this->getNRows( );
    double *localWork = new double[workSize];
    double *wr = new double[this->getNRows( )];
    double *wi = new double[this->getNRows( )];

    dhseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, wr, wi, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    for ( int i = 0; i < this->getNRows( ); i++ ) {
      eigs[i] = std::complex<double>( wr[i], wi[i] );
    }

    delete [] localWork;
    delete [] wr;
    delete [] wi;
  }

  template<>
  void FullMatrix<int, double>::eigs( double *eigenvectors, double *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    double *up_package = new double[this->nRows * this->nRows]; //(int) this->nRows*(this->nRows+1)/2
    double *workspace = new double[8 * this->nRows];
    int lwork = 8 * this->nRows;
    int *iwork = new int[5 * this->nRows];
    int *ifail = new int[this->nRows];
    int info;

    int idx = 0;
    int m;

    // pack matrix to upper form
    for ( int j = 0; j < this->nRows; j++ ) {
      for ( int i = 0; i < this->nRows; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    double eps = ( double ) EPS;
    dsyevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, eigenvalues, eigenvectors, &this->nRows, workspace, &lwork, iwork, ifail, &info );
    delete [] up_package;
    delete [] workspace;
    delete [] iwork;
    delete [] ifail;
  }

  template<>
  void FullMatrix<int, double>::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( int i = 0; i < this->nRows; i++ ) {
      for ( int j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }


  // <INT, SINGLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< int, float >::applyMIC(
    Vector< int, float > const & x,
    Vector< int, float > & y,
    bool transA,
    float alpha,
    float beta,
    int device
    ) {

    typedef int LO;
    typedef float SC;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SC * data = this->data;
    SC * xdata = x.data;
    SC * ydata = y.data;

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : transAChar, nRows, nCols, alpha, one, beta )
#endif
    {
      sgemv_( &transAChar, &nRows, &nCols, &alpha, data, &nRows, xdata, &one,
        &beta, ydata, &one );
    }
  }

  template<>
  void FullMatrix<int, float>::getCol( int idx, float *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (float ) );
  }

  template<>
  void FullMatrix<int, float>::getRow( int idx, float *outCol ) const {
    for ( int i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<int, float>::scale( float alpha ) {
    int dataSize = this->nRows * this->nCols;
    int incx = 1;
    sscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  float FullMatrix<int, float>::norm1( ) const {
    char norm = '1';
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  float FullMatrix<int, float>::normFro( ) {
    char norm = 'F';
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  float FullMatrix<int, float>::normI( ) {
    char norm = 'I';
    if ( work == nullptr ) {
      work = new float[this->nRows];
    }
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  void FullMatrix<int, float>::add( FullMatrix<int, float> &A, float alpha ) {
    int size = this->nRows * this->nCols;
    int inc = 1;
    saxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<int, float>::multiply( FullMatrix<int, float> &A, FullMatrix<int, float> &B, bool transA,
    bool transB, float alpha, float beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    sgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, float>::multiply(
    FullMatrix<int, float> &A,
    FullMatrix<int, float> &B,
    int ARows,
    int ACols,
    int BRows,
    int BCols,
    bool transA,
    bool transB,
    float alpha,
    float beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    sgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, float>::apply(
    Vector<int, float> const & x,
    Vector<int, float> & y,
    bool transA,
    float alpha,
    float beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    sgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows, x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, float>::applySubmatrix(
    Vector<int, float> const &x,
    Vector<int, float> &y,
    int ARows,
    int ACols,
    bool transA,
    float alpha,
    float beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    sgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, float>::LUSolve( Vector<int, float> &x, int nRhs ) {

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    int info;
    char trans = 'N';

    // d a matrix (LU)
    sgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    sgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, float>::LUSolve(
    FullMatrix<int, float> & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int info = 0;
    char trans = 'N';

    // factorize a matrix
    int *ipiv = new int[ipivLength];
    sgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    sgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, float>::CholeskiSolve(
    Vector<int, float> & x,
    int nRhs
    ) {

    int info = 0;
    char uplo = 'L';

    spotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    spotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, float>::CholeskiSolve(
    FullMatrix<int, float> & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int info = 0;
    char uplo = 'L';

    spotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    spotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, float>::backward( Vector<int, float> &x, int nRhs, int n ) {

    int ipivLength;
    int nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    for ( int i = 0; i < nr; i++ ) {
      ipiv[i] = i + 1;
    }
    int info;
    char trans = 'N';

    // solve a system
    sgetrs_( &trans, &nr, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, float>::schurOfHessenberg(
    FullMatrix<int, float> &Z,
    std::complex<float> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    int ilo = 1; //
    int ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( int i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    int info = 0;
    int workSize = 3 * this->getNRows( );
    float *localWork = new float[workSize];
    float *wr = new float[this->getNRows( )];
    float *wi = new float[this->getNRows( )];

    shseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, wr, wi, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    for ( int i = 0; i < this->getNRows( ); i++ ) {
      eigs[i] = std::complex<double>( wr[i], wi[i] );
    }

    delete [] localWork;
    delete [] wr;
    delete [] wi;
  }

  template<>
  void FullMatrix<int, float>::eigs( float *eigenvectors, float *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    float *up_package = new float[this->nRows * this->nRows];
    float *workspace = new float[8 * this->nRows];
    int lwork = 8 * this->nRows;
    int *iwork = new int[5 * this->nRows];
    int *ifail = new int[this->nRows];
    int info;

    int idx = 0;
    int m;

    // pack matrix to upper form
    for ( int j = 0; j < this->nRows; j++ ) {
      for ( int i = 0; i < this->nRows; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    float eps = ( float ) EPS;
    ssyevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, eigenvalues, eigenvectors, &this->nRows, workspace, &lwork, iwork, ifail, &info );
    delete [] up_package;
    delete [] workspace;
    delete [] iwork;
    delete [] ifail;
  }

  template<>
  void FullMatrix<int, float>::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( int i = 0; i < this->nRows; i++ ) {
      for ( int j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }


  // <INT, COMPLEX SINGLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< int, std::complex< float > >::applyMIC(
    Vector< int, std::complex< float > > const & x,
    Vector< int, std::complex< float > > & y,
    bool transA,
    std::complex< float > alpha,
    std::complex< float > beta,
    int device
    ) {

    typedef int LO;
    typedef std::complex< float > SC;
    typedef float SCVT;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SCVT * alpha_c = reinterpret_cast < SCVT * > ( &alpha );
    SCVT * beta_c = reinterpret_cast < SCVT * > ( &beta );
    SCVT * data_c = reinterpret_cast < SCVT * > ( this->data );
    SCVT * ydata_c = reinterpret_cast < SCVT * > ( y.data );
    SCVT * xdata_c = reinterpret_cast < SCVT * > ( x.data );

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : alpha_c[ 0 : 2 ], beta_c[ 0 : 2 ] ) \
map( to : transAChar, nRows, nCols, one )
#endif
    {
      cgemv_( &transAChar, &nRows, &nCols,
        reinterpret_cast < SC * > ( alpha_c ),
        reinterpret_cast < SC * > ( data_c ),
        &nRows,
        reinterpret_cast < SC * > ( xdata_c ),
        &one,
        reinterpret_cast < SC * > ( beta_c ),
        reinterpret_cast < SC * > ( ydata_c ),
        &one );
    }
  }

  template<>
  void FullMatrix<int, std::complex<float> >::getRow( int idx, std::complex<float> *outCol ) const {
    for ( int i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<int, std::complex<float> >::scale( std::complex<float> alpha ) {
    int dataSize = this->nRows * this->nCols;
    int incx = 1;
    cscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  std::complex<float> FullMatrix<int, std::complex<float> >::norm1( ) const {
    char norm = '1';
    float *work2 = new float[this->nRows];
    float norm1 = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( norm1, 0.0 );

  }

  template<>
  std::complex<float> FullMatrix<int, std::complex<float> >::normFro( ) {
    char norm = 'F';
    float *work2 = new float[this->nRows];
    float normF = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( normF, 0.0 );
  }

  template<>
  std::complex<float> FullMatrix<int, std::complex<float> >::normI( ) {
    char norm = 'I';
    float *work2 = new float[this->nRows];
    float normI = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( normI, 0.0 );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::add( FullMatrix<int, std::complex<float> > &A, std::complex<float> alpha ) {
    int size = this->nRows * this->nCols;
    int inc = 1;
    caxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::multiply( FullMatrix<int, std::complex<float> > &A, FullMatrix<int, std::complex<float> > &B, bool transA,
    bool transB, std::complex<float> alpha, std::complex<float> beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    cgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::multiply(
    FullMatrix<int, std::complex<float> > &A,
    FullMatrix<int, std::complex<float> > &B,
    int ARows,
    int ACols,
    int BRows,
    int BCols,
    bool transA,
    bool transB,
    std::complex<float> alpha,
    std::complex<float> beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    cgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::apply(
    Vector<int, std::complex<float> > const & x,
    Vector<int, std::complex<float> > & y,
    bool transA,
    std::complex<float> alpha,
    std::complex<float> beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    cgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::applySubmatrix(
    Vector<int, std::complex<float> > const &x,
    Vector<int, std::complex<float> > &y,
    int ARows,
    int ACols,
    bool transA,
    std::complex<float> alpha,
    std::complex<float> beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    cgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::eigs( std::complex<float> *eigenvectors, std::complex<float> *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    std::complex<float> *up_package = new std::complex<float>[this->nRows * this->nRows];
    std::complex<float> *workspace = new std::complex<float>[8 * this->nRows];
    int lwork = 8 * this->nRows;
    int *iwork = new int[5 * this->nRows];
    float *rwork = new float[7 * this->nRows];
    int *ifail = new int[this->nRows];
    int info;
    float *w = new float[this->nRows];

    int idx = 0;
    int m;

    // pack matrix to upper form
    // pack matrix to upper form
    for ( int j = 0; j < this->nRows; j++ ) {
      for ( int i = 0; i <= j; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    float eps = ( float ) EPS;
    cheevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, w, eigenvectors, &this->nRows, workspace, &lwork, rwork, iwork, ifail, &info );

    for ( int i = 0; i < this->nRows; i++ ) {
      eigenvalues[i] = w[i];
    }
    delete [] up_package;
    delete [] workspace;
    delete [] w;
    delete [] iwork;
    delete [] ifail;
    delete [] rwork;
  }

  template<>
  void FullMatrix<int, std::complex<float> >::LUSolve( Vector<int, std::complex<float> > &x, int nRhs ) {

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    int info;
    char trans = 'N';

    // factorize a matrix (LU)
    cgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    cgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<float> >::LUSolve(
    FullMatrix<int, std::complex<float> > & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int info = 0;
    char trans = 'N';

    // factorize a matrix
    int * ipiv = new int[ipivLength];
    cgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    cgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<float> >::CholeskiSolve(
    Vector<int, std::complex<float> > & x,
    int nRhs
    ) {

    int info = 0;
    char uplo = 'L';

    cpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    cpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, std::complex<float> >::CholeskiSolve(
    FullMatrix<int, std::complex<float> > & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int info = 0;
    char uplo = 'L';

    cpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    cpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, std::complex<float>>::backward( Vector<int, std::complex<float> > &x, int nRhs, int n ) {

    int ipivLength;
    int nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    for ( int i = 0; i < nr; i++ ) {
      ipiv[i] = i + 1;
    }
    int info;
    char trans = 'N';

    // solve a system
    cgetrs_( &trans, &nr, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<float> >::schurOfHessenberg(
    FullMatrix<int, std::complex<float> > &Z,
    std::complex<float> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    int ilo = 1; //
    int ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( int i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    int info = 0;
    int workSize = 3 * this->getNRows( );
    std::complex<float> *localWork = new std::complex<float>[workSize];

    chseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, eigs, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    delete [] localWork;
  }

  template<>
  void FullMatrix<int, std::complex<float> >::print( std::ostream &stream ) const {
    std::ios::fmtflags f( std::cout.flags( ) );

    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( int i = 0; i < this->nRows; i++ ) {
      for ( int j = 0; j < this->nCols; j++ ) {
        stream << std::real( get( i, j ) ) << std::showpos <<
          std::imag( get( i, j ) ) << "i ";
        std::cout.flags( f );
      }
      stream << "\n";
    }
    stream << "\n";
  }

  // <INT, COMPLEX DOUBLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< int, std::complex< double > >::applyMIC(
    Vector< int, std::complex< double > > const & x,
    Vector< int, std::complex< double > > & y,
    bool transA,
    std::complex< double > alpha,
    std::complex< double > beta,
    int device
    ) {

    typedef int LO;
    typedef std::complex< double > SC;
    typedef double SCVT;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SCVT * alpha_c = reinterpret_cast < SCVT * > ( &alpha );
    SCVT * beta_c = reinterpret_cast < SCVT * > ( &beta );
    SCVT * data_c = reinterpret_cast < SCVT * > ( this->data );
    SCVT * ydata_c = reinterpret_cast < SCVT * > ( y.data );
    SCVT * xdata_c = reinterpret_cast < SCVT * > ( x.data );

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : alpha_c[ 0 : 2 ], beta_c[ 0 : 2 ] ) \
map( to : transAChar, nRows, nCols, one )
#endif
    {
      zgemv_( &transAChar, &nRows, &nCols,
        reinterpret_cast < SC * > ( alpha_c ),
        reinterpret_cast < SC * > ( data_c ),
        &nRows,
        reinterpret_cast < SC * > ( xdata_c ),
        &one,
        reinterpret_cast < SC * > ( beta_c ),
        reinterpret_cast < SC * > ( ydata_c ),
        &one );
    }
  }

  template<>
  void FullMatrix<int, std::complex<double> >::getCol( int idx, std::complex<double> *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (std::complex<double> ) );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::getRow( int idx, std::complex<double> *outCol ) const {
    for ( int i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<int, std::complex<double> >::scale( std::complex<double> alpha ) {
    int dataSize = this->nRows * this->nCols;
    int incx = 1;
    zscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  std::complex<double> FullMatrix<int, std::complex<double> >::norm1( ) const {
    char norm = '1';
    double *work2 = new double[this->nRows];
    double norm1 = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( norm1, 0.0 );

  }

  template<>
  std::complex<double> FullMatrix<int, std::complex<double> >::normFro( ) {
    char norm = 'F';
    double *work2 = new double[this->nRows];
    double normF = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( normF, 0.0 );
  }

  template<>
  std::complex<double> FullMatrix<int, std::complex<double> >::normI( ) {
    char norm = 'I';
    double *work2 = new double[this->nRows];
    double normI = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( normI, 0.0 );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::add( FullMatrix<int, std::complex<double> > &A, std::complex<double> alpha ) {
    int size = this->nRows * this->nCols;
    int inc = 1;
    zaxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::multiply(
    FullMatrix<int, std::complex<double> > &A,
    FullMatrix<int, std::complex<double> > &B,
    bool transA,
    bool transB,
    std::complex<double> alpha,
    std::complex<double> beta
    ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    zgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::multiply(
    FullMatrix<int, std::complex<double> > &A,
    FullMatrix<int, std::complex<double> > &B,
    int ARows,
    int ACols,
    int BRows,
    int BCols,
    bool transA,
    bool transB,
    std::complex<double> alpha,
    std::complex<double> beta ) {

    char transAChar, transBChar;
    int lda, ldb;
    int nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    zgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::apply(
    Vector<int, std::complex<double> > const & x,
    Vector<int, std::complex<double> > & y,
    bool transA,
    std::complex<double> alpha,
    std::complex<double> beta

    ) {
    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    zgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::applySubmatrix(
    Vector<int, std::complex<double> > const &x,
    Vector<int, std::complex<double> > &y,
    int ARows,
    int ACols,
    bool transA,
    std::complex<double> alpha,
    std::complex<double> beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    int one = 1;

    zgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::LUSolve( Vector<int, std::complex<double> > &x, int nRhs ) {

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    int info;
    char trans = 'N';

    // factorize a matrix (LU)
    zgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    zgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<double> >::LUSolve(
    FullMatrix<int, std::complex<double> > & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int info = 0;
    char trans = 'N';

    // factorize a matrix
    int * ipiv = new int[ipivLength];
    zgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    zgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<double> >::CholeskiSolve(
    Vector<int, std::complex<double> > & x,
    int nRhs
    ) {

    int info = 0;
    char uplo = 'L';

    zpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    zpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::CholeskiSolve(
    FullMatrix<int, std::complex<double> > & x,
    int nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    int info = 0;
    char uplo = 'L';

    zpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    zpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<int, std::complex<double> >::backward( Vector<int, std::complex<double> > &x, int nRhs, int n ) {

    int ipivLength;
    int nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    int *ipiv = new int[ipivLength];
    for ( int i = 0; i < nr; i++ ) {
      ipiv[i] = i + 1;
    }
    int info;
    char trans = 'N';

    // solve a system
    zgetrs_( &trans, &nr, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<int, std::complex<double> >::schurOfHessenberg(
    FullMatrix<int, std::complex<double> > &Z,
    std::complex<double> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    int ilo = 1; //
    int ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( long i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    int info = 0;
    int workSize = 3 * this->getNRows( );
    std::complex<double> *localWork = new std::complex<double>[workSize];

    zhseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, eigs, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    delete [] localWork;
  }

  template<>
  void FullMatrix<int, std::complex<double> >::eigs( std::complex<double> *eigenvectors, std::complex<double> *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    std::complex<double> *up_package = new std::complex<double>[( int ) this->nRows * ( this->nRows + 1 ) / 2];
    std::complex<double> *workspace = new std::complex<double>[8 * this->nRows];
    int lwork = 8 * this->nRows;
    int *iwork = new int[5 * this->nRows];
    double *rwork = new double[7 * this->nRows];
    int *ifail = new int[this->nRows];
    int info;

    int idx = 0;
    int m;
    double *w = new double[this->nRows];


    // pack matrix to upper form
    for ( int j = 0; j < this->nRows; j++ ) {
      for ( int i = 0; i <= j; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }
    double eps = ( double ) EPS;
    zheevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, w, eigenvectors, &this->nRows, workspace, &lwork, rwork, iwork, ifail, &info );
    for ( int i = 0; i < this->nRows; i++ ) {
      eigenvalues[i] = w[i];
    }
    delete [] up_package;
    delete [] workspace;
    delete [] w;
    delete [] iwork;
    delete [] ifail;
    delete [] rwork;
  }

  template<>
  void FullMatrix<int, std::complex<double> >::print( std::ostream &stream ) const {
    std::ios::fmtflags f( std::cout.flags( ) );

    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( int i = 0; i < this->nRows; i++ ) {
      for ( int j = 0; j < this->nCols; j++ ) {
        stream << std::real( get( i, j ) ) << std::showpos <<
          std::imag( get( i, j ) ) << "i ";
        std::cout.flags( f );
      }
      stream << "\n";
    }
    stream << "\n";
  }

#elif defined(BLAS_LONG)

  // <LONG, DOUBLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< long, double >::applyMIC(
    Vector< long, double > const & x,
    Vector< long, double > & y,
    bool transA,
    double alpha,
    double beta,
    int device
    ) {

    typedef long LO;
    typedef double SC;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SC * data = this->data;
    SC * xdata = x.data;
    SC * ydata = y.data;

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : transAChar, nRows, nCols, alpha, one, beta )
#endif
    {
      dgemv_( &transAChar, &nRows, &nCols, &alpha, data, &nRows, xdata, &one,
        &beta, ydata, &one );
    }
  }

  template<>
  void FullMatrix<long, double>::getCol( long idx, double *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (double ) );
  }

  template<>
  void FullMatrix<long, double>::getRow( long idx, double *outCol ) const {
    for ( long i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<long, double>::scale( double alpha ) {
    long dataSize = this->nRows * this->nCols;
    long incx = 1;
    dscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  double FullMatrix<long, double>::norm1( ) const {
    char norm = '1';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  double FullMatrix<long, double>::normFro( ) {
    char norm = 'F';
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  double FullMatrix<long, double>::normI( ) {
    char norm = 'I';
    if ( work == nullptr ) {
      work = new double[this->nRows];
    }
    return dlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  void FullMatrix<long, double>::add( FullMatrix<long, double> &A, double alpha ) {
    long size = this->nRows * this->nCols;
    long inc = 1;
    daxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<long, double>::multiply( FullMatrix<long, double> &A, FullMatrix<long, double> &B, bool transA,
    bool transB, double alpha, double beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    dgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, double>::multiply(
    FullMatrix<long, double> &A,
    FullMatrix<long, double> &B,
    long ARows,
    long ACols,
    long BRows,
    long BCols,
    bool transA,
    bool transB,
    double alpha,
    double beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    dgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, double>::apply(
    Vector<long, double> const & x,
    Vector<long, double> & y,
    bool transA,
    double alpha,
    double beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    dgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows, x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, double>::applySubmatrix(
    Vector<long, double> const &x,
    Vector<long, double> &y,
    long ARows,
    long ACols,
    bool transA,
    double alpha,
    double beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    dgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, double>::LUSolve( Vector<long, double> &x, long nRhs ) {

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    long info;
    char trans = 'N';

    // factorize a matrix
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, double >::LUSolve(
    FullMatrix<long, double > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long info = 0;
    char trans = 'N';

    // factorize a matrix
    long * ipiv = new long[ipivLength];
    dgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    dgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix< long, double >::CholeskiSolve(
    Vector< long, double > & x,
    long nRhs
    ) {

    long info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix< long, double >::CholeskiSolve(
    FullMatrix< long, double > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long info = 0;
    char uplo = 'L';

    dpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    dpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<long, double>::backward(
    Vector<long, double> &x,
    long nRhs,
    long n
    ) {
    long ipivLength;
    long nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    for ( long i = 0; i < nr; i++ ) {
      ipiv[ i ] = i + 1;
    }
    long info;
    char trans = 'N';

    // solve a system
    dgetrs_( &trans, &nr, &nRhs, this->data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;

  }

  template<>
  void FullMatrix<long, double>::schurOfHessenberg(
    FullMatrix<long, double> &Z,
    std::complex<double> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    long ilo = 1; //
    long ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( long i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    long info = 0;
    long workSize = 3 * this->getNRows( );
    double *localWork = new double[workSize];
    double *wr = new double[this->getNRows( )];
    double *wi = new double[this->getNRows( )];

    dhseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, wr, wi, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    for ( long i = 0; i < this->getNRows( ); i++ ) {
      eigs[i] = std::complex<double>( wr[i], wi[i] );
    }

    delete [] localWork;
    delete [] wr;
    delete [] wi;
  }

  template<>
  void FullMatrix<long, double>::eigs( double *eigenvectors, double *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    double *up_package = new double[this->nRows * this->nRows]; //(long) this->nRows*(this->nRows+1)/2
    double *workspace = new double[8 * this->nRows];
    long lwork = 8 * this->nRows;
    long *iwork = new long[5 * this->nRows];
    long *ifail = new long[this->nRows];
    long info;

    long idx = 0;
    long m;

    // pack matrix to upper form
    for ( long j = 0; j < this->nRows; j++ ) {
      for ( long i = 0; i < this->nRows; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    double eps = ( double ) EPS;
    dsyevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, eigenvalues, eigenvectors, &this->nRows, workspace, &lwork, iwork, ifail, &info );
    delete [] up_package;
    delete [] workspace;
    delete [] iwork;
    delete [] ifail;
  }

  template<>
  void FullMatrix<long, double>::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( long i = 0; i < this->nRows; i++ ) {
      for ( long j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }


  // <LONG, SINGLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< long, float >::applyMIC(
    Vector< long, float > const & x,
    Vector< long, float > & y,
    bool transA,
    float alpha,
    float beta,
    int device
    ) {

    typedef long LO;
    typedef float SC;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SC * data = this->data;
    SC * xdata = x.data;
    SC * ydata = y.data;

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : transAChar, nRows, nCols, alpha, one, beta )
#endif
    {
      sgemv_( &transAChar, &nRows, &nCols, &alpha, data, &nRows, xdata, &one,
        &beta, ydata, &one );
    }
  }

  template<>
  void FullMatrix<long, float>::getCol( long idx, float *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (float ) );
  }

  template<>
  void FullMatrix<long, float>::getRow( long idx, float *outCol ) const {
    for ( long i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<long, float>::scale( float alpha ) {
    long dataSize = this->nRows * this->nCols;
    long incx = 1;
    sscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  float FullMatrix<long, float>::norm1( ) const {
    char norm = '1';
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  float FullMatrix<long, float>::normFro( ) {
    char norm = 'F';
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  float FullMatrix<long, float>::normI( ) {
    char norm = 'I';
    if ( work == nullptr ) {
      work = new float[this->nRows];
    }
    return slange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work );
  }

  template<>
  void FullMatrix<long, float>::add( FullMatrix<long, float> &A, float alpha ) {
    long size = this->nRows * this->nCols;
    long inc = 1;
    saxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<long, float>::multiply( FullMatrix<long, float> &A, FullMatrix<long, float> &B, bool transA,
    bool transB, float alpha, float beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    sgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, float>::multiply(
    FullMatrix<long, float> &A,
    FullMatrix<long, float> &B,
    long ARows,
    long ACols,
    long BRows,
    long BCols,
    bool transA,
    bool transB,
    float alpha,
    float beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    sgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, float>::apply(
    Vector<long, float> const & x,
    Vector<long, float> & y,
    bool transA,
    float alpha,
    float beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    sgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows, x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, float>::applySubmatrix(
    Vector<long, float> const &x,
    Vector<long, float> &y,
    long ARows,
    long ACols,
    bool transA,
    float alpha,
    float beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    sgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, float>::LUSolve( Vector<long, float> &x, long nRhs ) {

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    long info;
    char trans = 'N';

    // factorize a matrix (LU)
    sgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    sgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, float >::LUSolve(
    FullMatrix< long, float > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long info = 0;
    char trans = 'N';

    // factorize a matrix
    long * ipiv = new long[ipivLength];
    sgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    sgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix< long, float >::CholeskiSolve(
    Vector< long, float > & x,
    long nRhs
    ) {

    long info = 0;
    char uplo = 'L';

    spotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    spotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix< long, float >::CholeskiSolve(
    FullMatrix< long, float > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long info = 0;
    char uplo = 'L';

    spotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    spotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<long, float>::backward(
    Vector<long, float> &x,
    long nRhs,
    long n
    ) {

    long ipivLength;
    long nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    for ( long i = 0; i < nr; i++ ) {
      ipiv[ i ] = i + 1;
    }
    long info;
    char trans = 'N';

    // solve a system
    sgetrs_( &trans, &nr, &nRhs, this->data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, float>::schurOfHessenberg(
    FullMatrix<long, float> &Z,
    std::complex<float> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    long ilo = 1; //
    long ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( long i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    long info = 0;
    long workSize = 3 * this->getNRows( );
    float *localWork = new float[workSize];
    float *wr = new float[this->getNRows( )];
    float *wi = new float[this->getNRows( )];

    shseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, wr, wi, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    for ( long i = 0; i < this->getNRows( ); i++ ) {
      eigs[i] = std::complex<double>( wr[i], wi[i] );
    }

    delete [] localWork;
    delete [] wr;
    delete [] wi;
  }

  template<>
  void FullMatrix<long, float>::eigs( float *eigenvectors, float *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    float *up_package = new float[this->nRows * this->nRows];
    float *workspace = new float[8 * this->nRows];
    long lwork = 8 * this->nRows;
    long *iwork = new long[5 * this->nRows];
    long *ifail = new long[this->nRows];
    long info;

    long idx = 0;
    long m;

    // pack matrix to upper form
    for ( long j = 0; j < this->nRows; j++ ) {
      for ( long i = 0; i < this->nRows; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    float eps = ( float ) EPS;
    ssyevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, eigenvalues, eigenvectors, &this->nRows, workspace, &lwork, iwork, ifail, &info );
    delete [] up_package;
    delete [] workspace;
    delete [] iwork;
    delete [] ifail;
  }

  template<>
  void FullMatrix<long, float>::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( long i = 0; i < this->nRows; i++ ) {
      for ( long j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }


  // <LONG, COMPLEX SINGLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< long, std::complex< float > >::applyMIC(
    Vector< long, std::complex< float > > const & x,
    Vector< long, std::complex< float > > & y,
    bool transA,
    std::complex< float > alpha,
    std::complex< float > beta,
    int device
    ) {

    typedef long LO;
    typedef std::complex< float > SC;
    typedef float SCVT;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SCVT * alpha_c = reinterpret_cast < SCVT * > ( &alpha );
    SCVT * beta_c = reinterpret_cast < SCVT * > ( &beta );
    SCVT * data_c = reinterpret_cast < SCVT * > ( this->data );
    SCVT * ydata_c = reinterpret_cast < SCVT * > ( y.data );
    SCVT * xdata_c = reinterpret_cast < SCVT * > ( x.data );

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : alpha_c[ 0 : 2 ], beta_c[ 0 : 2 ] ) \
map( to : transAChar, nRows, nCols, one )
#endif
    {
      cgemv_( &transAChar, &nRows, &nCols,
        reinterpret_cast < SC * > ( alpha_c ),
        reinterpret_cast < SC * > ( data_c ),
        &nRows,
        reinterpret_cast < SC * > ( xdata_c ),
        &one,
        reinterpret_cast < SC * > ( beta_c ),
        reinterpret_cast < SC * > ( ydata_c ),
        &one );
    }
  }

  template<>
  void FullMatrix<long, std::complex<float> >::getCol( long idx, std::complex<float> *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (std::complex<float> ) );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::getRow( long idx, std::complex<float> *outCol ) const {
    for ( long i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<long, std::complex<float> >::scale( std::complex<float> alpha ) {
    long dataSize = this->nRows * this->nCols;
    long incx = 1;
    cscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  std::complex<float> FullMatrix<long, std::complex<float> >::norm1( ) const {
    char norm = '1';
    float *work2 = new float[this->nRows];
    float norm1 = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( norm1, 0.0 );

  }

  template<>
  std::complex<float> FullMatrix<long, std::complex<float> >::normFro( ) {
    char norm = 'F';
    float *work2 = new float[this->nRows];
    float normF = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( normF, 0.0 );
  }

  template<>
  std::complex<float> FullMatrix<long, std::complex<float> >::normI( ) {
    char norm = 'I';
    float *work2 = new float[this->nRows];
    float normI = clange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<float>( normI, 0.0 );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::add( FullMatrix<long, std::complex<float> > &A, std::complex<float> alpha ) {
    long size = this->nRows * this->nCols;
    long inc = 1;
    caxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::multiply( FullMatrix<long, std::complex<float> > &A, FullMatrix<long, std::complex<float> > &B, bool transA,
    bool transB, std::complex<float> alpha, std::complex<float> beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    cgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::multiply(
    FullMatrix<long, std::complex<float> > &A,
    FullMatrix<long, std::complex<float> > &B,
    long ARows,
    long ACols,
    long BRows,
    long BCols,
    bool transA,
    bool transB,
    std::complex<float> alpha,
    std::complex<float> beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    cgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::apply(
    Vector<long, std::complex<float> > const & x,
    Vector<long, std::complex<float> > & y,
    bool transA,
    std::complex<float> alpha,
    std::complex<float> beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    cgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::applySubmatrix(
    Vector<long, std::complex<float> > const &x,
    Vector<long, std::complex<float> > &y,
    long ARows,
    long ACols,
    bool transA,
    std::complex<float> alpha,
    std::complex<float> beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    cgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::eigs( std::complex<float> *eigenvectors, std::complex<float> *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    std::complex<float> *up_package = new std::complex<float>[this->nRows * this->nRows];
    std::complex<float> *workspace = new std::complex<float>[8 * this->nRows];
    long lwork = 8 * this->nRows;
    long *iwork = new long[5 * this->nRows];
    float *rwork = new float[7 * this->nRows];
    long *ifail = new long[this->nRows];
    long info;
    float *w = new float[this->nRows];

    long idx = 0;
    long m;

    // pack matrix to upper form
    // pack matrix to upper form
    for ( long j = 0; j < this->nRows; j++ ) {
      for ( long i = 0; i <= j; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }

    float eps = ( float ) EPS;
    cheevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, w, eigenvectors, &this->nRows, workspace, &lwork, rwork, iwork, ifail, &info );

    for ( long i = 0; i < this->nRows; i++ ) {
      eigenvalues[i] = w[i];
    }
    delete [] up_package;
    delete [] workspace;
    delete [] w;
    delete [] iwork;
    delete [] ifail;
    delete [] rwork;
  }

  template<>
  void FullMatrix<long, std::complex<float> >::LUSolve( Vector<long, std::complex<float> > &x, long nRhs ) {

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    long info;
    char trans = 'N';

    // factorize a matrix (LU)
    cgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    cgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, std::complex<float> >::LUSolve(
    FullMatrix<long, std::complex<float> > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long info = 0;
    char trans = 'N';

    // factorize a matrix
    long * ipiv = new long[ipivLength];
    cgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    cgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix< long, std::complex< float > >::CholeskiSolve(
    Vector< long, std::complex< float > > & x,
    long nRhs
    ) {

    long info = 0;
    char uplo = 'L';

    cpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    cpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix< long, std::complex< float > >::CholeskiSolve(
    FullMatrix< long, std::complex< float > > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long info = 0;
    char uplo = 'L';

    cpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    cpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<long, std::complex<float> >::backward(
    Vector<long, std::complex<float> > &x,
    long nRhs,
    long n
    ) {

    long ipivLength;
    long nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    for ( long i = 0; i < nr; i++ ) {
      ipiv[ i ] = i + 1;
    }
    long info;
    char trans = 'N';

    // solve a system
    cgetrs_( &trans, &nr, &nRhs, this->data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, std::complex<float> >::schurOfHessenberg(
    FullMatrix<long, std::complex<float> > &Z,
    std::complex<float> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    long ilo = 1; //
    long ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( long i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    long info = 0;
    long workSize = 3 * this->getNRows( );
    std::complex<float> *localWork = new std::complex<float>[workSize];

    chseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, eigs, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    delete [] localWork;
  }

  template<>
  void FullMatrix<long, std::complex<float> >::print( std::ostream &stream ) const {
    stream << this->nRows << "\n";
    stream << this->nCols << "\n\n";
    for ( long i = 0; i < this->nRows; i++ ) {
      for ( long j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }

  // <LONG, COMPLEX DOUBLE> PRECISION SPECIALIZATION //

  template<>
  void FullMatrix< long, std::complex< double > >::applyMIC(
    Vector< long, std::complex< double > > const & x,
    Vector< long, std::complex< double > > & y,
    bool transA,
    std::complex< double > alpha,
    std::complex< double > beta,
    int device
    ) {

    typedef long LO;
    typedef std::complex< double > SC;
    typedef double SCVT;

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    LO one = 1;
    LO nRows = this->nRows;
    LO nCols = this->nCols;

    SCVT * alpha_c = reinterpret_cast < SCVT * > ( &alpha );
    SCVT * beta_c = reinterpret_cast < SCVT * > ( &beta );
    SCVT * data_c = reinterpret_cast < SCVT * > ( this->data );
    SCVT * ydata_c = reinterpret_cast < SCVT * > ( y.data );
    SCVT * xdata_c = reinterpret_cast < SCVT * > ( x.data );

#if N_MIC > 0
#pragma omp target device( device ) \
map( to : alpha_c[ 0 : 2 ], beta_c[ 0 : 2 ] ) \
map( to : transAChar, nRows, nCols, one )
#endif
    {
      zgemv_( &transAChar, &nRows, &nCols,
        reinterpret_cast < SC * > ( alpha_c ),
        reinterpret_cast < SC * > ( data_c ),
        &nRows,
        reinterpret_cast < SC * > ( xdata_c ),
        &one,
        reinterpret_cast < SC * > ( beta_c ),
        reinterpret_cast < SC * > ( ydata_c ),
        &one );
    }
  }

  template<>
  void FullMatrix<long, std::complex<double> >::getCol( long idx, std::complex<double> *outCol ) const {
    memcpy( outCol, this->data + idx * this->nRows, this->nRows * sizeof (std::complex<double> ) );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::getRow( long idx, std::complex<double> *outCol ) const {
    for ( long i = 0; i < this->nCols; i++ ) {
      outCol[i] = data[i * this->nRows + idx];
    }
  }

  template<>
  void FullMatrix<long, std::complex<double> >::scale( std::complex<double> alpha ) {
    long dataSize = this->nRows * this->nCols;
    long incx = 1;
    zscal_( &dataSize, &alpha, data, &incx );
  }

  template<>
  std::complex<double> FullMatrix<long, std::complex<double> >::norm1( ) const {
    char norm = '1';
    double *work2 = new double[this->nRows];
    double norm1 = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( norm1, 0.0 );

  }

  template<>
  std::complex<double> FullMatrix<long, std::complex<double> >::normFro( ) {
    char norm = 'F';
    double *work2 = new double[this->nRows];
    double normF = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( normF, 0.0 );
  }

  template<>
  std::complex<double> FullMatrix<long, std::complex<double> >::normI( ) {
    char norm = 'I';
    double *work2 = new double[this->nRows];
    double normI = zlange_( &norm, &this->nRows, &this->nCols, data, &this->nRows, work2 );
    delete [] work2;
    return std::complex<double>( normI, 0.0 );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::add( FullMatrix<long, std::complex<double> > &A, std::complex<double> alpha ) {
    long size = this->nRows * this->nCols;
    long inc = 1;
    zaxpy_( &size, &alpha, A.data, &inc, data, &inc );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::multiply( FullMatrix<long, std::complex<double> > &A, FullMatrix<long, std::complex<double> > &B, bool transA,
    bool transB, std::complex<double> alpha, std::complex<double> beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nCols;

    lda = A.nRows;
    ldb = B.nRows;

    if ( transA ) {
      transAChar = 'T';
      nCols = A.nRows;
    } else {
      transAChar = 'N';
      nCols = A.nCols;
    }

    if ( transB ) {
      transBChar = 'T';
    } else {
      transBChar = 'N';
    }

    zgemm_( &transAChar, &transBChar, &this->nRows, &this->nCols, &nCols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::multiply(
    FullMatrix<long, std::complex<double> > &A,
    FullMatrix<long, std::complex<double> > &B,
    long ARows,
    long ACols,
    long BRows,
    long BCols,
    bool transA,
    bool transB,
    std::complex<double> alpha,
    std::complex<double> beta ) {

    char transAChar, transBChar;
    long lda, ldb;
    long nACols, nARows, nBCols;

    if ( transA ) {
      transAChar = 'T';
      lda = A.nCols;
      nACols = ARows;
      nARows = ACols;
    } else {
      transAChar = 'N';
      lda = A.nRows;
      nARows = ARows;
      nACols = ACols;
    }

    if ( transB ) {
      transBChar = 'T';
      ldb = B.nCols;
      nBCols = BRows;
    } else {
      transBChar = 'N';
      ldb = B.nRows;
      nBCols = BCols;
    }

    zgemm_( &transAChar, &transBChar, &nARows, &nBCols, &nACols, &alpha,
      A.data, &lda, B.data, &ldb, &beta, data, &this->nRows );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::apply(
    Vector<long, std::complex<double> > const & x,
    Vector<long, std::complex<double> > & y,
    bool transA,
    std::complex<double> alpha,
    std::complex<double> beta
    ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    zgemv_( &transAChar, &this->nRows, &this->nCols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::applySubmatrix(
    Vector<long, std::complex<double> > const &x,
    Vector<long, std::complex<double> > &y,
    long ARows,
    long ACols,
    bool transA,
    std::complex<double> alpha,
    std::complex<double> beta ) {

    char transAChar;

    if ( transA ) {
      transAChar = 'T';
    } else {
      transAChar = 'N';
    }

    long one = 1;

    zgemv_( &transAChar, &ARows, &ACols, &alpha, data, &this->nRows,
      x.data, &one, &beta, y.data, &one );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::LUSolve( Vector<long, std::complex<double> > &x, long nRhs ) {

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    long info;
    char trans = 'N';

    // factorize a matrix (LU)
    zgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    zgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, std::complex<double> >::LUSolve(
    FullMatrix<long, std::complex<double> > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long ipivLength;
    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long info = 0;
    char trans = 'N';

    // factorize a matrix
    long * ipiv = new long[ipivLength];
    zgetrf_( &this->nRows, &this->nCols, data, &this->nRows, ipiv, &info );

    // solve a system
    zgetrs_( &trans, &this->nRows, &nRhs, data, &this->nRows, ipiv, x.data,
      &this->nRows, &info );
    delete [] ipiv;
  }

  template<>
  void FullMatrix< long, std::complex< double > >::CholeskiSolve(
    Vector< long, std::complex< double > > & x,
    long nRhs
    ) {

    long info = 0;
    char uplo = 'L';

    zpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    zpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix< long, std::complex< double > >::CholeskiSolve(
    FullMatrix< long, std::complex< double > > & x,
    long nRhs
    ) {

    if ( nRhs <= 0 ) nRhs = x.getNCols( );

    long info = 0;
    char uplo = 'L';

    zpotrf_( &uplo, &this->nRows, this->data, &this->nRows, &info );
    zpotrs_( &uplo, &this->nRows, &nRhs, this->data, &this->nRows, x.data,
      &this->nRows, &info );
  }

  template<>
  void FullMatrix<long, std::complex<double> >::backward(
    Vector<long, std::complex<double> > &x,
    long nRhs,
    long n
    ) {

    long ipivLength;
    long nr = n;

    if ( n == 0 ) {
      nr = this->nRows;
    }

    if ( this->nRows < this->nCols ) {
      ipivLength = this->nRows;
    } else {
      ipivLength = this->nCols;
    }

    long *ipiv = new long[ipivLength];
    for ( long i = 0; i < nr; i++ ) {
      ipiv[ i ] = i + 1;
    }
    long info;
    char trans = 'N';

    // solve a system
    zgetrs_( &trans, &nr, &nRhs, this->data, &this->nRows, ipiv, x.data, &this->nRows, &info );

    delete [] ipiv;
  }

  template<>
  void FullMatrix<long, std::complex<double> >::schurOfHessenberg(
    FullMatrix<long, std::complex<double> > &Z,
    std::complex<double> * eigs
    ) {
    char job = 'E';
    char compz = 'I';

    long ilo = 1; //
    long ihi = this->getNRows( );
    ; //1;

    // if Z has not a good size, resize it
    if ( ( Z.getNRows( ) != this->getNRows( ) )
      && ( Z.getNCols( ) != this->getNCols( ) ) ) {
      Z.resize( this->getNRows( ), this->getNCols( ) );
    }
    Z.setAll( 0.0 );
    for ( long i = 0; i < this->getNRows( ); i++ ) {
      Z.set( i, i, 1.0 );
    }

    long info = 0;
    long workSize = 3 * this->getNRows( );
    std::complex<double> *localWork = new std::complex<double>[workSize];

    zhseqr_( &job, &compz, &this->nRows, &ilo, &ihi, this->data,
      &this->nRows, eigs, Z.data, &Z.nRows, localWork,
      &workSize, &info );

    delete [] localWork;
  }

  template<>
  void FullMatrix<long, std::complex<double> >::eigs( std::complex<double> *eigenvectors, std::complex<double> *eigenvalues ) {
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    std::complex<double> *up_package = new std::complex<double>[( long ) this->nRows * ( this->nRows + 1 ) / 2];
    std::complex<double> *workspace = new std::complex<double>[8 * this->nRows];
    long lwork = 8 * this->nRows;
    long *iwork = new long[5 * this->nRows];
    double *rwork = new double[7 * this->nRows];
    long *ifail = new long[this->nRows];
    long info;

    long idx = 0;
    long m;
    double *w = new double[this->nRows];


    // pack matrix to upper form
    for ( long j = 0; j < this->nRows; j++ ) {
      for ( long i = 0; i <= j; i++ ) {
        up_package[idx] = this->data[idx];
        idx++;
      }
    }
    double eps = ( double ) EPS;
    zheevx_( &jobz, &range, &uplo, &this->nRows, up_package, &this->nRows,
      nullptr, nullptr, nullptr, nullptr, &eps, &m, w, eigenvectors, &this->nRows, workspace, &lwork, rwork, iwork, ifail, &info );
    for ( long i = 0; i < this->nRows; i++ ) {
      eigenvalues[i] = w[i];
    }
    delete [] up_package;
    delete [] workspace;
    delete [] w;
    delete [] iwork;
    delete [] ifail;
    delete [] rwork;
  }

  template<>
  void FullMatrix<long, std::complex<double> >::print( std::ostream &stream ) const {
    stream << "\nFullMatrix<long, std::complex<double> >\n";
    stream << "Number of rows: " << this->nRows << "\n";
    stream << "Number of cols: " << this->nCols << "\n";
    for ( long i = 0; i < this->nRows; i++ ) {
      for ( long j = 0; j < this->nCols; j++ ) {
        stream << get( i, j ) << " ";
      }
      stream << "\n";
    }
    stream << "\n";
  }

#endif

}

#endif
