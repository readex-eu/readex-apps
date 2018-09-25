/*!
 * @file    IdentityOperator.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 18, 2013
 *
 */

#ifdef IDENTITYOPERATOR_H

namespace bem4i {

template<class LO, class SC>
IdentityOperator<LO, SC>::~IdentityOperator( ) {
}

template<class LO, class SC>
IdentityOperator<LO, SC>::IdentityOperator(
    BESpace<LO, SC>* space
    ) {

  this->space = space;

  basisType ansatz = space->getAnsatzFunctionType( );
  switch ( ansatz ) {
    case p0:
      this->dimDomain = space->getMesh( )->getNElements( );
      break;
    case p1:
      this->dimDomain = space->getMesh( )->getNNodes( );
      break;
    case p1dis:
      this->dimDomain = 3 * space->getMesh( )->getNElements( );
      break;
  }

  basisType test = space->getTestFunctionType( );
  switch ( test ) {
    case p0:
      this->dimRange = space->getMesh( )->getNElements( );
      break;
    case p1:
      this->dimRange = space->getMesh( )->getNNodes( );
      break;
    case p1dis:
      this->dimRange = 3 * space->getMesh( )->getNElements( );
      break;
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::apply(
    Vector<LO, SC> const & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  switch ( space->getTestFunctionType( ) ) {
    case p0:
      switch ( space->getAnsatzFunctionType( ) ) {
        case p0:
          this->applyP0P0( x, y, transA, alpha, beta );
          break;
        case p1:
          this->applyP0P1( x, y, transA, alpha, beta );
          break;
        default:
          std::cout << "Not yet implemented!" << std::endl;
          exit( 1 );
          break;
      }
      break;
    case p1:
      switch ( space->getAnsatzFunctionType( ) ) {
        case p0:
          this->applyP1P0( x, y, transA, alpha, beta );
          break;
        case p1:
          this->applyP1P1( x, y, transA, alpha, beta );
          break;
        default:
          std::cout << "Not yet implemented!" << std::endl;
          exit( 1 );
          break;
      }
      break;
    default:
      std::cout << "Not yet implemented!" << std::endl;
      exit( 1 );
      break;
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::assemble(
    SparseMatrix< LO, SC > & M
    ) const {

  switch ( space->getTestFunctionType( ) ) {
    case p0:
      switch ( space->getAnsatzFunctionType( ) ) {
        case p0:
          this->assembleP0P0( M );
          break;
        case p1:
          this->assembleP0P1( M );
          break;
        default:
          std::cout << "Not yet implemented!" << std::endl;
          exit( 1 );
          break;
      }
      break;
    case p1:
      switch ( space->getAnsatzFunctionType( ) ) {
        case p0:
          this->assembleP1P0( M );
          break;
        case p1:
          this->assembleP1P1( M );
          break;
        default:
          std::cout << "Not yet implemented!" << std::endl;
          exit( 1 );
          break;
      }
      break;
    default:
      std::cout << "Not yet implemented!" << std::endl;
      exit( 1 );
      break;
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::applyP0P0(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  for ( LO i = 0; i < y.getLength( ); i++ ) {
    y.data[i] = y.data[i] * beta +
        alpha * x.data[i] * space->getMesh( )->getElemArea( i );
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::applyP0P1(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  LO idx[3];
  SCVT zero = 0.0;
  if ( !transA ) {
    for ( LO i = 0; i < y.getLength( ); i++ ) {
      space->getMesh( )->getElement( i, idx );
      y.data[i] = y.data[i] * beta +
          alpha * space->getMesh( )->getElemArea( i ) / ( (SCVT) 3.0 ) *
          ( x.data[idx[0]] + x.data[idx[1]] + x.data[idx[2]] );
    }
  } else {
    SC scaledArea;
    if ( beta == zero ) {
      y.setAll( zero );
    } else {
      y.scale( beta );
    }
    for ( LO i = 0; i < x.getLength( ); i++ ) {
      scaledArea = alpha * space->getMesh( )->getElemArea( i ) / ( (SCVT) 3.0 );
      space->getMesh( )->getElement( i, idx );
      y.data[idx[0]] += scaledArea * x.data[i];
      y.data[idx[1]] += scaledArea * x.data[i];
      y.data[idx[2]] += scaledArea * x.data[i];
    }
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::applyP1P0(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta ) {

  this->applyP0P1( x, y, !transA, alpha, beta );
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::applyP1P1(
    Vector<LO, SC> const &x,
    Vector<LO, SC> &y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  LO nElems = this->space->getMesh( )->getNElements( );

  LO ind[ 3 ];
  SCVT area;
  SCVT aux[ 9 ];

  Vector< LO, SC > yCopy;

  y.scale( beta );

  for ( int j = 0; j < 9; ++j ) {
    aux[ j ] = 0.0;
  }

  for ( int q = 0; q < 3; ++q ) {
    aux[ 0 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
        quadPoints2[ 2 * q + 1 ] ) * ( 1.0 - quadPoints2[ 2 * q ] -
        quadPoints2[ 2 * q + 1 ] );
    aux[ 1 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
        quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q ];
    aux[ 2 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
        quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q + 1 ];
    aux[ 3 ] += quadWeights2[ q ] * quadPoints2[ 2 * q ] *
        quadPoints2[ 2 * q ];
    aux[ 4 ] += quadWeights2[ q ] * quadPoints2[ 2 * q ] *
        quadPoints2[ 2 * q + 1 ];
    aux[ 5 ] += quadWeights2[ q ] * quadPoints2[ 2 * q + 1 ] *
        quadPoints2[ 2 * q + 1 ];
  }

  for ( LO i = 0; i < nElems; ++i ) {
    this->space->getMesh( )->getElement( i, ind );
    area = this->space->getMesh( )->getElemArea( i );

    y.add( ind[ 0 ], x.get( ind[ 0 ] ) * aux[ 0 ] * area );
    y.add( ind[ 0 ], x.get( ind[ 1 ] ) * aux[ 1 ] * area );
    y.add( ind[ 0 ], x.get( ind[ 2 ] ) * aux[ 2 ] * area );
    y.add( ind[ 1 ], x.get( ind[ 0 ] ) * aux[ 1 ] * area );
    y.add( ind[ 1 ], x.get( ind[ 1 ] ) * aux[ 3 ] * area );
    y.add( ind[ 1 ], x.get( ind[ 2 ] ) * aux[ 4 ] * area );
    y.add( ind[ 2 ], x.get( ind[ 0 ] ) * aux[ 2 ] * area );
    y.add( ind[ 2 ], x.get( ind[ 1 ] ) * aux[ 4 ] * area );
    y.add( ind[ 2 ], x.get( ind[ 2 ] ) * aux[ 5 ] * area );

  }

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::assembleP0P0(
    SparseMatrix<LO, SC> & M
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );

  std::vector< LO > ind;
  ind.reserve( nElems );
  std::vector< SC > val;
  val.reserve( nElems );

  for ( LO i = 0; i < nElems; ++i ) {
    ind.push_back( i );
    val.push_back( this->space->getMesh( )->getElemArea( i ) );
  }

  M.setFromTriplets( nElems, nElems, ind, ind, val );
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::assembleP0P1(
    SparseMatrix<LO, SC> & M
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SC > val;
  val.reserve( 3 * nElems );

  LO elem[ 3 ];
  SCVT integral;

  for ( LO i = 0; i < nElems; ++i ) {
    this->space->getMesh( )->getElement( i, elem );
    integral = this->space->getMesh( )->getElemArea( i ) / 3.0;
    rowInd.push_back( i );
    rowInd.push_back( i );
    rowInd.push_back( i );
    colInd.push_back( elem[ 0 ] );
    colInd.push_back( elem[ 1 ] );
    colInd.push_back( elem[ 2 ] );
    val.push_back( integral );
    val.push_back( integral );
    val.push_back( integral );
  }

  M.setFromTriplets( nElems, nNodes, rowInd, colInd, val );
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::assembleP1P0(
    SparseMatrix<LO, SC> & M
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 3 * nElems );
  std::vector< LO > colInd;
  colInd.reserve( 3 * nElems );
  std::vector< SC > val;
  val.reserve( 3 * nElems );

  LO elem[ 3 ];
  SCVT integral;

  for ( LO i = 0; i < nElems; ++i ) {
    this->space->getMesh( )->getElement( i, elem );
    integral = this->space->getMesh( )->getElemArea( i ) / 3.0;
    colInd.push_back( i );
    colInd.push_back( i );
    colInd.push_back( i );
    rowInd.push_back( elem[ 0 ] );
    rowInd.push_back( elem[ 1 ] );
    rowInd.push_back( elem[ 2 ] );
    val.push_back( integral );
    val.push_back( integral );
    val.push_back( integral );
  }

  M.setFromTriplets( nNodes, nElems, rowInd, colInd, val );
}

template<class LO, class SC>
void IdentityOperator<LO, SC>::assembleP1P1(
    SparseMatrix< LO, SC > & M
    ) const {

  LO nNodes = this->space->getMesh( )->getNNodes( );
  LO nElems = this->space->getMesh( )->getNElements( );
  std::vector< LO > rowInd;
  std::vector< LO > colInd;
  std::vector< SC > values;
  rowInd.reserve( 9 * nElems );
  colInd.reserve( 9 * nElems );
  values.reserve( 9 * nElems );

  LO ind[ 3 ];
  SCVT area;
  SCVT aux[ 9 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  for ( LO i = 0; i < nElems; ++i ) {
    this->space->getMesh( )->getElement( i, ind );
    area = this->space->getMesh( )->getElemArea( i );
    for ( int left = 0; left < 3; ++left ) {
      for ( int right = 0; right < 3; ++right ) {
        rowInd.push_back( ind[ left ] );
        colInd.push_back( ind[ right ] );
      }
    }
    for ( int q = 0; q < 3; ++q ) {
      aux[ 0 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] ) * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] );
      aux[ 1 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q ];
      aux[ 2 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q + 1 ];
      aux[ 3 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q ];
      aux[ 4 ] += quadWeights2[ q ] * quadPoints2[ 2 * q ] *
          quadPoints2[ 2 * q ];
      aux[ 5 ] += quadWeights2[ q ] * quadPoints2[ 2 * q ] *
          quadPoints2[ 2 * q + 1 ];
      aux[ 6 ] += quadWeights2[ q ] * ( 1.0 - quadPoints2[ 2 * q ] -
          quadPoints2[ 2 * q + 1 ] ) * quadPoints2[ 2 * q + 1 ];
      aux[ 7 ] += quadWeights2[ q ] * quadPoints2[ 2 * q ] *
          quadPoints2[ 2 * q + 1 ];
      aux[ 8 ] += quadWeights2[ q ] * quadPoints2[ 2 * q + 1 ] *
          quadPoints2[ 2 * q + 1 ];
    }

    for ( int j = 0; j < 9; ++j ) {
      values.push_back( aux[ j ] * area );
      aux[ j ] = 0.0;
    }

  }

  M.setFromTriplets( nNodes, nNodes, rowInd, colInd, values );
}

}
#endif
