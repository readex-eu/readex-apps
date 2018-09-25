/*!
 * @file    LaplaceHypersingularOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta 
 * @date    November 21, 2014
 */

#ifdef LAPLACEHYPERSINGULAROPERATOR

namespace bem4i {

template<class LO, class SC>
LaplaceHypersingularOperator<LO, SC>::LaplaceHypersingularOperator( ) {
}

template<class LO, class SC>
LaplaceHypersingularOperator<LO, SC>::LaplaceHypersingularOperator(
    const LaplaceHypersingularOperator& orig
    ) {
}

template<class LO, class SC>
LaplaceHypersingularOperator<LO, SC>::LaplaceHypersingularOperator(
    BESpace<LO, SC> * space,
    Matrix<LO, SC> * V
    ) {

  this->space = space;
  this->V = V;
  this->T = nullptr;
  this->a = nullptr;
  this->regularized = false;
}

template<class LO, class SC>
LaplaceHypersingularOperator<LO, SC>::~LaplaceHypersingularOperator( ) {
  if ( this->T ) delete this->T;
  if ( this->a ) delete this->a;
}

template<class LO, class SC>
void LaplaceHypersingularOperator<LO, SC>::assembleT( ) {

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  std::vector< LO > rowInd;
  rowInd.reserve( 9 * nElems ); // 3 functions per element, curl has size 3
  std::vector< LO > colInd;
  colInd.reserve( 9 * nElems );
  std::vector< SC > values;
  values.reserve( 9 * nElems );
  LO elem[ 3 ];
  Vector< LO, SCVT > * curls = this->space->getMesh( )->getCurls( );

  for ( LO i = 0; i < nElems; ++i ) {

    this->space->getMesh( )->getElement( i, elem );

    for ( int node = 0; node < 3; ++node ) { // over basis functions
      for ( int ind = 0; ind < 3; ++ind ) { // over components of curl
        rowInd.push_back( ind * nElems + i );
        colInd.push_back( elem[ node ] );
        values.push_back( curls->get( 9 * i + 3 * node + ind ) );
      }
    }
  }

  if ( this->T ) delete T;
  this->T = new SparseMatrix< LO, SC >( 3 * nElems, nNodes, rowInd, colInd,
      values );
}

template<class LO, class SC>
void LaplaceHypersingularOperator<LO, SC>::assembleA( ) {

  if ( this->a ) delete this->a;

  this->a = new Vector< LO, SC >( this->space->getMesh( )->getNNodes( ), true );

  LO elem[ 3 ];
  SCVT integral;
  LO nElems = this->space->getMesh( )->getNElements( );

  for ( LO i = 0; i < nElems; ++i ) {
    this->space->getMesh( )->getElement( i, elem );
    integral = this->space->getMesh( )->getElemArea( i ) / 3.0;
    this->a->add( elem[ 0 ], integral );
    this->a->add( elem[ 1 ], integral );
    this->a->add( elem[ 2 ], integral );
  }
}

//template<class LO, class SC>
//void LaplaceHypersingularOperator<LO, SC>::getMatrix(
//    FullMatrix<LO, SC>& D
//    ) {
//
//  if ( typeid ( *this->V ) != typeid ( FullMatrix< LO, SC > ) ) {
//    std::cout << "Returning empty matrix, this->V is not a FullMatrix."
//        << std::endl;
//    return;
//  }
//
//}

template<class LO, class SC>
void LaplaceHypersingularOperator<LO, SC>::apply(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
) {

  if ( this->space->getAnsatzFunctionType( ) != p1 ||
      this->space->getTestFunctionType( ) != p1 ) {
    std::cout << "Not implemented!" << std::endl;
    return;
  }

  this->applyP1P1( x, y, transA, alpha, beta );

}

template<class LO, class SC>
void LaplaceHypersingularOperator<LO, SC>::applyP1P1(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  Vector< LO, SC > yCopy;

  if ( beta != 0.0 ) {
    y.copy( yCopy );
  }

  LO nElems = this->space->getMesh( )->getNElements( );

  // perform T*x
  Vector< LO, SC > firstRes( 3 * nElems );
  this->T->apply( x, firstRes, false, 1.0, 0.0 );

  // perform diag(V,V,V) * ( T * x )
  Vector< LO, SC > firstResPart( nElems );
  Vector< LO, SC > secondResPart( nElems );
  Vector< LO, SC > secondRes( 3 * nElems );
  for ( LO ind = 0; ind < 3; ++ind ) {

    for ( LO i = 0; i < nElems; ++i ) {
      firstResPart.set( i, firstRes.get( ind * nElems + i ) );
    }

    this->V->apply( firstResPart, secondResPart, transA, 1.0, 0.0 );

    for ( LO i = 0; i < nElems; ++i ) {
      secondRes.set( ind * nElems + i, secondResPart.get( i ) );
    }
  }

  // perform T^T * [ diag(V,V,V) * ( T * x ) ]
  this->T->apply( secondRes, y, true, 1.0, 0.0 );

  // stabilize by D*x + a*a'*x
  if ( this->regularized && this->a ) {
    SC dot = this->a->dot( x );
    y.add( *this->a, dot );
  }

  if ( alpha != 1.0 ) {
    y.scale( alpha );
  }

  if ( beta != 0 ) {
    y.add( yCopy, beta );
  }

}

}
#endif
