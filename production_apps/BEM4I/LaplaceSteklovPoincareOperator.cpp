/*!
 * @file    LaplaceSteklovPoincareOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta 
 * @date    November 24, 2014
 */

#ifdef LAPLACESTEKLOVPOINCAREOPERATOR_H

namespace bem4i {

template<class LO, class SC>
LaplaceSteklovPoincareOperator<LO, SC>::LaplaceSteklovPoincareOperator( ) {
}

template<class LO, class SC>
LaplaceSteklovPoincareOperator<LO, SC>::LaplaceSteklovPoincareOperator(
    const LaplaceSteklovPoincareOperator& orig
    ) {
}

template<class LO, class SC>
LaplaceSteklovPoincareOperator<LO, SC>::LaplaceSteklovPoincareOperator(
    BESpace<LO, SC> * space,
    LinearOperator<LO, SC> * V,
    LinearOperator<LO, SC> * K,
    LinearOperator<LO, SC> * D,
    LinearOperator<LO, SC> * M
    ) {

  this->space = space;
  this->V = V;
  this->K = K;
  this->D = D;
  this->M = M;
  this->epsSingleLayer = 1e-12;
  this->maxItSingleLayer = 1000;
}

template<class LO, class SC>
LaplaceSteklovPoincareOperator<LO, SC>::~LaplaceSteklovPoincareOperator( ) {
}

template<class LO, class SC>
void LaplaceSteklovPoincareOperator<LO, SC>::apply(
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
void LaplaceSteklovPoincareOperator<LO, SC>::applyP1P1(
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
  LO nNodes = this->space->getMesh( )->getNNodes( );

  this->D->apply( x, y, transA, 1.0, 0.0 );

  Vector< LO, SC > hidPlusKx( nElems );
  this->K->apply( x, hidPlusKx, transA, 1.0, 0.0 );
  Vector< LO, SC > Mx( nElems );
  this->M->apply( x, Mx, transA, 1.0, 0.0 );
  hidPlusKx.add( Mx, 0.5 );

  Vector< LO, SC > VinvHidPlusKx( nElems );
  IterativeSolver<LO, SC>::CGSolve( *this->V, hidPlusKx, VinvHidPlusKx,
      this->epsSingleLayer, this->maxItSingleLayer );

  Vector< LO, SC > HidtPlusKtVinvHidPlusKx( nNodes );
  this->K->apply( VinvHidPlusKx, HidtPlusKtVinvHidPlusKx, !transA, 1.0, 0.0 );
  Vector< LO, SC > MtVinvHidPlusKx( nNodes );
  this->M->apply( VinvHidPlusKx, MtVinvHidPlusKx, !transA, 1.0, 0.0 );
  HidtPlusKtVinvHidPlusKx.add( MtVinvHidPlusKx, 0.5 );

  y.add( HidtPlusKtVinvHidPlusKx, 1.0 );

  if ( alpha != 1.0 ) {
    y.scale( alpha );
  }

  if ( beta != 0 ) {
    y.add( yCopy, beta );
  }

}

}
#endif
