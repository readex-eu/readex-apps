/*!
 * @file    HelmholtzRegularizedExteriorNeumannOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta 
 * @date    April 7, 2015
 */

#ifdef HELMHOLTZREGULARIZEDEXTERIORNEUMANNOPERATOR_H

namespace bem4i {

template<class LO, class SC>
HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::
HelmholtzRegularizedExteriorNeumannOperator( ) {
}

template<class LO, class SC>
HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::
HelmholtzRegularizedExteriorNeumannOperator(
    const HelmholtzRegularizedExteriorNeumannOperator& orig
    ) {
}

template<class LO, class SC>
HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::
HelmholtzRegularizedExteriorNeumannOperator(
    BESpace<LO, SC> * space,
    LinearOperator<LO, SC> * Vlap,
    LinearOperator<LO, SC> * K,
    LinearOperator<LO, SC> * D,
    LinearOperator<LO, SC> * M
    ) {

  this->space = space;
  this->kappa = kappa;
  this->Vlap = Vlap;
  this->K = K;
  this->D = D;
  this->M = M;
  this->epsSingleLayer = 1e-18;
  this->maxItSingleLayer = 1000;
  this->eta = 1.0;
}

template<class LO, class SC>
HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::
~HelmholtzRegularizedExteriorNeumannOperator( ) {
}

template<class LO, class SC>
void HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::apply(
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
void HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::getW(
    const Vector< LO, SC > & density,
    Vector< LO, SC > & w
    ) const {

  if ( this->space->getAnsatzFunctionType( ) != p1 ||
      this->space->getTestFunctionType( ) != p1 ) {
    std::cout << "Not implemented!" << std::endl;
    return;
  }

  this->getWP1P1( density, w );
}
/*
template<class LO, class SC>
void HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::applyP1P1(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  if ( transA ) {
    std::cout << "Transposed application not implemented!" << std::endl;
    return;
  }

  Vector< LO, SC > yCopy;

  if ( std::abs( beta ) > 0.0 ) {
    y.copy( yCopy );
  }

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  this->D->apply( x, y, transA, 1.0, 0.0 );

  // -0.5I+K'
  Vector< LO, SC > hidPlusKx( nElems );
  Vector< LO, SC > xConj;
  x.copyToConjugate( xConj );
  this->K->apply( xConj, hidPlusKx, transA, 1.0, 0.0 );
  hidPlusKx.conjugate( );
  Vector< LO, SC > Mx( nElems );
  this->M->apply( x, Mx, transA, 1.0, 0.0 );
  hidPlusKx.add( Mx, -0.5 );

  Vector< LO, SC > VinvHidPlusKx( nElems );
  IterativeSolver<LO, SC>::CGSolve( *this->Vlap, hidPlusKx, VinvHidPlusKx,
      this->epsSingleLayer, this->maxItSingleLayer );

  Vector< LO, SC > HidtPlusKtVinvHidPlusKx( nNodes );
  this->K->apply( VinvHidPlusKx, HidtPlusKtVinvHidPlusKx, !transA, 1.0, 0.0 );
  Vector< LO, SC > MtVinvHidPlusKx( nNodes );
  this->M->apply( VinvHidPlusKx, MtVinvHidPlusKx, !transA, 1.0, 0.0 );
  HidtPlusKtVinvHidPlusKx.add( MtVinvHidPlusKx, -0.5 );

  SC iUnit( 0.0, 1.0 );
  y.add( HidtPlusKtVinvHidPlusKx, this->eta * iUnit );

  if ( alpha != 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > 0.0 ) {
    y.add( yCopy, beta );
  }

}
 */
///*

template<class LO, class SC>
void HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::applyP1P1(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  if ( transA ) {
    std::cout << "Transposed application not implemented!" << std::endl;
    return;
  }

  Vector< LO, SC > yCopy;

  if ( std::abs( beta ) > 0.0 ) {
    y.copy( yCopy );
  }

  LO nElems = this->space->getMesh( )->getNElements( );
  LO nNodes = this->space->getMesh( )->getNNodes( );

  SC iUnit( 0.0, 1.0 );
  Vector< LO, SC > xw( nElems, false );
  Vector< LO, SC > xv( nNodes, false );
  for ( LO i = 0; i < nElems; ++i ) {
    xw.set( i, x.get( i ) );
  }
  for ( LO i = 0; i < nNodes; ++i ) {
    xv.set( i, x.get( i + nElems ) );
  }

  Vector< LO, SC > Vw( nElems );
  Vector< LO, SC > KMw( nNodes );
  Vector< LO, SC > Mw( nNodes );

  this->Vlap->apply( xw, Vw );
  this->K->apply( xw, KMw, true, iUnit * this->eta );
  this->M->apply( xw, Mw, true, - (SCVT) 0.5 * iUnit * this->eta );
  KMw.add( Mw );

  Vector< LO, SC > Dv( nNodes );
  Vector< LO, SC > KMv( nElems );
  Vector< LO, SC > Mv( nElems );

  this->D->apply( xv, Dv );
  this->M->apply( xv, Mv, false, 0.5 );
  xv.conjugate( ); // now yv is invalid, but does not matter
  this->K->apply( xv, KMv, false, -1.0 );
  KMv.conjugate( );
  KMv.add( Mv );
  
  for ( LO i = 0; i < nElems; ++i ) {
    y.set( i, Vw.get( i ) + KMv.get( i ) );
  }
  for ( LO i = 0; i < nNodes; ++i ) {
    y.set( i + nElems, Dv.get( i ) + KMw.get( i ) );
  }

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > 0.0 ) {
    y.add( yCopy, beta );
  }

}
//*/

template<class LO, class SC>
void HelmholtzRegularizedExteriorNeumannOperator<LO, SC>::getWP1P1(
    const Vector<LO, SC> & density,
    Vector<LO, SC> & w
    ) const {

  LO nElems = this->space->getMesh( )->getNElements( );

  w.resize( nElems );

  Vector< LO, SC > hidPlusKx( nElems );
  Vector< LO, SC > densityConj;
  density.copyToConjugate( densityConj );
  this->K->apply( densityConj, hidPlusKx, false, 1.0, 0.0 );
  hidPlusKx.conjugate( );
  Vector< LO, SC > Mx( nElems );
  this->M->apply( density, Mx, false, 1.0, 0.0 );
  hidPlusKx.add( Mx, -0.5 );

  IterativeSolver<LO, SC>::CGSolve( *this->Vlap, hidPlusKx, w,
      this->epsSingleLayer, this->maxItSingleLayer );
}

}
#endif
