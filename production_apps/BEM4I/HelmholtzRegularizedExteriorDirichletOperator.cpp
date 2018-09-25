/*!
 * @file    HelmholtzRegularizedExteriorDirichletOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta 
 * @date    May 15, 2015
 */

#ifdef HELMHOLTZREGULARIZEDEXTERIORDIRICHLETOPERATOR_H

namespace bem4i {

template<class LO, class SC>
HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::
HelmholtzRegularizedExteriorDirichletOperator( ) {
}

template<class LO, class SC>
HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::
HelmholtzRegularizedExteriorDirichletOperator(
    const HelmholtzRegularizedExteriorDirichletOperator& orig
    ) {
}

template<class LO, class SC>
HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::
HelmholtzRegularizedExteriorDirichletOperator(
    BESpace<LO, SC> * space,
    LinearOperator<LO, SC> * V,
    LinearOperator<LO, SC> * K,
    LinearOperator<LO, SC> * Dlap,
    LinearOperator<LO, SC> * M
    ) {

  this->space = space;
  this->kappa = kappa;
  this->V = V;
  this->K = K;
  this->Dlap = Dlap;
  this->M = M;
  this->epsHypersingular = 1e-18;
  this->maxItHypersingular = 1000;
  this->eta = 1.0;
}

template<class LO, class SC>
HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::
~HelmholtzRegularizedExteriorDirichletOperator( ) {
}

template<class LO, class SC>
void HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::apply(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  if ( this->space->getAnsatzFunctionType( ) != p0 ||
      this->space->getTestFunctionType( ) != p0 ) {
    std::cout << "Not implemented!" << std::endl;
    return;
  }

  this->applyP0P0( x, y, transA, alpha, beta );

}

template<class LO, class SC>
void HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::getV(
    const Vector< LO, SC > & density,
    Vector< LO, SC > & v
    ) const {

  if ( this->space->getAnsatzFunctionType( ) != p0 ||
      this->space->getTestFunctionType( ) != p0 ) {
    std::cout << "Not implemented!" << std::endl;
    return;
  }

  this->getVP0P0( density, v );
}

template<class LO, class SC>
void HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::applyP0P0(
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

  this->V->apply( x, y, transA, 1.0, 0.0 );

  // 0.5I+K'
  Vector< LO, SC > hidPlusKx( nNodes );
  Vector< LO, SC > xConj;
  x.copyToConjugate( xConj );
  this->K->apply( xConj, hidPlusKx, !transA, 1.0, 0.0 );
  hidPlusKx.conjugate( );
  Vector< LO, SC > Mx( nNodes );
  this->M->apply( x, Mx, !transA, 1.0, 0.0 );
  hidPlusKx.add( Mx, 0.5 );

  Vector< LO, SC > DinvHidPlusKx( nNodes );
  IterativeSolver<LO, SC>::CGSolve( *this->Dlap, hidPlusKx, DinvHidPlusKx,
      this->epsHypersingular, this->maxItHypersingular );

//  Vector< LO, SC > DinvHidPlusKx( hidPlusKx );
//  FullMatrix<LO,SC> Dcopy(*((FullMatrix<LO,SC>*)this->Dlap));
//  Dcopy.LUSolve( DinvHidPlusKx );

  Vector< LO, SC > HidtPlusKtDinvHidPlusKx( nElems );
  this->K->apply( DinvHidPlusKx, HidtPlusKtDinvHidPlusKx, transA, 1.0, 0.0 );
  Vector< LO, SC > MtVinvHidPlusKx( nElems );
  this->M->apply( DinvHidPlusKx, MtVinvHidPlusKx, transA, 1.0, 0.0 );
  HidtPlusKtDinvHidPlusKx.add( MtVinvHidPlusKx, 0.5 );

  SC iUnit( 0.0, 1.0 );
  y.add( HidtPlusKtDinvHidPlusKx, this->eta * iUnit );

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > 0.0 ) {
    y.add( yCopy, beta );
  }

}

template<class LO, class SC>
void HelmholtzRegularizedExteriorDirichletOperator<LO, SC>::getVP0P0(
    const Vector<LO, SC> & density,
    Vector<LO, SC> & v
    ) const {

  LO nNodes = this->space->getMesh( )->getNNodes( );

  v.resize( nNodes );

  Vector< LO, SC > hidPlusKx( nNodes );
  Vector< LO, SC > densityConj;
  density.copyToConjugate( densityConj );
  this->K->apply( densityConj, hidPlusKx, true, 1.0, 0.0 );
  hidPlusKx.conjugate( );
  Vector< LO, SC > Mx( nNodes );
  this->M->apply( density, Mx, true, 1.0, 0.0 );
  hidPlusKx.add( Mx, 0.5 );

  IterativeSolver<LO, SC>::CGSolve( *this->Dlap, hidPlusKx, v,
      this->epsHypersingular, this->maxItHypersingular );
}

}
#endif
