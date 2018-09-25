/*!
 * @file    FixedTrackingSubproblem.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    August 12, 2014
 * 
 */

#ifdef FIXEDTRACKINGSUBPROBLEM_H

namespace bem4i {

template< class LO, class SC >
FixedTrackingSubproblem< LO, SC >::FixedTrackingSubproblem( ) {
}

template< class LO, class SC >
FixedTrackingSubproblem< LO, SC >::FixedTrackingSubproblem(
    SC dirichletDataFree,
    SC dirichletDataFixed,
    const Vector< LO, SC > & targetNeumannDataFixed
    ) {

  this->cost = 0.0;
  this->dirichletDataFree = dirichletDataFree;
  this->dirichletDataFixed = dirichletDataFixed;
  this->targetNeumannDataFixed =
      new Vector< LO, SC >( targetNeumannDataFixed );

  this->mesh = nullptr;
  this->primalDirichletData = nullptr;
  this->primalNeumannData = nullptr;
  this->dualDirichletData = nullptr;
  this->dualNeumannData = nullptr;
}

template< class LO, class SC >
FixedTrackingSubproblem< LO, SC >::~FixedTrackingSubproblem( ) {
  if ( this->primalDirichletData ) delete this->primalDirichletData;
  if ( this->primalNeumannData ) delete this->primalNeumannData;
  if ( this->dualDirichletData ) delete this->dualDirichletData;
  if ( this->dualNeumannData ) delete this->dualNeumannData;
  if ( this->targetNeumannDataFixed ) delete this->targetNeumannDataFixed;
  if ( this->mesh ) delete this->mesh;
}

template< class LO, class SC >
void FixedTrackingSubproblem<LO, SC>::printInfo( ) const {
  std::cout << "Class FixedTrackingSubproblem, " << std::endl;
  if ( mesh ) {
    std::cout << "\t";
    this->mesh->printInfo( );
    std::cout << "\tfree nodes: " << this->nNodesFree << ", fixed nodes: " <<
        this->nNodesFixed << "." << std::endl;
  }
  std::cout << "\tJ = " << this->cost << "." << std::endl;
}

template<class LO, class SC>
void FixedTrackingSubproblem<LO, SC>::printVtu(
    const string & filename
    ) const {

  std::vector< string > nodeNames;
  std::vector< string > elemNames;
  std::vector< Vector< LO, SC > * > nodalData;
  std::vector< Vector< LO, SC > * > elemData;

  if ( this->primalDirichletData ) {
    nodeNames.push_back( "Dirichlet" );
    nodalData.push_back( this->primalDirichletData );
  }

  if ( this->dualDirichletData ) {
    nodeNames.push_back( "Dirichlet_adjoint" );
    nodalData.push_back( this->dualDirichletData );
  }

  if ( this->primalNeumannData ) {
    elemNames.push_back( "Neumann" );
    elemData.push_back( this->primalNeumannData );
  }

  if ( this->dualNeumannData ) {
    elemNames.push_back( "Neumann_adjoint" );
    elemData.push_back( this->dualNeumannData );
  }

  Vector< LO, SC > targetNeumann;
  if ( this->targetNeumannDataFixed ) {
    elemNames.push_back( "Neumann_target" );
    targetNeumann.resize( this->nElemsFree + this->nElemsFixed, true );
    for ( LO i = 0; i < this->nElemsFixed; ++i ) {
      targetNeumann.set( this->nElemsFree + i,
          this->targetNeumannDataFixed->get( i ) );
    }
    elemData.push_back( &targetNeumann );
  }

  this->mesh->printParaviewVtu( filename, &nodeNames, &nodalData, &elemNames,
      &elemData );

}

template< class LO, class SC >
bool FixedTrackingSubproblem<LO, SC>::solve( ) {

  LO nElems = this->mesh->getNElements( );

  const SCVT CGeps = 1e-12;
  const LO maxIter = 5000;
  const int order = 3;
  int quadOrder [4] = { order, order, order, order };
  const quadratureType quadrature = SauterSchwab;
  const int maxElemPerCluster = 40;
  const SCVT eta = 1.0;
  const SCVT epsilon = 1e-5;

  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  this->mesh->nestedDissection( tree, maxElemPerCluster );
  this->mesh->nestedDissection( tree2, maxElemPerCluster );
  ProgressMonitor::step( );

  FastBESpace< LO, SC > bespace00( mesh, p0, p0, eta );
  bespace00.setEpsilonACA( epsilon );
  FastBESpace< LO, SC > bespace10( mesh, p1, p0, eta );
  bespace10.setEpsilonACA( epsilon );
  BESpace< LO, SC > bespace11( mesh, p1, p1 );

  std::cout << "Assembling V" << std::endl;
  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC > ( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespace00, quadOrder,
      quadrature );
  formV.assemble( *V );

  std::cout << "Assembling K" << std::endl;
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespace10, quadOrder,
      quadrature );
  formK.assemble( *K );

  IdentityOperator< LO, SC > id01( &bespace10 );

  IdentityOperator< LO, SC > id11( &bespace11 );
  ProgressMonitor::init( "Setting up M^11 identity" );
  id11.setUpIdentityMatrixP1P1( );
  ProgressMonitor::step( );

  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );

  ProgressMonitor::init( "Setting up the rhs for the primal problem" );
  id01.apply( *this->primalDirichletData, aux, false, 0.5, 0.0 );
  K->apply( *this->primalDirichletData, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the primal system" );
  if ( !V->CGSolve( rhs, *this->primalNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  //  ProgressMonitor::init( "Mapping primal Neumann data to nodes" );
  //  Vector< LO, SC > primalNeumannDataNodal( nNodes );
  //  SparseMatrix< LO, SC >* mass11 = id11.getIdentityMatrixP1P1( );
  //  id01.apply( *this->primalNeumannData, primalNeumannDataNodal, true, 1.0,
  //      0.0 );
  //  if ( !mass11->CGSolve( primalNeumannDataNodal, primalNeumannDataNodal, CGeps,
  //      maxIter ) )
  //    return false;
  //  ProgressMonitor::step( );

  ProgressMonitor::init( "Computing cost" );
  Vector< LO, SC > primalNeumannMinusTargetZeroed( nElems, true );
  for ( LO i = this->nElemsFree; i < nElems; ++i ) {
    primalNeumannMinusTargetZeroed.set( i, this->primalNeumannData->get( i ) -
        this->targetNeumannDataFixed->get( i - this->nElemsFree ) );
  }

  this->cost = 0.0;
  for ( LO i = this->nElemsFree; i < nElems; ++i ) {
    this->cost += primalNeumannMinusTargetZeroed.get( i ) *
        primalNeumannMinusTargetZeroed.get( i ) * this->mesh->getElemArea( i );
  }
  this->cost /= 2.0;
  ProgressMonitor::step( );

  ProgressMonitor::init(
      "Computing nodal Dirichlet BC [du/dn] for the adjoint" );
  //primalNeumannDataNodal.copy( *this->dualDirichletData );
  this->elem2NodalAreaWeighted( primalNeumannMinusTargetZeroed,
      *this->dualDirichletData );
  // potrebuju odecist uzlove hodnoty!!!
  //this->dualDirichletData->add( -Q );
  //  for ( LO i = this->nNodesFree; i < nNodes; ++i ) {
  //    this->dualDirichletData->set( i, 0.0 );
  //  }
  ProgressMonitor::step( );

  Vector< LO, SC > rhsAdjoint( nElems );

  ProgressMonitor::init( "Setting up the rhs for the adjoint problem" );
  id01.apply( *this->dualDirichletData, aux, false, 0.5, 0.0 );
  K->apply( *this->dualDirichletData, rhsAdjoint );
  rhsAdjoint.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the adjoint system" );
  if ( !V->CGSolve( rhsAdjoint, *this->dualNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  delete V;
  delete K;

  return true;
}

template< class LO, class SC >
void FixedTrackingSubproblem<LO, SC>::getCost( SCVT & cost ) const {
  cost = this->cost;
}

template< class LO, class SC >
void FixedTrackingSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbation,
    SCVT & dx1,
    SCVT & dx2,
    SCVT & dx3
    ) const {

  dx1 = dx2 = dx3 = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT elemPert;
  SC g;

  for ( LO i = 0; i < this->nElemsFree; ++i ) {
    this->mesh->getElement( i, elem );

    elemPert = ( perturbation.get( elem[ 0 ] ) +
        perturbation.get( elem[ 1 ] ) +
        perturbation.get( elem[ 2 ] ) ) / 3.0;

    if ( std::abs( elemPert ) < EPS ) continue;

    area = this->mesh->getElemArea( i );

    this->mesh->getNormal( i, n );

    g = -this->primalNeumannData->get( i ) * this->dualNeumannData->get( i );

    dx1 += area * g * elemPert * n[ 0 ];
    dx2 += area * g * elemPert * n[ 1 ];
    dx3 += area * g * elemPert * n[ 2 ];
  }

}

template< class LO, class SC >
void FixedTrackingSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbationX1,
    const Vector< LO, SCVT > & perturbationX2,
    const Vector< LO, SCVT > & perturbationX3,
    SCVT & dx
    ) const {

  dx = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT elemPertX1, elemPertX2, elemPertX3;
  SC g;

  for ( LO i = 0; i < this->nElemsFree; ++i ) {

    this->mesh->getElement( i, elem );

    elemPertX1 = ( perturbationX1.get( elem[ 0 ] ) +
        perturbationX1.get( elem[ 1 ] ) +
        perturbationX1.get( elem[ 2 ] ) ) / 3.0;
    elemPertX2 = ( perturbationX2.get( elem[ 0 ] ) +
        perturbationX2.get( elem[ 1 ] ) +
        perturbationX2.get( elem[ 2 ] ) ) / 3.0;
    elemPertX3 = ( perturbationX3.get( elem[ 0 ] ) +
        perturbationX3.get( elem[ 1 ] ) +
        perturbationX3.get( elem[ 2 ] ) ) / 3.0;

    if ( std::sqrt( elemPertX1 * elemPertX1 + elemPertX2 * elemPertX2 + elemPertX3 *
        elemPertX3 ) < EPS ) continue;

    area = this->mesh->getElemArea( i );

    this->mesh->getNormal( i, n );

    g = -this->primalNeumannData->get( i ) * this->dualNeumannData->get( i );

    dx += area * g * ( elemPertX1 * n[ 0 ] +
        elemPertX2 * n[ 1 ] + elemPertX3 * n[ 2 ] );
  }

}

template<class LO, class SC>
void FixedTrackingSubproblem<LO, SC>::getShapeGradient(
    Vector< LO, SC > & grad
    ) const {

  grad.resize( this->nNodesFree );
  SC g;

  Vector< LO, SC > primal, dual;
  this->elem2NodalAreaWeighted( *this->primalNeumannData, primal );
  this->elem2NodalAreaWeighted( *this->dualNeumannData, dual );

  for ( LO i = 0; i < this->nNodesFree; ++i ) {
    g = -primal.get( i ) * dual.get( i );

    grad.set( i, g );
  }
}

template<class LO, class SC>
void FixedTrackingSubproblem<LO, SC>::setProblemData(
    const SurfaceMesh3D< LO, SC > * freeMesh,
    const SurfaceMesh3D< LO, SC > * fixedMesh
    ) {

  if ( fixedMesh ) {
    this->nNodesFixed = fixedMesh->getNNodes( );
    this->nElemsFixed = fixedMesh->getNElements( );
  } else {
    this->nNodesFixed = 0;
    this->nElemsFixed = 0;
  }

  this->nNodesFree = freeMesh->getNNodes( );
  this->nElemsFree = freeMesh->getNElements( );

  LO nNodes = this->nNodesFree + this->nNodesFixed;
  LO nElems = this->nElemsFree + this->nElemsFixed;

  if ( this->mesh ) delete this->mesh;
  this->mesh = new SurfaceMesh3D< LO, SC >( *freeMesh );
  if ( fixedMesh ) {
    this->mesh->append( *fixedMesh );
  }

  if ( this->primalDirichletData ) delete this->primalDirichletData;
  this->primalDirichletData = new Vector< LO, SC >( nNodes );
  for ( LO i = 0; i < this->nNodesFree; ++i ) {
    this->primalDirichletData->set( i, this->dirichletDataFree );
  }
  for ( LO i = 0; i < this->nNodesFixed; ++i ) {
    this->primalDirichletData->set( i + this->nNodesFree,
        this->dirichletDataFixed );
  }

  if ( this->primalNeumannData ) delete this->primalNeumannData;
  this->primalNeumannData = new Vector< LO, SC >( nElems );

  if ( this->dualDirichletData ) delete this->dualDirichletData;
  this->dualDirichletData = new Vector< LO, SC >( nNodes );

  if ( this->dualNeumannData ) delete this->dualNeumannData;
  this->dualNeumannData = new Vector< LO, SC >( nElems );
}

//---PRIVATE-METHODS----------------------------------------------------------//
// todo: move to SurfaceMesh3d

template<class LO, class SC>
void FixedTrackingSubproblem<LO, SC>::elem2NodalAreaWeighted(
    const Vector<LO, SC> & elem,
    Vector<LO, SC> & nodal
    ) const {

  LO nNodes = this->nNodesFree + this->nNodesFixed;
  nodal.resize( nNodes, true );
  SCVT areaSum;
  SCVT sum;
  SCVT area;
  vector< LO > elems;

  for ( LO i = 0; i < nNodes; ++i ) {
    this->mesh->getElements( i, elems );
    areaSum = 0;
    sum = 0;

    for ( int j = 0; j < elems.size( ); ++j ) {
      area = this->mesh->getElemArea( elems[ j ] );
      sum += area * elem.get( elems[ j ] );
      areaSum += area;
    }

    nodal.set( i, sum / areaSum );
  }

}


} // end namespace bem4i

#endif /* FIXEDTRACKINGSUBPROBLEM_H */
