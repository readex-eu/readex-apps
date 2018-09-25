/*!
 * @file    BernoulliSubproblem.cpp
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    May 19, 2014
 *
 */

#ifdef BERNOULLISUBPROBLEM_H

namespace bem4i {

template< class LO, class SC >
BernoulliSubproblem< LO, SC >::BernoulliSubproblem( ) {
}

template< class LO, class SC >
BernoulliSubproblem< LO, SC >::BernoulliSubproblem(
    SC dirichletDataFree,
    const Vector< LO, SC > & dirichletDataFixed,
    SC Q,
    basisType dirBasis,
    basisType neuBasis
    ) {

  this->Q = Q;
  this->cost = 0.0;
  this->dirichletDataFree = dirichletDataFree;
  this->dirichletDataFixed = new Vector< LO, SC >( dirichletDataFixed );
  this->dirBasis = dirBasis;
  this->neuBasis = neuBasis;

  this->myCost = L2Norm;
  this->V = new FullMatrix< LO, SC >( 0, 0 );
  this->K = new FullMatrix< LO, SC >( 0, 0 );
  this->id = nullptr;
  this->mesh = nullptr;
  this->primalDirichletData = nullptr;
  this->primalNeumannData = nullptr;
  this->dualDirichletData = nullptr;
  this->dualNeumannData = nullptr;
}

template< class LO, class SC >
BernoulliSubproblem< LO, SC >::~BernoulliSubproblem( ) {
  if ( this->dirichletDataFixed ) delete this->dirichletDataFixed;
  if ( this->primalDirichletData ) delete this->primalDirichletData;
  if ( this->primalNeumannData ) delete this->primalNeumannData;
  if ( this->dualDirichletData ) delete this->dualDirichletData;
  if ( this->dualNeumannData ) delete this->dualNeumannData;
  if ( this->V ) delete this->V;
  if ( this->K ) delete this->K;
  if ( this->mesh ) delete this->mesh;
}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::printInfo( ) const {
  std::cout << "Class BernoulliSubproblem, " << std::endl;
  if ( mesh ) {
    std::cout << "\t";
    this->mesh->printInfo( );
    std::cout << "\tfree nodes: " << this->nNodesFree << ", fixed nodes: " <<
        this->nNodesFixed << "." << std::endl;
  }
  std::cout << "\tQ = " << this->Q << "." << std::endl;
  std::cout << "\tJ = " << this->cost << "." << std::endl;
}

template<class LO, class SC>
void BernoulliSubproblem<LO, SC>::printVtu(
    const string & filename
    ) const {

  std::vector< string > nodeNames;
  std::vector< string > elemNames;
  std::vector< Vector< LO, SC > * > nodalData;
  std::vector< Vector< LO, SC > * > elemData;

  if ( this->primalDirichletData ) {
    if ( this->dirBasis == p1 ) {
      nodeNames.push_back( "Dirichlet" );
      nodalData.push_back( this->primalDirichletData );
    } else if ( this->dirBasis == p0 ) {
      elemNames.push_back( "Dirichlet" );
      elemData.push_back( this->primalDirichletData );
    }
  }

  if ( this->dualDirichletData ) {
    if ( this->dirBasis == p1 ) {
      nodeNames.push_back( "Dirichlet_adjoint" );
      nodalData.push_back( this->dualDirichletData );
    } else if ( this->dirBasis == p0 ) {
      elemNames.push_back( "Dirichlet_adjoint" );
      elemData.push_back( this->dualDirichletData );
    }
  }

  Vector< LO, SC > curvature;
  if ( this->mesh->curvatureReady( ) ) {
    this->mesh->getCurvatureVector( curvature );
    nodeNames.push_back( "curvature" );
    nodalData.push_back( &curvature );
  }

  if ( this->primalNeumannData ) {
    if ( this->neuBasis == p0 ) {
      elemNames.push_back( "Neumann" );
      elemData.push_back( this->primalNeumannData );
    } else if ( this->neuBasis == p1 ) {
      nodeNames.push_back( "Neumann" );
      nodalData.push_back( this->primalNeumannData );
    }
  }

  if ( this->dualNeumannData ) {
    if ( this->neuBasis == p0 ) {
      elemNames.push_back( "Neumann_adjoint" );
      elemData.push_back( this->dualNeumannData );
    } else if ( this->neuBasis == p1 ) {
      nodeNames.push_back( "Neumann_adjoint" );
      nodalData.push_back( this->dualNeumannData );
    }
  }

  this->mesh->printParaviewVtu( filename, &nodeNames, &nodalData, &elemNames,
      &elemData );

}

template< class LO, class SC >
bool BernoulliSubproblem<LO, SC>::solve( ) {

  LO nNodes = this->mesh->getNNodes( );
  LO nElems = this->mesh->getNElements( );

  const SCVT CGeps = 1e-12;
  const LO maxIter = 5000;
  const int order = 3;
  int quadOrder [ 4 ] = { order, order, order, order };
  const int orderDisj = 4;
  int quadOrderDisj[ 2 ] = { orderDisj, orderDisj };
  quadratureType quadrature = SauterSchwab;

  BESpace< LO, SC > bespace00( mesh, p0, p0 );
  BESpace< LO, SC > bespace11( mesh, p1, p1 );
  BESpace< LO, SC > * bespaceAll = nullptr;
  LO sizeAll = 0;
  LO nFreeAll = 0;
  //LO nFixedAll = 0;
  if ( this->dirBasis == p0 && this->neuBasis == p0 ) {
    bespaceAll = &bespace00;
    sizeAll = nElems;
    nFreeAll = this->nElemsFree;
    //nFixedAll = this->nElemsFixed;
  } else if ( this->dirBasis == p1 && this->neuBasis == p1 ) {
    bespaceAll = &bespace11;
    sizeAll = nNodes;
    nFreeAll = this->nNodesFree;
    //nFixedAll = this->nNodesFixed;
  } else {
    std::cout << "Not implemented" << std::endl;
    exit( 0 );
  }

  BEBilinearFormLaplace1Layer< LO, SC > formV( bespaceAll, quadOrder,
      quadrature, quadOrderDisj );
  ProgressMonitor::init( "Assembling V" );
  formV.assemble( *this->V );
  ProgressMonitor::step( );

  BEBilinearFormLaplace2Layer< LO, SC > formK( bespaceAll, quadOrder,
      quadrature, quadOrderDisj );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( *this->K );
  ProgressMonitor::step( );

  if ( this->id ) delete this->id;
  this->id = new IdentityOperator< LO, SC >( bespaceAll );

  Vector< LO, SC > rhs( sizeAll );
  Vector< LO, SC > aux( sizeAll );

  ProgressMonitor::init( "Setting up the rhs for the primal problem" );
  this->id->apply( *this->primalDirichletData, aux, false, 0.5, 0.0 );
  this->K->apply( *this->primalDirichletData, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the primal system" );
  if ( !this->V->CGSolve( rhs, *this->primalNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  this->primalNeumannData->copy( *this->dualDirichletData );
  this->dualDirichletData->add( -Q );
  for ( LO i = nFreeAll; i < sizeAll; ++i ) {
    this->dualDirichletData->set( i, 0.0 );
  }

  ProgressMonitor::init( "Computing cost" );
  this->updateCost( );
  ProgressMonitor::step( );

  if ( this->myCost == VNorm ) {
    Vector< LO, SC > aux( sizeAll );
    this->V->apply( *this->dualDirichletData, aux );
    this->id->CGSolve( aux, *this->dualDirichletData, CGeps, maxIter );
    for ( LO i = nFreeAll; i < sizeAll; ++i ) {
      this->dualDirichletData->set( i, 0.0 );
    }
  }

  Vector< LO, SC > rhsAdjoint( sizeAll );

  ProgressMonitor::init( "Setting up the rhs for the adjoint problem" );
  this->id->apply( *this->dualDirichletData, aux, false, 0.5, 0.0 );
  this->K->apply( *this->dualDirichletData, rhsAdjoint );
  rhsAdjoint.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the adjoint system" );
  if ( !this->V->CGSolve( rhsAdjoint, *this->dualNeumannData, CGeps, maxIter ) )
    return false;
  ProgressMonitor::step( );

  return true;
}

template<class LO, class SC>
void BernoulliSubproblem<LO, SC>::updateCost( ) {

  if ( this->myCost == L2Norm ) {
    LO elem[ 3 ];
    SC elemData = 0.0;
    this->cost = 0.0;
    for ( LO i = 0; i < this->nElemsFree; ++i ) {
      if ( this->neuBasis == p0 ) {
        elemData = this->dualDirichletData->get( i );
      } else if ( this->neuBasis == p1 ) {
        this->mesh->getElement( i, elem );
        elemData = ( this->dualDirichletData->get( elem[ 0 ] ) +
            this->dualDirichletData->get( elem[ 1 ] ) +
            this->dualDirichletData->get( elem[ 2 ] ) ) / 3.0;
      }
      this->cost += elemData * elemData * this->mesh->getElemArea( i );
    }
    this->cost /= 2.0;

  } else if ( this->myCost == VNorm ) {
    LO size = 0;
    if ( this->neuBasis == p0 ) {
      size = this->mesh->getNElements( );
    } else if ( this->neuBasis == p1 ) {
      size = this->mesh->getNNodes( );
    }
    Vector< LO, SC > aux( size );
    this->V->apply( *this->dualDirichletData, aux );
    this->cost = 0.5 * this->dualDirichletData->dot( aux );

  } else {
    std::cout << "Not implemented" << std::endl;
    exit( 0 );
  }
}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::getCost( SCVT & cost ) const {
  cost = this->cost;
}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbation,
    SCVT & dx1,
    SCVT & dx2,
    SCVT & dx3
    ) const {

  std::cout << "Not used!" << std::endl;
  exit( 0 );

  dx1 = dx2 = dx3 = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT curv;
  SCVT elemPert;
  SCVT elemNeuPrimal = 0.0;
  SCVT elemNeuDual = 0.0;
  SC g;
  if ( !this->mesh->curvatureReady( ) )
    this->mesh->initAdditiveCurvature( );

  for ( LO i = 0; i < this->nElemsFree; ++i ) {
    this->mesh->getElement( i, elem );

    elemPert = ( perturbation.get( elem[ 0 ] ) +
        perturbation.get( elem[ 1 ] ) +
        perturbation.get( elem[ 2 ] ) ) / 3.0;

    if ( std::abs( elemPert ) < EPS ) continue;

    curv = ( this->mesh->getCurvature( elem[ 0 ] ) +
        this->mesh->getCurvature( elem[ 1 ] ) +
        this->mesh->getCurvature( elem[ 2 ] ) ) / 3.0;

    area = this->mesh->getElemArea( i );

    this->mesh->getNormal( i, n );

    if ( this->neuBasis == p0 ) {
      elemNeuPrimal = this->primalNeumannData->get( i );
      elemNeuDual = this->dualNeumannData->get( i );
    } else if ( this->neuBasis == p1 ) {
      elemNeuPrimal = ( this->primalNeumannData->get( elem[ 0 ] ) +
          this->primalNeumannData->get( elem[ 1 ] ) +
          this->primalNeumannData->get( elem[ 2 ] ) ) / 3.0;
      elemNeuDual = ( this->dualNeumannData->get( elem[ 0 ] ) +
          this->dualNeumannData->get( elem[ 1 ] ) +
          this->dualNeumannData->get( elem[ 2 ] ) ) / 3.0;
    }

    g = -elemNeuPrimal * elemNeuDual - ( curv / 2.0 ) *
        ( elemNeuPrimal * elemNeuPrimal - this->Q * this->Q );

    dx1 += area * g * elemPert * n[ 0 ];
    dx2 += area * g * elemPert * n[ 1 ];
    dx3 += area * g * elemPert * n[ 2 ];
  }

}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::getShapeGradient(
    const Vector< LO, SCVT > & perturbationX1,
    const Vector< LO, SCVT > & perturbationX2,
    const Vector< LO, SCVT > & perturbationX3,
    SCVT & dx
    ) const {

  switch ( this->myCost ) {
    case L2Norm:
      this->getShapeGradientL2Norm(
          perturbationX1, perturbationX2, perturbationX3, dx );
      break;
    case VNorm:
      this->getShapeGradientVNorm(
          perturbationX1, perturbationX2, perturbationX3, dx );
      break;
    default:
      std::cout << "Not implemented!" << std::endl;
      break;
  }
}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::getShapeGradientVNorm(
    const Vector< LO, SCVT > & perturbationX1,
    const Vector< LO, SCVT > & perturbationX2,
    const Vector< LO, SCVT > & perturbationX3,
    SCVT & dx
    ) const {

  dx = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT curv;
  SCVT elemPertX1, elemPertX2, elemPertX3;
  SCVT elemNeuPrimal = 0.0;
  SCVT elemNeuDual = 0.0;
  SCVT elemDirDual = 0.0;
  SC g;
  if ( !this->mesh->curvatureReady( ) )
    this->mesh->initAdditiveCurvature( );

  // (-du/dn dp/dn - pHQ/2)(V,n)
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

    if ( std::sqrt( elemPertX1 * elemPertX1 + elemPertX2 * elemPertX2 +
        elemPertX3 * elemPertX3 ) < EPS ) continue;

    curv = ( this->mesh->getCurvature( elem[ 0 ] ) +
        this->mesh->getCurvature( elem[ 1 ] ) +
        this->mesh->getCurvature( elem[ 2 ] ) ) / 3.0;

    area = this->mesh->getElemArea( i );

    this->mesh->getNormal( i, n );

    if ( this->neuBasis == p0 ) {
      elemNeuPrimal = this->primalNeumannData->get( i );
      elemNeuDual = this->dualNeumannData->get( i );
      elemDirDual = this->dualDirichletData->get( i );
    } else if ( this->neuBasis == p1 ) {
      elemNeuPrimal = ( this->primalNeumannData->get( elem[ 0 ] ) +
          this->primalNeumannData->get( elem[ 1 ] ) +
          this->primalNeumannData->get( elem[ 2 ] ) ) / 3.0;
      elemNeuDual = ( this->dualNeumannData->get( elem[ 0 ] ) +
          this->dualNeumannData->get( elem[ 1 ] ) +
          this->dualNeumannData->get( elem[ 2 ] ) ) / 3.0;
      elemDirDual = ( this->dualDirichletData->get( elem[ 0 ] ) +
          this->dualDirichletData->get( elem[ 1 ] ) +
          this->dualDirichletData->get( elem[ 2 ] ) ) / 3.0;
    }

    g = -elemNeuPrimal * elemNeuDual - ( curv * this->Q * elemDirDual ) / 2.0;

    dx += area * g * ( elemPertX1 * n[ 0 ] +
        elemPertX2 * n[ 1 ] + elemPertX3 * n[ 2 ] );
  }

  // 1/2 (du/dn-Q)(1/2+K*)(du/dn-Q)(V,n)
  LO sizeAll = 0;
  LO nFreeAll = 0;
  if ( this->dirBasis == p0 && this->neuBasis == p0 ) {
    sizeAll = this->mesh->getNElements( );
    nFreeAll = this->nElemsFree;
  } else {
    std::cout << "Not implemented" << std::endl;
    exit( 0 );
  }

  Vector< LO, SC > aux1( *this->primalNeumannData );
  Vector< LO, SC > aux2( sizeAll );
  aux1.add( -Q );
  for ( LO i = nFreeAll; i < sizeAll; ++i ) {
    aux1.set( i, 0.0 );
  }

  this->K->apply( aux1, aux2, true );
  this->id->apply( aux1, aux2, true, 0.5, 1.0 );
  for ( LO i = nFreeAll; i < sizeAll; ++i ) {
    aux2.set( i, 0.0 );
  }

  SCVT dot;
  for ( LO i = 0; i < this->nElemsFree; ++i ) {
    this->mesh->getElement( i, elem );
    this->mesh->getNormal( i, n );

    elemPertX1 = ( perturbationX1.get( elem[ 0 ] ) +
        perturbationX1.get( elem[ 1 ] ) +
        perturbationX1.get( elem[ 2 ] ) ) / 3.0;
    elemPertX2 = ( perturbationX2.get( elem[ 0 ] ) +
        perturbationX2.get( elem[ 1 ] ) +
        perturbationX2.get( elem[ 2 ] ) ) / 3.0;
    elemPertX3 = ( perturbationX3.get( elem[ 0 ] ) +
        perturbationX3.get( elem[ 1 ] ) +
        perturbationX3.get( elem[ 2 ] ) ) / 3.0;

    dot = n[ 0 ] * elemPertX1 + n[ 1 ] * elemPertX2 + n[ 2 ] * elemPertX3;
    aux1.set( i, aux1.get( i ) * dot );
  }

  dx += 0.5 * aux1.dot( aux2 );

}

template< class LO, class SC >
void BernoulliSubproblem<LO, SC>::getShapeGradientL2Norm(
    const Vector< LO, SCVT > & perturbationX1,
    const Vector< LO, SCVT > & perturbationX2,
    const Vector< LO, SCVT > & perturbationX3,
    SCVT & dx
    ) const {

  dx = 0.0;
  LO elem[ 3 ];
  SCVT n[ 3 ];
  SCVT area;
  SCVT curv;
  SCVT elemPertX1, elemPertX2, elemPertX3;
  SCVT elemNeuPrimal = 0.0;
  SCVT elemNeuDual = 0.0;
  SC g;
  if ( !this->mesh->curvatureReady( ) )
    this->mesh->initAdditiveCurvature( );

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

    if ( std::sqrt( elemPertX1 * elemPertX1 + elemPertX2 * elemPertX2 +
        elemPertX3 * elemPertX3 ) < EPS ) continue;

    curv = ( this->mesh->getCurvature( elem[ 0 ] ) +
        this->mesh->getCurvature( elem[ 1 ] ) +
        this->mesh->getCurvature( elem[ 2 ] ) ) / 3.0;

    area = this->mesh->getElemArea( i );

    this->mesh->getNormal( i, n );

    if ( this->neuBasis == p0 ) {
      elemNeuPrimal = this->primalNeumannData->get( i );
      elemNeuDual = this->dualNeumannData->get( i );
    } else if ( this->neuBasis == p1 ) {
      elemNeuPrimal = ( this->primalNeumannData->get( elem[ 0 ] ) +
          this->primalNeumannData->get( elem[ 1 ] ) +
          this->primalNeumannData->get( elem[ 2 ] ) ) / 3.0;
      elemNeuDual = ( this->dualNeumannData->get( elem[ 0 ] ) +
          this->dualNeumannData->get( elem[ 1 ] ) +
          this->dualNeumannData->get( elem[ 2 ] ) ) / 3.0;
    }

    g = -elemNeuPrimal * elemNeuDual - ( curv / 2.0 ) *
        ( elemNeuPrimal * elemNeuPrimal - this->Q * this->Q );

    dx += area * g * ( elemPertX1 * n[ 0 ] +
        elemPertX2 * n[ 1 ] + elemPertX3 * n[ 2 ] );
  }

}

template<class LO, class SC>
void BernoulliSubproblem<LO, SC>::getShapeGradient(
    Vector< LO, SC > & grad
    ) const {

  std::cout << "Not used!" << std::endl;
  exit( 0 );

  grad.resize( this->nNodesFree );
  SC g;
  SCVT curv;
  if ( !this->mesh->curvatureReady( ) )
    this->mesh->initAdditiveCurvature( );

  Vector< LO, SC > * primal = nullptr;
  Vector< LO, SC > * dual = nullptr;
  if ( this->neuBasis == p0 ) {
    primal = new Vector< LO, SC >( 0 );
    dual = new Vector< LO, SC >( 0 );
    this->elem2NodalAreaWeighted( *this->primalNeumannData, *primal );
    this->elem2NodalAreaWeighted( *this->dualNeumannData, *dual );
  } else if ( this->neuBasis == p1 ) {
    primal = this->primalNeumannData;
    dual = this->dualNeumannData;
  }

  for ( LO i = 0; i < this->nNodesFree; ++i ) {
    curv = this->mesh->getCurvature( i );

    g = -primal->get( i ) * dual->get( i ) - ( curv / 2.0 ) *
        ( primal->get( i ) * primal->get( i ) - this->Q * this->Q );

    grad.set( i, g );
  }

  if ( this->neuBasis == p0 ) {
    delete primal;
    delete dual;
  }
}

template<class LO, class SC>
void BernoulliSubproblem<LO, SC>::setProblemData(
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

  LO sizeAll = 0;
  LO nFreeAll = 0;
  LO nFixedAll = 0;
  if ( this->dirBasis == p0 && this->neuBasis == p0 ) {
    sizeAll = nElems;
    nFreeAll = this->nElemsFree;
    nFixedAll = this->nElemsFixed;
  } else if ( this->dirBasis == p1 && this->neuBasis == p1 ) {
    sizeAll = nNodes;
    nFreeAll = this->nNodesFree;
    nFixedAll = this->nNodesFixed;
  } else {
    std::cout << "Not implemented" << std::endl;
    exit( 0 );
  }

  if ( this->primalDirichletData ) delete this->primalDirichletData;
  this->primalDirichletData = new Vector< LO, SC >( sizeAll );
  for ( LO i = 0; i < nFreeAll; ++i ) {
    this->primalDirichletData->set( i, this->dirichletDataFree );
  }
  for ( LO i = 0; i < nFixedAll; ++i ) {
    this->primalDirichletData->set( i + nFreeAll,
        this->dirichletDataFixed->get( i ) );
  }

  if ( this->primalNeumannData ) delete this->primalNeumannData;
  this->primalNeumannData = new Vector< LO, SC >( sizeAll );

  if ( this->dualDirichletData ) delete this->dualDirichletData;
  this->dualDirichletData = new Vector< LO, SC >( sizeAll );

  if ( this->dualNeumannData ) delete this->dualNeumannData;
  this->dualNeumannData = new Vector< LO, SC >( sizeAll );
}

//---PRIVATE-METHODS----------------------------------------------------------//

template<class LO, class SC>
void BernoulliSubproblem<LO, SC>::elem2NodalAreaWeighted(
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

#endif /*BERNOULLISUBPROBLEM_H*/

//#ifdef BERNOULLISUBPROBLEM_H
//
////#define BERNOULLISUBPROBLEM_ACA
//
//namespace bem4i {
//
//template< class LO, class SC >
//BernoulliSubproblem< LO, SC >::BernoulliSubproblem( ) {
//}
//
//template< class LO, class SC >
//BernoulliSubproblem< LO, SC >::BernoulliSubproblem(
//    SC dirichletDataFree,
//    const Vector< LO, SC > & dirichletDataFixed,
//    SC Q
//    ) {
//
//  this->Q = Q;
//  this->cost = 0.0;
//  this->dirichletDataFree = dirichletDataFree;
//  this->dirichletDataFixed = new Vector< LO, SC >( dirichletDataFixed );
//
//  this->mesh = nullptr;
//  this->primalDirichletData = nullptr;
//  this->primalNeumannData = nullptr;
//  this->dualDirichletData = nullptr;
//  this->dualNeumannData = nullptr;
//}
//
//template< class LO, class SC >
//BernoulliSubproblem< LO, SC >::~BernoulliSubproblem( ) {
//  if ( this->dirichletDataFixed ) delete this->dirichletDataFixed;
//  if ( this->primalDirichletData ) delete this->primalDirichletData;
//  if ( this->primalNeumannData ) delete this->primalNeumannData;
//  if ( this->dualDirichletData ) delete this->dualDirichletData;
//  if ( this->dualNeumannData ) delete this->dualNeumannData;
//  if ( this->mesh ) delete this->mesh;
//}
//
//template< class LO, class SC >
//void BernoulliSubproblem<LO, SC>::printInfo( ) const {
//  std::cout << "Class BernoulliSubproblem, " << std::endl;
//  if ( mesh ) {
//    std::cout << "\t";
//    this->mesh->printInfo( );
//    std::cout << "\tfree nodes: " << this->nNodesFree << ", fixed nodes: " <<
//        this->nNodesFixed << "." << std::endl;
//  }
//  std::cout << "\tQ = " << this->Q << "." << std::endl;
//  std::cout << "\tJ = " << this->cost << "." << std::endl;
//}
//
//template<class LO, class SC>
//void BernoulliSubproblem<LO, SC>::printVtu(
//    const string & filename
//    ) const {
//
//  std::vector< string > nodeNames;
//  std::vector< string > elemNames;
//  std::vector< Vector< LO, SC > * > nodalData;
//  std::vector< Vector< LO, SC > * > elemData;
//
//  if ( this->primalDirichletData ) {
//    nodeNames.push_back( "Dirichlet" );
//    nodalData.push_back( this->primalDirichletData );
//  }
//
//  if ( this->dualDirichletData ) {
//    nodeNames.push_back( "Dirichlet_adjoint" );
//    nodalData.push_back( this->dualDirichletData );
//  }
//
//  Vector< LO, SC > curvature;
//  if ( this->mesh->curvatureReady( ) ) {
//    this->mesh->getCurvatureVector( curvature );
//    nodeNames.push_back( "curvature" );
//    nodalData.push_back( &curvature );
//  }
//
//  if ( this->primalNeumannData ) {
//    elemNames.push_back( "Neumann" );
//    elemData.push_back( this->primalNeumannData );
//  }
//
//  if ( this->dualNeumannData ) {
//    elemNames.push_back( "Neumann_adjoint" );
//    elemData.push_back( this->dualNeumannData );
//  }
//
//  this->mesh->printParaviewVtu( filename, &nodeNames, &nodalData, &elemNames,
//      &elemData );
//
//}
//
//template< class LO, class SC >
//bool BernoulliSubproblem<LO, SC>::solve( ) {
//
//  timeval start, stop;
//
//  LO nNodes = this->mesh->getNNodes( );
//  LO nElems = this->mesh->getNElements( );
//
//  const SCVT CGeps = 1e-12;
//  const LO maxIter = 5000;
//  const int order = 3;
//  int quadOrder [ 4 ] = { order, order, order, order };
//  const int orderDisj = 4;
//  int quadOrderDisj[ 2 ] = { orderDisj, orderDisj };
//  quadratureType quadrature = SauterSchwab;
//
//#ifdef BERNOULLISUBPROBLEM_ACA
//
//  LO ACAminClusterSize = 30;
//  SCVT ACAepsilon = 1e-4;
//  SCVT ACAeta = 1.2;
//  LO dummy = 0;
//
//  std::cout << "Using ACA, eps =  " << ACAepsilon << ", eta = " << ACAeta <<
//      " , min cluster size = " << ACAminClusterSize << "." << std::endl;
//
//  Tree< BECluster< LO, SC > *, LO > tree, tree2;
//  std::cout << "Nested dissection... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  mesh->nestedDissection( tree, ACAminClusterSize );
//  mesh->nestedDissection( tree2, ACAminClusterSize );
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) <<
//      " secs." << std::endl;
//
//  FastBESpace< LO, SC > bespace00( mesh, p0, p0, &tree, ACAeta, dummy, dummy,
//      dummy );
//  bespace00.setEpsilonACA( ACAepsilon );
//  FastBESpace< LO, SC > bespace10( mesh, p1, p0, &tree2, ACAeta, dummy, dummy,
//      dummy );
//  bespace10.setEpsilonACA( ACAepsilon );
//
//  ACAMatrix< LO, SC > * V = new ACAMatrix< LO, SC >( 0, 0 );
//  ACAMatrix< LO, SC > * K = new ACAMatrix< LO, SC >( 0, 0 );
//#else
//  std::cout << "Using full BEM." << std::endl;
//
//  BESpace< LO, SC > bespace00( mesh, p0, p0 );
//  BESpace< LO, SC > bespace10( mesh, p1, p0 );
//  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC > ( 0, 0 );
//  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
//#endif
//
//  BESpace< LO, SC > bespace11( mesh, p1, p1 );
//
//  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespace00, quadOrder,
//      quadrature, quadOrderDisj );
//  formV.assemble( *V );
//
//  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespace10, quadOrder,
//      quadrature, quadOrderDisj );
//  formK.assemble( *K );
//
//  IdentityOperator< LO, SC > id01( &bespace10 );
//
//  IdentityOperator< LO, SC > id11( &bespace11 );
//  std::cout << "Setting up M^11 identity ... ";
//  std::cout.flush( );
//  SparseMatrix< LO, SC > mass11;
//  gettimeofday( &start, nullptr );
//  id11.assemble( mass11 );
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  Vector< LO, SC > rhs( nElems );
//  Vector< LO, SC > aux( nElems );
//
//  std::cout << "Setting up the rhs for the primal problem ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  id01.apply( *this->primalDirichletData, aux, false, 0.5, 0.0 );
//  K->apply( *this->primalDirichletData, rhs );
//  rhs.add( aux );
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  std::cout << "Solving the primal system ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  if ( !V->CGSolve( rhs, *this->primalNeumannData, CGeps, maxIter ) )
//    return false;
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  std::cout << "Mapping primal Neumann data to nodes ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  Vector< LO, SC > primalNodal( nNodes );
//  id01.apply( *this->primalNeumannData, primalNodal, true, 1.0, 0.0 );
//  if ( !mass11.CGSolve( primalNodal, primalNodal, CGeps, maxIter ) )
//    return false;
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  std::cout << "Computing nodal Dirichlet BC [du/dn] for the adjoint ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  Vector< LO, SC > primalMinusQZeroed( *this->primalNeumannData );
//  primalMinusQZeroed.add( -Q );
//  for ( LO i = this->nElemsFree; i < nElems; ++i ) {
//    primalMinusQZeroed.set( i, 0.0 );
//  }
//
//  this->cost = 0.0;
//  for ( LO i = 0; i < this->nElemsFree; ++i ) {
//    this->cost += primalMinusQZeroed.get( i ) * primalMinusQZeroed.get( i ) *
//        this->mesh->getElemArea( i );
//  }
//  this->cost /= 2.0;
//
//  primalNodal.copy( *this->dualDirichletData );
//  this->dualDirichletData->add( -Q );
//  for ( LO i = this->nNodesFree; i < nNodes; ++i ) {
//    this->dualDirichletData->set( i, 0.0 );
//  }
//
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  Vector< LO, SC > rhsAdjoint( nElems );
//
//  std::cout << "Setting up the rhs for the adjoint problem ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  id01.apply( *this->dualDirichletData, aux, false, 0.5, 0.0 );
//  K->apply( *this->dualDirichletData, rhsAdjoint );
//  rhsAdjoint.add( aux );
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  std::cout << "Solving the adjoint system ... ";
//  std::cout.flush( );
//  gettimeofday( &start, nullptr );
//  if ( !V->CGSolve( rhsAdjoint, *this->dualNeumannData, CGeps, maxIter ) )
//    return false;
//  gettimeofday( &stop, nullptr );
//  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs."
//      << std::endl;
//
//  delete V;
//  delete K;
//
//  return true;
//}
//
//template< class LO, class SC >
//void BernoulliSubproblem<LO, SC>::getCost( SCVT & cost ) const {
//  cost = this->cost;
//}
//
//template< class LO, class SC >
//void BernoulliSubproblem<LO, SC>::getShapeGradient(
//    const Vector< LO, SCVT > & perturbation,
//    SCVT & dx1,
//    SCVT & dx2,
//    SCVT & dx3
//    ) const {
//
//  dx1 = dx2 = dx3 = 0.0;
//  LO elem[ 3 ];
//  SCVT n[ 3 ];
//  SCVT area;
//  SCVT curv;
//  SCVT elemPert;
//  SC g;
//  if ( !this->mesh->curvatureReady( ) )
//    this->mesh->initAdditiveCurvature( );
//
//  for ( LO i = 0; i < this->nElemsFree; ++i ) {
//    this->mesh->getElement( i, elem );
//
//    elemPert = ( perturbation.get( elem[ 0 ] ) +
//        perturbation.get( elem[ 1 ] ) +
//        perturbation.get( elem[ 2 ] ) ) / 3.0;
//
//    if ( std::abs( elemPert ) < EPS ) continue;
//
//    curv = ( this->mesh->getCurvature( elem[ 0 ] ) +
//        this->mesh->getCurvature( elem[ 1 ] ) +
//        this->mesh->getCurvature( elem[ 2 ] ) ) / 3.0;
//
//    area = this->mesh->getElemArea( i );
//
//    this->mesh->getNormal( i, n );
//
//    g = -this->primalNeumannData->get( i ) * this->dualNeumannData->get( i ) -
//        ( curv / 2.0 ) * ( this->primalNeumannData->get( i )
//        * this->primalNeumannData->get( i ) - this->Q * this->Q );
//
//    dx1 += area * g * elemPert * n[ 0 ];
//    dx2 += area * g * elemPert * n[ 1 ];
//    dx3 += area * g * elemPert * n[ 2 ];
//  }
//
//}
//
//template< class LO, class SC >
//void BernoulliSubproblem<LO, SC>::getShapeGradient(
//    const Vector< LO, SCVT > & perturbationX1,
//    const Vector< LO, SCVT > & perturbationX2,
//    const Vector< LO, SCVT > & perturbationX3,
//    SCVT & dx
//    ) const {
//
//  dx = 0.0;
//  LO elem[ 3 ];
//  SCVT n[ 3 ];
//  SCVT area;
//  SCVT curv;
//  SCVT elemPertX1, elemPertX2, elemPertX3;
//  SC g;
//  if ( !this->mesh->curvatureReady( ) )
//    this->mesh->initAdditiveCurvature( );
//
//  for ( LO i = 0; i < this->nElemsFree; ++i ) {
//
//    this->mesh->getElement( i, elem );
//
//    elemPertX1 = ( perturbationX1.get( elem[ 0 ] ) +
//        perturbationX1.get( elem[ 1 ] ) +
//        perturbationX1.get( elem[ 2 ] ) ) / 3.0;
//    elemPertX2 = ( perturbationX2.get( elem[ 0 ] ) +
//        perturbationX2.get( elem[ 1 ] ) +
//        perturbationX2.get( elem[ 2 ] ) ) / 3.0;
//    elemPertX3 = ( perturbationX3.get( elem[ 0 ] ) +
//        perturbationX3.get( elem[ 1 ] ) +
//        perturbationX3.get( elem[ 2 ] ) ) / 3.0;
//
//    if ( std::sqrt( elemPertX1 * elemPertX1 + elemPertX2 * elemPertX2 + elemPertX3 *
//        elemPertX3 ) < EPS ) continue;
//
//    curv = ( this->mesh->getCurvature( elem[ 0 ] ) +
//        this->mesh->getCurvature( elem[ 1 ] ) +
//        this->mesh->getCurvature( elem[ 2 ] ) ) / 3.0;
//
//    area = this->mesh->getElemArea( i );
//
//    this->mesh->getNormal( i, n );
//
//    g = -this->primalNeumannData->get( i ) * this->dualNeumannData->get( i ) -
//        ( curv / 2.0 ) * ( this->primalNeumannData->get( i )
//        * this->primalNeumannData->get( i ) - this->Q * this->Q );
//
//    dx += area * g * ( elemPertX1 * n[ 0 ] +
//        elemPertX2 * n[ 1 ] + elemPertX3 * n[ 2 ] );
//  }
//
//}
//
//template<class LO, class SC>
//void BernoulliSubproblem<LO, SC>::getShapeGradient(
//    Vector< LO, SC > & grad
//    ) const {
//
//  grad.resize( this->nNodesFree );
//  SC g;
//  SCVT curv;
//  if ( !this->mesh->curvatureReady( ) )
//    this->mesh->initAdditiveCurvature( );
//
//  Vector< LO, SC > primal, dual;
//  this->elem2NodalAreaWeighted( *this->primalNeumannData, primal );
//  this->elem2NodalAreaWeighted( *this->dualNeumannData, dual );
//
//  for ( LO i = 0; i < this->nNodesFree; ++i ) {
//    curv = this->mesh->getCurvature( i );
//
//    g = -primal.get( i ) * dual.get( i ) - ( curv / 2.0 ) * ( primal.get( i )
//        * primal.get( i ) - this->Q * this->Q );
//
//    grad.set( i, g );
//  }
//
//}
//
//template<class LO, class SC>
//void BernoulliSubproblem<LO, SC>::setProblemData(
//    const SurfaceMesh3D< LO, SC > * freeMesh,
//    const SurfaceMesh3D< LO, SC > * fixedMesh
//    ) {
//
//  if ( fixedMesh ) {
//    this->nNodesFixed = fixedMesh->getNNodes( );
//    this->nElemsFixed = fixedMesh->getNElements( );
//  } else {
//    this->nNodesFixed = 0;
//    this->nElemsFixed = 0;
//  }
//
//  this->nNodesFree = freeMesh->getNNodes( );
//  this->nElemsFree = freeMesh->getNElements( );
//
//  LO nNodes = this->nNodesFree + this->nNodesFixed;
//  LO nElems = this->nElemsFree + this->nElemsFixed;
//
//  if ( this->mesh ) delete this->mesh;
//  this->mesh = new SurfaceMesh3D< LO, SC >( *freeMesh );
//  if ( fixedMesh ) {
//    this->mesh->append( *fixedMesh );
//  }
//
//  if ( this->primalDirichletData ) delete this->primalDirichletData;
//  this->primalDirichletData = new Vector< LO, SC >( nNodes );
//  for ( LO i = 0; i < this->nNodesFree; ++i ) {
//    this->primalDirichletData->set( i, this->dirichletDataFree );
//  }
//  for ( LO i = 0; i < this->nNodesFixed; ++i ) {
//    this->primalDirichletData->set( i + this->nNodesFree,
//        this->dirichletDataFixed->get( i ) );
//  }
//
//  if ( this->primalNeumannData ) delete this->primalNeumannData;
//  this->primalNeumannData = new Vector< LO, SC >( nElems );
//
//  if ( this->dualDirichletData ) delete this->dualDirichletData;
//  this->dualDirichletData = new Vector< LO, SC >( nNodes );
//
//  if ( this->dualNeumannData ) delete this->dualNeumannData;
//  this->dualNeumannData = new Vector< LO, SC >( nElems );
//}
//
////---PRIVATE-METHODS----------------------------------------------------------//
//
//template<class LO, class SC>
//void BernoulliSubproblem<LO, SC>::elem2NodalAreaWeighted(
//    const Vector<LO, SC> & elem,
//    Vector<LO, SC> & nodal
//    ) const {
//
//  LO nNodes = this->nNodesFree + this->nNodesFixed;
//  nodal.resize( nNodes, true );
//  SCVT areaSum;
//  SCVT sum;
//  SCVT area;
//  vector< LO > elems;
//
//  for ( LO i = 0; i < nNodes; ++i ) {
//    this->mesh->getElements( i, elems );
//    areaSum = 0;
//    sum = 0;
//
//    for ( int j = 0; j < elems.size( ); ++j ) {
//      area = this->mesh->getElemArea( elems[ j ] );
//      sum += area * elem.get( elems[ j ] );
//      areaSum += area;
//    }
//
//    nodal.set( i, sum / areaSum );
//  }
//
//}
//
//
//} // end namespace bem4i
//
//#endif /*BERNOULLISUBPROBLEM_H*/
