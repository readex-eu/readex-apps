/*!
 * @file    HelmholtzNeumannRobinOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    July 30, 2015
 */

#ifdef HELMHOLTZNEUMANNROBINOPERATOR_H

namespace bem4i {

template<class LO, class SC>
HelmholtzNeumannRobinOperator<LO, SC>::
HelmholtzNeumannRobinOperator( ) {
}

template<class LO, class SC>
HelmholtzNeumannRobinOperator<LO, SC>::
HelmholtzNeumannRobinOperator(
    const HelmholtzNeumannRobinOperator & orig
    ) {
}

template<class LO, class SC>
HelmholtzNeumannRobinOperator<LO, SC>::HelmholtzNeumannRobinOperator(
    SurfaceMesh3D< LO, SC > * mesh,
    std::vector< LO > & inE,
    std::vector< LO > & outE,
    SC kappa,
    LinearOperator< LO, SC > * V,
    LinearOperator< LO, SC > * K,
    LinearOperator< LO, SC > * D,
    LinearOperator< LO, SC > * M,
    basisType type
    ) {

  this->type = type;
  this->ubar = 300;
  this->mesh = mesh;
  this->inE = inE;
  this->outE = outE;
  this->kappa = kappa;
  this->V = V;
  this->K = K;
  this->D = D;
  this->M = M;
  this->MTilde = nullptr;
  this->assembleMTilde( );
}

template<class LO, class SC>
HelmholtzNeumannRobinOperator<LO, SC>::
~HelmholtzNeumannRobinOperator( ) {
  if ( this->MTilde ) delete this->MTilde;
}

template<class LO, class SC>
bool HelmholtzNeumannRobinOperator<LO, SC>::LUSolve(
    Vector<LO, SC> & rhs
    ) {

  FullMatrix< LO, SC > * Vm =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->V );
  FullMatrix< LO, SC > * Km =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->K );
  FullMatrix< LO, SC > * Dm =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->D );
  SparseMatrix< LO, SC > * Mm =
      dynamic_cast<SparseMatrix< LO, SC > *> ( this->M );

  if ( !Vm || !Km || !Dm || !Mm ) {
    std::cout << "LUSolve only for full matrices!" << std::endl;
    return false;
  }
  
  LO size = 0;
  if ( this->type == p0 ) {
    size = this->mesh->getNElements( );
  } else if ( this->type == p1 ) {
    size = this->mesh->getNNodes( );
  }

  FullMatrix< LO, SC > * KId = new FullMatrix< LO, SC >( *Km );  
  KId->add( *Mm, 0.5 );
  FullMatrix< LO, SC > * Vinv = new FullMatrix< LO, SC >( *KId );
  FullMatrix< LO, SC > * Vcopy = new FullMatrix< LO, SC >( *Vm );
  Vcopy->LUSolve( *Vinv, size );
  delete Vcopy;
  FullMatrix< LO, SC > * S = new FullMatrix< LO, SC >( size, size );
  S->multiply( *KId, *Vinv, true, false, 1.0, 0.0 );
  delete KId;
  delete Vinv;
  
  S->add( *Dm );
  S->add( *this->MTilde );
  
  S->LUSolve( rhs );
  
  delete S;
  
  return true;
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::apply(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  SCVT prec = 1e-12;
  LO maxit = 2000;
  LO restarts = maxit;

  LO size = 0;
  if ( this->type == p0 ) {
    size = this->mesh->getNElements( );
  } else if ( this->type == p1 ) {
    size = this->mesh->getNNodes( );
  }

  Vector< LO, SC > yCopy;

  if ( std::abs( beta ) > 0.0 ) {
    y.copy( yCopy );
  }

  // Dx
  this->D->apply( x, y );
  // MTildex
  this->MTilde->apply( x, y, false, 1.0, 1.0 );

  // (0.5M+K)x
  Vector< LO, SC > aux1( size );
  this->M->apply( x, aux1, false, 0.5 );
  this->K->apply( x, aux1, false, 1.0, 1.0 );

  // (V^-1)aux1
  Vector< LO, SC > aux2( size );
  IterativeSolver<LO, SC>::GMRESSolve( *this->V, aux1, aux2, prec, maxit,
      restarts );

  this->M->apply( aux2, y, true, 0.5, 1.0 );
  this->K->apply( aux2, y, true, 1.0, 1.0 );

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > (SCVT) 0.0 ) {
    y.add( yCopy, beta );
  }

}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::assembleMTilde( ) {
  if ( this->type == p1 ) {
    this->assembleMTildeP1P1( );
  } else if
    ( this->type == p0 ) {
    this->assembleMTildeP0P0( );
  }
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::assembleMTildeP1P1( ) {

  if ( this->MTilde ) delete this->MTilde;
  this->MTilde = new SparseMatrix< LO, SC >( 0, 0 );

  LO nNodes = this->mesh->getNNodes( );

  LO size = 9 * ( this->inE.size( ) + this->outE.size( ) );
  std::vector< LO > rows;
  std::vector< LO > cols;
  std::vector< SC > vals;
  rows.reserve( size );
  cols.reserve( size );
  vals.reserve( size );

  SC f = SC( (SCVT) 0.0, (SCVT) 1.0 ) * this->kappa;
  LO ind[ 3 ];
  SCVT area;
  SCVT aux[ 9 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  for ( auto it = this->inE.begin( ); it != this->inE.end( ); ++it ) {
    this->mesh->getElement( *it, ind );
    area = this->mesh->getElemArea( *it );

    for ( int left = 0; left < 3; ++left ) {
      for ( int right = 0; right < 3; ++right ) {
        cols.push_back( ind[ left ] );
        rows.push_back( ind[ right ] );
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
      vals.push_back( aux[ j ] * area * f );
      aux[ j ] = 0.0;
    }
  }

  for ( auto it = this->outE.begin( ); it != this->outE.end( ); ++it ) {
    this->mesh->getElement( *it, ind );
    area = this->mesh->getElemArea( *it );

    for ( int left = 0; left < 3; ++left ) {
      for ( int right = 0; right < 3; ++right ) {
        cols.push_back( ind[ left ] );
        rows.push_back( ind[ right ] );
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
      vals.push_back( aux[ j ] * area * f );
      aux[ j ] = 0.0;
    }
  }

  this->MTilde->setFromTriplets( nNodes, nNodes, rows, cols, vals );
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::assembleMTildeP0P0( ) {

  if ( this->MTilde ) delete this->MTilde;
  this->MTilde = new SparseMatrix< LO, SC >( 0, 0 );

  LO nElems = this->mesh->getNElements( );

  LO size = this->inE.size( ) + this->outE.size( );
  std::vector< LO > rows;
  std::vector< LO > cols;
  std::vector< SC > vals;
  rows.reserve( size );
  cols.reserve( size );
  vals.reserve( size );

  SC f = SC( (SCVT) 0.0, (SCVT) 1.0 ) * this->kappa;
  SCVT area;

  for ( auto it = this->inE.begin( ); it != this->inE.end( ); ++it ) {
    area = this->mesh->getElemArea( *it );
    cols.push_back( *it );
    rows.push_back( *it );
    vals.push_back( area * f );
  }

  for ( auto it = this->outE.begin( ); it != this->outE.end( ); ++it ) {
    area = this->mesh->getElemArea( *it );
    cols.push_back( *it );
    rows.push_back( *it );
    vals.push_back( area * f );
  }

  this->MTilde->setFromTriplets( nElems, nElems, rows, cols, vals );
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::getRHS(
    Vector< LO, SC > & rhs
    ) const {

  SC g = (SCVT) 2.0 * SC( (SCVT) 0.0, (SCVT) 1.0 ) * this->kappa * this->ubar;

  if ( this->type == p1 ) {
    LO nNodes = this->mesh->getNNodes( );
    rhs.resize( nNodes, true );

    LO ind[ 3 ];
    SCVT areaSc;
    for ( auto it = this->inE.begin( ); it != this->inE.end( ); ++it ) {
      this->mesh->getElement( *it, ind );
      areaSc = this->mesh->getElemArea( *it ) / 3.0;
      rhs.add( ind[ 0 ], g * areaSc );
      rhs.add( ind[ 1 ], g * areaSc );
      rhs.add( ind[ 2 ], g * areaSc );
    }

  } else if ( this->type == p0 ) {
    LO nElems = this->mesh->getNElements( );
    rhs.resize( nElems );

    SCVT area;
    for ( auto it = this->inE.begin( ); it != this->inE.end( ); ++it ) {
      area = this->mesh->getElemArea( *it );
      rhs.set( *it, g * area );
    }
  }
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::solveDirichletProblem(
    const Vector< LO, SC > & u,
    Vector< LO, SC > & dudn
    ) const {

  SCVT prec = 1e-8;
  LO maxit = 2000;
  LO restarts = maxit;

  LO size = 0;
  if ( this->type == p0 ) {
    size = this->mesh->getNElements( );
  } else if ( this->type == p1 ) {
    size = this->mesh->getNNodes( );
  }

  Vector< LO, SC > rhs( size );
  this->M->apply( u, rhs, false, 0.5 );
  this->K->apply( u, rhs, false, 1.0, 1.0 );

  dudn.resize( size );
  IterativeSolver<LO, SC>::GMRESSolve( *this->V, rhs, dudn, prec, maxit,
      restarts );
}

template<class LO, class SC>
void HelmholtzNeumannRobinOperator<LO, SC>::evalInside(
    SurfaceMesh3D<LO, SC>& grid,
    Vector<LO, SC>& u,
    Vector<LO, SC>& dudn,
    Vector<LO, SC>& res
    ) const {

  int order = 4;

  LO nPoints = grid.getNNodes( );
  res.resize( nPoints );
  BESpace< LO, SC > bespace( this->mesh, this->type, this->type );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespace, &u, &dudn,
      this->kappa, order );
  formula.evaluate( grid, true, res );
}

}
#endif
