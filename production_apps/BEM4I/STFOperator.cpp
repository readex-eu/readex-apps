/*!
 * @file    STFOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    June 26, 2015
 */

#ifdef STFOPERATOR_H

namespace bem4i {

template<class LO, class SC>
STFOperator<LO, SC>::
STFOperator( ) {
}

template<class LO, class SC>
STFOperator<LO, SC>::
STFOperator(
    const STFOperator & orig
    ) {
}

template<class LO, class SC>
STFOperator<LO, SC>::STFOperator(
    SurfaceMesh3D< LO, SC > * mesh,
    SC ki,
    SC ke,
    LinearOperator< LO, SC > * Vi,
    LinearOperator< LO, SC > * Ve,
    LinearOperator< LO, SC > * Ki,
    LinearOperator< LO, SC > * Ke,
    LinearOperator< LO, SC > * Di,
    LinearOperator< LO, SC > * De,
    LinearOperator< LO, SC > * M, 
    bool p0p0
    ) {

  this->mesh = mesh;
  this->ki = ki;
  this->ke = ke;
  this->Vi = Vi;
  this->Ve = Ve;
  this->Ki = Ki;
  this->Ke = Ke;
  this->Di = Di;
  this->De = De;
  this->M = M;
  this->p0p0 = p0p0;

}

template<class LO, class SC>
STFOperator<LO, SC>::
~STFOperator( ) {
}

template<class LO, class SC>
bool STFOperator<LO, SC>::LUSolve(
    Vector<LO, SC>& rhs
    ) {

  FullMatrix< LO, SC > * Vmi =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->Vi );
  FullMatrix< LO, SC > * Vme =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->Ve );
  FullMatrix< LO, SC > * Kmi =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->Ki );
  FullMatrix< LO, SC > * Kme =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->Ke );
  FullMatrix< LO, SC > * Dmi =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->Di );
  FullMatrix< LO, SC > * Dme =
      dynamic_cast<FullMatrix< LO, SC > *> ( this->De );

  if ( !Vmi || !Vme || !Kmi || !Kme || !Dmi || !Dme ) {
    std::cout << "LUSolve only for full matrices!" << std::endl;
    return false;
  }

  LO nElems = this->mesh->getNElements( );
  LO nNodes = this->mesh->getNNodes( );

  FullMatrix< LO, SC > wtf( nElems + nNodes, nElems + nNodes );

  SC val;
  
  for ( LO j = 0; j < nElems; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = Vmi->get( i, j );
      val += Vme->get( i, j );
      wtf.set( i, j, val );
    }
  }
  
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nNodes; ++i ) {
      val = Dmi->get( i, j );
      val += Dme->get( i, j );
      wtf.set( i + nElems, j + nElems, val );
    }
  }
  
  for ( LO j = 0; j < nNodes; ++j ) {
    for ( LO i = 0; i < nElems; ++i ) {
      val = Kmi->get( i, j );
      val += Kme->get( i, j );
      wtf.set( i, j + nElems, -val );
      wtf.set( j + nElems, i, val );
    }
  }

  wtf.LUSolve( rhs );
  return true;
}

//template<class LO, class SC>
//void STFOperator<LO, SC>::apply(
//    const Vector<LO, SC> & x,
//    Vector<LO, SC> & y,
//    bool transA,
//    SC alpha,
//    SC beta
//    ) {
//
//  if ( transA ) {
//    std::cout << "Transposed application not implemented!" << std::endl;
//    return;
//  }
//
//  Vector< LO, SC > yCopy;
//
//  if ( std::abs( beta ) > 0.0 ) {
//    y.copy( yCopy );
//  }
//
//  LO nElems = this->mesh->getNElements( );
//  LO nNodes = this->mesh->getNNodes( );
//
//  Vector< LO, SC > t( nElems );
//  memcpy( t.getData( ), x.getData( ), nElems * sizeof ( SC ) );
//  Vector< LO, SC > u( nNodes );
//  memcpy( u.getData( ), x.getData( ) + nElems, nNodes * sizeof ( SC ) );
//
//  Vector< LO, SC > vkie( nElems );
//  this->Vi->apply( t, vkie );
//  this->Ve->apply( t, vkie, false, 1.0, 1.0 );
//
//  this->Ki->apply( u, vkie, false, -1.0, 1.0 );
//  this->Ke->apply( u, vkie, false, -1.0, 1.0 );
//
//  Vector< LO, SC > ktdie( nNodes );
//  this->Ki->apply( t, ktdie, true );
//  this->Ke->apply( t, ktdie, true, 1.0, 1.0 );
//
//  this->Di->apply( u, ktdie, false, 1.0, 1.0 );
//  this->De->apply( u, ktdie, false, 1.0, 1.0 );
//
//  memcpy( y.getData( ), vkie.getData( ), nElems * sizeof ( SC ) );
//  memcpy( y.getData( ) + nElems, ktdie.getData( ), nNodes * sizeof ( SC ) );
//
//  if ( alpha != 1.0 ) {
//    y.scale( alpha );
//  }
//
//  if ( std::abs( beta ) > 0.0 ) {
//    y.add( yCopy, beta );
//  }
//
//}
//
//template<class LO, class SC>
//void STFOperator<LO, SC>::getRHS(
//    const Vector<LO, SC> & g,
//    const Vector<LO, SC> & f,
//    Vector< LO, SC > & rhs
//    ) const {
//
//  LO nElems = this->mesh->getNElements( );
//  LO nNodes = this->mesh->getNNodes( );
//
//  rhs.resize( nNodes + nElems );
//
//  Vector< LO, SC > v1( nElems );
//  Vector< LO, SC > v2( nNodes );
//
//  this->Ve->apply( g, v1 );
//  this->M->apply( f, v1, false, 0.5, 1.0 );
//  this->Ke->apply( f, v1, false, -1.0, 1.0 );
//
//  this->M->apply( g, v2, true );
//  this->Ke->apply( g, v2, true, 1.0, 1.0 );
//  this->De->apply( f, v2, false, 1.0, 1.0 );
//
//  memcpy( rhs.getData( ), v1.getData( ), nElems * sizeof ( SC ) );
//  memcpy( rhs.getData( ) + nElems, v2.getData( ), nNodes * sizeof ( SC ) );
//}
//
//template<class LO, class SC>
//void STFOperator<LO, SC>::getTraces(
//    const Vector< LO, SC > & res,
//    const Vector< LO, SC > & f,
//    const Vector< LO, SC > & g,
//    Vector< LO, SC > & ui,
//    Vector< LO, SC > & ue,
//    Vector< LO, SC > & ti,
//    Vector< LO, SC > & te
//    ) const {
//
//  LO nElems = this->mesh->getNElements( );
//  LO nNodes = this->mesh->getNNodes( );
//
//  ui.resize( nNodes );
//  ue.resize( nNodes );
//  ti.resize( nElems );
//  te.resize( nElems );
//
//  memcpy( ti.getData( ), res.getData( ), nElems * sizeof ( SC ) );
//  memcpy( ui.getData( ), res.getData( ) + nElems, nNodes * sizeof ( SC ) );
//
//  ti.add( g, te, -1.0 );
//  ui.add( f, ue, -1.0 );
//}

template<class LO, class SC>
void STFOperator<LO, SC>::apply(
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

  LO size = 0;
  if (!p0p0) {
    size = this->mesh->getNNodes( );
  } else {
    size = this->mesh->getNElements( );
  }
  

  Vector< LO, SC > t( size );
  memcpy( t.getData( ), x.getData( ), size * sizeof ( SC ) );
  Vector< LO, SC > u( size );
  memcpy( u.getData( ), x.getData( ) + size, size * sizeof ( SC ) );

  Vector< LO, SC > vkie( size );
  this->Vi->apply( t, vkie );
  this->Ve->apply( t, vkie, false, 1.0, 1.0 );

  this->Ki->apply( u, vkie, false, -1.0, 1.0 );
  this->Ke->apply( u, vkie, false, -1.0, 1.0 );

  Vector< LO, SC > ktdie( size );
  this->Ki->apply( t, ktdie, true );
  this->Ke->apply( t, ktdie, true, 1.0, 1.0 );

  this->Di->apply( u, ktdie, false, 1.0, 1.0 );
  this->De->apply( u, ktdie, false, 1.0, 1.0 );

  memcpy( y.getData( ), vkie.getData( ), size * sizeof ( SC ) );
  memcpy( y.getData( ) + size, ktdie.getData( ), size * sizeof ( SC ) );

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > 0.0 ) {
    y.add( yCopy, beta );
  }

}

template<class LO, class SC>
void STFOperator<LO, SC>::getRHS(
    const Vector<LO, SC> & g,
    const Vector<LO, SC> & f,
    Vector< LO, SC > & rhs
    ) const {

  LO size = 0;
  if (!p0p0) {
    size = this->mesh->getNNodes( );
  } else {
    size = this->mesh->getNElements( );
  }

  rhs.resize( size + size );

  Vector< LO, SC > v1( size );
  Vector< LO, SC > v2( size );

  this->Ve->apply( g, v1 );
  this->M->apply( f, v1, false, 0.5, 1.0 );
  this->Ke->apply( f, v1, false, -1.0, 1.0 );

  this->M->apply( g, v2, true );
  this->Ke->apply( g, v2, true, 1.0, 1.0 );
  this->De->apply( f, v2, false, 1.0, 1.0 );

  memcpy( rhs.getData( ), v1.getData( ), size * sizeof ( SC ) );
  memcpy( rhs.getData( ) + size, v2.getData( ), size * sizeof ( SC ) );
}

template<class LO, class SC>
void STFOperator<LO, SC>::getTraces(
    const Vector< LO, SC > & res,
    const Vector< LO, SC > & f,
    const Vector< LO, SC > & g,
    Vector< LO, SC > & ui,
    Vector< LO, SC > & ue,
    Vector< LO, SC > & ti,
    Vector< LO, SC > & te
    ) const {

  LO size = 0;
  if (!p0p0) {
    size = this->mesh->getNNodes( );
  } else {
    size = this->mesh->getNElements( );
  }

  ui.resize( size );
  ue.resize( size );
  ti.resize( size );
  te.resize( size );

  memcpy( ti.getData( ), res.getData( ), size * sizeof ( SC ) );
  memcpy( ui.getData( ), res.getData( ) + size, size * sizeof ( SC ) );

  ti.add( g, te, -1.0 );
  ui.add( f, ue, -1.0 );
}

}
#endif
