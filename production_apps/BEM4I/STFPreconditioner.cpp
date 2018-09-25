/*!
 * @file    STFOperator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    June 26, 2015
 */

#ifdef STFPRECONDITIONER_H

namespace bem4i {

template<class LO, class SC>
STFPreconditioner<LO, SC>::STFPreconditioner( ) {

}

template<class LO, class SC>
STFPreconditioner<LO, SC>::
STFPreconditioner(
    const STFPreconditioner & orig
    ) {
}

template<class LO, class SC>
STFPreconditioner<LO, SC>::STFPreconditioner(
    SurfaceMesh3D<LO, SC> * mesh,
    LinearOperator< LO, SC > * Vi,
    LinearOperator< LO, SC > * Ve,
    LinearOperator< LO, SC > * Ki,
    LinearOperator< LO, SC > * Ke,
    LinearOperator< LO, SC > * Di,
    LinearOperator< LO, SC > * De,
    SparseMatrix< LO, SC > * M,
    bool p0p0
    ) {
  this->mesh = mesh;
  this->Vi = Vi;
  this->Ve = Ve;
  this->Ki = Ki;
  this->Ke = Ke;
  this->Di = Di;
  this->De = De;
  this->M = M;
  this->p0p0 = p0p0;
}


//template<class LO, class SC>
//void STFPreconditioner<LO, SC>::apply(
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
//  Vector< LO, SC > yCopy;
//
//  if ( std::abs( beta ) > 0.0 ) {
//    y.copy( yCopy );
//  }
//
//  LO nElems = this->mesh->getNElements( );
//  LO nNodes = this->mesh->getNNodes( );
//
//  Vector< LO, SC > x1( nNodes );
//  memcpy( x1.getData( ), x.getData( ), nNodes * sizeof ( SC ) );
//  Vector< LO, SC > x2( nElems );
//  memcpy( x2.getData( ), x.getData( ) + nNodes, nElems * sizeof ( SC ) );
//
//  Vector< LO, SC > dktie( nNodes );
//  this->Di->apply( x1, dktie );
//  this->De->apply( x1, dktie, false, 1.0, 1.0 );
//
//  this->Ki->apply( x2, dktie, true, 1.0, 1.0 );
//  this->Ke->apply( x2, dktie, true, 1.0, 1.0 );
//
//  Vector< LO, SC > kvie( nElems );
//  this->Ki->apply( x1, kvie, false );
//  this->Ke->apply( x1, kvie, false, 1.0, 1.0 );
//
//  this->Vi->apply( x2, kvie, false, -1.0, 1.0 );
//  this->Ve->apply( x2, kvie, false, -1.0, 1.0 );
//
//  memcpy( y.getData( ), dktie.getData( ), nNodes * sizeof ( SC ) );
//  memcpy( y.getData( ) + nNodes, kvie.getData( ), nElems * sizeof ( SC ) );
//
//  if ( alpha != 1.0 ) {
//    y.scale( alpha );
//  }
//
//  if ( std::abs( beta ) > 0.0 ) {
//    y.add( yCopy, beta );
//  }
//
//
//}

template<class LO, class SC>
void STFPreconditioner<LO, SC>::apply(
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

  Vector< LO, SC > x1aux( size );
  memcpy( x1aux.getData( ), x.getData( ), size * sizeof ( SC ) );
  Vector< LO, SC > x2aux( size );
  memcpy( x2aux.getData( ), x.getData( ) + size, size * sizeof ( SC ) );
  
  Vector< LO, SC > x1( size );
  Vector< LO, SC > x2( size );
  
  this->M->CGSolve( x1aux, x1, 1e-8, 1000 );
  this->M->CGSolve( x2aux, x2, 1e-8, 1000 );

  Vector< LO, SC > dktie( size );
  this->Di->apply( x1, dktie );
  this->De->apply( x1, dktie, false, 1.0, 1.0 );

  this->Ki->apply( x2, dktie, true, 1.0, 1.0 );
  this->Ke->apply( x2, dktie, true, 1.0, 1.0 );

  Vector< LO, SC > kvie( size );
  this->Ki->apply( x1, kvie, false, -1.0 );
  this->Ke->apply( x1, kvie, false, -1.0, 1.0 );

  this->Vi->apply( x2, kvie, false, 1.0, 1.0 );
  this->Ve->apply( x2, kvie, false, 1.0, 1.0 );
  
  Vector< LO, SC > mkvie( size );
  Vector< LO, SC > mdktie( size );

  this->M->CGSolve( kvie, mkvie, 1e-8, 1000 );
  this->M->CGSolve( dktie, mdktie, 1e-8, 1000 );  

  memcpy( y.getData( ), mdktie.getData( ), size * sizeof ( SC ) );
  memcpy( y.getData( ) + size, mkvie.getData( ), size * sizeof ( SC ) );

  if ( alpha != (SCVT) 1.0 ) {
    y.scale( alpha );
  }

  if ( std::abs( beta ) > 0.0 ) {
    y.add( yCopy, beta );
  }

}

}

#endif
