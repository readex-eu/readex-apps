#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorLaplace.h"
#include "../../BEBilinearFormLaplace1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"
#include "../../IdentityOperator.h"

using namespace std;
using namespace bem4i;

void testLaplaceDirichletP0P0(
  string const &filename,
  int refine4,
  int refine9,
  bool mapToUnitBall,
  int quadTypeInt,
  int order,
  int nPoints
  );

int main(
  int argc,
  char** argv
  ) {
  intro( );

  string filename = "input/icosahedron.txt";
  testLaplaceDirichletP0P0( filename, 3, 0, true, 0, 5, 1 );

  return 0;
}

void testLaplaceDirichletP0P0(
  string const &filename,
  int refine4,
  int refine9,
  bool mapToUnitBall,
  int quadTypeInt,
  int order,
  int nPoints
  ) {

  typedef double SC;
  typedef double SCVT;
  typedef long LO;

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  ProgressMonitor::init( "Initializing mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  mesh.printInfo( );
  ProgressMonitor::step( );

  SC CGeps = 1e-12;

  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order " << order << "."
      << std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order " << order << "."
      << std::endl;
  }

  LO nElems = mesh.getNElements( );

  //int quadFar [] = { 4, 4 };
  int * quadFar = nullptr;

  ///*
  BESpace< LO, SC > bespaceV( &mesh, p0, p0 );
  FullMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespaceV, quadOrder, quadType,
    quadFar, false );
  ProgressMonitor::init( "V" );
  formV.assemble( V );
  ProgressMonitor::step( );
  //V.print( );
  //return;

  BESpace< LO, SC > bespaceK( &mesh, p0, p0 );
  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespaceK, quadOrder, quadType,
    quadFar );
  ProgressMonitor::init( "K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print( );
  //return;

  IdentityOperator< LO, SC > id( &bespaceK );

  Vector< LO, SC > dir( nElems );
  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );
  Vector< LO, SC > neu( nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  int rhsOrder = 5;
  int qSize = quadSizes[ rhsOrder ];
  SCVT * quadNodes = new SCVT[ 3 * qSize ];
  SCVT * ya;
  SC val;

  for ( LO i = 0; i < nElems; ++i ) {
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getQuadratureNodes( x1, x2, x3, quadNodes, rhsOrder );
    val = 0.0;
    for ( LO j = 0; j < qSize; ++j ) {
      ya = quadNodes + 3 * j;
      val += (SCVT) quadWeights[ rhsOrder ][ j ] * ( 1.0 + ya[0] ) *
        std::exp( 2.0 * M_PI * ya[1] ) * std::cos( 2.0 * M_PI * ya[2] );
    }
    dir.set( i, val );
  }
  delete [] quadNodes;

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    mesh.getNormal( i, n );

    neu.set( i, std::exp( 2.0 * M_PI * y[1] )
      * ( n[0] * std::cos( 2.0 * M_PI * y[2] )
      + 2.0 * M_PI * ( 1.0 + y[0] ) * n[1] * std::cos( 2.0 * M_PI * y[2] )
      - 2.0 * M_PI * ( 1.0 + y[0] ) * n[2] * std::sin( 2.0 * M_PI * y[2] ) ) );
  }

  ProgressMonitor::init( "Setting up the rhs" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  K.apply( dir, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> rhs2( rhs );
  rhs.setAll( 0.0 );
  V.CGSolve( rhs2, rhs, CGeps, 2000 );
  ProgressMonitor::step( );

  std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( rhs )
    << "." << std::endl;

}