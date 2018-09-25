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
#include "../../RepresentationFormulaLaplace.h"

using namespace std;
using namespace bem4i;

void testLaplaceDirichletP1P1(
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
  testLaplaceDirichletP1P1( filename, 3, 0, true, 0, 5, 1 );

  return 0;
}

void testLaplaceDirichletP1P1(
  string const &filename,
  int refine4,
  int refine9,
  bool mapToUnitBall,
  int quadTypeInt,
  int order,
  int nPoints ) {

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

  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );

  //int quadFar [] = { 4, 4 };
  int * quadFar = nullptr;

  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  FullMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormLaplace1Layer< LO, SC > formV( &bespace11, quadOrder, quadType,
    quadFar, false );
  ProgressMonitor::init( "V" );
  formV.assemble( V );
  ProgressMonitor::step( );

  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespace11, quadOrder, quadType,
    quadFar );
  ProgressMonitor::init( "K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print( );

  IdentityOperator< LO, SC > id( &bespace11 );

  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nNodes );
  Vector< LO, SC > aux( nNodes );
  Vector< LO, SC > neu( nNodes );
  SC y[ 3 ], n[ 3 ];

  IdentityOperator< LO, SC > id11( &bespace11 );
  SparseMatrix< LO, SC > M11;
  id11.assemble( M11 );
  Vector< LO, SC > dirRhs( nNodes );
  int rhsOrder = 5;
  int qSize = quadSizes[ rhsOrder ];
  SCVT * quadNodes = new SCVT[ 3 * qSize ];
  SCVT * ya;
  SC val, val1, val2, val3;
  SCVT area;
  LO elem[ 3 ];
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];
  for ( LO i = 0; i < nElems; ++i ) {
    val1 = val2 = val3 = 0.0;
    mesh.getElement( i, elem );
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getQuadratureNodes( x1, x2, x3, quadNodes, rhsOrder );
    area = mesh.getElemArea( i );
    for ( LO j = 0; j < qSize; ++j ) {
      ya = quadNodes + 3 * j;
      val = (SCVT) quadWeights[ rhsOrder ][ j ] * ( 1.0 + ya[0] ) *
        std::exp( 2.0 * M_PI * ya[1] ) * std::cos( 2.0 * M_PI * ya[2] );
      val1 += val * ( (SCVT) 1.0 - (SCVT) quadPoints[ rhsOrder ][ 2 * j ] -
        (SCVT) quadPoints[ rhsOrder ][ 2 * j + 1 ] );
      val2 += val * (SCVT) quadPoints[ rhsOrder ][ 2 * j ];
      val3 += val * (SCVT) quadPoints[ rhsOrder ][ 2 * j + 1 ];
    }
    val1 *= area;
    val2 *= area;
    val3 *= area;
    dirRhs.add( elem[ 0 ], val1 );
    dirRhs.add( elem[ 1 ], val2 );
    dirRhs.add( elem[ 2 ], val3 );
  }
  delete [] quadNodes;
  M11.CGSolve( dirRhs, dir, 1e-9 );

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, y );
    mesh.getNormalNodal( i, n );
    neu.set( i, std::exp( 2.0 * M_PI * y[1] )
      * ( n[0] * std::cos( 2.0 * M_PI * y[2] )
      + 2.0 * M_PI * ( 1.0 + y[0] ) * n[1] * std::cos( 2.0 * M_PI * y[2] )
      - 2.0 * M_PI * ( 1.0 + y[0] ) * n[2] * std::sin( 2.0 * M_PI * y[2] ) ) );
  }

  ProgressMonitor::init( "Setting up the rhs" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  K.apply( dir, rhs );
  // EXTERIOR ONLY
  //rhs.scale( -1.0 );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> rhs2( rhs );
  rhs.setAll( 0.0 );
  V.CGSolve( rhs2, rhs, CGeps );
  ProgressMonitor::step( );

  std::cout << "L2 relative error: " << mesh.l2RelativeErrorLin( rhs )
    << "." << std::endl;

  string nodeNames[] = { "Dirichlet_anal", "Neumann_anal", "Neumann_comp",
    "Neumann_err" };
  Vector< LO, SC > err( rhs.getLength( ) );
  rhs.add( neu, err, -1.0 );
  Vector< LO, SC >* nodalData[] = { &dir, &neu, &rhs, &err };
  mesh.printParaviewVtu( "output/output.vtu", 4, nodeNames, nodalData, 0,
    nullptr, nullptr );

  SC * evalPoint = new SC[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.250685;
    evalPoint[3 * i + 1] = 0.417808;
    evalPoint[3 * i + 2] = 0.584932;
  }
  SC exact = ( 1.0 + evalPoint[ 0 ] )
    * std::exp( 2.0 * M_PI * evalPoint[ 1 ] )
    * std::cos( 2.0 * M_PI * evalPoint[ 2 ] );
  Vector< LO, SC > res( nPoints );

  RepresentationFormulaLaplace<LO, SC> formula( &bespace11, &dir, &rhs );
  ProgressMonitor::init( "Evaluating representation formula" );
  formula.evaluate( evalPoint, nPoints, true, res );
  ProgressMonitor::step( );
  std::cout << "Point evaluation in ["
    << evalPoint[ 0 ] << ", "
    << evalPoint[ 1 ] << ", "
    << evalPoint[ 2 ] << "]: "
    << res.get( 0 ) << "." << std::endl;

  std::cout << "Absolute error: "
    << std::abs( res.get( 0 ) - exact ) << std::endl;

  std::cout << "Relative error: "
    << std::abs( ( res.get( 0 ) - exact ) / exact ) << std::endl;

  delete [] evalPoint;
}