#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzDirichletP1P1(
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
  testHelmholtzDirichletP1P1( filename, 3, 0, true, 0, 5, 1 );

  return 0;
}

void testHelmholtzDirichletP1P1(
  string const &filename,
  int refine4,
  int refine9,
  bool mapToUnitBall,
  int quadTypeInt,
  int order,
  int nPoints
  ) {

  typedef complex<double> SC;
  typedef long LO;
  typedef typename SC::value_type SCVT;

  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  ProgressMonitor::init( "Reading mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order " << order << "." <<
      std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order " << order << "." <<
      std::endl;
  }

  //int quadDisj [] = { 4, 4 };
  int * quadDisj = nullptr;

  //SC kappa(2.0, 1.0);
  SC kappa = 2.0;
  SC b = (SCVT) 2.0 * kappa;
  SC a = -std::sqrt( (SCVT) 3.0 ) * kappa;
  ///*

  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  FullMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace11, quadOrder, kappa,
    quadType, quadDisj );
  ProgressMonitor::init( "Assembling V" );
  formV.assemble( V );
  ProgressMonitor::step( );
  //V.print( );
  //return;

  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace11, quadOrder, kappa,
    quadType, quadDisj );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print( );
  //return;

  SC iUnit( 0.0, 1.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nNodes );
  Vector< LO, SC > aux( nNodes );
  Vector< LO, SC > neu( nNodes );
  SCVT y[ 3 ], n[ 3 ];

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
      val = (SCVT) quadWeights[ rhsOrder ][ j ] * b * ya[ 0 ] *
        std::exp( -a * ya[ 1 ] + iUnit * b * ya[ 2 ] );
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
    neu.set( i, b * std::exp( -a * y[ 1 ] + iUnit * b * y[ 2 ] ) * ( n[ 0 ] -
      a * y[ 0 ] * n[ 1 ] + iUnit * b * y[ 0 ] * n[ 2 ] ) );
  }

  ProgressMonitor::init( "Assembling rhs" );
  M11.apply( dir, aux, false, 0.5, 0.0 );
  K.apply( dir, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the system" );
  Vector< LO, SC > sol( nElems );
  V.GMRESSolve( rhs, sol, 1e-12, 1000, 1000 );
  ProgressMonitor::step( );

  std::cout << "L2 relative error: " <<
    mesh.l2RelativeErrorLin( sol, 5, kappa ) << "." << std::endl;

  string nodeNames[] = { "Dirichlet_anal", "Neumann_anal", "Neumann_comp",
    "Neumann_err" };
  Vector< LO, SC > err( sol );
  sol.add( neu, err, -1.0 );
  Vector< LO, SC >* nodalData[] = { &dir, &neu, &sol, &err };
  mesh.printParaviewVtu( "output/output.vtu", 4, nodeNames, nodalData, 0,
    nullptr, nullptr );

  SCVT * evalPoint = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.285910;
    evalPoint[3 * i + 1] = 0.476517;
    evalPoint[3 * i + 2] = 0.667123;
  }
  SC exact = b * evalPoint[ 0 ] * std::exp( -a * evalPoint[ 1 ] + iUnit * b *
    evalPoint[ 2 ] );

  Vector< LO, SC > res( nPoints );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespace11, &dir, &sol,
    kappa, order );
  ProgressMonitor::init( "Representation formula" );
  formula.evaluate( evalPoint, nPoints, true, res );
  ProgressMonitor::step( );

  std::cout << "Point evaluation in ["
    << evalPoint[ 0 ] << ", "
    << evalPoint[ 1 ] << ", "
    << evalPoint[ 2 ] << "]: "
    << res.get( 0 ) << "." << std::endl;

  if ( nPoints > 1 ) {
    std::cout << "Point evaluation in ["
      << evalPoint[ 0 ] << ", "
      << evalPoint[ 1 ] << ", "
      << evalPoint[ 2 ] << "]: "
      << res.get( nPoints - 1 ) << "." << std::endl;
  }

  std::cout << "Absolute error: "
    << std::abs( res.get( 0 ) - exact ) << std::endl;

  std::cout << "Relative error: "
    << std::abs( ( res.get( 0 ) - exact ) / exact ) << std::endl;

  delete [] evalPoint;
}
