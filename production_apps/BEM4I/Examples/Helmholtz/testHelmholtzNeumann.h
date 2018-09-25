#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzNeumann(
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
  intro();
  
  string filename = "input/icosahedron.txt";
  testHelmholtzNeumann( filename, 3, 0, true, 0, 5, 1 );
  
  return 0;
}

void testHelmholtzNeumann(
    string const & filename,
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

  //int quadDisjoint[] = { 4, 4 };
  int * quadDisjoint = nullptr;

  //SC kappa( 2.0, -0.2 );
  //SC kappa = 4.493;
  //SC kappa = 4.5;
  SC kappa = 2.0;
  SC b = (SCVT) 2.0 * kappa;
  SC a = -std::sqrt( (SCVT) 3.0 ) * kappa;

  BESpace< LO, SC > bespaceK( &mesh, p1, p0 );
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, quadOrder, kappa,
      quadType, quadDisjoint );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( *K );
  ProgressMonitor::step( );
  //K->print();

  IdentityOperator< LO, SC > id( &bespaceK );

  SC iUnit( 0.0, 1.0 );
  //SCVT mult = 2.0 / std::sqrt( 3.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nNodes );
  Vector< LO, SC > aux( nNodes );
  Vector< LO, SC > neu( nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ], n[ 3 ];

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    dir.set( i, b * x1[ 0 ] * std::exp( -a * x1[ 1 ] + iUnit * b * x1[ 2 ] ) );
    //dir.set( i,x1[ 0 ] * x1[ 1 ] * x1[ 2 ] );
  }

  int rhsOrder = 5;
  int qSize = quadSizes[ rhsOrder ];
  SCVT * quadNodes = new SCVT[ 3 * qSize ];
  SC val = 0.0;
  SCVT * ya;
  for ( LO i = 0; i < nElems; ++i ) {
    mesh.getNodes( i, x1, x2, x3 );
    mesh.getNormal( i, n );
    mesh.getQuadratureNodes( x1, x2, x3, quadNodes, rhsOrder );
    val = 0.0;
    for ( LO j = 0; j < qSize; ++j ) {
      ya = quadNodes + 3 * j;
      val += (SCVT) quadWeights[ rhsOrder ][ j ] *
          ( b * std::exp( -a * ya[ 1 ] + iUnit * b * ya[ 2 ] ) *
          ( n[ 0 ] - a * ya[ 0 ] * n[ 1 ] + iUnit * b * ya[ 0 ] * n[ 2 ] ) );
    }
    neu.set( i, val );
  }
  delete [] quadNodes;

  ProgressMonitor::init( "Setting up rhs" );
  id.apply( neu, rhs, true, 0.5, 0.0 );
  K->apply( neu, aux, true );
  rhs.add( aux, -1.0 );
  ProgressMonitor::step( );

  delete K;

  BESpace< LO, SC > bespaceD( &mesh, p1, p1 );
  FullMatrix< LO, SC > D( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespaceD, quadOrder,
      kappa, quadType, quadDisjoint );
  ProgressMonitor::init( "Assembling D" );
  formD.assemble( D );
  //formD.assemble( D, V );
  ProgressMonitor::step( );
  //D.print( );
  //return;

  ProgressMonitor::init( "Solving the system" );
  Vector< LO, SC > sol( nNodes );
  D.GMRESSolve( rhs, sol, 1e-12, 1000, 1000 );
  ProgressMonitor::step( );

  std::cout << "L2 relative error: " <<
      mesh.l2RelativeErrorLin( sol, 5, kappa ) << "." << std::endl;

  string nodeNames[] = { "Dirichlet_anal", "Dirichlet_comp", "Dirichlet_err" };
  string elemNames[] = { "Neumann_anal" };
  Vector< LO, SC > err( sol );
  sol.add( dir, err, -1.0 );
  Vector< LO, SC >* nodalData[] = { &dir, &sol, &err };
  Vector< LO, SC >* elemData[] = { &neu };
  mesh.printParaviewVtu( "output/output.vtu", 3, nodeNames, nodalData, 1,
      elemNames, elemData );


  SCVT * evalPoint = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {
    evalPoint[3 * i] = 0.285910;
    evalPoint[3 * i + 1] = 0.476517;
    evalPoint[3 * i + 2] = 0.667123;
  }
  SC exact = b * evalPoint[ 0 ] * std::exp( -a * evalPoint[ 1 ] + iUnit * b *
      evalPoint[ 2 ] );
  Vector< LO, SC > res( nPoints );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &sol,
      &neu, kappa, order );
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

  std::cout << "Exact solution: " << exact << std::endl;
  std::cout << "Absolute error: "
      << std::abs( res.get( 0 ) - exact ) << std::endl;

  std::cout << "Relative error: "
      << std::abs( ( res.get( 0 ) - exact ) / exact ) << std::endl;

  delete [] evalPoint;
}
