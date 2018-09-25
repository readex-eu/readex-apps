#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../BESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../HelmholtzHypersingularP0P0Operator.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzNeumannP0P0(
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

  string filename = "input/cube_12.txt";
  testHelmholtzNeumannP0P0( filename, 2, 0, true, 1, 3, 1 );

  return 0;
}

void testHelmholtzNeumannP0P0(
    string const & filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int nPoints
    ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  timeval start, stop;
  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  gettimeofday( &start, nullptr );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." <<
      std::endl;
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
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

  int quadDisj[2] = { 4, 4 };


  //SC kappa( 2.0, -0.2 );
  //SC kappa = 4.493;
  SC kappa = 1.0;
  SC b = (SCVT) 2.0 * kappa;
  SC a = -std::sqrt( (SCVT) 3.0 ) * kappa;

  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  std::cout << "nested dissection... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  mesh.nestedDissection( tree, 30 );
  mesh.nestedDissection( tree2, 30 );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  FastBESpace< LO, SC > bespaceK( &mesh, p0, p0, &tree2, 1.2, 8, 5, 0 );
  bespaceK.setEpsilonACA( 1e-4 );

  //ACAMatrix< LO, SC >* K = new ACAMatrix< LO, SC >( 0, 0 );
  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, quadOrder, kappa,
      quadType, quadDisj );
  std::cout << "Assembling K ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formK.assemble( *K );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." <<
      std::endl;
  //K->print(); delete K; return;

  IdentityOperator< LO, SC > id( &bespaceK );

  SC iUnit( 0.0, 1.0 );
  //SCVT mult = 2.0 / std::sqrt( 3.0 );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nElems );
  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );
  Vector< LO, SC > neu( nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    dir.set( i, b * y[ 0 ] * std::exp( -a * y[ 1 ] + iUnit * b * y[ 2 ] ) );
    //dir.set( i,y[ 0 ] * y[ 1 ] * y[ 2 ] );

  }

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    mesh.getNormal( i, n );
    neu.set( i, b * std::exp( -a * y[ 1 ] + iUnit * b * y[ 2 ] ) * ( n[ 0 ] -
        a * y[ 0 ] * n[ 1 ] + iUnit * b * y[ 0 ] * n[ 2 ] ) );
    //neu.set( i, n[ 0 ] * y[ 1 ] * y[ 2 ] + y[ 0 ] * n[ 1 ] * y[ 2 ] +
    //  y[ 0 ] * y[ 1 ] * n[ 2 ] );

  }

  std::cout << "Setting up the rhs ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  id.apply( neu, rhs, true, 0.5, 0.0 );
  K->apply( neu, aux, true );
  rhs.add( aux, -1.0 );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." <<
      std::endl;

  delete K;

  FastBESpace< LO, SC > bespaceD( &mesh, p0, p0, &tree, 1.2, 8, 5, 0 );
  bespaceK.setEpsilonACA( 1e-4 );
  ACAMatrix< LO, SC > H1( 0, 0 );
  ACAMatrix< LO, SC > H2( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespaceD, quadOrder,
      kappa, quadType, quadDisj );
  std::cout << "Assembling D ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formD.assembleH1P0P0( H1 );
  formD.assembleH2P0P0( H2 );
  //H2.print();
  //formD.assemble( D, V );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." <<
      std::endl;

  HelmholtzHypersingularP0P0Operator< LO, SC > D( &H1, &H2, mesh );


  std::cout << "Solving the system ... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  Vector< LO, SC > sol( nElems );
  D.FGMRESSolve( rhs, sol, 1e-12, 1000, 1000 );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." <<
      std::endl;


  std::cout << "L2 relative error: " <<
      mesh.l2RelativeErrorConst( sol, 5, kappa ) << "." << std::endl;

  //string nodeNames[] = { };
  string elemNames[] = { "Neumann_anal", "Dirichlet_anal", "Dirichlet_comp", "Dirichlet_err" };
  Vector< LO, SC > err( sol );
  sol.add( dir, err, -1.0 );

  //  Vector< LO, SC >* nodalData[] = { , };
  Vector< LO, SC >* elemData[] = { &neu, &dir, &sol, &err };
  mesh.printParaviewVtu( "output/output.vtu", 0, nullptr, nullptr, 4,
      elemNames, elemData );


  std::cout << "Evaluating representation formula ... ";
  std::cout.flush( );
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
      &neu, kappa, 4 );
  gettimeofday( &start, nullptr );
  formula.evaluate( evalPoint, nPoints, true, res );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
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