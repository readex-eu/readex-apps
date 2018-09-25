#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzNeumannACA(
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
  testHelmholtzNeumannACA( filename, 1, 0, true, 1, 3, 1 );

  return 0;
}

void testHelmholtzNeumannACA(
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

  int quadDisj [] = { 4, 4 };

  SCVT two = 2.0;
  SCVT three = 3.0;
  SCVT four = 4.0;

  SC kappa = 2.0;

  SCVT ACAeps = 1e-4;
  SCVT ACAscale = 1e-4;
  SCVT ACAeta = 1.2;
  LO ACAgroupSize = 100;
  LO dummy = 0;
  LO maxElems = 10;

  Tree< BECluster< LO, SC > *, LO > tree, tree2;
  ProgressMonitor::init( "Nested dissection" );
  mesh.nestedDissection( tree, maxElems );
  mesh.nestedDissection( tree2, maxElems );
  ProgressMonitor::step( );

  FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree2, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );
  bespaceK.setEpsilonACA( ACAeps );
  bespaceK.setScaleACA( ACAscale );
  ACAMatrix< LO, SC > * K = new ACAMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, quadOrder, kappa,
      quadType, quadDisj );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( *K );
  ProgressMonitor::step( );

  std::cout << K->getCompressionRatio( ) << std::endl;

  IdentityOperator< LO, SC > id( &bespaceK );

  SC iUnit( 0.0, 1.0 );
  //SCVT mult = 2.0 / std::sqrt( 3.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nNodes );
  Vector< LO, SC > sol( nNodes );
  Vector< LO, SC > aux( nNodes );
  Vector< LO, SC > neu( nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    //dir.set( i, std::exp( iUnit * mult * ( x1[0]+x1[1]+x1[2] ) ) );
    dir.set( i, four * x1[0] * std::exp(
        two * std::sqrt( three ) * x1[1] + four * iUnit * x1[2] ) );
  }
  //dir.print();

  for ( LO i = 0; i < nElems; i++ ) {

    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    //neu.set( i, iUnit * mult * std::exp( iUnit * mult * ( y[0]+y[1]+y[2] ) ) * ( y[0]+y[1]+y[2] ) );
    //neu.set( i, std::exp(2.0*sqrt(3.0)*y[1]+4.0*iUnit*y[2]) * 4.0 * y[0] * ( 1.0 + 2.0*sqrt(3.0)*y[1] + 4.0*iUnit*y[2] ) );
    mesh.getNormal( i, n );
    neu.set( i, std::exp( two * std::sqrt( three ) * y[1] +
        four * iUnit * y[2] ) * four * ( n[0] + two * std::sqrt( three ) *
        n[1] * y[0] + four * iUnit * y[0] * n[2] ) );
  }
  //neu.print();

  ProgressMonitor::init( "Setting up RHS" );
  id.apply( neu, rhs, true, 0.5, 0.0 );
  //rhs.print();
  K->apply( neu, aux, true );
  //aux.print();
  rhs.add( aux, -1.0 );
  //rhs.print();
  ProgressMonitor::step( );

  //delete K;

  FastBESpace< LO, SC > bespaceD( &mesh, p1, p1, &tree, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );
  bespaceD.setEpsilonACA( ACAeps );
  bespaceD.setScaleACA( ACAscale );
  ACAMatrix< LO, SC > D( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespaceD, quadOrder,
      kappa, quadType, quadDisj );
  ProgressMonitor::init( "Assembling D" );
  formD.assemble( D );
  ProgressMonitor::step( );

  std::cout << D.getCompressionRatio( ) << std::endl;

  ProgressMonitor::init( "Solving the system" );
  SCVT prec = 1e-10;
  LO maxIter = 3000;
  D.GMRESSolve( rhs, sol, prec, maxIter, maxIter );
  ProgressMonitor::step( );

  std::cout << "L2 relative error: " << mesh.l2RelativeErrorLin( sol ) << "." <<
      std::endl;

  string nodeNames[] = { "Dirichlet_anal", "Dirichlet_comp", "Dirichlet_err" };
  string elemNames[] = { "Neumann_anal" };
  Vector< LO, SC > err( sol );
  sol.add( dir, err, -1.0 );
  Vector< LO, SC >* nodalData[] = { &dir, &sol, &err };
  Vector< LO, SC >* elemData[] = { &neu };
  mesh.printParaviewVtu( "output/output.vtu", 3, nodeNames, nodalData, 1,
      elemNames, elemData );

  ProgressMonitor::init( "Representation formula" );
  SCVT * evalPoint = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.285910;
    evalPoint[3 * i + 1] = 0.476517;
    evalPoint[3 * i + 2] = 0.667123;
  }
  SC exact = four * evalPoint[0] * std::exp( two * std::sqrt( three ) *
      evalPoint[1] + four * iUnit * evalPoint[2] );
  Vector< LO, SC > res( nPoints );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &sol, &neu, kappa,
      order );
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