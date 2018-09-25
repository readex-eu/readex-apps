#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"

using namespace std;
using namespace bem4i;

void testHelmholtzDirichletACA(
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
  testHelmholtzDirichletACA( filename, 3, 0, true, 1, 3, 1 );

  return 0;
}

void testHelmholtzDirichletACA(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int nPoints ) {

  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  std::cout.precision( 8 );
  std::cout.setf( std::ios::scientific );

  ProgressMonitor::init( "Loading mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  //mesh.refine(4);
  if ( mapToUnitBall ) mesh.mapToUnitBall( );
  ProgressMonitor::step( );
  mesh.printInfo( );

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

  int quadDisj [] = { 4, 4 };

  SCVT ACAeps = 1e-6;
  SCVT ACAscale = 1e-3;
  SCVT ACAeta = 1.0;
  //  bool ACAgroup = true;
  LO ACAgroupSize = 500;
  LO dummy = 0;
  LO maxElems = 50;

  ProgressMonitor::init( "Nested dissection" );
  Tree<BECluster<LO, SC>*, LO> tree, tree2;
  mesh.nestedDissection( tree, maxElems );
  mesh.nestedDissection( tree2, maxElems );
  ProgressMonitor::step( );

  SC kappa = 2.0;
  ///*

  FastBESpace< LO, SC > bespaceV( &mesh, p0, p0, &tree, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );
  bespaceV.setEpsilonACA( ACAeps );
  bespaceV.setScaleACA( ACAscale );

  ACAMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespaceV, quadOrder, kappa,
      quadType, quadDisj );
  //ProgressMonitor::init( "Assembling V" );
  formV.assemble( V );
  //ProgressMonitor::step( );
  std::cout << "Compression rate of V: " <<
      V.getCompressionRatio( ) << std::endl;

  FastBESpace< LO, SC > bespaceK( &mesh, p1, p0, &tree2, ACAeta, dummy, dummy,
      dummy, ACAgroupSize );
  bespaceK.setEpsilonACA( ACAeps );
  bespaceK.setScaleACA( ACAscale );
  ACAMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespaceK, quadOrder, kappa,
      quadType, quadDisj );
  //ProgressMonitor::init( "Assembling K" );
  formK.assemble( K );
  //ProgressMonitor::step( );
  std::cout << "Compression rate of K: " <<
      K.getCompressionRatio( ) << std::endl;

  IdentityOperator< LO, SC > id( &bespaceK );
  //*/
  SC iUnit( 0.0, 1.0 );
  //SCVT mult = 2.0 / std::sqrt( 3.0 );
  LO nNodes = mesh.getNNodes( );
  LO nElems = mesh.getNElements( );
  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );
  Vector< LO, SC > neu( nElems );
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];
  SCVT four = 4.0;
  SCVT three = 3.0;
  SCVT two = 2.0;

  //SC y1[] = { 0.15, 0.24, 0.7 };
  //SC y2[] = { -0.125, 0.23, 0.74 };
  //SC norm1, norm2;
  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    //dir.set( i, std::exp( iUnit * mult * ( x1[0]+x1[1]+x1[2] ) ) );
    dir.set( i, four * x1[0] * std::exp( two * std::sqrt( three ) * x1[1] +
        four * iUnit * x1[2] ) );
    //norm1 = std::sqrt( (x1[0]-y1[0])*(x1[0]-y1[0]) + (x1[1]-y1[1])*(x1[1]-y1[1]) + (x1[2]-y1[2])*(x1[2]-y1[2]) );
    //norm2 = std::sqrt( (x1[0]-y2[0])*(x1[0]-y2[0]) + (x1[1]-y2[1])*(x1[1]-y2[1]) + (x1[2]-y2[2])*(x1[2]-y2[2]) );
    //dir.set( i, 1.0 / (4.0*M_PI) * ( std::exp( iUnit * 2.0 * norm1 ) / norm1 + std::exp( iUnit * 2.0 * norm2 ) / norm2 ) );
  }
  //dir.print();
  //string nodeNames[] = { "Dirichlet_anal" };
  //Vector< LO, SC >* nodalData[] = { &dir };
  //mesh.printParaviewVtu( "output/output.vtu", 1, nodeNames, nodalData, 0, nullptr, nullptr );
  //return;
  ///*
  //LO ind[3];
  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    //neu.set( i, iUnit * mult * std::exp( iUnit * mult * ( y[0]+y[1]+y[2] ) ) * ( y[0]+y[1]+y[2] ) );
    //neu.set( i, std::exp(2.0*sqrt(3.0)*y[1]+4.0*iUnit*y[2]) * 4.0 * y[0] * ( 1.0 + 2.0*sqrt(3.0)*y[1] + 4.0*iUnit*y[2] ) );
    mesh.getNormal( i, n );
    neu.set( i,
        std::exp( two * std::sqrt( three ) * y[1] + four * iUnit * y[2] ) *
        four * ( n[0] + two * std::sqrt( three ) * n[1] * y[0] +
        four * iUnit * y[0] * n[2] ) );

    //mesh.getElement( i, ind );
    //neu.set( i, 1.0 / 3.0 * iUnit * mult * ( dir.get(ind[0])*(x1[0]+x1[1]+x1[2]) + dir.get(ind[1])*(x2[0]+x2[1]+x2[2]) + dir.get(ind[2])*(x2[0]+x2[1]+x2[2]) ) );
  }
  //neu.print();
  ///*
  ProgressMonitor::init( "Setting up RHS" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  //aux.print();
  K.apply( dir, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> sol( rhs.getLength( ) );
  V.GMRESSolve( rhs, sol, 1e-10, 1000, 100 );
  ProgressMonitor::step( );

  //std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( neu, rhs ) << "." << std::endl;
  std::cout << "L2 relative error: " <<
      mesh.l2RelativeErrorConst( sol ) << "." << std::endl;

  string nodeNames[] = { "Dirichlet_anal" };
  string elemNames[] = { "Neumann_anal", "Neumann_comp", "Neumann_err" };
  Vector< LO, SC >* nodalData[] = { &dir };
  Vector< LO, SC > err( sol );
  sol.add( neu, err, -1.0 );
  Vector< LO, SC >* elemData[] = { &neu, &sol, &err };
  mesh.printParaviewVtu( "output/output.vtu", 1, nodeNames, nodalData, 3,
      elemNames, elemData );
  //*/

  ProgressMonitor::init( "Representation formula" );
  SCVT * evalPoint = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.285910;
    evalPoint[3 * i + 1] = 0.476517;
    evalPoint[3 * i + 2] = 0.667123;
  }
  SC exact = four * evalPoint[0] * std::exp(
      two * std::sqrt( three ) * evalPoint[1] + four * iUnit * evalPoint[2] );
  Vector< LO, SC > res( nPoints );

  BESpace< LO, SC > bespaceRep( &mesh, p1, p0 );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceRep, &dir, &sol,
      kappa, 5 );
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
  /*
    const int aCount = 30;
    const int bCount = aCount;
    const int cCount = aCount;
    SCVT aSize = 1.98;
    SCVT bSize = aSize;
    SCVT cSize = aSize;
    SCVT aStep = aSize / (SCVT) aCount;
    SCVT bStep = bSize / (SCVT) bCount;
    SCVT cStep = cSize / (SCVT) cCount;
    const int nPoints = ( aCount + 1 ) * ( bCount + 1 ) * ( cCount + 1 );
    SCVT pointCloud[ 3 * nPoints ];
    int counter = 0;
    for ( int a = 0; a <= aCount; a++ ) {
      for ( int b = 0; b <= bCount; b++ ) {
        for ( int c = 0; c <= cCount; c++ ) {
          pointCloud[ counter++ ] = a * aStep - aSize / 2.0;
          pointCloud[ counter++ ] = b * bStep - bSize / 2.0;
          pointCloud[ counter++ ] = c * cStep - cSize / 2.0;
        }
      }
    }

    Vector< LO, SC > pointSolution( nPoints );
    std::cout << "Evaluating representation formula ... ";
    std::cout.flush( );
    RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceK, &dir, &rhs, kappa, order );
    gettimeofday( &start, nullptr );
    formula.evaluate( pointCloud, nPoints, true, pointSolution );
    gettimeofday( &stop, nullptr );
    Vector< LO, SCVT > pointSolutionReal( nPoints );
    Vector< LO, SCVT > pointSolutionImag( nPoints );
    for ( LO i = 0; i < pointSolution.getLength( ); i++ ) {
      pointSolutionReal.set( i, pointSolution.get( i ).real( ) );
      pointSolutionImag.set( i, pointSolution.get( i ).imag( ) );
    }
    std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
    string nodeNamesSol[] = { "solution_real", "solution_imag" };
    Vector< LO, SCVT >* nodalDataSol[] = { &pointSolutionReal, &pointSolutionImag };
    formula.printParaviewVtu( "output/solution_point_cloud.vtp", pointCloud, nPoints, 2, nodeNamesSol, nodalDataSol );
   */
}