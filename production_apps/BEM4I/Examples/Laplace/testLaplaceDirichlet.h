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
#include "../../RepresentationFormulaLaplace.h"
#include "../../IdentityOperator.h"
#include "../../Laplace1LayerP0P0MultilvlPrecond.h"

#include "readex.h"

using namespace std;
using namespace bem4i;

void testLaplaceDirichlet(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int nPoints
    );

template< class LO >
bool readData(
    string const & filename,
    std::vector< LO > & data
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/icosahedron.txt";
  testLaplaceDirichlet( filename, 3, 1, true, 0, 5, 1 );

  return 0;
}

void testLaplaceDirichlet(
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

  READEX_INIT()
  READEX_REGION_DEFINE(main)
  READEX_PHASE_START(main,"main",SCOREP_USER_REGION_TYPE_COMMON);

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

  READEX_REGION_DEFINE(assemble_v)
  READEX_REGION_START(assemble_v,"assemble_v",SCOREP_USER_REGION_TYPE_COMMON);
  
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
  
  READEX_REGION_STOP(assemble_v) 

  READEX_REGION_DEFINE(assemble_k)
  READEX_REGION_START(assemble_k,"assemble_k",SCOREP_USER_REGION_TYPE_COMMON);
  
  BESpace< LO, SC > bespaceK( &mesh, p1, p0 );
  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormLaplace2Layer< LO, SC > formK( &bespaceK, quadOrder, quadType,
      quadFar );
  ProgressMonitor::init( "K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print( );
  //return;
  
  READEX_REGION_STOP(assemble_k) 

  /*
  ProgressMonitor::init( "V, K" );
  formV.assembleWith2Layer( V, K, bespaceK );
  ProgressMonitor::step();
   */
  IdentityOperator< LO, SC > id( &bespaceK );
  //*/

  Vector< LO, SC > dir( nNodes );
  Vector< LO, SC > rhs( nElems );
  Vector< LO, SC > aux( nElems );
  Vector< LO, SC > neu( nElems );
  SC x1[ 3 ], x2[ 3 ], x3[ 3 ], y[ 3 ], n[ 3 ];

  /*
  for ( LO i = 0; i < nNodes; i++ ) {
    mesh.getNode( i, x1 );
    dir.set( i, x1[ 0 ] + x1[ 1 ] + x1[ 2 ] );
  }
   */

  ///*
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );
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
  //*/
  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    mesh.getNormal( i, n );
    ///*
    neu.set( i, std::exp( 2.0 * M_PI * y[1] ) * ( n[0] * std::cos( 2.0 * M_PI * y[2] )
        + 2.0 * M_PI * ( 1.0 + y[0] ) * n[1] * std::cos( 2.0 * M_PI * y[2] )
        - 2.0 * M_PI * ( 1.0 + y[0] ) * n[2] * std::sin( 2.0 * M_PI * y[2] ) ) );
    //*/
    //neu.set( i, n[ 0 ] + n[ 1 ] + n[ 2 ] );
  }

  ///*
  ProgressMonitor::init( "Setting up the rhs" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  //aux.print();
  K.apply( dir, rhs );
  // EXTERIOR ONLY
  //rhs.scale( -1.0 );
  rhs.add( aux );
  //rhs.print();
  //return;
  ProgressMonitor::step( );
  //*/

  ProgressMonitor::init( "Setting up the preconditioner" );
  /*
  / Create the class Laplace1LayerP0P0MultilvlPrecond 
  /   - containing tree<BECluster> and mesh
  / Init the preconditioner and apply within the CG method
   */
  //Laplace1LayerP0P0MultilvlPrecond<LO, SC> precond( &mesh, 15 );
  ProgressMonitor::step( );

  //V.print( );

  READEX_REGION_DEFINE(cg_solve)
  READEX_REGION_START(cg_solve,"cg_solve",SCOREP_USER_REGION_TYPE_COMMON);
  
  ProgressMonitor::init( "Solving the system" );
  Vector<LO, SC> rhs2( rhs );
  rhs.setAll( 0.0 );
  V.CGSolve( rhs2, rhs, CGeps, 2000 );
  //V.CGSolve( rhs2, rhs, CGeps, 2000, &precond, "V solve" );
  ProgressMonitor::step( );

  READEX_REGION_STOP(cg_solve)

  //std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( neu, rhs )
  //  << "." << std::endl;
  std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( rhs )
      << "." << std::endl;
  /*
  Vector< LO, SC > dal;
  readData< LO, SC >( "input/dalibor2.txt", dal );
  std::cout << "L2 relative error: " << mesh.l2RelativeErrorConst( dal )
      << "." << std::endl;
   */
  
  READEX_REGION_DEFINE(print_vtu)
  READEX_REGION_START(print_vtu,"print_vtu",SCOREP_USER_REGION_TYPE_COMMON);
  
  string nodeNames[] = { "Dirichlet_anal" };
  string elemNames[] = { "Neumann_anal", "Neumann_comp", "Neumann_err" };
  Vector< LO, SC >* nodalData[] = { &dir };
  Vector< LO, SC > err( rhs.getLength( ) );
  //Vector< LO, SC > errRel( rhs );
  rhs.add( neu, err, -1.0 );
  //errRel.errorRelative( neu );
  Vector< LO, SC >* elemData[] = { &neu, &rhs, &err };
  mesh.printParaviewVtu( "output/output.vtu", 1, nodeNames, nodalData, 3,
      elemNames, elemData );
  //*/
  
  READEX_REGION_STOP(print_vtu)

  SC * evalPoint = new SC[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {
    evalPoint[3 * i] = 0.250685;
    evalPoint[3 * i + 1] = 0.417808;
    evalPoint[3 * i + 2] = 0.584932;
  }
  ///*
  SC exact = ( 1.0 + evalPoint[ 0 ] )
      * std::exp( 2.0 * M_PI * evalPoint[ 1 ] )
      * std::cos( 2.0 * M_PI * evalPoint[ 2 ] );
  //*/
  //SC exact = evalPoint[ 0 ] + evalPoint[ 1 ] + evalPoint[ 2 ];
  Vector< LO, SC > res( nPoints );

  BESpace< LO, SC > bespaceRep( &mesh, p1, p0 );
  RepresentationFormulaLaplace<LO, SC> formula( &bespaceRep, &dir, &rhs );
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

  READEX_PHASE_STOP(main)
  READEX_CLOSE()
}

template< class LO >
bool readData(
    string const &filename,
    std::vector< LO > & data
    ) {

  std::cout << "Reading file " << filename << " ... ";
  std::ifstream file( filename.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "not found." << std::endl;
    return false;
  }

  LO value;
  LO length;

  file >> length;
  data.clear( );
  data.reserve( length );

  for ( LO i = 0; i < length; ++i ) {
    file >> value;
    data.push_back( value );
  }

  file.close( );
  std::cout << "done." << std::endl;

  return true;
}
