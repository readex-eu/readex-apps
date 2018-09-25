#include "../../Settings.h"

#include <iostream>
#include <complex>

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

#include "readex.h"

using namespace std;
using namespace bem4i;

void testHelmholtzDirichlet(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int orderFar,
    int nPoints
    );

int main(
    int argc,
    char** argv
    ) {

  std::string filename = "input/cube_12.txt";
  int refine = 5;
  int quadTypeInt = 1;
  int order = 4;
  int orderFar = 4;

  if( argc == 6 ){
    filename = argv[ 1 ];
    refine = atoi( argv[ 2 ] );
    quadTypeInt = atoi( argv[ 3 ] );
    order = atoi( argv[ 4 ] );
    orderFar = atoi( argv[ 5 ] );
  }

  MPI_Init( &argc, &argv );
  READEX_INIT();
  
  READEX_PHASE_REGION_DEFINE(main_phase);
  READEX_PHASE_START(main_phase,"main_phase",SCOREP_USER_REGION_TYPE_COMMON);
  
  intro();
  
  testHelmholtzDirichlet( filename, refine, 0, false, quadTypeInt, order, orderFar, 1 );
  
  READEX_PHASE_STOP(main_phase);

  READEX_CLOSE();
  MPI_Finalize( );

  return 0;
}

void testHelmholtzDirichlet(
    string const &filename,
    int refine4,
    int refine9,
    bool mapToUnitBall,
    int quadTypeInt,
    int order,
    int orderFar,
    int nPoints ) {

  typedef complex<double> SC;
  typedef long LO;
  typedef typename SC::value_type SCVT;

  READEX_SIGNIFICANT_REGION_DEFINE(helmholtz_dirichlet);
  READEX_REGION_START(helmholtz_dirichlet,"helmholtz_dirichlet",SCOREP_USER_REGION_TYPE_COMMON);

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

  int * quadDisj = nullptr;
  if( orderFar > 0 ){  
    quadDisj = new int[ 2 ];
    quadDisj[ 0 ] = quadDisj[ 1 ] = orderFar;
  }
  
  int quadOrder[ 4 ];
  quadratureType quadType;
  if ( quadTypeInt == 0 ) {
    quadType = Steinbach;
    quadOrder[ 0 ] = quadOrder[ 1 ] = order;
    std::cout << "Using Steinbach quadrature, order near-field " << order;
    if( quadDisj )
      std::cout << ", order far-field " << orderFar;
    std::cout << "." << std::endl;
  } else {
    quadType = SauterSchwab;
    quadOrder[ 0 ] = quadOrder[ 1 ] = quadOrder[ 2 ] = quadOrder[ 3 ] = order;
    std::cout << "Using Sauter-Schwab quadrature, order near-field " << order;
    if( quadDisj )
      std::cout << ", order far-field " << orderFar;
    std::cout << "." << std::endl;
  }

  //SC kappa(2.0, 1.0);
  SC kappa = 2.0;
  SC b = (SCVT) 2.0 * kappa;
  SC a = -std::sqrt( (SCVT) 3.0 ) * kappa;
  ///*
  BESpace< LO, SC > bespace00( &mesh, p0, p0 );
  BESpace< LO, SC > bespace10( &mesh, p1, p0 );
  BESpace< LO, SC > bespace11( &mesh, p1, p1 );

  READEX_SIGNIFICANT_REGION_DEFINE(assemble_v);
  READEX_REGION_START(assemble_v,"assemble_v",SCOREP_USER_REGION_TYPE_COMMON);
  
  FullMatrix< LO, SC > V( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace00, quadOrder, kappa,
      quadType, quadDisj );
  ProgressMonitor::init( "Assembling V" );
  formV.assemble( V );
  ProgressMonitor::step( );
  //V.print( );
  //return;
  
  READEX_REGION_STOP(assemble_v); 

  READEX_SIGNIFICANT_REGION_DEFINE(assemble_k);
  READEX_REGION_START(assemble_k,"assemble_k",SCOREP_USER_REGION_TYPE_COMMON);
  
  FullMatrix< LO, SC > K( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace10, quadOrder, kappa,
      quadType, quadDisj );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( K );
  ProgressMonitor::step( );
  //K.print( );
  //return;
  
  READEX_REGION_STOP(assemble_k); 

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

  for ( LO i = 0; i < nElems; i++ ) {
    mesh.getNodes( i, x1, x2, x3 );
    y[0] = 1.0 / 3.0 * ( x1[0] + x2[0] + x3[0] );
    y[1] = 1.0 / 3.0 * ( x1[1] + x2[1] + x3[1] );
    y[2] = 1.0 / 3.0 * ( x1[2] + x2[2] + x3[2] );
    mesh.getNormal( i, n );
    neu.set( i, b * std::exp( -a * y[ 1 ] + iUnit * b * y[ 2 ] ) * ( n[ 0 ] -
        a * y[ 0 ] * n[ 1 ] + iUnit * b * y[ 0 ] * n[ 2 ] ) );
  }

  IdentityOperator< LO, SC > id( &bespace10 );
  ProgressMonitor::init( "Assembling rhs" );
  id.apply( dir, aux, false, 0.5, 0.0 );
  K.apply( dir, rhs );
  rhs.add( aux );
  ProgressMonitor::step( );

  READEX_SIGNIFICANT_REGION_DEFINE(gmres_solve);
  READEX_REGION_START(gmres_solve,"gmres_solve",SCOREP_USER_REGION_TYPE_COMMON);
  
  ProgressMonitor::init( "Solving the system" );
  Vector< LO, SC > sol( nElems );
  V.GMRESSolve( rhs, sol, 1e-12, 1000, 1000 );
  ProgressMonitor::step( );
  
  READEX_REGION_STOP(gmres_solve);

  std::cout << "L2 relative error: " <<
      mesh.l2RelativeErrorConst( sol, 5, kappa ) << "." << std::endl;

  Vector< LO, SC > err( sol );
  sol.add( neu, err, -1.0 );
  system("mkdir -p /tmp/zapletal");

  READEX_SIGNIFICANT_REGION_DEFINE(print_vtu);
  READEX_REGION_START(print_vtu,"print_vtu",SCOREP_USER_REGION_TYPE_COMMON);
  
  ProgressMonitor::init( "Printing VTU" );
  
  std::string nodeNames[ ] = { "Dirichlet_anal" };
  std::string elemNames[ ] = { "Neumann_anal", "Neumann_comp", "Neumann_err" };
  Vector< LO, SC > * nodalData[ ] = { &dir };
  Vector< LO, SC > * elemData[ ] = { &neu, &sol, &err };
  std::stringstream name;
  name << "/tmp/output";
  name << ".vtu";
  mesh.printParaviewVtu( name.str( ), 1, nodeNames, nodalData, 3,
      elemNames, elemData );

  SCVT c = 340;
  SCVT omega = c * std::real( kappa );
  SCVT T = 2.0 * M_PI / omega;
  int n_steps = 24;

#pragma omp parallel for schedule( static, 1 )
for( int i = 0; i < n_steps; ++i){
  std::string my_nodeNames[ ] = { "Dirichlet_anal" };
  std::string my_elemNames[ ] = { "Neumann_anal", "Neumann_comp" };
  Vector< LO, SC > my_dir( dir );
  Vector< LO, SC > my_neu( neu );
  Vector< LO, SC > my_sol( sol );

  SCVT time = i * T / n_steps;
  SC time_fact = std::exp( -iUnit * omega * time );
  my_dir.scale( time_fact );
  my_neu.scale( time_fact );
  my_sol.scale( time_fact );

  Vector< LO, SC > * my_nodalData[] = { &my_dir };
  Vector< LO, SC >* my_elemData[] = { &my_neu, &my_sol };
  std::stringstream my_name;
  my_name << "/tmp/output_";
  my_name << i;
  my_name << ".vtu";
  mesh.printParaviewVtu( my_name.str( ), 1, my_nodeNames, my_nodalData, 2,
      my_elemNames, my_elemData );
}
  ProgressMonitor::step( );
  
  READEX_REGION_STOP(print_vtu);

  /*
  SCVT * evalPoint = new SCVT[ 3 * nPoints ];
  for ( LO i = 0; i < nPoints; i++ ) {

    evalPoint[3 * i] = 0.285910;
    evalPoint[3 * i + 1] = 0.476517;
    evalPoint[3 * i + 2] = 0.667123;
  }
  SC exact = b * evalPoint[ 0 ] * std::exp( -a * evalPoint[ 1 ] + iUnit * b *
      evalPoint[ 2 ] );

  Vector< LO, SC > res( nPoints );
  BESpace< LO, SC > bespaceRep( &mesh, p1, p0 );
  RepresentationFormulaHelmholtz<LO, SC> formula( &bespaceRep, &dir, &sol,
      kappa, 5 );
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
  */

  if( quadDisj ) delete [ ] quadDisj;

  READEX_REGION_STOP(helmholtz_dirichlet);

}
