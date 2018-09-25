#include "../../Settings.h"

#include <iostream>
#include <vector>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../SurfaceMesh3D.h"
#include "../../FullMatrix.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"

using namespace std;
using namespace bem4i;

void testRohan( );

template< class LO, class SC >
bool readData(
    string const & filename,
    Vector< LO, SC > & data,
    LO begin = 0,
    LO end = 0
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

  testRohan( );

  return 0;
}

void testRohan( ) {

  typedef double SCVT;
  typedef std::complex< SCVT > SC;
  typedef int LO;
  /*
  std::string gg = "input/rohan_28k.txt";
  SurfaceMesh3D< LO, SCVT > ggg( gg );
  //ggg.refine( 3 );
  ggg.printInfo( );

  string bla = "input/scatter_grid_y.txt";
  SurfaceMesh3D< LO, SCVT > sit;
  sit.load( bla.c_str( ) );
  sit.scale( 0.5, 1.0, 0.5 );
  sit.move( 0.5, 0.0, 0.5 );
  sit.scale( 0.7, 1.0, 0.2 );
  sit.move( -0.2, 0.0, -0.1 );
  sit.refine( 7 );

  for ( LO i = 1; i < 20; ++i ) {
    sit.move( 0.0, 0.01, 0.0 );
    SurfaceMesh3D< LO, SCVT > mi;
    ProgressMonitor::init( "Trimming" );
    ggg.trimEvaluationGrid( sit, mi, false );
    ProgressMonitor::step( );
    mi.printInfo( );
    std::stringstream file;
    file << "output/output_i_" << i << ".vtu";
    mi.printParaviewVtu( file.str( ) );
    std::stringstream file2;
    file2 << "output/output_i_" << i << ".txt";
    mi.print( file2.str( ) );
  }
  return;
   */

  std::string outputPath = "output/rohan/ac3d_wg01_surf";
  const std::string createDir = "mkdir -p ./" + outputPath;
  if ( system( createDir.c_str( ) ) ) exit( -1 );

  std::string inputPathPrefix = "input/rohan/ac3d_wg01_surf";
  std::stringstream inputFilePath;
  inputFilePath << inputPathPrefix << ".txt";

  SurfaceMesh3D< LO, SC > mesh( inputFilePath.str( ) );
  mesh.printInfo( );

  std::vector< LO > inE;
  std::vector< LO > outE;

  basisType type = p1;

  LO nElems = mesh.getNElements( );

  /*
  SCVT x1[ 3 ], x2[ 3 ], x3[ 3 ];
  LO nNodes = mesh.getNNodes( );
  for ( LO i = 0; i < nElems; ++i ) {
    mesh.getNodes( i, x1, x2, x3 );
    if ( ( x1[ 0 ] == -1 ) && ( x2[ 0 ] == -1 ) && ( x3[ 0 ] == -1 ) ) {
      inE.push_back( i );
    } else if ( ( x1[ 0 ] == 1 ) && ( x2[ 0 ] == 1 ) && ( x3[ 0 ] == 1 ) ) {
      outE.push_back( i );
    }
  }
   */
  inputFilePath.clear( );
  inputFilePath.str( "" );
  inputFilePath << inputPathPrefix << "_in.txt";
  readData< LO >( inputFilePath.str( ), inE );
  inputFilePath.clear( );
  inputFilePath.str( "" );
  inputFilePath << inputPathPrefix << "_out.txt";
  readData< LO >( inputFilePath.str( ), outE );

  int order = 3;
  int quadOrder[ 4 ] = { order, order, order, order };
  quadratureType quadType = SauterSchwab;

  int quadDisjoint[] = { 4, 4 };

  SC kappa = 20.0;

  BESpace< LO, SC > bespace( &mesh, type, type );

  FullMatrix< LO, SC > * V = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz1Layer< LO, SC > formV( &bespace, quadOrder, kappa,
      quadType, quadDisjoint );
  ProgressMonitor::init( "Assembling V" );
  formV.assemble( *V );
  ProgressMonitor::step( );

  FullMatrix< LO, SC > * K = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtz2Layer< LO, SC > formK( &bespace, quadOrder, kappa,
      quadType, quadDisjoint );
  ProgressMonitor::init( "Assembling K" );
  formK.assemble( *K );
  ProgressMonitor::step( );

  FullMatrix< LO, SC > * D = new FullMatrix< LO, SC >( 0, 0 );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formD( &bespace, quadOrder,
      kappa, quadType, quadDisjoint );
  ProgressMonitor::init( "Assembling D" );
  if ( type == p0 ) {
    formD.assembleP0P0( *D );
  } else if ( type == p1 ) {
    formD.assemble( *D );
  }
  ProgressMonitor::step( );

  SparseMatrix< LO, SC > M;
  IdentityOperator< LO, SC > id( &bespace );
  id.assemble( M );

  HelmholtzNeumannRobinOperator< LO, SC > op( &mesh, inE, outE, kappa, V, K, D,
      &M, type );

  Vector< LO, SC > rhs;
  op.getRHS( rhs );

  ProgressMonitor::init( "Solving the system" );
  //  SCVT prec = 1e-8;
  //  LO maxIt = 2000;
  //  LO restarts = maxIt;
  //  LO size = 0;
  //  if ( type == p0 ) {
  //    size = nElems;
  //  } else if ( type == p1 ) {
  //    size = nNodes;
  //  } else {
  //    return;
  //  }
  //  Vector< LO, SC > u( size );
  //  op.GMRESSolve( rhs, u, prec, maxIt, restarts );
  Vector< LO, SC > u( rhs );
  op.LUSolve( u );
  ProgressMonitor::step( );

  ProgressMonitor::init( "Solving Dirichlet problem" );
  Vector< LO, SC > dudn;
  op.solveDirichletProblem( u, dudn );
  ProgressMonitor::step( );

  Vector< LO, SC > bc( nElems, true );
  for ( auto it = inE.begin( ); it != inE.end( ); ++it ) {
    bc.set( *it, 1.0 );
  }
  for ( auto it = outE.begin( ); it != outE.end( ); ++it ) {
    bc.set( *it, 2.0 );
  }

  std::stringstream outputFilePath;
  outputFilePath << outputPath << "/volume.vtu";

  if ( type == p1 ) {
    string nodeNames[] = { "dirichlet", "neumann" };
    string elemNames[] = { "bc" };
    Vector< LO, SC >* nodalData[] = { &u, &dudn };
    Vector< LO, SC >* elemData[] = { &bc };
    mesh.printParaviewVtu( outputFilePath.str( ), 2, nodeNames,
        nodalData, 1, elemNames, elemData );
  } else if ( type == p0 ) {
    string elemNames[] = { "dirichlet", "neumann", "bc" };
    Vector< LO, SC >* elemData[] = { &u, &dudn, &bc };
    mesh.printParaviewVtu( outputFilePath.str( ), 0, nullptr, nullptr, 3,
        elemNames, elemData );
  }

  delete V;
  delete K;
  delete D;

  /*
  for ( LO i = 1; i < 20; ++i ) {

    inputFilePath.clear( );
    inputFilePath.str( "" );
    inputFilePath << inputPathPrefix << "_grid_" << i << ".txt";
    SurfaceMesh3D< LO, SC > grid( inputFilePath.str( ) );

    ProgressMonitor::init( "Evaluating inside" );
    Vector< LO, SC > res( grid.getNNodes( ) );
    op.evalInside( grid, u, dudn, res );
    ProgressMonitor::step( );

    string nodeNames[] = { "solution" };
    Vector< LO, SC >* nodalData[] = { &res };
    outputFilePath.clear( );
    outputFilePath.str( "" );
    outputFilePath << outputPath << "/grid_i_" << i << ".vtu";
    grid.printParaviewVtu( outputFilePath.str( ), 1, nodeNames, nodalData, 0,
        nullptr, nullptr );
  }
   */
}

template< class LO, class SC >
bool readData(
    string const &filename,
    Vector< LO, SC > &data,
    LO begin,
    LO end
    ) {

  std::cout << "Reading file " << filename << " ... ";
  std::ifstream file( filename.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "not found." << std::endl;
    return false;
  }

  LO length;
  //int idx;
  SC value;

  file >> length;
  if ( end == 0 ) end = length - 1;

  if ( end + 1 > length ) {
    data.resize( end + 1, true );
  } else {
    data.resize( length, true );
  }

  for ( LO i = begin; i <= end; ++i ) {
    //file >> idx;
    file >> value;
    //data.set( idx, value );
    data.set( i, value );
  }

  file.close( );
  std::cout << "done." << std::endl;

  return true;
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