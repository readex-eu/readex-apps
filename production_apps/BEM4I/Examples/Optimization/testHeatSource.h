#include "../../Settings.h"

#include <iostream>
#include <vector>

#include "../auxiliary.h"
#include "../../SurfaceMesh3D.h"
#include "../../HeatSourceSubproblem.h"
#include "../../MultiresolutionOptimizer.h"
#include "../../Vector.h"

using namespace std;
using namespace bem4i;

void testHeatSource( );

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

  testHeatSource( );

  return 0;
}

void testHeatSource( ) {

  typedef double SC;
  typedef int LO;

  SurfaceMesh3D< LO, SC > sensor( "input/icosahedron.txt" );
  sensor.refine( 3 );
  sensor.mapToUnitBall( );
  //sensor.scale( 1.1 );

  SurfaceMesh3D< LO, SC > source( "input/icosahedron.txt" );
  source.flipNormals( );

  /*
  SurfaceMesh3D< LO, SC > source( "input/duck2.txt" );
  source.printParaviewVtu( "output/duck2.vtu" );
  return;
  //source.refine( 3 );
  //source.mapToUnitBall( );
  //source.scale( 0.8, 1.0, 1.2 );
  Vector< LO, SC > dummy( sensor.getNElements( ), true );
  HeatSourceSubproblem< LO, SC > forwardProblem( dummy );
  forwardProblem.setProblemData( &source, &sensor );
  forwardProblem.solve( );
  forwardProblem.printVtu( "output/forward_heat_ellipsoid_320.vtu" );
  //source.printParaviewVtu( "output/ellipsoid.vtu" );
  return;
   */

  Vector< LO, SC > targetNeumann( sensor.getNNodes( ) );
  readData< LO, SC >( "input/target_sphere_heat_1280.txt", targetNeumann,
      0, sensor.getNElements( ) - 1 );
  HeatSourceSubproblem< LO, SC > bvp( targetNeumann );
  bvp.setRegularizationParameter( 1e-6 );
  //bvp.setCostMultiplier( 10.0 );

  LO analLevel = 3;
  string tagFile = "";
  string bndFile = "input/icosa_alpha.bnd";
  MultiresolutionOptimizer< LO, SC > optimizer( fixedRef );
  optimizer.setFixedPart( sensor );
  optimizer.setFreePart( source, analLevel, tagFile, bndFile );
  optimizer.setMaxOptimLevel( analLevel );
  optimizer.setBVP( bvp );
  optimizer.setOutputPath( "output/optim/test_ipopt_3" );
#ifdef IPOPT
  optimizer.setIpoptTol( 1e-3 );
  optimizer.setIpoptAcceptableTol( 1e-2 );
  optimizer.setIpoptAcceptableObjChangeTol( 1e-2 );
  optimizer.setIpoptAcceptableIter( 2 );
  optimizer.setIpoptMaxIter( 20 );
#else
  optimizer.setNloptFRel( 1e-2 );
  optimizer.setNloptJumpTol( 2 );
#endif
  optimizer.setNormalizeCost( false );
  optimizer.setNormalizeGrad( false );
  optimizer.setNormalizeOnEachLevel( false );

  optimizer.optimize( );
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