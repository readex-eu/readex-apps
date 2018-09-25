#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../SurfaceMesh3D.h"
#include "../../FullMatrix.h"
#include "../../HomogenizationProblem.h"
#include "../../HelmholtzNeumannRobinOperator.h"

using namespace std;
using namespace bem4i;

void testHomogenizationProblem( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  testHomogenizationProblem( );

  return 0;
}

void testHomogenizationProblem( ) {

  typedef double SC;
  typedef int LO;

  FullMatrix< LO, SC > homoMatrix;

  HomogenizationProblem< LO, SC > problem;

  //  std::string meshFile = "input/icosahedron.txt";
  //  SurfaceMesh3D< LO, SC > meshInclusion( meshFile );
  //  meshInclusion.flipNormals( );
  //  meshInclusion.refine( 1 );
  //  meshInclusion.mapToUnitBall( );
  //  meshInclusion.scale( 0.15, 0.25, 0.35 );
  //  meshInclusion.move( 0.5, 0.5, 0.5 );

  //  std::string meshFile = "input/cube_12.txt";
  //  SurfaceMesh3D< LO, SC > meshInclusion( meshFile );
  //  meshInclusion.flipNormals( );
  //  meshInclusion.refine( 3 );
  //  meshInclusion.scale( 0.315 );
  //  meshInclusion.move( 0.5, 0.5, 0.5 );

  //  std::string meshFile = "input/homo/homocube_1.txt";
  //  SurfaceMesh3D< LO, SC > meshInclusion( meshFile );
  //  meshInclusion.flipNormals( );
  //  meshInclusion.refine( 5 );
  //  meshInclusion.scaleAroundCentroid( 0.5 );

  std::string meshFile = "input/fichera.txt";
  SurfaceMesh3D< LO, SC > meshInclusion( meshFile );
  meshInclusion.flipNormals( );
  meshInclusion.refine( 0 );
  meshInclusion.scale( 0.25 );
  SC u[ 3 ] = { 1.0, 0.5, 0.25 };
  SC alpha = M_PI / 3;
  meshInclusion.rotate( u, alpha );
  meshInclusion.move( 0.5, 0.5, 0.5 );

  problem.setMeshInclusion( meshInclusion );
  problem.setNSegsPerEdge( 4 );
  problem.solve( );
  problem.getHomoMatrix( homoMatrix );
  problem.printVtu( "output/homo/test.vtu" );

  homoMatrix.print( );
}

