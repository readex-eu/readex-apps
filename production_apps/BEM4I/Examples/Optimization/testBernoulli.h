#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../SurfaceMesh3D.h"
#include "../../BernoulliSubproblem.h"
#include "../../MultiresolutionOptimizer.h"
#include "../../Vector.h"

using namespace std;
using namespace bem4i;

void testBernoulli( );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  testBernoulli( );

  return 0;
}

void testBernoulli( ) {

  typedef double SC;
  typedef int LO;

  basisType basis = p0;

  SurfaceMesh3D< LO, SC > fixedMesh( "input/lshape_28.txt" );
  //SurfaceMesh3D< LO, SC > fixedMesh( "input/icosahedron.txt" );
  fixedMesh.flipNormals( );
  fixedMesh.refine( 2 );
  //fixedMesh.mapToUnitBall( );
  //fixedMesh.scale( 2.0 );

  SurfaceMesh3D< LO, SC > freeMesh( "input/icosahedron.txt" );
  //SurfaceMesh3D< LO, SC > freeMesh( "input/cube_12.txt" );
  //freeMesh.scale( 1.5, 2.0, 2.0 );
  //freeMesh.refine( 2 );
  freeMesh.scale( 2.5 );
  //SurfaceMesh3D< LO, SC > freeMesh( "input/lshape_envelope_2_0.txt" );

  SC Q = -1.5;
  SC dirFree = 0.0;
  LO sizeFixed = 0;
  if ( basis == p0 ) {
    sizeFixed = fixedMesh.getNElements( );
  } else if ( basis == p1 ) {
    sizeFixed = fixedMesh.getNNodes( );
  }
  Vector< LO, SC > dirFixed( sizeFixed );
  dirFixed.setAll( 1.0 );
  BernoulliSubproblem< LO, SC > bvp( dirFree, dirFixed, Q, basis, basis );

  bvp.setCostType( VNorm );

  LO analLevel = 2;
  string tagFile = "";
  //string tagFile = "input/lshape_envelope.tag";
  string bndFile = "";
  //string bndFile = "input/lshape_bounds2.bnd";

  MultiresolutionOptimizer< LO, SC > optimizer( freeForm );

  optimizer.setFixedPart( fixedMesh );
  optimizer.setFreePart( freeMesh, analLevel, tagFile, bndFile );
  optimizer.setMaxOptimLevel( analLevel );
  optimizer.setBVP( bvp );
  optimizer.setOutputPath( "output/optim/bernoulli_nlopt" );
  optimizer.setNormalizeCost( true );
  optimizer.setNormalizeGrad( true );
  optimizer.setNormalizeOnEachLevel( true );
  optimizer.setNloptFRel( 1e-2 );
  optimizer.setNloptJumpTol( 1 );

  optimizer.optimize( );
}

//void testBernoulli( ) {
//
//  typedef double SC;
//  typedef int LO;
///*
//  OpenMeshWrapper< LO, SC > bb( "input/icosahedron.txt" );
//  bb.printParaviewVtu( "output/ico20.vtu" );
//  bb.subdivideLoop( 2 );
//  bb.printParaviewVtu( "output/ico320.vtu" );
//  bb.subdivideLoop( 3 );
//  bb.printParaviewVtu( "output/ico20480.vtu" );
//
//  bb.load( "input/icosahedron.txt" );
//  SC xx[ 3 ], nn[ 3 ];
//  bb.getNode( 0, xx );
//  bb.getNodalNormal( 0, nn );
//  xx[ 0 ] += nn[ 0 ];
//  xx[ 1 ] += nn[ 1 ];
//  xx[ 2 ] += nn[ 2 ];
//  bb.setNode( 0, xx );
//  bb.printParaviewVtu( "output/ico20pert.vtu" );
//  bb.subdivideLoop( 2 );
//  bb.printParaviewVtu( "output/ico320pert.vtu" );
//  bb.subdivideLoop( 3 );
//  bb.printParaviewVtu( "output/ico20480pert.vtu" );
//
//  return;
//*/
//  SurfaceMesh3D< LO, SC > fixedMesh( "input/lshape_28.txt" );
//  //SurfaceMesh3D< LO, SC > fixedMesh( "input/icosahedron.txt" );
//  fixedMesh.flipNormals( );
//  fixedMesh.refine( 3 );
//  //fixedMesh.mapToUnitBall( );
//  //fixedMesh.scale( 2.0 );
//
//  SurfaceMesh3D< LO, SC > freeMesh( "input/icosahedron.txt" );
//  //SurfaceMesh3D< LO, SC > freeMesh( "input/cube_12.txt" );
//  //freeMesh.scale( 1.5, 2.0, 2.0 );
//  //freeMesh.refine( 2 );
//  freeMesh.scale( 2.5 );
//  //SurfaceMesh3D< LO, SC > freeMesh( "input/lshape_envelope_2_0.txt" );
//
//  SC Q = -1.5;
//  SC dirFree = 0.0;
//  Vector< LO, SC > dirFixed( fixedMesh.getNNodes( ) );
//  dirFixed.setAll( 1.0 );
//  BernoulliSubproblem< LO, SC > bvp( dirFree, dirFixed, Q );
//
//  LO analLevel = 3;
//  string tagFile = "";
//  //string tagFile = "input/lshape_envelope.tag";
//  string bndFile = "";
//  //string bndFile = "input/lshape_bounds2.bnd";
//
//  MultiresolutionOptimizer< LO, SC > optimizer( fixedRef );
//
//  optimizer.setFixedPart( fixedMesh );
//  optimizer.setFreePart( freeMesh, analLevel, tagFile, bndFile );
//  optimizer.setMaxOptimLevel( analLevel );
//  optimizer.setBVP( bvp );
//  optimizer.setOutputPath( "output/optim/bernoulli_nlopt" );
//  optimizer.setNormalizeCost( true );
//  optimizer.setNormalizeGrad( true );
//  optimizer.setNormalizeOnEachLevel( true );
//#ifdef IPOPT
//  optimizer.setIpoptTol( 1e-2 );
//  optimizer.setIpoptAcceptableTol( 1e-2 );
//  optimizer.setIpoptAcceptableObjChangeTol( 1e-2 );
//  optimizer.setIpoptAcceptableIter( 2 );
//  optimizer.setIpoptMaxIter( 20 );
//#else
//  optimizer.setNloptFRel( 1e-3 );
//  optimizer.setNloptJumpTol( 1 );
//#endif
//  optimizer.optimize( );
//}