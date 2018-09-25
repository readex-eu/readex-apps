#include "../../Settings.h"

#include <iostream>

#include "../auxiliary.h"
#include "../../FullMatrix.h"
#include "../../FMMMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorLaplace.h"
#include "../../BEBilinearFormLaplace1Layer.h"
#include "../../BEBilinearFormLaplace2Layer.h"

using namespace std;
using namespace bem4i;

void testFastBEM(
    string const &filename,
    int printPrecision
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = "input/cube_12.txt";
  int printPrecision = 15;
  testFastBEM( filename, printPrecision );

  return 0;
}

void testFastBEM( string const &filename, int printPrecision ) {
  typedef double SC;
  typedef int LO;
  timeval start, stop;
  std::cout.setf( std::ios::showpoint | std::ios::fixed );
  std::cout.precision( printPrecision );

  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );
  mesh.refine( 3 );
  mesh.printInfo( );

  Tree<BECluster<LO, SC>*, LO> tree;
  std::cout << "nested dissection... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  mesh.nestedDissection( tree, 20 );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "creating block cluster tree... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  FastBESpace< LO, SC > bespaceV( &mesh, p0, p0, &tree );
  bespaceV.setEpsilonACA( 1e-4 );


  //FastBESpace< LO, SC > bespaceV( &mesh, p0, p0);
  //bespaceV.setFastBEMParams(1.0, 5, 6, 0);//
  //bespaceV.createBlockClusterTree();

  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;
  BEBilinearFormLaplace1Layer<LO, SC> formV( &bespaceV );
  BEBilinearFormLaplace2Layer<LO, SC> formK( &bespaceV );

  //FMMMatrix<LO, SC> matrixV;
  ACAMatrix<LO, SC> matrixV;


  FMMMatrix<LO, SC> matrixK;
  FullMatrix<LO, SC> fullMatrixV;
  FullMatrix<LO, SC> fullMatrixK;
  gettimeofday( &start, nullptr );
  std::cout << "assembling V... ";
  std::cout.flush( );
  formV.assemble( matrixV );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  //  Vector<LO, SC> test( mesh.getNElements( ) );
  //  test.setAll( 1.0 );
  //  Vector<LO, SC> test2( mesh.getNElements( ) );
  //  test2.setAll( 1.0 );
  //  matrixV.apply( test, test2 );
  //
  //  test2.print( );


  //std::cout << "compression rate: " << matrixV.getCompressionRatio( ) * 100 << " %" << std::endl;

  std::cout << "assembling K... ";
  std::cout.flush( );
  formK.assemble( matrixK );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "compression rate: " << matrixK.getCompressionRatio( ) * 100 << " %" << std::endl;

  // TESTING MULTIPLICATION BY MATRIX V

  Vector<LO, SC> v( matrixV.getNCols( ) );
  v.setAll( 1.0 );
  Vector<LO, SC> v2( matrixV.getNRows( ) );
  v2.setAll( 0.0 );

  std::cout << "multiplying by V... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  matrixV.apply( v, v2 );
  //v2.setAll(0.0);
  //matrixV.apply(v,v2);
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "assembling full V... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formV.assemble( fullMatrixV );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  Vector<LO, SC> v3( matrixV.getNRows( ) );
  v3.setAll( 0.0 );
  v.setAll( 1.0 );
  fullMatrixV.apply( v, v3 );
  for ( int i = 0; i < matrixV.getNRows( ); i++ ) {

    std::cout << v2.get( i ) << ", " << v3.get( i ) << std::endl;
  }
  v3.add( v2, -1.0 );
  std::cout << "Norm of difference for multiplication by V: " << v3.norm2( ) / v2.norm2( ) << std::endl;

  // TESTING MULTIPLICATION BY MATRIX K

  Vector<LO, SC> vk( matrixK.getNCols( ) );
  vk.setAll( 1.0 );
  Vector<LO, SC> vk2( matrixK.getNRows( ) );
  vk2.setAll( 0.0 );

  std::cout << "multiplying by K... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  matrixK.apply( vk, vk2 );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  std::cout << "assembling full K... ";
  std::cout.flush( );
  gettimeofday( &start, nullptr );
  formK.assemble( fullMatrixK );
  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  Vector<LO, SC> vk3( matrixK.getNRows( ) );
  vk3.setAll( 0.0 );
  vk.setAll( 1.0 );
  fullMatrixK.apply( vk, vk3 );

  vk3.add( vk2, -1.0 );

  std::cout << "Norm of difference for multiplication by K: " << vk3.norm2( ) / vk2.norm2( ) << std::endl;

  //fullMatrixK.print(std::cout);
  //for (int i=0; i < v.getLength(); i++) {
  //  //std::cout << v2.get(i) - v3.get(i) << std::endl;
  //  std::cout << vk3.get(i) << ", " << vk2.get(i) << std::endl;
  //}

  //  BESpace< LO, SC > bespaceV( &mesh, p0, p0 );
  //  FullMatrix< LO, SC > V( 0, 0 );
  //  BEBilinearFormLaplace1Layer< LO, SC > formV( bespaceV );
  //  std::cout << "Assembling V ... ";
  //  std::cout.flush();
  //  gettimeofday(&start, nullptr);
  //  formV.assemble( V );
  //  gettimeofday(&stop, nullptr);
  //  std::cout << "done in " << (double) ( timeDiff(start, stop) )  << " secs." << std::endl;
  //V.print();

}