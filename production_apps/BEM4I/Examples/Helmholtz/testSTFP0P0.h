#include "../../Settings.h"


#include <iostream>
#include <mpi.h>

#include "../auxiliary.h"
#include "../../ProgressMonitor.h"
#include "../../FullMatrix.h"
#include "../../SurfaceMesh3D.h"
#include "../../FastBESpace.h"
#include "../../Vector.h"
#include "../../BEIntegratorHelmholtz.h"
#include "../../BEBilinearFormHelmholtzHypersingular.h"
#include "../../BEBilinearFormHelmholtz1Layer.h"
#include "../../BEBilinearFormHelmholtz2Layer.h"
#include "../../HelmholtzHypersingularP0P0Operator.h"
#include "../../IdentityOperator.h"
#include "../../RepresentationFormulaHelmholtz.h"
#include "../../STFOperator.h"
#include "../../STFPreconditioner.h"

using namespace std;
using namespace bem4i;

void testSTFP0P0(
    string const &filename,
    int refine4,
    int refine9,
    int nSources,
    double * sources,
    int nDevices,
    double * devices,
    double * leftBottom,
    double * rightTop,
    int nSlices,
    double startZ,
    double endZ,
    double waveLength,
    string const &outputFolder
    );

int main(
    int argc,
    char** argv
    ) {
  intro( );

  string filename = string( argv[ 1 ] );
  int refine4 = atoi( argv[ 2 ] );
  int refine9 = atoi( argv[ 3 ] );
  string sourceDevFilename = string( argv[4] );
  double leftBottom[2];
  double rightTop[2];
  leftBottom[0] = atof( argv[ 5 ] );
  leftBottom[1] = atof( argv[ 6 ] );
  rightTop[0] = atof( argv[ 7 ] );
  rightTop[1] = atof( argv[ 8 ] );
  int nSlices = atoi( argv[ 9 ] );
  double startZ = atof( argv[ 10 ] );
  double endZ = atof( argv[ 11 ] );
  double waveLength = atof( argv[ 12 ] );
  string outputFolder = string( argv[ 13 ] );

  std::ifstream file( sourceDevFilename.c_str( ) );

  int nSources = 0;
  int nDevices = 0;
  if ( !file.is_open( ) ) {
    std::cout << "Device file could not be opened!" << std::endl;
  }

  file >> nSources;
  double * sources = new double[nSources * 3];
  for ( int i = 0; i < nSources * 3; i++ ) {
    file >> sources[i];
  }
  file >> nDevices;
  double * devices = new double[nDevices * 3];
  for ( int i = 0; i < nDevices * 3; i++ ) {
    file >> devices[i];
  }
  testSTFP0P0( filename, refine4, refine9,
      nSources, sources, nDevices, devices,
      leftBottom, rightTop, nSlices, startZ, endZ, waveLength, outputFolder );
  delete [] devices;
  delete [] sources;

  return 0;
}

void testSTFP0P0(
    string const &filename,
    int refine4,
    int refine9,
    int nSources,
    double * sources,
    int nDevices,
    double * devices,
    double * leftBottom,
    double * rightTop,
    int nSlices,
    double startZ,
    double endZ,
    double waveLength,
    string const &outputFolder
    ) {
  MPI_Init( nullptr, nullptr );
  typedef complex<double> SC;
  typedef int LO;
  typedef typename SC::value_type SCVT;
  timeval start, stop;

  SCVT gmresprec = 1e-4;
  LO gmresmaxit = 200;

  ProgressMonitor::init( "Initializing mesh" );
  SurfaceMesh3D< LO, SC > mesh;
  mesh.load( filename.c_str( ) );

  mesh.refine( refine4, 2 );
  mesh.refine( refine9, 3 );
  ProgressMonitor::step( );
  mesh.printInfo( );

  quadratureType quadType = SauterSchwab;
  int orderNear = 3;
  int orderFar = 3;
  int quadOrder[] = { orderNear, orderNear, orderNear, orderNear };
  int quadDisjoint[] = { orderFar, orderFar };

  SC iUnit( 0.0, 1.0 );
  SC etai = 1.0;
  SC etae( 5.0, -0.05 );
  SCVT lambda = (SCVT) waveLength;
  SCVT kappa = 2.0 * M_PI / lambda;
  SC refrIndex( std::sqrt( ( std::abs( etae ) + etae.real( ) ) / 2.0 ),
      std::sqrt( ( std::abs( etae ) - etae.real( ) ) / 2.0 ) );
  SC ki = kappa * std::sqrt( etai );
  SC ke = kappa * refrIndex;

  std::cout << "kappa_i = " << ki << std::endl;
  std::cout << "kappa_e = " << ke << std::endl;

  Tree<BECluster<LO, SC>*, LO> tree;
  mesh.nestedDissection( tree, 120 );
  FastBESpace< LO, SC > bespace( &mesh, p0, p0, &tree, 0.9, 8, 5, 0 );
  bespace.setScaleACA( 1e-12 );
  bespace.setEpsilonACA( 1e-2 );

  // matrices V (interior, exterior)
  gettimeofday( &start, nullptr );
  MPIACAMatrix< LO, SC > * Vi = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtz1Layer< LO, SC > formVi( &bespace, quadOrder,
      ki, quadType, quadDisjoint );
  formVi.assemble( *Vi );
  MPI_Barrier( MPI_COMM_WORLD );

  MPIACAMatrix< LO, SC > * Ve = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtz1Layer< LO, SC > formVe( &bespace, quadOrder,
      ke, quadType, quadDisjoint );
  formVe.assemble( *Ve );
  MPI_Barrier( MPI_COMM_WORLD );

  // matrices K (interior, exterior)
  MPIACAMatrix< LO, SC > * Ki = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtz2Layer< LO, SC > formKi( &bespace, quadOrder,
      ki, quadType, quadDisjoint );
  formKi.assemble( *Ki );
  MPI_Barrier( MPI_COMM_WORLD );

  MPIACAMatrix< LO, SC > * Ke = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtz2Layer< LO, SC > formKe( &bespace, quadOrder,
      ke, quadType, quadDisjoint );
  formKe.assemble( *Ke );
  MPI_Barrier( MPI_COMM_WORLD );

  // matrices D (interior, exterior)
  MPIACAMatrix< LO, SC > * H1i = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  MPIACAMatrix< LO, SC > * H2i = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formDi( &bespace, quadOrder,
      ki, quadType, quadDisjoint );
  formDi.assembleH1P0P0( *H1i );
  MPI_Barrier( MPI_COMM_WORLD );
  formDi.assembleH2P0P0( *H2i );
  HelmholtzHypersingularP0P0Operator<LO, SC> * Di =
      new HelmholtzHypersingularP0P0Operator<LO, SC>( H1i, H2i, mesh );
  delete H1i;
  MPI_Barrier( MPI_COMM_WORLD );

  MPIACAMatrix< LO, SC > * H1e = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  MPIACAMatrix< LO, SC > * H2e = new MPIACAMatrix< LO, SC >( 0, 0,
      MPI_COMM_WORLD );
  BEBilinearFormHelmholtzHypersingular< LO, SC > formDe( &bespace, quadOrder,
      ke, quadType, quadDisjoint );
  formDe.assembleH1P0P0( *H1e );
  MPI_Barrier( MPI_COMM_WORLD );
  formDe.assembleH2P0P0( *H2e );
  HelmholtzHypersingularP0P0Operator<LO, SC> * De =
      new HelmholtzHypersingularP0P0Operator<LO, SC>( H1e, H2e, mesh );
  delete H1e;
  MPI_Barrier( MPI_COMM_WORLD );

  // identity operator
  SparseMatrix<LO, SC> M;
  IdentityOperator< LO, SC > id( &bespace );
  id.assemble( M );

  gettimeofday( &stop, nullptr );
  std::cout << "done in " << (double) ( timeDiff( start, stop ) ) << " secs." << std::endl;

  STFOperator< LO, SC > stf( &mesh, ki, ke, Vi, Ve, Ki, Ke, Di, De, &M, true );

  LO nElems = mesh.getNElements( );

  Vector< LO, SC > f( nElems );
  Vector< LO, SC > g( nElems );

  string file_scatter = "input/scatter_grid_z_big.txt";
  SurfaceMesh3D< LO, SC > ms;
  SurfaceMesh3D< LO, SC > mi;
  SurfaceMesh3D< LO, SC > me;
  ms.load( file_scatter.c_str( ) );
  //ms.scale(5 ); //orig mesh
  ms.refine( 7 ); // orig mesh
  //ms.move(4.2, 3.5, 1.2); //orig mesh

  ms.move( 1.0, 1.0, startZ ); //orig mesh
  ms.scale( ( rightTop[0] - leftBottom[0] ) / 2.0,
      ( rightTop[1] - leftBottom[1] ) / 2.0, 1.0 ); //orig mesh
  ms.move( leftBottom[0], leftBottom[1], 0.0 );
  //ms.refine( 9 ); // orig mesh

  ms.printInfo( );
  mesh.trimEvaluationGrid( ms, mi, me );

  // loop over source positions
  for ( int s = 0; s < nSources; s++ ) {
    const std::string createDir = "mkdir -p ./" + outputFolder
        + "/router_" + std::to_string( (long long) s );
    if ( system( createDir.c_str( ) ) ) exit( -1 );
    //SCVT y[ 3 ] = { 5.0, 1.0, 1.0 }; // initial mesh
    SCVT * y = &( sources[3 * s] ); //{ 8.0, 10.0, 1.0 }; // flat2 mesh
    SCVT x[ 3 ], x1[3], x2[3], x3[3], n[ 3 ];
    SCVT norm, dot;

    int rhsOrder = 5;
    int qSize = quadSizes[ rhsOrder ];
    SCVT * quadNodes = new SCVT[ 3 * qSize ];
    SCVT * ya;
    SC val;

    for ( LO i = 0; i < nElems; ++i ) {
      val = 0.0;
      mesh.getNodes( i, x1, x2, x3 );
      mesh.getQuadratureNodes( x1, x2, x3, quadNodes, rhsOrder );
      for ( LO j = 0; j < qSize; j++ ) {
        ya = quadNodes + 3 * j;
        norm = DIST3( ya, y );
        val += -PI_FACT * std::exp( iUnit * ki * norm ) / norm *
            quadWeights[ rhsOrder ][ j ];
      }
      f.set( i, val );
    }

    for ( LO i = 0; i < nElems; ++i ) {
      mesh.getNormal( i, n );
      mesh.getNodes( i, x1, x2, x3 );
      val = 0.0;

      mesh.getQuadratureNodes( x1, x2, x3, quadNodes, rhsOrder );
      for ( LO j = 0; j < qSize; j++ ) {
        ya = quadNodes + 3 * j;
        norm = DIST3( ya, y );
        dot = n[ 0 ] * ( ya[ 0 ] - y[ 0 ] ) + n[ 1 ] * ( ya[ 1 ] - y[ 1 ] ) +
            n[ 2 ] * ( ya[ 2 ] - y[ 2 ] );
        val += -( iUnit * ki * norm - 1.0 ) *
            std::exp( iUnit * ki * norm ) * dot / ( norm * norm * norm ) *
            quadWeights[ rhsOrder ][ j ];
      }
      g.set( i, val );
    }
    delete [] quadNodes;
    Vector< LO, SC > rhs;
    stf.getRHS( g, f, rhs );

    //rhs.print();

    ProgressMonitor::init( "Solving the system" );
    Vector< LO, SC > res( nElems + nElems );
    STFPreconditioner< LO, SC > stfPrec( &mesh, Vi, Ve, Ki, Ke, Di, De, &M, true );
    stf.FGMRESSolve( rhs, res, gmresprec, gmresmaxit, gmresmaxit, &stfPrec );
    ProgressMonitor::step( );

    //res.print(); return;

    Vector< LO, SC > ui, ue, ti, te;
    stf.getTraces( res, f, g, ui, ue, ti, te );

    std::stringstream file;
    file.fill( '0' );
    file.width( 4 );

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    double step = ( endZ - startZ ) / ( (double) nSlices );

    Vector<LO, SC> * timeEval = nullptr;
    RepresentationFormulaHelmholtz<LO, SC> fi( &bespace, &ui, &ti, ki );
    RepresentationFormulaHelmholtz<LO, SC> fe( &bespace, &ue, &te, ke );
    for ( int i = 0; i < nSlices; i++ ) {

      ProgressMonitor::init( "Evaluating inside" );
      LO npi = mi.getNNodes( );
      Vector< LO, SC > resi( npi );

      fi.evaluate( mi, true, resi );
      for ( LO j = 0; j < mi.getNNodes( ); ++j ) {
        mi.getNode( j, x );
        norm = DIST3( x, y );
        resi.add( j, PI_FACT * std::exp( iUnit * ki * norm ) / norm );
      }
      ProgressMonitor::step( );

      file.str( std::string( ) );
      file.clear( );

      file << outputFolder << "/router_" << s << "/interior_slice_" << i << ".vtu";

      string namesi[] = { "ui" };
      Vector< LO, SC >* datai[] = { &resi };
      if ( rank == 0 ) {
        mi.printParaviewVtu( file.str( ).c_str( ), 1, namesi, datai, 0, nullptr,
            nullptr );
      }
      MPI_Barrier( MPI_COMM_WORLD );

      if ( i == nSlices / 2 ) {
        timeEval = new Vector<LO, SC>( resi );
      }
      /*
      ProgressMonitor::init( "Evaluating outside" );
      LO npe = me.getNNodes( );
      Vector< LO, SC > rese( npe );

      fe.evaluate( me, false, rese );
      ProgressMonitor::step( );

      file.str( std::string( ) );
      file.clear( );
      file << outputFolder << "/exterior_slice_" << i << ".vtu";

      string namese[] = { "ue" };
      Vector< LO, SC >* datae[] = { &rese };
      if (rank == 0) {
        me.printParaviewVtu( file.str().c_str(), 1, namese, datae, 0, nullptr,
            nullptr );
      }
      me.move(0.0, 0.0, step);
       */
      mi.move( 0.0, 0.0, step );


    }

    SCVT velocity = 350.0;
    SCVT endTime = 2.0 * M_PI / ( velocity * ki.real( ) );
    int nSteps = 40;
    SCVT timeStep = endTime / ( (double) nSteps );
    Vector< LO, SC > *pointSolutionTotalTime;
    std::stringstream file2;
    file2.fill( '0' );
    file2.width( 4 );
    Vector< LO, SC >* nodalDataSol[1];
    string nodeNamesSol[1];
    for ( int i = 0; i < nSteps; i++ ) {
      pointSolutionTotalTime = new Vector< LO, SC >( *timeEval );
      pointSolutionTotalTime->scale(
          std::exp( -iUnit * velocity * ki * (SCVT) i * timeStep ) );
      nodalDataSol[ 0 ] = pointSolutionTotalTime;
      nodeNamesSol[ 0 ] = "ui";
      file2.str( std::string( ) );
      file2.clear( );
      file2 << outputFolder << "/router_" << s << "/solution_time_" << i
          << ".vtu";
      if ( rank == 0 ) {
        mi.printParaviewVtu( file2.str( ).c_str( ), 1, nodeNamesSol,
            nodalDataSol, 0, nullptr, nullptr );
      }
      delete pointSolutionTotalTime;
    }


    ui.add( f, -1.0 );
    ti.add( g, -1.0 );

    string elemNames[] = { "ui", "ue", "ti", "te" };
    //string elemNames[] = { "ti", "te" };
    Vector< LO, SC >* elemData[] = { &ui, &ue, &ti, &te };
    //Vector< LO, SC >* elemData[] = { &ti, &te };
    if ( rank == 0 ) {
      mesh.printParaviewVtu( "output/surface.vtu", 0, nullptr, nullptr, 4,
          elemNames, elemData );
    }

    // get the strength of signal in given points
    SCVT devX[3];
    for ( int i = 0; i < nDevices; i++ ) {
      devX[0] = devices[3 * i];
      devX[1] = devices[3 * i + 1];
      devX[2] = devices[3 * i + 2];
      SCVT startX = devX[0] - 0.05;
      SCVT endX = devX[0] + 0.05;
      SCVT startY = devX[1] - 0.05;
      //SCVT endY = devX[1] + 0.05;
      nSteps = 10;
      SCVT step = ( endX - startX ) / ( (SCVT) nSteps );
      SCVT integral = 0.0;
      for ( LO j = 0; j < nSteps; j++ ) {
        for ( LO k = 0; k < nSteps; k++ ) {
          SCVT center[3] = { startX + step / 2.0 + (SCVT) step * j,
            startY + step / 2.0 + (SCVT) step * k, devX[2] };
          Vector<LO, SC> resQuad( 1 );
          fi.evaluate( center, 1, true, resQuad );
          SC val = resQuad.get( 0 );
          integral += val.real( ) * val.real( ) * step*step;
        }
      }
      std::cout << "Device " << i << ": " << integral << std::endl;
      ;
    }

    delete timeEval;
  }
  delete Vi;
  delete Ve;
  delete Ki;
  delete Ke;
  delete H2i;
  delete H2e;
  delete Di;
  delete De;
  MPI_Finalize( );
}