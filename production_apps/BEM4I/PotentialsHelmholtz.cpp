/*!
 * @file    PotentialsHelmholtz.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    April 17, 2015
 * 
 */

#ifdef POTENTIALSHELMHOLTZ_H

namespace bem4i {

template<class LO, class SC>
PotentialsHelmholtz<LO, SC>::PotentialsHelmholtz(
    BESpace<LO, SC> * space,
    SC kappa,
    Vector<LO, SC> * density,
    int quadOrder
    ) {

  this->space = space;
  this->kappa = kappa;
  this->density = density;
  this->quadOrder = quadOrder;
}

template<class LO, class SC>
PotentialsHelmholtz<LO, SC>::~PotentialsHelmholtz( ) {
}

template<class LO, class SC>
void PotentialsHelmholtz<LO, SC>::singleLayerPotential(
    const SCVT * x,
    LO n,
    Vector<LO, SC> & values,
    SCVT t
    ) const {

  int qOrder[ 2 ] = { this->quadOrder, this->quadOrder };

#pragma omp parallel
  {
    int rank = omp_get_thread_num( );
    LO nLocalElems = n / omp_get_num_threads( );
    LO myElems = nLocalElems;
    if ( rank == omp_get_num_threads( ) - 1 ) {
      LO remainder = n % omp_get_num_threads( );
      myElems = nLocalElems + remainder;
    }

    BEIntegratorHelmholtz<LO, SC> integrator( this->space, qOrder,
        this->kappa );
    Vector<LO, SC> localValues( myElems );
    integrator.singleLayerPotential( x + 3 * omp_get_thread_num( )
        * nLocalElems, myElems, *this->density, localValues );
    for ( LO i = 0; i < myElems; ++i ) {
      values.set( rank * nLocalElems + i, localValues.get( i ) );
    }
  }
}

template<class LO, class SC>
void PotentialsHelmholtz<LO, SC>::singleLayerPotential(
    Mesh<LO, SC> & mesh,
    Vector<LO, SC> & values,
    SCVT t ) const {

  SCVT * nodes = mesh.getNodes( )->data( );
  LO nNodes = mesh.getNNodes( );

  this->singleLayerPotential( nodes, nNodes, values, t );
}

template<class LO, class SC>
void PotentialsHelmholtz<LO, SC>::doubleLayerPotential(
    const SCVT * x,
    LO n,
    Vector<LO, SC> & values,
    SCVT t
    ) const {

  int qOrder[ 2 ] = { this->quadOrder, this->quadOrder };

  int mpiSize = 1;
  int mpiRank = 0;
  int mpiInitialized = 0;
  MPI_Initialized( &mpiInitialized );
  if ( mpiInitialized != 0 ) {
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
  }

  int *numElems = new int[mpiSize];
  int *displs = new int[mpiSize];
  displs[0] = 0;
  int nElemsPerMPI = n / mpiSize;
  for ( int i = 0; i < mpiSize - 1; i++ ) {
    numElems[i] = nElemsPerMPI;
    displs[i + 1] = displs[i] + numElems[i];
  }
  numElems[mpiSize - 1] = nElemsPerMPI + ( n % mpiSize );

  SC *valuesPerMPI = new SC[numElems[mpiRank]];

  #pragma omp parallel
  {
    int ompSize = omp_get_num_threads( );
    int ompRank = omp_get_thread_num( );

    LO nLocalElems = numElems[mpiRank] / ompSize;
    LO myElems = nLocalElems;
    if ( ompRank == ompSize - 1 ) {
      myElems = nLocalElems + numElems[mpiRank] % ompSize;
    }

    BEIntegratorHelmholtz<LO, SC> integrator( this->space, qOrder,
        this->kappa );
    Vector<LO, SC> localValues( myElems );
    integrator.doubleLayerPotential( x + 3 * ( displs[mpiRank] +
        ompRank * nLocalElems ), myElems, *this->density, localValues );

    for ( LO i = 0; i < myElems; ++i ) {
      valuesPerMPI[ompRank * nLocalElems + i] = localValues.get( i );
    }
  }

  if ( mpiInitialized != 0 ) {
    MPI_Datatype MPI_SC = GetType<LO, SC>::MPI_SC( );
    MPI_Allgatherv( valuesPerMPI, numElems[mpiRank],
        MPI_SC, values.getData( ), numElems, displs,
        MPI_SC, MPI_COMM_WORLD );
  } else {
    memcpy( values.getData( ), valuesPerMPI, sizeof (SC ) * numElems[0] );
  }

  delete [] numElems;
  delete [] displs;
  delete [] valuesPerMPI;
}

template<class LO, class SC>
void PotentialsHelmholtz<LO, SC>::doubleLayerPotential(
    Mesh<LO, SC> & mesh,
    Vector<LO, SC> & values,
    SCVT t
    ) const {

  SCVT * nodes = mesh.getNodes( )->data( );
  LO nNodes = mesh.getNNodes( );

  this->doubleLayerPotential( nodes, nNodes, values, t );
}

}

#endif
