/* 
 * File:   PotentialsWave.cpp
 * Author: mer126
 * 
 * Created on May 20, 2014, 2:58 PM
 */


#ifdef POTENTIALSWAVE_H

namespace bem4i {

template<class LO, class SC>
PotentialsWave<LO, SC>::PotentialsWave( const PotentialsWave& orig ) {
  this->space = orig.space;
}

template<class LO, class SC>
PotentialsWave<LO,SC>::PotentialsWave(
    BESpace<LO, SC>* space,
    Vector<LO, SCVT> *density
    ) {
  this->space = space;
  this->density = density;
}

template<class LO, class SC>
PotentialsWave<LO, SC>::~PotentialsWave( ) {
}

template<class LO, class SC>
void PotentialsWave<LO, SC>::doubleLayerPotential(
    const SCVT *x,
    LO n,
    Vector<LO, SC> & values,
    SCVT t
    ) const {

#pragma omp parallel
  {
    int rank = omp_get_thread_num( );
    LO nLocalElems = n / omp_get_num_threads( );
    LO myElems = nLocalElems;
    if ( rank == omp_get_num_threads( ) - 1 ) {
      LO remainder = n % omp_get_num_threads( );
      myElems = nLocalElems + remainder;
    }

    BEIntegratorWave<LO, SC> integrator( this->space,
        nullptr, 9 );
    Vector<LO, SC> localValues( myElems );
    integrator.doubleLayerPotential( x + 3 * omp_get_thread_num( )
        * nLocalElems, myElems, t, *( this->density ), localValues );
    for ( int i = 0; i < myElems; i++ ) {
      values.set( rank * nLocalElems + i, localValues.get( i ) );
    }
  }
}

template<class LO, class SC>
void PotentialsWave<LO, SC>::doubleLayerPotential(
    Mesh<LO, SC> &mesh,
    Vector<LO, SC> & values,
    SCVT t
    ) const {

  SCVT *nodes = &( ( *( mesh.getNodes( ) ) )[0] );
  LO nNodes = mesh.getNNodes( );

  this->doubleLayerPotential( nodes, nNodes, values, t );
}

}

#endif
