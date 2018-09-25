/*!
 * @file    PotentialsLaplace.cpp
 * @author  Jan Zapletal
 * @date    February 17, 2015
 * 
 */

#ifdef POTENTIALSLAPLACE_H

namespace bem4i {

template<class LO, class SC>
PotentialsLaplace<LO, SC>::PotentialsLaplace(
    BESpace<LO, SC> * space,
    Vector<LO, SCVT> * density
    ) {

  this->space = space;
  this->density = density;
}

template<class LO, class SC>
PotentialsLaplace<LO, SC>::~PotentialsLaplace( ) {
}

template<class LO, class SC>
void PotentialsLaplace<LO, SC>::doubleLayerPotential(
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
    int qOrder[4] = {3, 3, 3, 3};
    if ( rank == omp_get_num_threads( ) - 1 ) {
      LO remainder = n % omp_get_num_threads( );
      myElems = nLocalElems + remainder;
    }

    BEIntegratorLaplace<LO, SC> integrator( this->space, qOrder );
    Vector<LO, SC> localValues( myElems );
    integrator.doubleLayerPotential( x + 3 * omp_get_thread_num( )
        * nLocalElems, myElems, *this->density, localValues );
    for ( int i = 0; i < myElems; i++ ) {
      values.set( rank * nLocalElems + i, localValues.get( i ) );
    }
  }
  
}

template<class LO, class SC>
void PotentialsLaplace<LO, SC>::doubleLayerPotential(
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
