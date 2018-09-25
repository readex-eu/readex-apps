/*!
 * @file    Lame1LayerP0P0MultilvlPrecond.cpp
 * @author  Lukas Maly
 * @date    June 20, 2016
 */

//#include "Lame1LayerP0P0MultilvlPrecond.h"

#ifdef LAME1LAYERP0P0MULTILVLPRECOND_H

namespace bem4i {

template<class LO, class SC>
Lame1LayerP0P0MultilvlPrecond<LO, SC>::
Lame1LayerP0P0MultilvlPrecond(
    SurfaceMesh3D<LO, SC> * mesh,
    LO maxElems
    ) {
  
  this->Vprecond = new Laplace1LayerP0P0MultilvlPrecond<LO, SC>( mesh, maxElems );

  this->init( );
  
}

template<class LO, class SC>
Lame1LayerP0P0MultilvlPrecond<LO, SC>::~Lame1LayerP0P0MultilvlPrecond( ) {

    if ( this->Vprecond ) delete this->Vprecond;
    if ( this->LameBlockPrecond ) delete this->LameBlockPrecond;
    if ( this->O ) delete this->O;
}

template<class LO, class SC>
void Lame1LayerP0P0MultilvlPrecond<LO, SC>::init( ) 
{

  LO VDimDomain = this->Vprecond->getDimDomain( );
  LO VDimRange = this->Vprecond->getDimRange( );

  this->dimDomain = 3 * VDimDomain;
  this->dimRange = 3 * VDimRange;
  
  this->O = new Zero<LO, SC>( VDimDomain , VDimRange );

  this->LameBlockPrecond = new BlockLinearOperator< LO, SC >( 3, 3 );
  this->LameBlockPrecond->setBlock( 0, 0, this->Vprecond );
  this->LameBlockPrecond->setBlock( 0, 1, this->O );
  this->LameBlockPrecond->setBlock( 0, 2, this->O );
  this->LameBlockPrecond->setBlock( 1, 0, this->O );
  this->LameBlockPrecond->setBlock( 1, 1, this->Vprecond );
  this->LameBlockPrecond->setBlock( 1, 2, this->O );
  this->LameBlockPrecond->setBlock( 2, 0, this->O );
  this->LameBlockPrecond->setBlock( 2, 1, this->O );
  this->LameBlockPrecond->setBlock( 2, 2, this->Vprecond );
  if ( !this->LameBlockPrecond->isValid( ) ) 
    std::cout << "LameBlockPreconditioner invalid!" << std::endl;

}

template<class LO, class SC> 
void Lame1LayerP0P0MultilvlPrecond< LO, SC >::apply(
      Vector< LO, SC > const & x,
      Vector< LO, SC > & y,
      bool transA,
      SC alpha,
      SC beta
      ) {

  this->LameBlockPrecond->apply( x, y, transA, alpha, beta );
}

/*
#ifdef BLAS_INT
template class Lame1LayerP0P0MultilvlPrecond<int, double>;
template class Lame1LayerP0P0MultilvlPrecond<int, float>;
#endif

#ifdef BLAS_LONG
template class Lame1LayerP0P0MultilvlPrecond<long, double>;
template class Lame1LayerP0P0MultilvlPrecond<long, float>;
#endif
*/

} // end of namespace

#endif
