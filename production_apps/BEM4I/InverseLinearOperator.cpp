/*!
 * @file    InverseLinearOperator.cpp
 * @author  Jan Zapletal
 * @date    April 26, 2016
 * 
 */


#ifdef INVERSELINEAROPERATOR_H

namespace bem4i {

template< class LO, class SC >
void InverseLinearOperator< LO, SC >::apply(
    Vector< LO, SC > const & x,
    Vector< LO, SC > & y,
    bool transA,
    SC alpha,
    SC beta
    ) {

  Vector< LO, SC > z( this->dimRange, true );
  //z.setAll( 1.0 );

  switch ( this->params->solver ) {
    case( SolverParameters< LO, SC >::CG ):
      this->op->CGSolve( x, z, this->params->precision, this->params->maxIter, 
      this->precond, this->params->msg );
      break;
    case( SolverParameters< LO, SC >::GMRES ):
      this->op->GMRESSolve( x, z, this->params->precision,
          this->params->maxIter, this->params->maxIter, this->precond,
          this->params->msg );
      break;
    case( SolverParameters< LO, SC >::DGMRES ):
      std::cout << "Not implemented" << std::endl;
      break;
    case( SolverParameters< LO, SC >::FGMRES ):
      this->FGMRESSolve( x, z, this->params->precision,
          this->params->maxIter, this->params->maxIter, this->precond,
          this->params->msg );
      break;
  }

  y.scale( beta );
  y.add( z, alpha );
}

}
#endif
