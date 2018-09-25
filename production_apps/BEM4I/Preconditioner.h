/*!
 * @file    Preconditioner.h
 * @author  Michal Merta 
 * @date    March 12, 2014
 * @brief   Header file for abstract class Preconditioner
 * 
 */

#ifndef PRECONDITIONER_H
#define	PRECONDITIONER_H

#include "LinearOperator.h"

namespace bem4i {

template <class LO, class SC>
class LinearOperator;

/*! 
 * Abstract class representing a preconditioner
 * 
 * the class applies the inverse of a matrix spectrally equivalent
 * to a given system matrix
 */
template<class LO, class SC>
class Preconditioner : public LinearOperator<LO, SC> {

public:

private:

};

}

#endif	/* PRECONDITIONER_H */
