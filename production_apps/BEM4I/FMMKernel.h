/*!
 * @file    FMMKernel.h
 * @author  Michal Merta 
 * @date    August 19, 2013
 * @brief   Header file for pure virtual class FMMKernel
 * 
 */

#ifndef FMMKERNEL_H
#define	FMMKERNEL_H

#include "Vector.h"
#include "Tree.h"
#include "FastBESpace.h"

namespace bem4i {

/*! 
 * Class representing the basis of all FMM kernels
 * 
 */
template<class LO, class SC>
class FMMKernel {
   
public:

  //! destructor
  virtual ~FMMKernel();
  
  //! method performs upward pass of the given tree
  virtual void upward() = 0;
  
  //! method performs downward pass of the given tree
  virtual void downward() = 0;
  
  //! setter for the multiplication vector
  inline void setVector(Vector<LO, SC> const * vector) {
    this->multVector = vector;
  }
  
  //! method returns left element cluster tree
  inline Tree<BECluster<LO, SC>*, LO >* getLeftClusterTree() const {
    return leftClusterTree;
  }
  
  //! method returns right element cluster tree
  inline Tree<BECluster<LO, SC>*, LO >* getRightClusterTree() const {
    return rightClusterTree;
  }
  
  //! returns whether the object is ready for matrix-vector multiplication 
  virtual bool isReady() const = 0;  
  
  virtual SC computeApproximation(BECluster<LO, SC>* cluster, LO element) = 0;
  
protected:
  
  //! vector by which the matrix is multiplied
  Vector<LO, SC> const * multVector;
  
  //! cluster tree associated with rows of matrix
  Tree<BECluster<LO, SC>*, LO >* leftClusterTree;
  
  //! cluster tree associated with columns of matrix
  Tree<BECluster<LO, SC>*, LO >* rightClusterTree;
  
  //! fast BE space
  FastBESpace<LO, SC>* beSpace;

  
};

}

// include .cpp file to overcome linking problems due to templates
#include "FMMKernel.cpp"

#endif	/* FMMKERNEL_H */

