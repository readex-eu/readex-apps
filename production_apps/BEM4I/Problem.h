/*!
 * @file    Problem.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    June 10, 2014
 * @brief   Header file for class pure virtual class Problem
 * 
 */

#ifndef PROBLEM_H
#define	PROBLEM_H

namespace bem4i {

/*!
 * Pure virtual class representing a problem
 * 
 */
template< class LO, class SC >
class Problem {

  typedef typename GetType<LO, SC>::SCVT SCVT;

public:

  //! destructor

  virtual ~Problem( ) {
  };

  /*!
   * Solves the problem 
   */
  virtual bool solve( ) = 0;
  
  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      SCVT value
      ) {}


  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      const string & value
      ) {}

  // todo: deprecated, use templated version instead
  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      void* value
      ) {}

  /*!
   * Sets parameters of a solver
   */
  virtual void set(
      int name,
      int value
      ) {}  

  protected:
    const int maxParameters = 50;
  
  private:   
};

}


#endif	/* PROBLEM_H */

