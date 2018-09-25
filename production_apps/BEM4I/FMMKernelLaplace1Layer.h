/*!
 * @file    FMMKernelLaplace1Layer.h
 * @author  Michal Merta 
 * @date    August 19, 2013
 * @brief   Header file for class FMMKernelLaplace1Layer
 * 
 */

#ifndef FMMKERNELLAPLACE1LAYER_H
#define	FMMKERNELLAPLACE1LAYER_H

#include "FMMKernel.h"
#include "Quadratures.h"
#include <complex>

namespace bem4i {

  template<class LO, class SC>
  class FMMKernelLaplace1Layer : public FMMKernel<LO, SC> {
  public:

    //! default constructor
    FMMKernelLaplace1Layer();

    //! copy constructor
    FMMKernelLaplace1Layer(const FMMKernelLaplace1Layer& orig);

    /*! constructor taking FastBESpace as argument
     *
     * @param[in]       outerElem index of the outer integration triangle
     * @param[in]       innerElem index of the inner integration triangle
     * @param[in,out]   element matrix
     */
    FMMKernelLaplace1Layer(FastBESpace<LO, SC>* space, int nMax,
            int startLevel = 0, int quadOrder = 5);

    //! destructor 
    virtual ~FMMKernelLaplace1Layer();

    //! performs upward pass in tree
    virtual void upward();

    //! performs downward pass in tree
    virtual void downward();

    virtual bool isReady() const {
      return (up && down);
    }

    //! computes approximation
    virtual SC computeApproximation(BECluster<LO, SC>* cluster, LO element);

  private:

    //! order of the multipole expansion
    int nMax;

    //! order of quadrature for moments evaluation
    int quadOrder;

    //! starting level for tree passing
    int startLevel;

    //! mapping from index i to indices n,m 
    int **ind2nm;

    //! whether the upward and downward pass were performed
    bool up, down;

    inline void cart2sph(SC *cart, SC *sph) {
      // r, phi, theta
      sph[0] = std::sqrt(cart[0] * cart[0] + cart[1] * cart[1] + 
          cart[2] * cart[2]);
      sph[1] = atan2(cart[1], cart[0]);
      sph[2] = acos(cart[2] / sph[0]);
    }

    /*!
     * 
     * returns array of values of spherical harmonic R^n_m
     * use function getRnm to obtain specific value
     * 
     * values are stored in the 1D array in the form
     * [R^0_0, R^1_0, ..., R^n_0, R^1_1, R^2_1, ..., R^n_1, R^2_2, ...]
     * @param[in]     nMax    order of approximation
     * @param[in]     x       cartesian coordinate of the point
     * @param[in,out] Rarray  spherical harmonic R^n_m
     */
    void computeR(int nMax, SC *x, std::complex<SC> *Rarray);

    /*!
     * 
     * returns array of values of spherical harmonic S^n_m
     * use function getRnm to obtain specific value
     * 
     * values are stored in the 1D array in the form
     * [S^0_0, S^1_0, ..., S^n_0, S^1_1, S^2_1, ..., S^n_1, S^2_2, ...]
     * @param[in]     nMax    order of approximation
     * @param[in]     x       cartesian coordinate of the point
     * @param[in,out] Sarray  spherical harmonic R^n_m
     */
    void computeS(int nMax, SC *x, std::complex<SC> *Sarray);

    //! method returns R^n_m from given array

    inline std::complex<SC> getRnm(std::complex<SC> *R, int nMax, int n, int m) {
      if ((abs(m) > abs(n)) || (n > nMax)) {
        return 0.0;
      }
      if (m >= 0) {
        return R[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ];
      }
      m = -m;
      return ((SC) pow(-1.0, m)) * conj(R[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ]);
    }

    //! method returns S^n_m from given array

    inline std::complex<SC> getSnm(std::complex<SC> *S, int nMax, int n, int m) {
      if ((abs(m) > abs(n)) || (n > nMax)) {
        return 0.0;
      }
      if (m >= 0) {
        return S[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ];
      }
      m = -m;
      return ((SC) pow(-1.0, m)) * conj(S[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ]);
    }

    //! method returns multipole moments (n,m) from given arru

    inline std::complex<SC> getMoment(std::complex<SC> *M, int nMax, int n, int m) {
      if ((abs(m) > abs(n)) || (n > nMax)) {
        return 0.0;
      }
      if (m >= 0) {
        return M[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ];
      }
      m = -m;
      return ((SC) pow(-1.0, m)) * conj(M[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ]);
    }

    //! sets value of array M

    inline void setMnm(std::complex<SC> *M, int nMax, int n, int m, std::complex<SC> val) {
      if (m >= 0) {
        M[ (int) ((2 * nMax + 3 - m) * m / 2) + (n - m) ] = val;
      }
    }

    //! perform recursively upward pass in the tree
    void doUpwardPass(TreeMember<BECluster<LO, SC>* >* node);

    //! perform recursively downward pass in the tree
    void doDownwardPass(TreeMember<BECluster<LO, SC>* >* node);

    //! computes multipole moments of given cluster
    void computeMoments(TreeMember<BECluster<LO, SC>* >* node);

    //! computes local expansion coefficients of given cluster
    void computeLocExp(TreeMember<BECluster<LO, SC>* >* node);

    // include variables dependent on parallelization version
#ifdef FMM_SCALAR
    //! temporary store for R coefficients computations
    std::complex<SC> *tempR;

    //! temporary store for M coefficients computations
    std::complex<SC> *tmpM;
#elif defined(FMM_OMP)
    //! temporary store for R coefficients computations
    std::complex<SC> **tempR;

    //! temporary store for M coefficients computations
    std::complex<SC> **tmpM;
#endif

  };
}

// include .cpp file to overcome linking problems due to templates
#include "FMMKernelLaplace1Layer.cpp"

#endif	/* FMMKERNELLAPLACE1LAYER_H */

