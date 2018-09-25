/*!
 * @file    PotentialsHelmholtz.h
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    April 17, 2015
 * @brief   Header file for class PotentialsHelmholtz
 * 
 */

#ifndef POTENTIALSHELMHOLTZ_H
#define	POTENTIALSHELMHOLTZ_H

#include "Vector.h"
#include "BESpace.h"
#include "Macros.h"
#include "Mesh.h"
#include "Potentials.h"
#include "BEIntegratorHelmholtz.h"

namespace bem4i {

    /*! 
     * Class representing single and double layer potentials for Laplace equation
     * 
     */
    template<class LO, class SC>
    class PotentialsHelmholtz : public Potentials<LO, SC> {
        typedef typename GetType<LO, SC>::SCVT SCVT;

    public:

        /*!
         * constructor taking Dirichlet and Neumann data as arguments
         * 
         * @param[in]     space
         * @param[in]     dirichlet - vector of dirichlet data
         * @param[in]     neumann - vector of neumann data
         * 
         */
        PotentialsHelmholtz(
                BESpace<LO, SC> * space,
                SC kappa,
                Vector<LO, SC> * density,
                int quadOrder = 5
                );

        //! destructor
        virtual ~PotentialsHelmholtz();

        /*!
         * function evaluates single layer potential in given array of points
         *  
         * @param[in]       x
         * @param[in,out]   result
         */
        virtual void singleLayerPotential(
                const SCVT * x,
                LO n,
                Vector<LO, SC> & values,
                SCVT t = 0.0
                ) const;

        /*!
         * function evaluates double layer potential in given array of points
         *  
         * @param[in]       x
         * @param[in,out]   result
         */
        virtual void doubleLayerPotential(
                const SCVT * x,
                LO n,
                Vector<LO, SC> & values,
                SCVT t = 0.0
                ) const;

        /*!
         * function evaluates single layer potential in nodes of given mesh
         *  
         * @param[in]       mesh input mesh
         * @param[in,out]   values
         */
        virtual void singleLayerPotential(
                Mesh<LO, SC> & mesh,
                Vector<LO, SC> & values,
                SCVT t = 0.0
                ) const;

        /*!
         * function evaluates double layer potential in nodes of given mesh
         *  
         * @param[in]       mesh input mesh
         * @param[in,out]   values
         */
        virtual void doubleLayerPotential(
                Mesh<LO, SC> & mesh,
                Vector<LO, SC> & values,
                SCVT t = 0.0
                ) const;

    private:

        PotentialsHelmholtz() {
        };

        //! copy constructor
        PotentialsHelmholtz(
                const PotentialsHelmholtz & orig
                );

        SC kappa;
        
        int quadOrder;

    };

} // end of namespace bem4i

// include .cpp file to overcome linking problems due to templates
#include "PotentialsHelmholtz.cpp"

#endif	/* POTENTIALSHELMHOLTZ_H */
