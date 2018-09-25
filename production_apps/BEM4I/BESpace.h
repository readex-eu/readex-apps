/*!
 * @file    BESpace.h
 * @author  Michal Merta 
 * @date    July 9, 2013
 * @brief   Header file for class BESpace
 * 
 */

#ifndef BESPACE_H
#define	BESPACE_H

#include <vector>
//#include "SurfaceMesh3D.h"

namespace bem4i {

template< class LO, class SC >
class SurfaceMesh3D;

using std::vector;

enum basisType {
  p0,
  p1,
  p1dis

};

// mapping from basis type to number of basis functions
const int basisType2functions[3] = { 1, 3, 3 };

/*!
 * Class representing a boundary element space over particular mesh
 * Only for scalar quantities
 */
template<class LO, class SC>
class BESpace {
public:

  //! default constructor
  BESpace( );

  //! copy constructor
  BESpace(
      const BESpace & orig
      );

  /*!
   * constructor taking a mesh and ansatzs types
   * 
   * @param[in]   mesh
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (constant/linear)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (constant/linear)
   */
  BESpace(
      SurfaceMesh3D<LO, SC> * mesh,
      basisType ansatzFunctionType,
      basisType testFunctionType
      );

  /*
   * constructor taking a mesh and ansatz types
   * 
   * @param[in]   mesh
   * @param[in]   innerElems        indices of elements over which the inner integral is evaluated
   * @param[in]   outerElems        indices of elements over which the outer integral is evaluated
   * @param[in]   innerAnsatzType   ansatz type of the inner integral (p0/p1)
   * @param[in]   outerAnsatzType   ansatz type of the outer integral (p0/p1)
   */
  // BESpace(SurfaceMesh3D<LO>& mesh, vector<LO> innerElems, vector<LO> outerElems,
  //         basisType innerAnsatzType, basisType outerAnsatzType);

  //! destructor
  virtual ~BESpace( );

  //! returns the associated mesh

  inline SurfaceMesh3D<LO, SC>* getMesh( ) const {
    return mesh;
  }

  //! returns the associated left (row) mesh

  inline SurfaceMesh3D<LO, SC>* getLeftMesh( ) const {
    return leftMesh;
  }

  //! returns the associated right (column) mesh

  inline SurfaceMesh3D<LO, SC>* getRightMesh( ) const {
    return rightMesh;
  }

  //! returns the ansatz type of the inner integral

  inline basisType getAnsatzFunctionType( ) const {
    return ansatzFunctionType;
  }

  //! returns the ansatz type of the outer integral

  inline basisType getTestFunctionType( ) const {
    return testFunctionType;
  }

  //! returns a vector of indices of the elements over which the inner integral is evaluated

  inline vector<LO> const & getInnerElems( ) const {
    return innerElems;
  }

  //! returns a vector of indices of the elements over which the outer integral is evaluated

  inline vector<LO> const & getOuterElems( ) const {
    return outerElems;
  }

  inline LO getNInnerElems( ) const {
    return nInnerElems;
  }

  inline LO getNOuterElems( ) const {
    return nOuterElems;
  }

  //! returns degrees of freedom associated to given outer element
  void getOuterElemDOFs( LO elemIdx, LO *DOFs ) const;

  //! returns degrees of freedom associated to given inner element
  void getInnerElemDOFs( LO elemIdx, LO *DOFs ) const;

  //! returns local degrees of freedom associated to given outer element in a BECluster
  void getClusterOuterElemDOFs( BECluster<LO, SC> const *cluster, LO elemIdx, LO *DOFs ) const;

  //! returns local degrees of freedom associated to given inner element in a BECluster
  void getClusterInnerElemDOFs( BECluster<LO, SC> const *cluster, LO elemIdx, LO *DOFs ) const;

  //! returns a total number of DOFs per inner element

  inline LO getDOFsPerInnerElem( ) const {
    return totalDOFsPerInnerElem;
  }

  //! returns a total number of DOFs per outer element

  inline LO getDOFsPerOuterElem( ) const {
    return totalDOFsPerOuterElem;
  }

  inline LO getOuterDOFs( ) const {
    switch ( testFunctionType ) {
      case p0:
        return mesh->getNElements( );
      case p1:
        return mesh->getNNodes( );
      case p1dis:
        return 3 * mesh->getNElements( );
      default:
        return 0;
    }
  }

  inline LO getInnerDOFs( ) const {
    switch ( ansatzFunctionType ) {
      case p0:
        return mesh->getNElements( );
      case p1:
        return mesh->getNNodes( );
      case p1dis:
        return 3 * mesh->getNElements( );
      default:
        return 0;
    }
  }

  inline LO getClusterOuterDOFs( BECluster<LO, SC> const *cluster ) const {
    switch ( testFunctionType ) {
      case p0:
        return cluster->nelems;
      case p1:
        return cluster->nnodes;
      case p1dis:
        return 3 * cluster->nelems;
      default:
        return 0;
    }
  }

  inline LO getClusterInnerDOFs( BECluster<LO, SC> const *cluster ) const {
    switch ( ansatzFunctionType ) {
      case p0:
        return cluster->nelems;
      case p1:
        return cluster->nnodes;
      case p1dis:
        return 3 * cluster->nelems;
      default:
        return 0;
    }
  }

  //! returns support of idx-th ansatz function

  inline void getInnerSupport(
      LO idx,
      vector< LO > & elems
      ) {

    elems.clear( );

    switch ( ansatzFunctionType ) {
      case p0:
        elems.reserve( 1 );
        elems.push_back( idx );
        break;
      case p1:
        this->mesh->getElements( idx, elems );
        break;
      case p1dis:
        elems.reserve( 1 );
        elems.push_back( idx / 3 );
        break;
    }
  }

  //! returns support of idx-th test function 
  //  (only intersected with cluster for ACA!)

  inline void getInnerSupport(
      LO idx,
      std::vector< LO > & elems,
      const std::vector< LO > & cluster
      ) {

    std::vector< LO > aux;
    elems.clear( );

    switch ( ansatzFunctionType ) {
      case p0:
        //aux.reserve( 1 );
        //aux.push_back( idx );
        elems.reserve( 1 );
        elems.push_back( idx );
        return;
      case p1:
        this->mesh->getElements( idx, aux );
        break;
      case p1dis:
        elems.reserve( 1 );
        elems.push_back( idx / 3 );
        return;
    }


    elems.resize( aux.size( ) );
    std::vector< LO > clusterCopy = cluster;
    std::sort( aux.begin( ), aux.end( ) );
    std::sort( clusterCopy.begin( ), clusterCopy.end( ) );

    auto it =
        std::set_intersection( aux.begin( ), aux.end( ), clusterCopy.begin( ),
        clusterCopy.end( ), elems.begin( ) );

    elems.resize( it - elems.begin( ) );
  }

  //! returns support of idx-th test function

  inline void getOuterSupport(
      LO idx,
      vector< LO > & elems
      ) {
    elems.clear( );

    switch ( testFunctionType ) {
      case p0:
        elems.reserve( 1 );
        elems.push_back( idx );
        break;
      case p1:
        this->mesh->getElements( idx, elems );
        break;
      case p1dis:
        elems.reserve( 1 );
        elems.push_back( idx / 3 );
        break;
    }
  }

  //! returns support of idx-th test function 
  //  (only intersected with cluster for ACA!)

  inline void getOuterSupport(
      LO idx,
      std::vector< LO > & elems,
      const std::vector< LO > & cluster
      ) {

    std::vector< LO > aux;

    elems.clear( );

    switch ( testFunctionType ) {
      case p0:
        //aux.reserve( 1 );
        //aux.push_back( idx );
        elems.reserve( 1 );
        elems.push_back( idx );
        return;
        break;
      case p1:
        this->mesh->getElements( idx, aux );
        break;
      case p1dis:
        elems.reserve( 1 );
        elems.push_back( idx / 3 );
        return;
    }


    elems.resize( aux.size( ) );
    std::vector< LO > clusterCopy = cluster;
    std::sort( aux.begin( ), aux.end( ) );
    std::sort( clusterCopy.begin( ), clusterCopy.end( ) );

    auto it =
        std::set_intersection( aux.begin( ), aux.end( ), clusterCopy.begin( ),
        clusterCopy.end( ), elems.begin( ) );

    elems.resize( it - elems.begin( ) );
  }

  void setAnsatzFunctionType( basisType type ) {
    this->ansatzFunctionType = type;
  }

  void setTestFunctionType( basisType type ) {
    this->testFunctionType = type;
  }



protected:

  //! mesh over which the space is defined
  SurfaceMesh3D<LO, SC>* mesh;

  //! left (row) mesh 
  SurfaceMesh3D<LO, SC>* leftMesh;

  //! right (column) mesh
  SurfaceMesh3D<LO, SC>* rightMesh;

  //! ansatz type of the inner collocation integral 
  basisType ansatzFunctionType;

  //! ansatz type of the outer Galerkin integral 
  basisType testFunctionType;

  //! subset of mesh elements
  vector<LO> innerElems;

  //! subset of mesh elements
  vector<LO> outerElems;

  //! number of inner elements
  LO nInnerElems;

  //! number of outer elements
  LO nOuterElems;

  //! total degrees of freedom per inner element
  LO totalDOFsPerInnerElem;

  //! total degrees of freedom per inner element
  LO totalDOFsPerOuterElem;

  //void defineMappingToDOFs( );

private:

  //! array of mapping of inner element index to global degrees of freedom
  //vector<LO> innerIdx2DOFs;

  //! array of mapping of outer element index to global degrees of freedom
  //vector<LO> outerIdx2DOFs;




};

}

// include .cpp file to overcome linking problems due to templates
#include "BESpace.cpp"


#endif	/* BESPACE_H */

