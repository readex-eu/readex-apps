/*!
 * @file    Vector.h
 * @author  Michal Merta 
 * @author  Jan Zapletal 
 * @date    July 7, 2013
 * @brief   Header file for class Vector
 * 
 */

#ifndef VECTOR_H2
#define VECTOR_H2

#include <cstring>
#include <complex>
#include <iostream>
#include <fstream>
#include <type_traits>
#include "Macros.h"
#include "Offloadable.h"

#include "BLAS_LAPACK_wrapper.h"

namespace bem4i {

  // forward declaration
  template<class LO2, class SC2>
  class FullMatrix;
  template<class LO3, class SC3>
  class IdentityOperator;

  /*! 
   * Class representing a full vector
   * 
   * Vector operations are performed by calling BLAS/LAPACK routines.
   * 
   */
  template<class LO, class SC>
  class Vector : public Offloadable {
    typedef typename GetType<LO, SC>::SCVT SCVT;

    //template<class LO2, class SCVT>
    friend class FullMatrix<LO, SC>;
    //template<class LO3, class SC3> 
    friend class IdentityOperator<LO, SC>;

  public:

    //! default constructor
    Vector( );

    //! copy constructor
    Vector(
      const Vector& orig
      );

    //! constructor allocating a vector of length length
    Vector(
      LO length,
      bool zeroOut = true
      );

    //! constructor taking a user preallocated array to store vector entries
    Vector(
      LO length,
      SC * data,
      bool copy = false
      );

    //! destructor
    virtual ~Vector( );

    //! returns a length of a vector

    inline LO getLength( ) const {
      return length;
    }

    void copy(
      Vector<LO, SC> & copy
      ) const;

    void copyToConjugate(
      Vector<LO, SC> & copy
      ) const;

    //! returns an element idx of a vector

    inline SC get(
      LO idx
      ) const {
      return data[idx];
    }

    //! returns pointer to data

    inline SC* getData( ) const {
      return data;
    }

    inline void setDeleteData(
      bool deleteData
      ) {
      this->deleteData = deleteData;
    }

    void setData(
      LO length,
      SC * data,
      bool copy = false
      ) {

      this->length = length;
      if ( copy ) {
        if ( this->data && this->deleteData ) {
          delete [ ] this->data;

        }
        this->data = new SC[ length ];
        memcpy( this->data, data, length * sizeof ( SC ) );
        deleteData = true;

      } else {
        this->data = data;
        deleteData = false;
      }
    }

    void append(
      Vector< LO, SC > & x
      ) {

      SC * tmp = new SC[ this->length + x.length ];
      memcpy( tmp, this->data, this->length * sizeof ( SC ) );
      memcpy( tmp + this->length, x.data, x.length * sizeof ( SC ) );

      if ( this->data && this->deleteData ) {
        delete [ ] this->data;
      }
      this->data = tmp;
      this->length += x.length;
      this->deleteData = true;
    }

    //! returns real part of the vector
    void real(
      Vector< LO, SCVT > & real
      ) const;

    //! returns imag part of the vector
    void imag(
      Vector< LO, SCVT > & imag
      ) const;

    //! sets an element idx of a vector to value val

    inline void set(
      LO idx,
      SC val
      ) {
      data[idx] = val;
    }

    //! sets all elements of a vector to a value val

    inline void setAll( SC val ) {
      SC *p = data, *last = data + length;
      while ( p != last ) *( p++ ) = val;
    }

    inline SC sum( ) {
      SC sum = 0.0;
      for ( LO i = 0; i < this->length; i++ ) {
        sum += this->data[ i ];
      }
      return sum;
    }

    //! scales all values of a vector by alpha
    void scale( SC alpha );

    //! performs in-place conjugation
    void conjugate( );

    //! elementwise relative error
    //void errorRelative( Vector<LO,SC> const & sol );

    void resize( LO length, bool zero = true );

    //! computes 2-norm of a vector
    SCVT norm2( ) const;

    LO findAbsMax( ) const;

    //! computes a dot product of vector this and vector a
    SC dot(
      Vector<LO, SC> const & a
      ) const;

    /*!
     * @brief computes cross product this \times b
     * 
     * Computes a sum this = this + alpha*x
     * @param b second operand, size = nB * sizeof this 
     * @param ret preallocated return vector
     * @param nB optional argument specifying number of inputs stored in b
     */
    void cross(
      Vector<LO, SC> const & b,
      Vector<LO, SC> & ret,
      LO nB = 1
      ) const;

    /*!
     * @brief Vector-vector addition
     * 
     * Computes a sum this = this + alpha*x
     * @param x 
     * @param alpha
     */
    void add(
      Vector<LO, SC> const & x,
      SC alpha = 1.0
      );

    /*!
     * @brief Vector-vector addition
     * 
     * Computes a sum ret = this + alpha*x
     * @param x not overwritten
     * @param ret preallocated vector
     * @param alpha
     */
    void add(
      Vector<LO, SC> const &x,
      Vector<LO, SC> &ret,
      SC alpha = 1.0
      );

    void add( LO i, SC val ) {
      this->data[i] += val;
    }

    void addAtomic(
      LO i,
      SC val
      ) {

#pragma omp atomic update
      this->data[ i ] += val;
    }

    void add( SC val ) {
      for ( LO i = 0; i < this->length; ++i ) {
        this->data[i] += val;
      }
    }

    //! prints the vector
    void print( std::ostream &stream = std::cout ) const;

    void load(
      const std::string & fileName
      );

    void xferToMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

#pragma omp target enter data device( device ) \
map( to : data[ 0 : mult * length ] )

#endif
    }

    void xferToHost(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

#pragma omp target update device( device ) \
from( data[ 0 : mult * length ] )

#endif
    }

    void updateMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );
      int mult = 1;
      if ( !std::is_same< SC, SCVT >::value ) mult = 2;

#pragma omp target update device( device ) \
to( data[ 0 : mult * length ] )

#endif
    }

    void deleteMIC(
      int device = 0
      ) const {

#if N_MIC > 0

      const SCVT * data = reinterpret_cast < SCVT * > ( this->data );

#pragma omp target exit data device( device ) \
map( delete : data )

#endif
    }

  private:

    //! number of vector elements
    LO length;

    //! array of vector data
    SC *data;

    //! determines whether to destructor deletes data
    bool deleteData;

  };

}
// include .cpp file to overcome linking problems due to templates
#include "Vector.cpp"


#endif /* VECTOR_H2 */

