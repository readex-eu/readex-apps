/*!
 * @file    Vector.cpp
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    July 7, 2013
 * 
 */

#ifdef VECTOR_H2

namespace bem4i {

template<class LO, class SC>
Vector<LO, SC>::Vector( ) {
  this->data = nullptr;
  this->length = 0;
  this->deleteData = true;
}

template<class LO, class SC>
Vector<LO, SC>::Vector(
  const Vector & orig
  ) {

  this->length = orig.length;
  this->data = new SC[ this->length ];
  if ( orig.data ) {
    memcpy( this->data, orig.data, this->length * sizeof ( SC ) );

  } else {
    this->data = nullptr;
  }
  this->deleteData = true;
}

template<class LO, class SC>
Vector<LO, SC>::Vector(
  LO length,
  bool zeroOut
  ) {

  this->length = length;
  this->data = new SC[ length ];
  if ( zeroOut ) this->setAll( 0.0 );
  this->deleteData = true;
}

template<class LO, class SC>
Vector<LO, SC>::Vector(
  LO length,
  SC * data,
  bool copy
  ) {

  this->length = length;
  if ( copy ) {
    this->data = new SC[ length ];
    memcpy( this->data, data, length * sizeof ( SC ) );
    this->deleteData = true;
  } else {
    this->data = data;
    this->deleteData = false;
  }
}

template<class LO, class SC>
Vector<LO, SC>::~Vector( ) {
  if ( this->deleteData && this->data ) {
    delete [ ] this->data;
  }
}

// GENERIC IMPLEMENTATION //

template<class LO, class SC>
void Vector<LO, SC>::copy(
  Vector<LO, SC> & copy
  ) const {

  copy.length = this->length;
  if ( copy.data && copy.deleteData ) {
    delete [ ] copy.data;
  }

  if ( this->data ) {
    copy.data = new SC[ this->length ];
    memcpy( copy.data, this->data, this->length * sizeof ( SC ) );
    copy.deleteData = true;
  } else {
    copy.data = nullptr;
    copy.deleteData = false;
  }


}

template<class LO, class SC>
void Vector<LO, SC>::copyToConjugate(
  Vector<LO, SC> & copy
  ) const {

  copy.length = this->length;
  if ( copy.data && copy.deleteData ) {
    delete [ ] copy.data;
  }

  if ( this->data ) {
    copy.data = new SC[ this->length ];
    SC val;

    for ( LO i = 0; i < this->length; ++i ) {
      val = (SC) std::conj( this->data[ i ] );
      copy.data[ i ] = val;
    }
    copy.deleteData = true;
  } else {
    copy.data = nullptr;
    copy.deleteData = false;
  }


}

template<class LO, class SC>
void Vector<LO, SC>::resize(
  LO length,
  bool zero
  ) {

  this->length = length;

  if ( this->data && this->deleteData ) {
    delete [ ] this->data;
  }

  this->data = new SC[ length ];
  this->deleteData = true;

  if ( zero ) {
    this->setAll( 0.0 );
  }
}

template<class LO, class SC>
void Vector<LO, SC>::real( Vector< LO, SCVT > & real ) const {
  real.resize( this->length );
  for ( LO i = 0; i < this->length; i++ ) {
    real.data[ i ] = this->data[ i ];
  }
}

template<class LO, class SC>
void Vector<LO, SC>::imag( Vector< LO, SCVT > & imag ) const {
  imag.resize( this->length );
  imag.setAll( 0.0 );
}

template<class LO, class SC>
typename bem4i::Vector<LO, SC>::SCVT Vector<LO, SC>::norm2( ) const {
  LO one = 1;
  return dnrm2_( &length, data, &one );
}

template<class LO, class SC>
LO Vector<LO, SC>::findAbsMax( ) const {
  SCVT max = 0.0;
  LO maxInd = 0;
  SCVT data;

  for ( LO i = 0; i < this->length; ++i ) {
    data = std::abs( this->data[ i ] );
    if ( data > max ) {
      max = data;
      maxInd = i;
    }
  }

  return maxInd;
}

template<class LO, class SC>
SC Vector<LO, SC>::dot( Vector<LO, SC> const &a ) const {
  LO one = 1;
  return ddot_( &length, &data, &one, a.data, &one );
}

template<class LO, class SC>
void Vector<LO, SC>::cross( Vector<LO, SC> const &b, Vector<LO, SC> &ret, LO nB ) const {
  for ( LO i = 0; i < nB; i++ ) {
    ret.data[ 3 * i + 0 ] = this->data[1] * b.data[3 * i + 2] - this->data[2] * b.data[3 * i + 1];
    ret.data[ 3 * i + 1 ] = this->data[2] * b.data[3 * i + 0] - this->data[0] * b.data[3 * i + 2];
    ret.data[ 3 * i + 2 ] = this->data[0] * b.data[3 * i + 1] - this->data[1] * b.data[3 * i + 0];
  }
}

template<class LO, class SC>
void Vector<LO, SC>::scale( SC alpha ) {
  int incx = 1;
  dscal_( &length, &alpha, data, &incx );
}

template<class LO, class SC>
void Vector<LO, SC>::conjugate( ) {

  SC val;

#pragma omp parallel for private( val )
  for ( LO i = 0; i < this->length; ++i ) {
    val = (SC) std::conj( this->data[ i ] );
    this->data[ i ] = val;
  }
}

//template<class LO, class SC>
//void Vector<LO, SC>::errorRelative( Vector<LO,SC> const & sol ) {
//  SC solData;
//  for( LO i = 0; i < this->length; i++ ){
//    solData = sol.data[ i ];
//    this->data[ i ] = std::abs( this->data[ i ] - solData ) / std::abs( solData );
//  }
//}

template<class LO, class SC>
void Vector<LO, SC>::add( Vector<LO, SC> const &x, SC alpha ) {
  LO inc = 1;
  daxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<class LO, class SC>
void Vector<LO, SC>::add( Vector<LO, SC> const &x, Vector<LO, SC> &ret, SC alpha ) {
  for ( LO i = 0; i < length; i++ ) {
    ret.data[ i ] = data[ i ] + alpha * x.data[ i ];
  }
}

template<>
void Vector< int, std::complex< double > >::addAtomic(
  int i,
  std::complex< double > val
  ) {

#pragma omp atomic update
  reinterpret_cast<double *> ( this->data )[ 2 * i ] += std::real( val );
#pragma omp atomic update
  reinterpret_cast<double *> ( this->data )[ 2 * i + 1 ] += std::imag( val );
}

template<>
void Vector< long, std::complex< double > >::addAtomic(
  long i,
  std::complex< double > val
  ) {
#pragma omp atomic update
  reinterpret_cast<double *> ( this->data )[ 2 * i ] += std::real( val );
#pragma omp atomic update
  reinterpret_cast<double *> ( this->data )[ 2 * i + 1 ] += std::imag( val );
}

template<>
void Vector< int, std::complex< float > >::addAtomic(
  int i,
  std::complex< float > val
  ) {
#pragma omp atomic update
  reinterpret_cast<float *> ( this->data )[ 2 * i ] += std::real( val );
#pragma omp atomic update
  reinterpret_cast<float *> ( this->data )[ 2 * i + 1 ] += std::imag( val );
}

template<>
void Vector< long, std::complex< float > >::addAtomic(
  long i,
  std::complex< float > val
  ) {
#pragma omp atomic update
  reinterpret_cast<float *> ( this->data )[ 2 * i ] += std::real( val );
#pragma omp atomic update
  reinterpret_cast<float *> ( this->data )[ 2 * i + 1 ] += std::imag( val );
}

template<class LO, class SC>
void Vector<LO, SC>::print( std::ostream &stream ) const {
  stream << "\nVector\n";
  stream << "Length: " << length << "\n";
  for ( int i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";

}

template<class LO, class SC>
void Vector<LO, SC>::load(
  const std::string& fileName
  ) {

  std::ifstream file( fileName.c_str( ) );

  if ( !file.is_open( ) ) {
    std::cout << "File could not be opened!" << std::endl;
    return;
  }

  file >> this->length;

  if ( this->data ) delete this->data;
  this->data = new SC[ this->length ];
  this->deleteData = true;

  for ( int k = 0; k < this->length; ++k ) {
    file >> this->data[ k ];
  }

  file.close( );
}

#ifdef BLAS_INT
// <INT, DOUBLE> PRECISION SPECIALIZATION //

template<>
double Vector<int, double>::norm2( ) const {
  int one = 1;
  return dnrm2_( &length, data, &one );
}

template<>
double Vector<int, double>::dot( Vector<int, double> const &a ) const {
  int one = 1;
  return ddot_( &length, data, &one, a.data, &one );
}

template<>
void Vector<int, double>::scale( double alpha ) {
  int incx = 1;
  dscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<int, double>::add( Vector<int, double> const &x, double alpha ) {
  int inc = 1;
  daxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<int, double>::print( std::ostream &stream ) const {
  stream << "\nVector<int, double>\n";
  stream << "Length: " << length << "\n";
  for ( int i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";

}

// <INT, SINGLE> PRECISION SPECIALIZATION //

template<>
float Vector<int, float>::norm2( ) const {
  int one = 1;
  return snrm2_( &length, data, &one );
}

template<>
float Vector<int, float>::dot( Vector<int, float> const &a ) const {
  int one = 1;
  return sdot_( &length, data, &one, a.data, &one );
}

template<>
void Vector<int, float>::scale( float alpha ) {
  int incx = 1;
  sscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<int, float>::add( Vector<int, float> const &x, float alpha ) {
  int inc = 1;
  saxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<int, float>::print( std::ostream &stream ) const {
  stream << "\nVector<int, float>\n";
  stream << "Length: " << length << "\n";
  for ( int i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";

}

// <INT, COMPLEX SINGLE> PRECISION SPECIALIZATION //

template<>
float Vector<int, std::complex<float> >::norm2( ) const {
  int one = 1;
  float norm = scnrm2_( &length, data, &one );
  return norm;
}

template<>
std::complex<float> Vector<int, std::complex<float> >::dot( Vector<int, std::complex<float> > const &a ) const {
  //  int one = 1;
  //  return cdotc_( &length, data, &one, a.data, &one );
  std::complex<float> ret( 0.0, 0.0 );

  for ( int i = 0; i < length; ++i ) {
    ret += std::conj( data[ i ] ) * a.data[ i ];
  }
  return ret;
}

template<>
void Vector<int, std::complex<float> >::scale( std::complex<float> alpha ) {
  int incx = 1;
  cscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<int, std::complex<float> >::add( Vector<int, std::complex<float> > const &x, std::complex<float> alpha ) {
  int inc = 1;
  caxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<int, std::complex<float> >::print( std::ostream &stream ) const {
  std::ios::fmtflags f( std::cout.flags( ) );

  stream << "\nVector<int, std::complex<float> >\n";
  stream << "Length: " << length << "\n";
  for ( int i = 0; i < length; i++ ) {
    stream << std::real( get( i ) ) << std::showpos << std::imag( get( i ) ) <<
      "i ";
    std::cout.flags( f );
  }
  stream << "\n";
}

template<>
void Vector<int, std::complex<float> >::real( Vector<int, float> & real ) const {
  real.resize( this->length );
  for ( int i = 0; i < this->length; i++ ) {
    real.set( i, this->data[ i ].real( ) );
  }
}

template<>
void Vector<int, std::complex<float> >::imag( Vector<int, float > & imag ) const {
  imag.resize( this->length );
  for ( int i = 0; i < this->length; i++ ) {
    imag.set( i, this->data[ i ].imag( ) );
  }
}

// <INT, COMPLEX DOUBLE> PRECISION SPECIALIZATION //

template<>
double Vector<int, std::complex<double> >::norm2( ) const {
  int one = 1;
  double norm = dznrm2_( &length, data, &one );
  return norm;
}

template<>
std::complex<double> Vector<int, std::complex<double> >::dot( Vector<int, std::complex<double> > const &a ) const {
  //int one = 1;
  //return zdotc_( &length, data, &one, a.data, &one );
  std::complex<double> ret( 0.0, 0.0 );

  for ( int i = 0; i < length; ++i ) {
    ret += std::conj( data[ i ] ) * a.data[ i ];
  }
  return ret;
}

template<>
void Vector<int, std::complex<double> >::scale( std::complex<double> alpha ) {
  int incx = 1;
  zscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<int, std::complex<double> >::add( Vector<int, std::complex<double> > const &x, std::complex<double> alpha ) {
  int inc = 1;
  zaxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<int, std::complex<double> >::print( std::ostream &stream ) const {
  std::ios::fmtflags f( std::cout.flags( ) );

  stream << "\nVector<int, std::complex<double> >\n";
  stream << "Length: " << length << "\n";
  for ( int i = 0; i < length; i++ ) {
    stream << std::real( get( i ) ) << std::showpos << std::imag( get( i ) ) <<
      "i ";
    std::cout.flags( f );
  }
  stream << "\n";
}

template<>
void Vector<int, std::complex<double> >::real( Vector<int, double> & real ) const {
  real.resize( this->length );
  for ( int i = 0; i < this->length; i++ ) {
    real.set( i, this->data[ i ].real( ) );
  }
}

template<>
void Vector<int, std::complex<double> >::imag( Vector<int, double> & imag ) const {
  imag.resize( this->length );
  for ( int i = 0; i < this->length; i++ ) {
    imag.set( i, this->data[ i ].imag( ) );
  }
}

#elif defined(BLAS_LONG)

// <LONG, DOUBLE> PRECISION SPECIALIZATION //

template<>
double Vector<long, double>::norm2( ) const {
  long one = 1;
  return dnrm2_( &length, data, &one );
}

template<>
double Vector<long, double>::dot( Vector<long, double> const &a ) const {
  long one = 1;
  return ddot_( &length, data, &one, a.data, &one );
}

template<>
void Vector<long, double>::scale( double alpha ) {
  long incx = 1;
  dscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<long, double>::add( Vector<long, double> const &x, double alpha ) {
  long inc = 1;
  daxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<long, double>::print(
  std::ostream &stream
  ) const {
  stream << "\nVector<long, double>\n";
  stream << "Length: " << length << "\n";
  for ( long i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";

}

// <LONG, SINGLE> PRECISION SPECIALIZATION //

template<>
float Vector<long, float>::norm2( ) const {
  long one = 1;
  return snrm2_( &length, data, &one );
}

template<>
float Vector<long, float>::dot( Vector<long, float> const &a ) const {
  long one = 1;
  return sdot_( &length, data, &one, a.data, &one );
}

template<>
void Vector<long, float>::scale( float alpha ) {
  long incx = 1;
  sscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<long, float>::add( Vector<long, float> const &x, float alpha ) {
  long inc = 1;
  saxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<long, float>::print(
  std::ostream &stream
  ) const {
  stream << "\nVector<long, float>\n";
  stream << "Length: " << length << "\n";
  for ( long i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";

}


// <LONG, COMPLEX SINGLE> PRECISION SPECIALIZATION //

template<>
float Vector<long, std::complex<float> >::norm2( ) const {
  long one = 1;
  float norm = scnrm2_( &length, data, &one );
  return norm;
}

template<>
std::complex<float> Vector<long, std::complex<float> >::dot( Vector<long, std::complex<float> > const &a ) const {
  long one = 1;
  return cdotc_( &length, data, &one, a.data, &one );
}

template<>
void Vector<long, std::complex<float> >::scale( std::complex<float> alpha ) {
  long incx = 1;
  cscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<long, std::complex<float> >::add(
  Vector<long, std::complex<float> > const &x,
  std::complex<float> alpha
  ) {
  long inc = 1;
  caxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<long, std::complex<float> >::print(
  std::ostream &stream
  ) const {
  stream << "\nVector<long, std::complex<float> >\n";
  stream << "Length: " << length << "\n";
  for ( long i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";
}

template<>
void Vector<long, std::complex<float> >::real( Vector<long, float> & real ) const {
  real.resize( this->length );
  for ( long i = 0; i < this->length; i++ ) {
    real.set( i, this->data[ i ].real( ) );
  }
}

template<>
void Vector<long, std::complex<float> >::imag( Vector<long, float> & imag ) const {
  imag.resize( this->length );
  for ( long i = 0; i < this->length; i++ ) {
    imag.set( i, this->data[ i ].imag( ) );
  }
}

// <LONG, COMPLEX DOUBLE> PRECISION SPECIALIZATION //

template<>
double Vector<long, std::complex<double> >::norm2( ) const {
  long one = 1;
  double norm = dznrm2_( &length, data, &one );
  return norm;
}

template<>
void Vector<long, std::complex<double> >::scale( std::complex<double> alpha ) {
  long incx = 1;
  zscal_( &length, &alpha, data, &incx );
}

template<>
void Vector<long, std::complex<double> >::add( Vector<long, std::complex<double> > const &x, std::complex<double> alpha ) {
  long inc = 1;
  zaxpy_( &length, &alpha, x.data, &inc, data, &inc );
}

template<>
void Vector<long, std::complex<double> >::print(
  std::ostream &stream
  ) const {
  stream << "\nVector<long, std::complex<double> >\n";
  stream << "Length: " << length << "\n";
  for ( long i = 0; i < length; i++ ) {
    stream << get( i ) << " ";
  }
  stream << "\n";
}

template<>
std::complex<double> Vector<long, std::complex<double> >::dot( Vector<long, std::complex<double> > const &a ) const {

  std::complex<double> result = 0.0;
  std::complex<double> b;
  for ( long i = 0; i < length; i++ ) {
    b = std::conj( data[i] );
    result += b * a.data[i];
  }
  return result;

  // for some reason the BLAS dot product does not work with long
  //long one = 1;
  //return zdotc_( &length, data, &one, a.data, &one );
}

template<>
void Vector<long, std::complex<double> >::real( Vector<long, double> & real ) const {
  real.resize( this->length );
  for ( long i = 0; i < this->length; i++ ) {
    real.set( i, this->data[ i ].real( ) );
  }
}

template<>
void Vector<long, std::complex<double> >::imag( Vector<long, double> & imag ) const {
  imag.resize( this->length );
  for ( long i = 0; i < this->length; i++ ) {
    imag.set( i, this->data[ i ].imag( ) );
  }
}

#endif

}

#endif
