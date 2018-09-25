/*!
 * @file    FMMKernelLaplace2LayerOMP.cpp
 * @author  Michal Merta 
 * @date    September 25, 2013
 * 
 */
#ifdef FMMKERNELLAPLACE2LAYER_H

#include <omp.h>

// we need this hack because openmp (mainly Intel) is not able to see class member variables
static void* this_class;

namespace bem4i {

template<class LO, class SC>
FMMKernelLaplace2Layer<LO, SC>::FMMKernelLaplace2Layer( ) {
  tempR = nullptr;
  tmpM = nullptr;
  ind2nm = nullptr;

  up = false;
  down = false;
}

template<class LO, class SC>
FMMKernelLaplace2Layer<LO, SC>::FMMKernelLaplace2Layer( const FMMKernelLaplace2Layer& orig ) {
}

template<class LO, class SC>
FMMKernelLaplace2Layer<LO, SC>::FMMKernelLaplace2Layer( FastBESpace<LO, SC>* space, int nMax,
    int startLevel, int quadOrder ) {

  this->beSpace = space;
  this->nMax = nMax;
  this->startLevel = startLevel;
  this->quadOrder = quadOrder;

  this->leftClusterTree = space->getLeftClusterTree( );
  this->rightClusterTree = space->getRightClusterTree( );

  down = false;
  up = false;
  int numThreads;

  // hackity hack
  this_class = this;

  // create the temporary stores
  //get number of running threads and allocate adequate number of temporary arrays
#pragma omp parallel 
  {
#pragma omp critical
    numThreads = omp_get_num_threads( );
  }
  //std::cout << numThreads <<std::endl;                                   
  tempR = new std::complex<SC>*[numThreads];
  tmpM = new std::complex<SC>*[numThreads];


#pragma omp parallel 
  {
#pragma omp critical
    {
      ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR[omp_get_thread_num( )] = new std::complex<SC>[ (int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 )];
      ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tmpM[omp_get_thread_num( )] = new std::complex<SC>[ (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 )];
    }
  }

  // create the mapping from index i to coefficients n,m

  ind2nm = new int*[(int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 )];
  for ( int i = 0; i < (int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 ); i++ ) {
    ind2nm[i] = new int[2];
  }
  for ( int n = 0; n <= nMax; n++ ) {
    for ( int m = 0; m <= n; m++ ) {
      ind2nm[ (int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m ) ][0] = n;
      ind2nm[ (int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m ) ][1] = m;
    }
  }
}

template<class LO, class SC>
FMMKernelLaplace2Layer<LO, SC>::~FMMKernelLaplace2Layer( ) {

#pragma omp parallel
  {
    if ( ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR[omp_get_thread_num( )] )
      delete [] ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR[omp_get_thread_num( )];
    if ( ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tmpM[omp_get_thread_num( )] )
      delete [] ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tmpM[omp_get_thread_num( )];
  }

  if ( tempR )
    delete [] tempR;
  if ( tmpM )
    delete [] tmpM;
  if ( ind2nm ) {

    for ( int i = 0; i < (int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 ); i++ ) {
      if ( ind2nm[i] )
        delete ind2nm[i];
    }
  }

}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::computeR( int nMax, SC *x, std::complex<SC> *Rarray ) {
  int m, n, idx, p = nMax;
  SC r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  std::complex<SC> Rii;

  if ( Rarray == nullptr ) {
    Rarray = new std::complex<SC>[(int) ( p + 2 )*( p + 1 ) / 2];
    //totalSize += (int) (p + 2)*(p + 1) / 2;
  }

  for ( idx = 0; idx < ( p + 2 )*( p + 1 ) / 2; idx++ ) {
    Rarray[idx] = 0;
  }

  for ( m = 0, idx = 0; m <= p; m++ ) {
    for ( n = m; n <= p; n++, idx++ ) {
      if ( m == n && n == 0 ) {
        Rarray[idx] = 1;
        Rii = Rarray[idx];
      } else if ( m == n ) {
        Rarray[idx] = (SC) ( 0.5 / n ) * Rii * std::complex<SC>( x[0], x[1] );
        Rii = Rarray[idx];
      } else if ( n - m == 1 ) {
        Rarray[idx] = ( ( 2 * n - 1 ) * x[2] / ( n + m ) / ( n - m ) ) * Rarray[idx - 1];
      } else
        Rarray[idx] = (SC) ( 1.0 / ( n + m ) / ( n - m ) ) * ( ( 2 * n - 1 ) * x[2] * Rarray[idx - 1] - r2 * Rarray[idx - 2] );
    }
  }
}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::computeS( int nMax, SC *x, std::complex<SC> *Sarray ) {
  int m, n, idx, p = nMax;
  SC r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  std::complex<SC> Sii;

  if ( Sarray == nullptr ) {
    Sarray = new std::complex<SC>[(int) ( p + 2 )*( p + 1 ) / 2];
    //totalSize += (int) (p + 2)*(p + 1) / 2;
  }
  for ( idx = 0; idx < ( p + 2 )*( p + 1 ) / 2; idx++ ) {
    Sarray[idx] = 0;
  }
  for ( m = 0, idx = 0; m <= p; m++ ) {
    for ( n = m; n <= p; n++, idx++ ) {
      if ( m == n && n == 0 ) {
        Sarray[idx] = 1.0 / std::sqrt( r2 );
        Sii = Sarray[idx];
      } else if ( m == n ) {
        Sarray[idx] = ( ( 2 * n - 1 ) / r2 ) * Sii * std::complex<SC>( x[0], x[1] );
        Sii = Sarray[idx];
      } else if ( n - m == 1 )
        Sarray[idx] = ( ( 2 * n - 1 ) * x[2] / r2 ) * Sarray[idx - 1];

      else
        Sarray[idx] = (SC) ( 1.0 / r2 ) * ( ( (SC) ( 2 * n - 1 ) * x[2] ) * Sarray[idx - 1]-( (SC) ( n - 1 + m )*( n - 1 - m ) ) * Sarray[idx - 2] );
    }
  }
}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::upward( ) {

  TreeMember<BECluster<LO, SC>* >* node = this->rightClusterTree->getRootNode( );

  //recursively travel through tree
#pragma omp parallel
  {

#pragma omp single
    doUpwardPass( node );
  }
  this->up = true;

}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::doUpwardPass( TreeMember<BECluster<LO, SC>* >* node ) {

  // traverse the tree from the bottom
  for ( int i = 0; i < ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->rightClusterTree->getNSons( node ); i++ ) {
    node = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->rightClusterTree->getSonNode( node, i );
#pragma omp task
    doUpwardPass( node );
    node = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->rightClusterTree->getParentNode( node );

  }

  // compute multipole moments
#pragma omp taskwait
  if ( ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->rightClusterTree->getNodeLevel( node ) >= ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->startLevel ) {
    computeMoments( node );
  }
}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::downward( ) {



  //recursively travel through tree
#pragma omp parallel 
  {
    TreeMember<BECluster<LO, SC>* >* node = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->leftClusterTree->getRootNode( );
#pragma omp single
    doDownwardPass( node );
  }

  this->down = true;

}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::doDownwardPass( TreeMember<BECluster<LO, SC>* >* node ) {

  // compute coefficients of local expansion
  if ( ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->leftClusterTree->getNodeLevel( node ) >= ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->startLevel ) {
    computeLocExp( node );
  }


  // traverse the tree from the top
  for ( int i = 0; i < ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->leftClusterTree->getNSons( node ); i++ ) {

    node = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->leftClusterTree->getSonNode( node, i );
#pragma omp task
    doDownwardPass( node );

    node = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->leftClusterTree->getParentNode( node );
  }

}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::computeMoments( TreeMember<BECluster<LO, SC>* >* node ) {

  // get pointer to *this* class
  FMMKernelLaplace2Layer<LO, SC>* thisClass = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class );
  int nMax = thisClass->nMax;
  int quadOrder = thisClass->quadOrder;
  std::complex<SC> **tempR = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR;

  int momentArrSize = (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 );

  // variables for gaussian quadrature 
  int quad_size = quadSizes[quadOrder];
  double *quad_points = quadPoints[quadOrder];
  double *quad_weight = quadWeights[quadOrder];

  // some aux. variables
  LO elemIdx;
  LO *vectorIdx = new LO[thisClass->beSpace->getDOFsPerInnerElem( )];
  SC *val = new SC[thisClass->beSpace->getDOFsPerInnerElem( )];
  SC x1[3], x2[3], x3[3]; // element nodes
  SC R_xi[3], yMinYc[3];
  SC normal[3];
  SC half = (SC) 0.5;

  SC intPx, intPy;

  std::complex<SC> Rdx1, Rdx2, Rdx3;
  std::complex<SC> imagUnit( 0.0, 1.0 );

  // current cluster
  BECluster<LO, SC>* cluster = thisClass->rightClusterTree->get( node );
  SC *yc = cluster->centroid;

  // preallocate memory for moments and initialize it
  if ( cluster->moments == nullptr ) {
    cluster->moments = new MultipoleMoments<LO, SC>( );
  }

  if ( cluster->moments->momentsDL == nullptr ) {
    cluster->moments->momentsDL = new std::complex<SC>[momentArrSize];
    //totalSize += momentArrSize;
  }
  for ( int i = 0; i < momentArrSize; i++ ) {
    cluster->moments->momentsDL[i] = 0.0;
  }

  if ( thisClass->rightClusterTree->getNSons( node ) == 0 ) {

    // if leaf, compute moments directly

    for ( LO i = 0; i < cluster->nelems; i++ ) {

      elemIdx = ( *cluster->elems )[i];
      thisClass->beSpace->getInnerElemDOFs( elemIdx, vectorIdx );
      for ( int j = 0; j < thisClass->beSpace->getDOFsPerInnerElem( ); j++ ) {
        val[j] = thisClass->multVector->get( vectorIdx[j] );
      }

      // get Jacobian
      SC Jacob = thisClass->beSpace->getMesh( )->getElemArea( elemIdx );
      thisClass->beSpace->getMesh( )->getNodes( elemIdx, x1, x2, x3 );

      thisClass->beSpace->getMesh( )->getNormal( elemIdx, normal );

      // now perform actual quadrature

      for ( int j = 0; j < quad_size; j++ ) {
        intPx = quad_points[ 2 * j ];
        intPy = quad_points[ 2 * j + 1 ];
        //xi = quad_points[j];
        R_xi[0] = x1[0] + ( x2[0] - x1[0] ) * intPx + ( x3[0] - x1[0] ) * intPy;
        R_xi[1] = x1[1] + ( x2[1] - x1[1] ) * intPx + ( x3[1] - x1[1] ) * intPy;
        R_xi[2] = x1[2] + ( x2[2] - x1[2] ) * intPx + ( x3[2] - x1[2] ) * intPy;
        yMinYc[0] = R_xi[0] - yc[0];
        yMinYc[1] = R_xi[1] - yc[1];
        yMinYc[2] = R_xi[2] - yc[2];

        computeR( nMax, yMinYc, tempR[omp_get_thread_num( )] );

        for ( int n = 0; n <= nMax; n++ ) {
          for ( int m = 0; m <= n; m++ ) {
            if ( m == 0 ) {
              Rdx1 = half * ( getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, -1 ) - getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, 1 ) );
              Rdx2 = imagUnit * half * ( getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, -1 ) + getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, 1 ) );
              Rdx3 = getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, 0 );
            } else {
              Rdx1 = half * ( getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, m - 1 ) - getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, m + 1 ) );
              Rdx2 = imagUnit * half * ( getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, m - 1 ) + getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, m + 1 ) );
              Rdx3 = getRnm( tempR[omp_get_thread_num( )], nMax, n - 1, m );
            }
            if ( m == 0 && n == 0 ) {
              Rdx1 = 0.0;
              Rdx2 = 0.0;
              Rdx3 = 0.0;
            }
            if ( thisClass->beSpace->getAnsatzFunctionType( ) == p0 ) {
              cluster->moments->momentsDL[ (int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m ) ] +=
                  ( val[0] * ( (SC) quad_weight[j] ) * Jacob * ( Rdx1 * normal[0] + Rdx2 * normal[1] + Rdx3 * normal[2] ) );
            } else if ( thisClass->beSpace->getAnsatzFunctionType( ) == p1 ) {

              SC phi1 = 1 - intPx - intPy;
              SC phi2 = intPx;
              SC phi3 = intPy;
              cluster->moments->momentsDL[ (int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m ) ] +=
                  ( val[0] * phi1 * ( (SC) quad_weight[j] ) * Jacob * ( Rdx1 * normal[0] + Rdx2 * normal[1] + Rdx3 * normal[2] ) +
                  val[1] * phi2 * ( (SC) quad_weight[j] ) * Jacob * ( Rdx1 * normal[0] + Rdx2 * normal[1] + Rdx3 * normal[2] ) +
                  val[2] * phi3 * ( (SC) quad_weight[j] ) * Jacob * ( Rdx1 * normal[0] + Rdx2 * normal[1] + Rdx3 * normal[2] )
                  );
            }
          }
        }
      }

    }
  } else {

    // if not leaf, use M2M translation

    for ( int i = 0; i < thisClass->rightClusterTree->getNSons( node ); i++ ) {
      node = thisClass->rightClusterTree->getSonNode( node, i );
      //this->rightClusterTree->get(node);
      SC* ycSon = thisClass->rightClusterTree->get( node )->centroid;
      SC y[3];
      y[0] = ycSon[0] - yc[0];
      y[1] = ycSon[1] - yc[1];
      y[2] = ycSon[2] - yc[2];

      computeR( nMax, y, tempR[omp_get_thread_num( )] );
      std::complex<SC> *sonMDL = thisClass->rightClusterTree->get( node )->moments->momentsDL;
      std::complex<SC> result;

      // outer two loops - setting values of moments depending on n, m
      for ( int n = 0; n <= nMax; n++ ) {
        for ( int m = 0; m <= n; m++ ) {
          // inner two loops - sums in the formula
          for ( int k = 0; k <= n; k++ ) {
            for ( int l = -k; l <= k; l++ ) {
              result = getRnm( tempR[omp_get_thread_num( )], nMax, k, l ) *
                  getMoment( sonMDL, nMax, n - k, m - l );
              cluster->moments->momentsDL[ (int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m ) ] += result;
            }
          }
        }
      }
      node = thisClass->rightClusterTree->getParentNode( node );
    }

  }
}

template<class LO, class SC>
void FMMKernelLaplace2Layer<LO, SC>::computeLocExp( TreeMember<BECluster<LO, SC>* >* node ) {

  // get pointer to *this* class
  FMMKernelLaplace2Layer<LO, SC>* thisClass = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class );
  int nMax = thisClass->nMax;
  //int quadOrder = thisClass->quadOrder;
  int **ind2nm = thisClass->ind2nm;
  std::complex<SC> **tempR = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR;

  int momentArrSize = (int) ( ( nMax ) + 2 )*( nMax + 1 ) / 2;

  // get current cluster
  BECluster<LO, SC>* cluster = thisClass->leftClusterTree->get( node );

  SC *xL = cluster->centroid;
  SC y[3];
  SC *yc;

  std::vector<BECluster<LO, SC>* >* admClusters = cluster->admClusters;

  int nAdmissible;
  if ( admClusters ) {
    nAdmissible = admClusters->size( );
  } else {
    nAdmissible = 0;
  }

  BECluster<LO, SC>* distClust;

  // create the array to store coefficients of local exp.
  // preallocate memory for moments and initialize it

  if ( cluster->moments == nullptr ) {
    cluster->moments = new MultipoleMoments<LO, SC>( );
  }

  if ( cluster->M2L == nullptr && nAdmissible > 0 ) {
    cluster->M2L = new std::complex<SC>*[nAdmissible];
    for ( int i = 0; i < nAdmissible; i++ ) {
      cluster->M2L[i] = nullptr;
    }
  }

  // create the array to store coefficients of local expansion

  if ( cluster->moments->localExpCoefDL == nullptr ) {
    cluster->moments->localExpCoefDL = new std::complex<SC>[momentArrSize];

  }


  for ( int i = 0; i < momentArrSize; i++ ) {
    cluster->moments->localExpCoefDL[i] = std::complex<SC>( 0.0, 0.0 );
  }

  // get the parent cluster of the current cluster
  BECluster<LO, SC> *parent;
  if ( thisClass->leftClusterTree->getNodeLevel( node ) > 0 ) {
    // get BE cluster from the parent node
    parent = thisClass->leftClusterTree->get( thisClass->leftClusterTree->getParentNode( node ) );
  } else {
    parent = nullptr;
  }

  // at first, perform moment-to-local translation 

  for ( int i = 0; i < nAdmissible; i++ ) {

    distClust = ( *admClusters )[i];

    yc = distClust->centroid;

    y[0] = xL[0] - yc[0];
    y[1] = xL[1] - yc[1];
    y[2] = xL[2] - yc[2];

    // if the first run, precompute M2L coefficients
    if ( cluster->M2L[i] == nullptr ) {
      cluster->M2L[i] = new std::complex<SC>[(int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 )];
      computeS( 2 * nMax, y, cluster->M2L[i] );
      for ( int j = 0; j < (int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 ); j++ ) {
        cluster->M2L[i][j] = conj( cluster->M2L[i][j] );
      }
    }

    for ( int j = 0; j < (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 ); j++ ) {
      std::complex<SC> tmp1( 0, 0 );

      int n = ind2nm[j][0];
      int m = ind2nm[j][1];

      SC power = pow( -1.0, n );

      // inner two loops - sums in the formula
      for ( int k = 0; k <= nMax; k++ ) { // !! -n !!
        for ( int l = 0; l <= k; l++ ) {
          tmp1 += ( getSnm( cluster->M2L[i], 2 * nMax, n + k, m + l ) ) *
              getMoment( distClust->moments->momentsDL, nMax, k, l ) * power;
          if ( l != 0 )
            tmp1 += ( getSnm( cluster->M2L[i], 2 * nMax, n + k, m - l ) ) *
            getMoment( distClust->moments->momentsDL, nMax, k, -l ) * power;
        }
      }
      cluster->moments->localExpCoefDL[ j ] += tmp1;
      tmp1 = std::complex<SC>( 0, 0 );
    }
  }


  // now perform local-to-local translation
  if ( thisClass->leftClusterTree->getNodeLevel( node ) != 0 ) {
    if ( parent->moments ) {
      SC *xLp = parent->centroid;
      y[0] = xL[0] - xLp[0];
      y[1] = xL[1] - xLp[1];
      y[2] = xL[2] - xLp[2];

      computeR( nMax, y, tempR[omp_get_thread_num( )] );

      for ( int j = 0; j < (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 ); j++ ) {
        int n = ind2nm[j][0];
        int m = ind2nm[j][1];

        // inner two loops - sums in the formula
        for ( int k = n; k <= nMax; k++ ) {
          for ( int l = 0; l <= k; l++ ) {
            if ( parent->moments->localExpCoefDL ) {
              cluster->moments->localExpCoefDL[ j ] += ( getRnm( tempR[omp_get_thread_num( )], nMax, k - n, l - m ) ) *
                  getMoment( parent->moments->localExpCoefDL, nMax, k, l );
              if ( l != 0 ) {
                cluster->moments->localExpCoefDL[ j ] += ( getRnm( tempR[omp_get_thread_num( )], nMax, k - n, -l - m ) ) *
                    getMoment( parent->moments->localExpCoefDL, nMax, k, -l );
              }
            }
          }
        }
      }
    }

  }
}

template<class LO, class SC>
SC FMMKernelLaplace2Layer<LO, SC>::computeApproximation( BECluster<LO, SC>* cluster, LO element ) {

  // get pointer to *this* class
  FMMKernelLaplace2Layer<LO, SC>* thisClass = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class );
  int nMax = thisClass->nMax;
  int quadOrder = thisClass->quadOrder;
  std::complex<SC> **tempR = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tempR;
  std::complex<SC> **tmpM = ( ( FMMKernelLaplace2Layer<LO, SC>* )this_class )->tmpM;

  int momentArrSize = (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 );
  std::complex<SC> result( 0.0, 0.0 );

  for ( int i = 0; i < momentArrSize; i++ ) {
    tmpM[omp_get_thread_num( )][i] = std::complex<SC>( 0.0, 0.0 );
  }

  // variables for gaussian quadrature 
  int quad_size = quadSizes[quadOrder];
  double *quad_points = quadPoints[quadOrder];
  double *quad_weight = quadWeights[quadOrder];

  // element nodes
  SC x1[3], x2[3], x3[3];
  SC R_xi[3], yMinYc[3];
  SC intPx, intPy;

  // cluster centroid
  SC *yc = cluster->centroid;

  // aux. arr. - nodes in element
  //int nodes[3];
  SC jacob = thisClass->beSpace->getMesh( )->getElemArea( ( *cluster->elems )[element] );

  thisClass->beSpace->getMesh( )->getNodes( ( *cluster->elems )[element], x1, x2, x3 );

  for ( int j = 0; j < quad_size; j++ ) {
    intPx = quad_points[ 2 * j ];
    intPy = quad_points[ 2 * j + 1 ];
    R_xi[0] = x1[0] + ( x2[0] - x1[0] ) * intPx + ( x3[0] - x1[0] ) * intPy;
    R_xi[1] = x1[1] + ( x2[1] - x1[1] ) * intPx + ( x3[1] - x1[1] ) * intPy;
    R_xi[2] = x1[2] + ( x2[2] - x1[2] ) * intPx + ( x3[2] - x1[2] ) * intPy;
    yMinYc[0] = R_xi[0] - yc[0];
    yMinYc[1] = R_xi[1] - yc[1];
    yMinYc[2] = R_xi[2] - yc[2];

    computeR( nMax, yMinYc, tempR[omp_get_thread_num( )] );

    for ( int k = 0; k < momentArrSize; k++ ) {
      tmpM[omp_get_thread_num( )][k] += ( tempR[omp_get_thread_num( )][k] * (SC) quad_weight[j] * jacob );

    }


  }
  SC piRep = 1.0 / ( (SC) 4 * M_PI );
  if ( cluster->moments->localExpCoefDL ) {
    for ( int i = 0; i <= nMax; i++ ) {
      for ( int j = -i; j <= i; j++ ) {
        result += ( getMoment( tmpM[omp_get_thread_num( )], nMax, i, j ) *
            getMoment( cluster->moments->localExpCoefDL, nMax, i, j ) * piRep );

      }
    }
  }
  return real( result );
}
}
#endif
