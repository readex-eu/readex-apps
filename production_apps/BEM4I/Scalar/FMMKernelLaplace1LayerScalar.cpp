/*!
 * @file    FMMKernelLaplace1LayerScalar.cpp
 * @author  Michal Merta 
 * @date    August 19, 2013
 * 
 */

#ifdef FMMKERNELLAPLACE1LAYER_H

namespace bem4i {

template<class LO, class SC>
FMMKernelLaplace1Layer<LO, SC>::FMMKernelLaplace1Layer( ) {
  tempR = nullptr;
  tmpM = nullptr;
  ind2nm = nullptr;

  up = false;
  down = false;
}

template<class LO, class SC>
FMMKernelLaplace1Layer<LO, SC>::FMMKernelLaplace1Layer( const FMMKernelLaplace1Layer& orig ) {
}

template<class LO, class SC>
FMMKernelLaplace1Layer<LO, SC>::FMMKernelLaplace1Layer( FastBESpace<LO, SC>* space, int nMax,
    int startLevel, int quadOrder ) {

  this->beSpace = space;
  this->nMax = nMax;
  this->startLevel = startLevel;
  this->quadOrder = quadOrder;

  this->leftClusterTree = space->getLeftClusterTree( );
  this->rightClusterTree = space->getRightClusterTree( );

  down = false;
  up = false;

  // create the temporary stores
  tempR = new std::complex<SC>[ (int) ( ( 2 * nMax + 2 )*( 2 * nMax + 1 ) / 2 )];
  tmpM = new std::complex<SC>[ (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 )];

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
FMMKernelLaplace1Layer<LO, SC>::~FMMKernelLaplace1Layer( ) {
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
void FMMKernelLaplace1Layer<LO, SC>::computeR( int nMax, SC *x, std::complex<SC> *Rarray ) {
  int m, n, idx, p = nMax;
  SC r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  std::complex<SC> Rii;

  if ( Rarray == nullptr ) {
    Rarray = new std::complex<SC>[(int) ( p + 2 )*( p + 1 ) / 2];
    //totalSize += (int) (p + 2)*(p + 1) / 2;
  }

  for ( idx = 0; idx < ( p + 2 )*( p + 1 ) / 2; idx++ ) {
    Rarray[idx] = 0.0;
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
void FMMKernelLaplace1Layer<LO, SC>::computeS( int nMax, SC *x, std::complex<SC> *Sarray ) {
  int m, n, idx, p = nMax;
  SC r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  std::complex<SC> Sii;

  if ( Sarray == nullptr ) {
    Sarray = new std::complex<SC>[(int) ( p + 2 )*( p + 1 ) / 2];
    //totalSize += (int) (p + 2)*(p + 1) / 2;
  }
  for ( idx = 0; idx < ( p + 2 )*( p + 1 ) / 2; idx++ ) {
    Sarray[idx] = 0.0;
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
void FMMKernelLaplace1Layer<LO, SC>::upward( ) {

  TreeMember<BECluster<LO, SC>* >* node = this->rightClusterTree->getRootNode( );

  //recursively travel through tree
  doUpwardPass( node );

  this->up = true;

}

template<class LO, class SC>
void FMMKernelLaplace1Layer<LO, SC>::doUpwardPass( TreeMember<BECluster<LO, SC>* >* node ) {

  // traverse the tree from the bottom
  for ( int i = 0; i < this->rightClusterTree->getNSons( node ); i++ ) {
    node = this->rightClusterTree->getSonNode( node, i );
    doUpwardPass( node );
    node = this->rightClusterTree->getParentNode( node );
  }

  // compute multipole moments
  if ( this->rightClusterTree->getNodeLevel( node ) >= startLevel ) {

    computeMoments( node );
  }
}

template<class LO, class SC>
void FMMKernelLaplace1Layer<LO, SC>::downward( ) {

  TreeMember<BECluster<LO, SC>* >* node = this->leftClusterTree->getRootNode( );

  //recursively travel through tree
  doDownwardPass( node );

  this->down = true;
}

template<class LO, class SC>
void FMMKernelLaplace1Layer<LO, SC>::doDownwardPass( TreeMember<BECluster<LO, SC>* >* node ) {

  // compute coefficients of local expansion
  if ( this->leftClusterTree->getNodeLevel( node ) >= startLevel ) {
    computeLocExp( node );
  }

  // traverse the tree from the top
  for ( int i = 0; i < this->leftClusterTree->getNSons( node ); i++ ) {

    node = this->leftClusterTree->getSonNode( node, i );
    doDownwardPass( node );
    node = this->leftClusterTree->getParentNode( node );
  }
}

template<class LO, class SC>
void FMMKernelLaplace1Layer<LO, SC>::computeMoments( TreeMember<BECluster<LO, SC>* >* node ) {

  int momentArrSize = (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 );

  // variables for gaussian quadrature 
  int quad_size = quadSizes[quadOrder];
  double *quad_points = quadPoints[quadOrder];
  double *quad_weight = quadWeights[quadOrder];

  // some aux. variables
  LO elemIdx;
  LO *vectorIdx = new LO[this->beSpace->getDOFsPerInnerElem( )];
  SC *val = new SC[this->beSpace->getDOFsPerInnerElem( )];
  SC x1[3], x2[3], x3[3]; // element nodes
  SC R_xi[3], yMinYc[3];

  SC intPx, intPy;
  //Point2D xi;

  BECluster<LO, SC>* cluster = this->rightClusterTree->get( node );
  SC *yc = cluster->centroid;

  // preallocate memory for moments and initialize it
  if ( cluster->moments == nullptr ) {
    cluster->moments = new MultipoleMoments<LO, SC>( );
  }
  if ( cluster->moments->momentsSL == nullptr ) {
    cluster->moments->momentsSL = new std::complex<SC>[momentArrSize];
    //totalSize += momentArrSize;
  }
  for ( int i = 0; i < momentArrSize; i++ ) {
    cluster->moments->momentsSL[i] = std::complex<SC>( 0, 0 );
  }

  if ( this->rightClusterTree->getNSons( node ) == 0 ) {

    // if leaf, compute moments directly

    for ( LO i = 0; i < cluster->nelems; i++ ) {
      elemIdx = ( *cluster->elems )[i];
      this->beSpace->getInnerElemDOFs( elemIdx, vectorIdx );
      for ( int j = 0; j < this->beSpace->getDOFsPerInnerElem( ); j++ ) {
        val[j] = this->multVector->get( vectorIdx[j] );
      }


      // get Jacobian
      SC Jacob = this->beSpace->getMesh( )->getElemArea( elemIdx );
      this->beSpace->getMesh( )->getNodes( elemIdx, x1, x2, x3 );

      // now perform actual quadrature
      for ( int j = 0; j < quad_size; j++ ) {
        intPx = (SC) quad_points[ 2 * j ];
        intPy = (SC) quad_points[ 2 * j + 1 ];
        //xi = quad_points[j];
        R_xi[0] = x1[0] + ( x2[0] - x1[0] ) * intPx + ( x3[0] - x1[0] ) * intPy;
        R_xi[1] = x1[1] + ( x2[1] - x1[1] ) * intPx + ( x3[1] - x1[1] ) * intPy;
        R_xi[2] = x1[2] + ( x2[2] - x1[2] ) * intPx + ( x3[2] - x1[2] ) * intPy;
        yMinYc[0] = R_xi[0] - yc[0];
        yMinYc[1] = R_xi[1] - yc[1];
        yMinYc[2] = R_xi[2] - yc[2];

        computeR( nMax, yMinYc, tempR );

        if ( this->beSpace->getAnsatzFunctionType( ) == p0 ) {
          for ( int k = 0; k < momentArrSize; k++ ) {
            cluster->moments->momentsSL[k] += ( tempR[k] * val[0] * ( (SC) quad_weight[j] ) * Jacob );
          }
        } else if ( this->beSpace->getAnsatzFunctionType( ) == p1 ) {
          std::cout << "Not implemented!" << std::endl;
          /*
          SC phi1 = intPx; 
          SC phi2 = intPy;
          SC phi3 = 1 - intPx - intPy;
          for (int k = 0; k < momentArrSize; k++) {
            cluster->moments->momentsSL[k] += ((tempR[k] * phi1 * val[0] * ((SC) quad_weight[j]) * Jacob) +
                                               (tempR[k] * phi2 * val[1] * ((SC) quad_weight[j]) * Jacob) +
                                               (tempR[k] * phi3 * val[2] * ((SC) quad_weight[j]) * Jacob));
          }
           */
        }
      }

    }
  } else {

    // if not leaf, use M2M translation

    for ( int i = 0; i < this->rightClusterTree->getNSons( node ); i++ ) {
      node = this->rightClusterTree->getSonNode( node, i );
      //this->rightClusterTree->get(node);
      SC* ycSon = this->rightClusterTree->get( node )->centroid;
      SC y[3];
      y[0] = ycSon[0] - yc[0];
      y[1] = ycSon[1] - yc[1];
      y[2] = ycSon[2] - yc[2];

      computeR( nMax, y, tempR );
      std::complex<SC> *sonMSL = this->rightClusterTree->get( node )->moments->momentsSL;
      std::complex<SC> result;

      // outer two loops - setting values of moments depending on n, m
      for ( int n = 0; n <= nMax; n++ ) {
        for ( int m = 0; m <= n; m++ ) {
          // inner two loops - sums in the formula
          for ( int k = 0; k <= n; k++ ) {
            for ( int l = -k; l <= k; l++ ) {

              result = getRnm( tempR, nMax, k, l ) *
                  getMoment( sonMSL, nMax, n - k, m - l );
              cluster->moments->momentsSL[(int) ( ( 2 * nMax + 3 - m ) * m / 2 ) + ( n - m )] += result;
            }
          }
        }
      }
      node = this->rightClusterTree->getParentNode( node );
    }

  }
  delete [] val;
  delete [] vectorIdx;
}

template<class LO, class SC>
void FMMKernelLaplace1Layer<LO, SC>::computeLocExp( TreeMember<BECluster<LO, SC>* >* node ) {

  int momentArrSize = (int) ( ( nMax ) + 2 )*( nMax + 1 ) / 2;

  BECluster<LO, SC>* cluster = this->leftClusterTree->get( node );

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

  if ( cluster->M2L == nullptr ) {
    cluster->M2L = new std::complex<SC>*[nAdmissible];
    for ( int i = 0; i < nAdmissible; i++ ) {
      cluster->M2L[i] = nullptr;
    }
  }

  // create the array to store coefficients of local expansion
  if ( cluster->moments->localExpCoefSL == nullptr ) {
    cluster->moments->localExpCoefSL = new std::complex<SC>[momentArrSize];
  }
  for ( int i = 0; i < momentArrSize; i++ ) {
    cluster->moments->localExpCoefSL[i] = std::complex<SC>( 0.0, 0.0 );
  }

  // get the parent cluster of the current cluster
  BECluster<LO, SC> *parent;
  if ( this->leftClusterTree->getNodeLevel( node ) > 0 ) {
    // get BE cluster from the parent node
    parent = this->leftClusterTree->get( this->leftClusterTree->getParentNode( node ) );
  } else {
    parent = nullptr;
  }

  // at first, perform moment-to-local translation 
  //#pragma omp parallel for
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
          if ( l != 0 ) {
            tmp1 += ( ( getSnm( cluster->M2L[i], 2 * nMax, n + k, m + l ) ) *
                getMoment( distClust->moments->momentsSL, nMax, k, l ) +
                ( getSnm( cluster->M2L[i], 2 * nMax, n + k, m - l ) ) *
                getMoment( distClust->moments->momentsSL, nMax, k, -l ) ) * power;
          } else {
            tmp1 += ( getSnm( cluster->M2L[i], 2 * nMax, n + k, m + l ) ) *
                getMoment( distClust->moments->momentsSL, nMax, k, l ) * power;
          }
        }
      }
      cluster->moments->localExpCoefSL[ j ] += tmp1;
      tmp1 = std::complex<SC>( 0, 0 );
    }
  }


  // now perform local-to-local translation
  if ( this->leftClusterTree->getNodeLevel( node ) != 0 ) {
    if ( parent->moments ) {
      SC *xLp = parent->centroid;
      y[0] = xL[0] - xLp[0];
      y[1] = xL[1] - xLp[1];
      y[2] = xL[2] - xLp[2];

      computeR( nMax, y, tempR );

      for ( int j = 0; j < (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 ); j++ ) {
        int n = ind2nm[j][0];
        int m = ind2nm[j][1];

        // inner two loops - sums in the formula
        for ( int k = n; k <= nMax; k++ ) {
          for ( int l = 0; l <= k; l++ ) {
            if ( parent->moments->localExpCoefSL ) {
              cluster->moments->localExpCoefSL[ j ] += ( getRnm( tempR, nMax, k - n, l - m ) ) *
                  getMoment( parent->moments->localExpCoefSL, nMax, k, l );
              if ( l != 0 ) {
                cluster->moments->localExpCoefSL[ j ] += ( getRnm( tempR, nMax, k - n, -l - m ) ) *
                    getMoment( parent->moments->localExpCoefSL, nMax, k, -l );

              }
            }
          }
        }
      }
    }

  }
}

template<class LO, class SC>
SC FMMKernelLaplace1Layer<LO, SC>::computeApproximation( BECluster<LO, SC>* cluster, LO element ) {

  int momentArrSize = (int) ( ( nMax + 2 )*( nMax + 1 ) / 2 );
  std::complex<SC> result( 0.0, 0.0 );

  for ( int i = 0; i < momentArrSize; i++ ) {
    tmpM[i] = std::complex<SC>( 0.0, 0.0 );
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
  SC jacob = this->beSpace->getMesh( )->getElemArea( ( *cluster->elems )[element] );

  this->beSpace->getMesh( )->getNodes( ( *cluster->elems )[element], x1, x2, x3 );

  for ( int j = 0; j < quad_size; j++ ) {
    intPx = quad_points[ 2 * j ];
    intPy = quad_points[ 2 * j + 1 ];
    R_xi[0] = x1[0] + ( x2[0] - x1[0] ) * intPx + ( x3[0] - x1[0] ) * intPy;
    R_xi[1] = x1[1] + ( x2[1] - x1[1] ) * intPx + ( x3[1] - x1[1] ) * intPy;
    R_xi[2] = x1[2] + ( x2[2] - x1[2] ) * intPx + ( x3[2] - x1[2] ) * intPy;
    yMinYc[0] = R_xi[0] - yc[0];
    yMinYc[1] = R_xi[1] - yc[1];
    yMinYc[2] = R_xi[2] - yc[2];

    computeR( nMax, yMinYc, tempR );

    for ( int k = 0; k < momentArrSize; k++ ) {
      tmpM[k] += ( tempR[k] * (SC) quad_weight[j] * jacob );
    }


  }
  SC piRep = 1.0 / ( (SC) 4.0 * M_PI );
  if ( cluster->moments->localExpCoefSL ) {
    for ( int i = 0; i <= nMax; i++ ) {
      for ( int j = -i; j <= i; j++ ) {
        //std::cout << getMoment(tmpM, nMax, i, j) << ", " << getMoment(cluster->moments->localExpCoefSL, nMax, i, j) << std::endl;
        result += ( getMoment( tmpM, nMax, i, j ) *
            getMoment( cluster->moments->localExpCoefSL, nMax, i, j ) * piRep );
      }
    }
  }
  return real( result );
}
}
#endif
