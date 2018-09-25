/*!
 * @file    IterativeSolver.h
 * @author  Michal Merta 
 * @date    April 16, 2014
 * @brief   Servant IterativeSolver
 * 
 */


#ifndef ITERATIVESOLVER_H
#define	ITERATIVESOLVER_H

#include "Vector.h"
#include "LeftPreconditioner.h"
#include "FullMatrix.h"

#include <vector>
#include <algorithm>
#include <typeinfo>
#include <mpi.h>

namespace bem4i {

//!  global variable stores info about level of an iterative solver
int solver_level(0);

/*! 
 * Abstract class representing arbitrary linear operator
 * 
 */
template<class LO, class SC>
class IterativeSolver {
public:

  typedef typename GetType<LO, SC>::SCVT SCVT;

  //todo: cannot work for complex matrices!

  /*!
   * @brief CG solver for symmetric, positive definite matrices
   * 
   * @param rhs           right-hand side vector
   * @param sol           solution vector
   * @param prec[in,out]  on input: solver relative precision, on output:
   *                      real obtained precision
   * @param maxIt[in,out] on input: max. number of iterations, on output:
   *                      real number of iterations
   */
  static bool CGSolve(
      LinearOperator<LO, SC> & A,
      const Vector<LO, SC> & rhs,
      Vector<LO, SC> & sol,
      SCVT prec,
      LO maxIt,
      LinearOperator<LO, SC> * M = nullptr,
      const std::string & msg = ""
      ) {

    ++solver_level;

    bool isMsg = ( msg.compare( "" ) != 0 );
    std::string msgLvl;

    if( isMsg ){
#ifdef VERBOSE
      for (LO i = 0; i < solver_level; ++i)
        msgLvl.append( "   " );
      msgLvl.append( msg );
      std::cout << msgLvl << ": ";
#endif
    }

#ifdef VERBOSE
    if ( M ) std::cout << "P";
    std::cout << "CG started" << std::endl;
#endif

    //LO n = A.getNRows( );
    LO n = rhs.getLength( );
    Vector<LO, SC> r( n );
    Vector<LO, SC> s( n );
    Vector<LO, SC> v( n );

    SC rho0, rho, rhoNext;
    SC sigma, alpha, beta;
    LO iter = 0;

    A.apply( sol, r );
    r.add( rhs, -1.0 );
    if ( M ) {
      M->apply( r, v );
    } else {
      r.copy( v );
    }
    Vector<LO, SC> p( v );
    rho = rho0 = r.dot( v );
    SCVT res0 = r.norm2( ); // only for stopping

    if ( rhs.norm2( ) == 0.0 ) {
      sol.setAll( 0.0 );
      std::cout << "CG iterations: " << iter << ". " <<
          "Relative error: " << std::sqrt( std::abs( rho0 ) ) << 
          "." << std::endl;
      --solver_level;
      return true;
    }
    if ( std::sqrt( std::abs( rho0 ) ) / rhs.norm2( ) < EPS ) {
#ifdef VERBOSE
      if( isMsg ){
        std::cout << msgLvl << ": ";
      }
      if ( M ) std::cout << "P";
      std::cout << "CG iterations: " << iter << ". " <<
          "Relative error: " << std::sqrt( std::abs( rho0 ) ) / rhs.norm2( ) << 
          "." << std::endl;
#endif
      --solver_level;
      return true;
    }

    while ( std::abs( rho / rho0 ) > prec * prec && iter < maxIt ) {
      //while ( r.norm2( ) / res0 > prec && iter < maxIt ) {
      A.apply( p, s );
      sigma = p.dot( s );
      alpha = rho / sigma;
      sol.add( p, -alpha );
      r.add( s, -alpha );
      if ( M ) {
        //r.print();
        M->apply( r, v );
        //v.print();
      } else {
        r.copy( v );
      }
      rhoNext = r.dot( v );
      beta = rhoNext / rho;
      p.scale( beta );
      p.add( v );
      rho = rhoNext;
      ++iter;
      //if ( M ) std::cout << std::sqrt( std::abs( rho / rho0 ) ) << std::endl;
      //if ( M ) std::cout << r.norm2( ) / res0 << std::endl;
#ifdef MORE_VERBOSE
      if( isMsg ){
        std::cout << msgLvl << ": ";
      }
      if ( M ) std::cout << "P";
      std::cout << "CG iteration: " << iter << ". " <<
          "Relative error: " << std::sqrt( rho / rho0 ) << "." << std::endl;
#endif
    }

    bool ret = true;
    // if convergence, return true
    if ( iter == maxIt ) ret = false;

#ifdef VERBOSE
    A.apply( sol, r );
    r.add( rhs, -1.0 );
    if( isMsg ){
      std::cout << msgLvl << ": ";
    }
    if ( M ) std::cout << "P";
    std::cout << "CG iterations: " << iter << ". " <<
        "CG relative error: " << r.norm2( ) / res0 << "." << std::endl;
    if ( M ) {
      M->apply( r, v );
      if( isMsg ){
        std::cout << msgLvl << ": ";
      }
      std::cout << "  , " << "PCG relative error: " <<
          std::sqrt( std::abs( r.dot( v ) / rho0 ) ) << "." << std::endl;
    }
#endif
    --solver_level;
    return ret;
  }

  /*!
   * @brief Solves a system with general square matrix using GMRES
   * 
   * Computes solution of Ax = b where A is square, non-symmetric matrix
   * using Restarted Preconditioned Generalized Minimal Residual Method
   * 
   * @param M                 preconditioner
   * @param sol               solution vector 
   * @param rhs               right-hand side
   * @param precision[in,out] on input: solver relative precision, on output:
   *                          real obtained precision
   * @param maxIt[in,out]     on input: max. n. of iteration, on output:
   *                          real number of iterations
   */
  static bool GMRESSolve(
      LinearOperator<LO, SC> & A,
      Vector<LO, SC> const & rhs,
      Vector<LO, SC> & sol,
      SCVT & precision,
      LO maxIt,
      LO restarts,
      LinearOperator<LO, SC> * M = nullptr,
      const std::string & msg = ""
      ) {

    ++solver_level;
    bool isMsg = ( msg.compare( "" ) != 0 );
    std::string msgLvl;

    if( isMsg ){
#ifdef VERBOSE
      for (LO i = 0; i < solver_level; ++i)
        msgLvl.append( "   " );
      msgLvl.append( msg );
      std::cout << msgLvl << ": ";
#endif
    }

#ifdef VERBOSE
    std::cout << "GMRES started" << std::endl;
#endif

    // initialization
    restarts = std::min( restarts, rhs.getLength( ) );

    LO n = rhs.getLength( );
    LO m = restarts;
    LO itCount = 0;

    Vector<LO, SC> rhs2( rhs.getLength( ) );
    if ( M ) {
      M->apply( rhs, rhs2 );
    } else {
      rhs.copy( rhs2 );
    }
    SCVT bNorm = rhs.norm2( );//rhs2.norm2( );


    SCVT rNorm = 0.0;
    SC wVi;
    SC temp;

    if ( bNorm < EPS ) {
      bNorm = 1.0;
    }

    Vector<LO, SC> r0( rhs );

    Vector<LO, SC> r( n );
    Vector<LO, SC> s( m + 1 );

    //r0.setAll( 0.0 );
    A.apply( sol, r0, false, -1.0, 1.0 );
    if ( M ) {
      M->apply( r0, r );
    } else {
      r0.copy( r );
    }

    SCVT error = r.norm2( ) / bNorm;
    if ( error < precision ) {
      return true;
    }

    // initialize vectors and variables
    FullMatrix<LO, SC> V( n, m + 1 );
    FullMatrix<LO, SC> H( m + 1, m );
    Vector<LO, SC> cs( m );
    Vector<LO, SC> sn( m );
    Vector<LO, SC> e1( n );
    Vector<LO, SC> vi( n );
    Vector<LO, SC> Avi( n );
    Vector<LO, SC> w( n );
    cs.setAll( 0.0 );
    sn.setAll( 0.0 );
    e1.setAll( 0.0 );
    e1.set( 0, 1.0 );
    Avi.setAll( 0.0 );
    w.setAll( 0.0 );

    // current number of columns of matrix V
    LO currVSize = 0;

    while ( error > precision && itCount < maxIt ) {

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );
      if ( M ) {
        M->apply( r0, r );
      } else {
        r0.copy( r );
      }
      //r0.copy( r );
      rNorm = r.norm2( );

      for ( LO i = 0; i < n; i++ ) {
        V.set( i, 0, r.get( i ) / rNorm );
      }
      s.setAll( 0.0 );
      s.set( 0, rNorm );

      for ( LO i = 0; i < m; i++ ) {
        if ( itCount == maxIt ) {
          currVSize = i;
          break;
        }

        // construct orthonormal basis using Gram-Schmidt
        V.getCol( i, vi.getData( ) );
        A.apply( vi, Avi );
        if ( M ) {
          M->apply( Avi, w );
        } else {
          Avi.copy( w );
        }

        //Avi.copy( w );
        for ( LO k = 0; k <= i; k++ ) {
          V.getCol( k, vi.getData( ) );
          wVi = w.dot( vi );
          H.set( k, i, wVi );
          w.add( vi, -wVi );
        }

        H.set( i + 1, i, w.norm2( ) );

        for ( LO j = 0; j < n; j++ ) {
          V.set( j, i + 1, w.get( j ) / H.get( i + 1, i ) );
        }

        // apply Givens rotation
        for ( LO j = 0; j < i; j++ ) {
          temp = cs.get( j ) * H.get( j, i ) + sn.get( j ) * H.get( j + 1, i );
          H.set( j + 1, i, -sn.get( j ) * H.get( j, i ) + cs.get( j ) *
              H.get( j + 1, i ) );
          H.set( j, i, temp );
        }
        // form i-th rotation matrix
        SC a = H.get( i, i );
        SC b = H.get( i + 1, i );

        if ( b == (SCVT) 0.0 ) {
          cs.set( i, 1.0 );
          sn.set( i, 0.0 );
        } else if ( std::abs( b ) > std::abs( a ) ) {
          SC tmp = a / b;
          sn.set( i, (SCVT) 1.0 / std::sqrt( (SCVT) 1.0 + tmp * tmp ) );
          cs.set( i, tmp * sn.get( i ) );
        } else {
          SC tmp = b / a;
          cs.set( i, (SCVT) 1.0 / std::sqrt( (SCVT) 1.0 + tmp * tmp ) );
          sn.set( i, tmp * cs.get( i ) );
        }
        temp = cs.get( i ) * s.get( i );
        s.set( i + 1, -sn.get( i ) * s.get( i ) );
        s.set( i, temp );

        H.set( i, i, cs.get( i ) * H.get( i, i ) + sn.get( i ) *
            H.get( i + 1, i ) );
        H.set( i + 1, i, 0.0 );
        error = std::abs( s.get( i + 1 ) ) / bNorm;

        if ( error < precision ) {
          H.backward( s, 1, i + 1 );

          for ( LO j = 0; j < i + 1; j++ ) {
            V.getCol( j, vi.getData( ) );
            sol.add( vi, s.get( j ) );
            //std::cout << s.get(j) << std::endl;;
          }
          break;
        }

        ++currVSize;
        ++itCount;
#ifdef MORE_VERBOSE
        if( isMsg ){
          std::cout << msgLvl << ": ";
        }
        std::cout << "GMRES iteration: " << itCount << ". " <<
            "Relative error: " << error << std::endl;
#endif
      }
      if ( error < precision ) {
        break;
      }
      H.backward( s, 1, currVSize );
      for ( LO j = 0; j < currVSize; j++ ) {
        V.getCol( j, vi.getData( ) );
        sol.add( vi, s.get( j ) );
      }

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );
      if ( M ) {
        M->apply( r0, r );
      } else {
        r0.copy( r );
      }

      //r0.copy( r );
      s.set( m, r.norm2( ) );

      error = std::abs( s.get( m ) / bNorm );
      currVSize = 0;
    }


    rhs.copy( r0 );
    A.apply( sol, r0, false, -1.0, 1.0 );
    //Vector<LO, SC> bb( r0.getLength( ) );
    //M->apply( r0, bb );

#ifdef VERBOSE
    if( isMsg ){
      std::cout << msgLvl << ": ";
    }
    std::cout << "GMRES finished in " << itCount <<
        " iterations with relative precision " << r0.norm2( ) / bNorm <<
        "." << std::endl;
#endif

    // if convergence, return true
    --solver_level;
    return ( error <= precision );
  }

  /*!
   * @brief Solves a system with general square matrix using GMRES
   * 
   * Computes solution of Ax = b where A is square, non-symmetric matrix
   * using Flexible Restarted Preconditioned Generalized Minimal Residual Method
   * The method uses a right preconditioning and allows changing preconditioner 
   * in every iteration.
   * 
   * @param M                 right preconditioner
   * @param sol               solution vector 
   * @param rhs               right-hand side
   * @param precision[in,out] on input: solver relative precision, on output:
   *                          real obtained precision
   * @param maxIt[in,out]     on input: max. n. of iteration, on output:
   *                          real number of iterations
   */
  static bool FGMRESSolve(
      LinearOperator<LO, SC> &A,
      Vector<LO, SC> const &rhs,
      Vector<LO, SC> &sol,
      SCVT &precision,
      LO maxIt,
      LO restarts,
      LinearOperator<LO, SC> *M = nullptr,
      const std::string & msg = ""
      ) {

    ++solver_level;
    bool isMsg = ( msg.compare( "" ) != 0 );\
    std::string msgLvl;

#ifdef VERBOSE
    if( isMsg ){
      for (LO i = 0; i < solver_level; ++i)
        msgLvl.append( "   " );
      msgLvl.append( msg );
      std::cout << msgLvl << ": ";
    }
    std::cout << "FGMRES started" << std::endl;
#endif
    // initialization
    restarts = std::min( restarts, rhs.getLength( ) );

    LO n = rhs.getLength( );
    LO m = restarts;
    LO itCount = 0;

    SCVT bNorm = rhs.norm2( );

    SCVT rNorm = 0.0;
    SC wVi;
    SC temp;

    if ( bNorm < EPS ) {
      bNorm = 1.0;
    }

    Vector<LO, SC> r0( rhs );

    Vector<LO, SC> r( n );
    Vector<LO, SC> s( m + 1 );

    //r0.setAll( 0.0 );
    A.apply( sol, r0, false, -1.0, 1.0 );

    SCVT error = r0.norm2( ) / bNorm;
    if ( error < precision ) {
      --solver_level;
      return true;
    }

    // initialize vectors and variables
    FullMatrix<LO, SC> V( n, m + 1 );
    FullMatrix<LO, SC> Z( n, m + 1 );
    FullMatrix<LO, SC> H( m + 1, m );
    Vector<LO, SC> cs( m );
    Vector<LO, SC> sn( m );
    Vector<LO, SC> e1( n );
    Vector<LO, SC> vi( n );
    Vector<LO, SC> Mvi( n );
    Vector<LO, SC> w( n );
    Vector<LO, SC> partSol( n );
    Vector<LO, SC> MPartSol( n );
    cs.setAll( 0.0 );
    sn.setAll( 0.0 );
    e1.setAll( 0.0 );
    e1.set( 0, 1.0 );
    Mvi.setAll( 0.0 );
    w.setAll( 0.0 );

    // current number of columns of matrix V
    LO currVSize = 0;

    while ( error > precision && itCount < maxIt ) {

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );
      //r0.copy( r );
      rNorm = r0.norm2( );

      for ( LO i = 0; i < n; i++ ) {
        V.set( i, 0, r0.get( i ) / rNorm );
      }
      s.setAll( 0.0 );
      s.set( 0, rNorm );

      for ( LO i = 0; i < m; i++ ) {
        if ( itCount == maxIt ) {
          currVSize = i;
          break;
        }

        // construct orthonormal basis using Gram-Schmidt
        V.getCol( i, vi.getData( ) );
        if ( M ) {
          M->apply( vi, Mvi );
        } else {
          vi.copy( Mvi );
        }
        A.apply( Mvi, w );

        for ( LO k = 0; k < n; k++ ) {
          Z.set( k, i, Mvi.get( k ) );
        }

        //Avi.copy( w );
        for ( LO k = 0; k <= i; k++ ) {
          V.getCol( k, vi.getData( ) );
          wVi = w.dot( vi );
          H.set( k, i, wVi );
          w.add( vi, -wVi );
        }

        H.set( i + 1, i, w.norm2( ) );

        for ( LO j = 0; j < n; j++ ) {
          V.set( j, i + 1, w.get( j ) / H.get( i + 1, i ) );
        }

        // apply Givens rotation
        for ( LO j = 0; j < i; j++ ) {
          temp = cs.get( j ) * H.get( j, i ) + sn.get( j ) * H.get( j + 1, i );
          H.set( j + 1, i, -sn.get( j ) * H.get( j, i ) + cs.get( j ) * H.get( j + 1, i ) );
          H.set( j, i, temp );
        }
        // form i-th rotation matrix
        SC a = H.get( i, i );
        SC b = H.get( i + 1, i );
        if ( b == (SCVT) 0.0 ) {
          cs.set( i, 1.0 );
          sn.set( i, 0.0 );
        } else if ( std::abs( b ) > std::abs( a ) ) {
          SC tmp = a / b;
          sn.set( i, (SCVT) 1.0 / std::sqrt( (SCVT) 1.0 + tmp * tmp ) );
          cs.set( i, tmp * sn.get( i ) );
        } else {
          SC tmp = b / a;
          cs.set( i, (SCVT) 1.0 / std::sqrt( (SCVT) 1.0 + tmp * tmp ) );
          sn.set( i, tmp * cs.get( i ) );
        }
        temp = cs.get( i ) * s.get( i );
        s.set( i + 1, -sn.get( i ) * s.get( i ) );
        s.set( i, temp );

        H.set( i, i, cs.get( i ) * H.get( i, i ) + sn.get( i ) * H.get( i + 1, i ) );
        H.set( i + 1, i, 0.0 );
        error = std::abs( s.get( i + 1 ) ) / bNorm;

        if ( error < precision ) {

          H.backward( s, 1, i + 1 );

          for ( LO j = 0; j < i + 1; j++ ) {
            Z.getCol( j, vi.getData( ) );
            sol.add( vi, s.get( j ) );
          }
          break;
        }

        currVSize++;
        itCount++;
#ifdef MORE_VERBOSE
        if( isMsg ){
          std::cout << msgLvl << ": ";
        }
        std::cout << "FGMRES iteration: " << itCount << ". " <<
            "Relative error: " << error << std::endl;
#endif
      }
      if ( error < precision ) {
        break;
      }
      H.backward( s, 1, currVSize );
      for ( LO j = 0; j < currVSize; j++ ) {
        Z.getCol( j, vi.getData( ) );
        sol.add( vi, s.get( j ) );
      }

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );
      //M->apply( r0, r );

      //r0.copy( r );
      s.set( m, r0.norm2( ) );

      error = std::abs( s.get( m ) / bNorm );
      currVSize = 0;
    }


    rhs.copy( r0 );
    A.apply( sol, r0, false, -1.0, 1.0 );
    //Vector<LO, SC> bb( r0.getLength( ) );
    //M->apply( r0, bb );

#ifdef VERBOSE
    if( isMsg ){
      std::cout << msgLvl << ": ";
    }
    std::cout << "FGMRES finished in " << itCount << " iterations with relative " <<
        "precision " << r0.norm2( ) / bNorm << "." << std::endl;
#endif

    // if convegence, return true
    --solver_level;
    return ( error <= precision );
  }

  /*!
   * @brief Solves a system with general square matrix using deflated GMRES 
   * 
   * Computes solution of Ax = b where A is square, non-symmetric matrix
   * using Restarted Preconditioned Generalized Minimal Residual Method
   * 
   * @param sol       solution vector 
   * @param rhs       right-hand side
   * @param M         left preconditioner
   * @param precision[in,out] on input: solver relative precision, on output:
   *                          real obtained precision
   * @param maxIt[in,out]     on input: max. n. of iteration, on output:
   *                          real number of iterations
   */
  static bool DGMRESSolve(
      LinearOperator<LO, SC> &A,
      Vector<LO, SC> const &rhs,
      Vector<LO, SC> &sol,
      SCVT precision,
      LO maxIt,
      LO restarts,
      LO maxInvSpaceDim,
      LO maxEigsPerRestart,
      LinearOperator<LO, SC> *M,
      const std::string & msg = ""
      ) {

    ++solver_level;
    bool isMsg = ( msg.compare( "" ) != 0 );
    std::string msgLvl;

#ifdef VERBOSE
    if( isMsg ){
      for (LO i = 0; i < solver_level; ++i)
        msgLvl.append( "   " );
      msgLvl.append( msg );
      std::cout << msgLvl << ": ";
    }
    std::cout << "DGMRES started" << std::endl;
#endif

    // initialization
    LO n = rhs.getLength( );
    LO m = restarts;
    LO itCount = 0;

    SCVT bNorm = std::abs( rhs.norm2( ) );
    SCVT rNorm = 0.0;
    SC wVi;
    SC temp;

    if ( bNorm < EPS ) {
      bNorm = 1.0;
    }

    Vector<LO, SC> r0( rhs );

    Vector<LO, SC> r( rhs.getLength( ) );
    Vector<LO, SC> s( m + 1 );

    //r0.setAll( 0.0 );
    A.apply( sol, r0, false, -1.0, 1.0 );
    M->apply( r0, r );

    //r0.copy( r );


    SCVT error = r.norm2( ) / bNorm;
    if ( error < precision ) {
      --solver_level;
      return true;
    }

    // initialize vectors and variables
    FullMatrix<LO, SC> V( n, m + 1 );
    FullMatrix<LO, SC> H( m + 1, m );
    Vector<LO, SC> cs( m );
    Vector<LO, SC> sn( m );
    Vector<LO, SC> e1( n );
    Vector<LO, SC> vi( n );
    Vector<LO, SC> Avi( n );
    Vector<LO, SC> w( n );
    Vector<LO, SC> sVi( n );
    Vector<LO, SC> deflSVi( n );
    Vector<LO, SC> viDefl( n );

    // data necessary for deflation
    FullMatrix<LO, SC> deflH( m, m );
    FullMatrix<LO, SC> U( n, maxInvSpaceDim );
    FullMatrix<LO, SC> AU( n, maxInvSpaceDim );
    FullMatrix<LO, SC> T( maxInvSpaceDim, maxInvSpaceDim );
    FullMatrix<LO, SC> locT;
    //todo: allow reusie of factorization again
    //locT.setFactReuseOn( ); // keep the factorization
    LO currentInvSpaceDim = 0;
    SCVT lambdaN = 0.0;
    bool deflInitialized = false;

    cs.setAll( 0.0 );
    sn.setAll( 0.0 );
    e1.setAll( 0.0 );
    e1.set( 0, 1.0 );
    Avi.setAll( 0.0 );
    w.setAll( 0.0 );

    // current number of columns of matrix V
    LO currVSize = 0;

    while ( error > precision && itCount < maxIt ) {

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );
      M->apply( r0, r );
      //r0.copy( r );
      rNorm = r.norm2( );

      for ( LO i = 0; i < n; i++ ) {
        V.set( i, 0, r.get( i ) / rNorm );
      }
      s.setAll( 0.0 );
      s.set( 0, rNorm );

      for ( LO i = 0; i < m; i++ ) {

        if ( itCount == maxIt ) {
          currVSize = i;
          break;
        }

        // construct orthonormal basis using Gram-Schmidt
        V.getCol( i, vi.getData( ) );
        if ( deflInitialized ) {
          applyDeflation( U, locT, currentInvSpaceDim, vi, viDefl, lambdaN );
        } else {
          vi.copy( viDefl );
        }
        A.apply( viDefl, w );
        //    M->apply( Avi, w );

        //Avi.copy( w );

        // orthogonalization by modified Gram-Schmidt
        for ( LO k = 0; k <= i; k++ ) {
          V.getCol( k, vi.getData( ) );
          wVi = w.dot( vi );
          H.set( k, i, wVi );
          deflH.set( k, i, wVi );
          w.add( vi, -wVi );
        }

        H.set( i + 1, i, w.norm2( ) );
        if ( i < m - 1 ) {
          deflH.set( i + 1, i, w.norm2( ) );
        }

        for ( LO j = 0; j < n; j++ ) {
          V.set( j, i + 1, w.get( j ) / H.get( i + 1, i ) );
        }

        // apply Givens rotation
        for ( LO j = 0; j < i; j++ ) {
          temp = cs.get( j ) * H.get( j, i ) + sn.get( j ) * H.get( j + 1, i );
          H.set( j + 1, i, -sn.get( j ) * H.get( j, i ) + cs.get( j ) * H.get( j + 1, i ) );
          H.set( j, i, temp );
        }
        // form i-th rotation matrix
        SC a = H.get( i, i );
        SC b = H.get( i + 1, i );
        if ( b == 0.0 ) {
          cs.set( i, 1.0 );
          sn.set( i, 0.0 );
        } else if ( std::abs( b ) > std::abs( a ) ) {
          SC tmp = a / b;
          sn.set( i, 1.0 / std::sqrt( 1.0 + tmp * tmp ) );
          cs.set( i, tmp * sn.get( i ) );
        } else {
          SC tmp = b / a;
          cs.set( i, 1.0 / std::sqrt( 1.0 + tmp * tmp ) );
          sn.set( i, tmp * cs.get( i ) );
        }
        temp = cs.get( i ) * s.get( i );
        s.set( i + 1, -sn.get( i ) * s.get( i ) );
        s.set( i, temp );

        H.set( i, i, cs.get( i ) * H.get( i, i ) +
            sn.get( i ) * H.get( i + 1, i ) );
        H.set( i + 1, i, 0.0 );
        error = std::abs( s.get( i + 1 ) ) / bNorm;

#ifdef MORE_VERBOSE
        if( isMsg ){
          std::cout << msgLvl << ": ";
        }
        std::cout << "DGMRES iteration: " << itCount << ". " <<
            "Relative error: " << error << ". " <<
            "Invariant subspace dimension: " << currentInvSpaceDim << std::endl;
#endif

        if ( error <= precision ) {
          H.backward( s, 1, i + 1 );
          sVi.setAll( 0.0 );
          for ( LO j = 0; j < i + 1; j++ ) {
            V.getCol( j, vi.getData( ) );
            sVi.add( vi, s.get( j ) );
            //sol.add( vi, s.get( j ) );
          }
          if ( deflInitialized ) {
            applyDeflation( U, locT, currentInvSpaceDim, sVi, deflSVi, lambdaN );
          } else {
            sVi.copy( deflSVi );
          }
          sol.add( deflSVi );
          break;
        }
        itCount++;
        currVSize++;

      }
      if ( error <= precision ) {
        break;
      }

      H.backward( s, 1, currVSize );
      sVi.setAll( 0.0 );
      for ( LO j = 0; j < currVSize; j++ ) {
        V.getCol( j, vi.getData( ) );
        sVi.add( vi, s.get( j ) );
      }

      if ( deflInitialized ) {
        applyDeflation( U, locT, currentInvSpaceDim, sVi, deflSVi, lambdaN );
      } else {
        sVi.copy( deflSVi );
      }

      sol.add( deflSVi );

      rhs.copy( r0 ); // probably not very efficient, watch out!
      A.apply( sol, r0, false, -1.0, 1.0 );

      M->apply( r0, r );

      computeDeflationData( A, V, deflH, U, AU, T, locT, lambdaN,
          currentInvSpaceDim, maxInvSpaceDim, maxEigsPerRestart );
      deflH.setAll( 0.0 );

      deflInitialized = true;


      //r0.copy( r );
      s.set( m, r.norm2( ) );
      H.setAll( 0.0 );
      error = std::abs( s.get( m ) / bNorm );

      currVSize = 0;

    }

    // check the error (remove later)
    rhs.copy( r0 );
    A.apply( sol, r0, false, -1.0, 1.0 );

#ifdef VERBOSE
    if( isMsg ){
      std::cout << msgLvl << ": ";
    }
    std::cout << "DGMRES finished in " << itCount << " iterations with relative " <<
        "precision " << r0.norm2( ) << ".";
#endif
    --solver_level;
    return ( error <= precision );
  }

protected:

  static void computeDeflationData(
      LinearOperator<LO, SC> &mat,
      FullMatrix<LO, SC> &V,
      FullMatrix<LO, SC> &deflH,
      FullMatrix<LO, SC> &U,
      FullMatrix<LO, SC> &AU,
      FullMatrix<LO, SC> &T,
      FullMatrix<LO, SC> &locT,
      SCVT &lambdaN,
      LO &currentEigCount,
      LO maxNEigs,
      LO maxEigsToAdd ) {

    if ( currentEigCount >= maxNEigs ) {
      return;
    }

    LO m = deflH.getNRows( );
    LO n = V.getNRows( );

    // decompose deflH into a Hessenberg form
    std::complex<SCVT> *eigs = new std::complex<SCVT>[m];

    FullMatrix<LO, SC> schurU( m, deflH.getNCols( ) );

    deflH.schurOfHessenberg( schurU, eigs );


    //    for (int i = 0; i < m ; i++) {
    //      std::cout << eigs[i].real() << "+"<< eigs[i].imag()<<"j" << std::endl;
    //    }

    // count the number of eigenvalues
    LO nEigs = 0;

    while ( ( nEigs + currentEigCount < maxNEigs ) && ( nEigs < maxEigsToAdd ) ) {
      //if ( eigs[ m - nEigs - 1].imag( ) == 0.0 ) {
      nEigs++;
      //} else {
      //  nEigs += 1;
      //}
    }


    std::vector<SC> vEigs;
    std::vector<LO> perm;
    for ( LO i = 0; i < m; i++ ) {
      vEigs.push_back( std::abs( eigs[i] ) );
      perm.push_back( i );
    }
    // a pretty little lambda function inserted to sort (we need to sort 
    // eigenvalues and get the permutation vector)
    std::sort( perm.begin( ), perm.end( ), [&vEigs]( LO a, LO b ) -> bool {
      return vEigs[a] > vEigs[b];
    } );

    if ( std::abs( eigs[perm[0]] ) > lambdaN ) {
      lambdaN = std::abs( eigs[perm[0]] );
    }

    // save Schur vectors corresponding to the smallest eigenvalues
    FullMatrix<LO, SC> S( m, nEigs );
    for ( LO i = 0; i < nEigs; i++ ) {
      for ( LO j = 0; j < m; j++ ) {
        S.set( j, i, schurU.get( j, perm[m - i - 1] ) );
      }
    }

    // orthogonalize VS against U
    FullMatrix<LO, SC> X( n, nEigs );
    X.multiply( V, S, n, m, S.getNRows( ), S.getNCols( ) );

    SC* tmpXColData = new SC[n];
    SC* tmpUColData = new SC[n];

    if ( currentEigCount > 0 ) {
      for ( LO i = 0; i < nEigs; i++ ) {
        for ( LO j = 0; j < currentEigCount; j++ ) {
          X.getCol( i, tmpXColData );
          Vector<LO, SC> tmpXCol( n, tmpXColData, false );
          U.getCol( j, tmpUColData );
          Vector<LO, SC> tmpUCol( n, tmpUColData, false );
          SC dot = tmpXCol.dot( tmpUCol );
          for ( int k = 0; k < n; k++ ) {
            X.set( k, i, X.get( k, i ) - dot * U.get( k, j ) );
          }
        }
      }
    }

    int mpiInitialized;
    MPI_Initialized( &mpiInitialized );
    if ( mpiInitialized ) {
      if ( typeid (SC ) == typeid (float ) ) {
        MPI_Bcast( X.getData( ), nEigs*n, MPI_FLOAT, 0, MPI_COMM_WORLD );
      } else if ( typeid (SC ) == typeid (double ) ) {
        MPI_Bcast( X.getData( ), nEigs*n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      } else {
        std::cout << "Not yet implemented for complex matrices." << std::endl;
      }
    }

    FullMatrix<LO, SC> AX( n, nEigs );
    Vector<LO, SC> Ax( n );
    for ( LO i = 0; i < nEigs; i++ ) {
      X.getCol( i, tmpXColData );
      Vector<LO, SC> tmpXCol( n, tmpXColData, false );
      mat.apply( tmpXCol, Ax );
      for ( LO j = 0; j < n; j++ ) {
        AX.set( j, i, Ax.get( j ) );
      }
    }

    // update matrix T = [ UAU, UAX; XAU, XAX ]
    FullMatrix<LO, SC> XAX( nEigs, nEigs );
    XAX.multiply( X, AX, true, false );

    for ( LO i = currentEigCount; i < currentEigCount + nEigs; i++ ) {
      for ( LO j = currentEigCount; j < currentEigCount + nEigs; j++ ) {
        T.set( i, j, XAX.get( i - currentEigCount, j - currentEigCount ) );
      }
    }

    if ( currentEigCount > 0 ) {
      FullMatrix<LO, SC> UAX( currentEigCount, nEigs );
      FullMatrix<LO, SC> XAU( nEigs, currentEigCount );
      UAX.multiply( U, AX, n, currentEigCount, n, nEigs, true, false );
      XAU.multiply( X, AU, n, nEigs, n, currentEigCount, true, false );

      for ( LO i = 0; i < currentEigCount; i++ ) {
        for ( LO j = currentEigCount; j < currentEigCount + nEigs; j++ ) {
          T.set( i, j, UAX.get( i, j - currentEigCount ) );
        }
      }
      for ( LO i = currentEigCount; i < currentEigCount + nEigs; i++ ) {
        for ( LO j = 0; j < currentEigCount; j++ ) {
          T.set( i, j, XAU.get( i - currentEigCount, j ) );
        }
      }
    }

    for ( LO i = 0; i < nEigs; i++ ) {
      for ( LO j = 0; j < n; j++ ) {
        U.set( j, currentEigCount + i, X.get( j, i ) );
        AU.set( j, currentEigCount + i, AX.get( j, i ) );
      }
    }

    currentEigCount += nEigs;

    locT.resize( currentEigCount, currentEigCount );
    for ( LO i = 0; i < currentEigCount; i++ ) {
      for ( LO j = 0; j < currentEigCount; j++ ) {
        locT.set( i, j, T.get( i, j ) );
      }
    }

    delete [] eigs;
    delete [] tmpXColData;
    delete [] tmpUColData;
  }

  static void applyDeflation(
      FullMatrix<LO, SC> &U,
      FullMatrix<LO, SC> &locT,
      LO currentInvSpaceDim,
      Vector<LO, SC> &x,
      Vector<LO, SC> &y,
      SCVT &lambdaN
      ) {

    Vector<LO, SC> y1( currentInvSpaceDim );
    //    FullMatrix<LO, SC> locT( currentInvSpaceDim, currentInvSpaceDim );

    //    for ( LO i = 0; i < currentInvSpaceDim; i++ ) {
    //      for ( LO j = 0; j < currentInvSpaceDim; j++ ) {
    //        locT.set( i, j, T.get( i, j ) );
    //      }
    //    }

    U.applySubmatrix( x, y1, U.getNRows( ), currentInvSpaceDim, true );

    Vector<LO, SC> y2( y1 );
    locT.LUSolve( y2 );

    y2.scale( lambdaN );

    y2.add( y1, -1.0 );
    x.copy( y );

    U.applySubmatrix( y2, y, U.getNRows( ), currentInvSpaceDim, false, 1.0, 1.0 );

  }


private:

};
//
//template<>
//void IterativeSolver<int, std::complex<double> >::GMRESSolve(
//    LinearOperator<int, std::complex<double> > &A,
//    Vector<int, std::complex<double> > const &rhs,
//    Vector<int, std::complex<double> > &sol,
//    double precision,
//    int restarts,
//    int maxIt,
//    LeftPreconditioner<int, std::complex<double> > *M
//    ) {
//
//  // initialization
//  int n = rhs.getLength( );
//  int m = restarts;
//  int itCount = 0;
//
//  double bNorm = std::abs( rhs.norm2( ) );
//  double rNorm = 0.0;
//  std::complex<double> wVi;
//  std::complex<double> temp;
//
//  if ( bNorm < EPS ) {
//    bNorm = 1.0;
//  }
//
//  Vector<int, std::complex<double> > r0( rhs );
//
//  Vector<int, std::complex<double> > r( rhs.getLength( ) );
//  Vector<int, std::complex<double> > s( rhs.getLength( ) );
//
//  //r0.setAll( 0.0 );
//  A.apply( sol, r0, false, -1.0, 1.0 );
//  M->apply( r0, r );
//
//  //r0.copy( r );
//
//
//  double error = std::abs( r.norm2( ) / bNorm );
//  if ( error < precision ) {
//    return;
//  }
//
//  // initialize vectors and variables
//  FullMatrix<int, std::complex<double> > V( n, m + 1 );
//  FullMatrix<int, std::complex<double> > H( m + 1, m );
//  Vector<int, std::complex<double> > cs( m );
//  Vector<int, std::complex<double> > sn( m );
//  Vector<int, std::complex<double> > e1( n );
//  Vector<int, std::complex<double> > vi( n );
//  Vector<int, std::complex<double> > Avi( n );
//  Vector<int, std::complex<double> > w( n );
//  cs.setAll( 0.0 );
//  sn.setAll( 0.0 );
//  e1.setAll( 0.0 );
//  e1.set( 0, 1.0 );
//  Avi.setAll( 0.0 );
//  w.setAll( 0.0 );
//
//  while ( error > precision && itCount < maxIt ) {
//
//    rhs.copy( r0 ); // probably not very efficient, watch out!
//    A.apply( sol, r0, false, -1.0, 1.0 );
//    M->apply( r0, r );
//    //r0.copy( r );
//    rNorm = std::abs( r.norm2( ) );
//
//    for ( int i = 0; i < n; i++ ) {
//      V.set( i, 0, r.get( i ) / rNorm );
//    }
//    s.setAll( 0.0 );
//    s.set( 0, rNorm );
//
//    for ( int i = 0; i < m; i++ ) {
//      if ( maxIt == 3000 )
//        std::cout << "Inner error: " << error << "iteration: " << itCount << std::endl;
//
//      itCount++;
//      // construct orthonormal basis using Gram-Schmidt
//      V.getCol( i, vi.getData( ) );
//      A.apply( vi, Avi );
//      M->apply( Avi, w );
//
//      //Avi.copy( w );
//      for ( int k = 0; k <= i; k++ ) {
//        V.getCol( k, vi.getData( ) );
//        wVi = w.dot( vi );
//        H.set( k, i, wVi );
//        w.add( vi, -wVi );
//      }
//
//      H.set( i + 1, i, w.norm2( ) );
//
//      for ( int j = 0; j < n; j++ ) {
//        V.set( j, i + 1, w.get( j ) / H.get( i + 1, i ) );
//      }
//
//      // apply Givens rotation
//      for ( int j = 0; j < i; j++ ) {
//        temp = cs.get( j ) * H.get( j, i ) + sn.get( j ) * H.get( j + 1, i );
//        H.set( j + 1, i, -sn.get( j ) * H.get( j, i ) + cs.get( j ) * H.get( j + 1, i ) );
//        H.set( j, i, temp );
//      }
//      // form i-th rotation matrix
//      std::complex<double> a = H.get( i, i );
//      std::complex<double> b = H.get( i + 1, i );
//      if ( b == 0.0 ) {
//        cs.set( i, 1.0 );
//        sn.set( i, 0.0 );
//      } else if ( std::abs( b ) > std::abs( a ) ) {
//        std::complex<double> tmp = a / b;
//        sn.set( i, 1.0 / std::sqrt( 1.0 + tmp * tmp ) );
//        cs.set( i, tmp * sn.get( i ) );
//      } else {
//        std::complex<double> tmp = b / a;
//        cs.set( i, 1.0 / std::sqrt( 1.0 + tmp * tmp ) );
//        sn.set( i, tmp * cs.get( i ) );
//      }
//      temp = cs.get( i ) * s.get( i );
//      s.set( i + 1, -sn.get( i ) * s.get( i ) );
//      s.set( i, temp );
//
//      H.set( i, i, cs.get( i ) * H.get( i, i ) + sn.get( i ) * H.get( i + 1, i ) );
//      H.set( i + 1, i, 0.0 );
//      error = std::abs( s.get( i + 1 ) ) / bNorm;
//      if ( error < precision ) {
//
//        H.backward( s, 1, i + 1 );
//
//        for ( int j = 0; j < i; j++ ) {
//          V.getCol( j, vi.getData( ) );
//          sol.add( vi, s.get( j ) );
//        }
//        break;
//      }
//    }
//    if ( error <= precision ) {
//      break;
//    }
//    H.backward( s, 1, m );
//    for ( int j = 0; j < m; j++ ) {
//      V.getCol( j, vi.getData( ) );
//      sol.add( vi, s.get( j ) );
//    }
//
//    rhs.copy( r0 ); // probably not very efficient, watch out!
//    A.apply( sol, r0, false, -1.0, 1.0 );
//    //std::cout << "RESIDUAL" << r0.norm2()  << std::endl;
//    M->apply( r0, r );
//
//    //r0.copy( r );
//    s.set( m + 1, r.norm2( ) );
//    error = std::abs( s.get( m + 1 ) / bNorm );
//  }
//
//  if ( maxIt == 3000 ) {
//    std::cout << "Number of iterations: " << itCount << ", error: " << error << std::endl;
//  }
//  if ( maxIt == 200 ) {
//    //  std::cout << "  First prec. solve. Number of iterations: " << itCount << ", error: " << error << ", n: " << rhs.getLength() << std::endl;
//  }
//  if ( maxIt == 201 ) {
//    //  std::cout << "  Secnd prec. solve. Number of iterations: " << itCount << ", error: " << error <<  ", n: " << rhs.getLength() <<std::endl;
//  }
//
//}

}

#endif	/* LINEAROPERATOR_H */

