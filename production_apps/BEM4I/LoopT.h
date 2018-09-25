/*===========================================================================*\
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2014 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------* 
 *  This file is part of OpenMesh.                                           *
 *                                                                           *
 *  OpenMesh is free software: you can redistribute it and/or modify         * 
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenMesh is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenMesh.  If not,                                    *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *             
 *   $Revision: 995 $                                                        *
 *   $Date: 2014-02-12 15:29:15 +0100 (Mi, 12 Feb 2014) $                    *
 *   Modified for bem4i                                                      *
 *                                                                           *
\*===========================================================================*/

/** \file LoopT.h

 */

//=============================================================================
//
//  CLASS LoopT
//
//=============================================================================

#ifndef OPENMESH_SUBDIVIDER_UNIFORM_LOOPT_HH
#define OPENMESH_SUBDIVIDER_UNIFORM_LOOPT_HH


//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/Property.hh>
// -------------------- STL
#include <vector>
#if defined(OM_CC_MIPS)
#include <math.h>
#else
#include <cmath>
#endif

#include "SubdividerT.h"
#include "SparseMatrix.h"
#include "MeshTags.h"

//== NAMESPACE ================================================================

namespace OpenMesh { // BEGIN_NS_OPENMESH
namespace Subdivider { // BEGIN_NS_SUBDIVIDER
namespace Uniform { // BEGIN_NS_UNIFORM


//== CLASS DEFINITION =========================================================

/** %Uniform Loop subdivision algorithm.
 *
 *  Implementation as described in
 *
 *  C. T. Loop, "Smooth Subdivision Surfaces Based on Triangles",
 *  M.S. Thesis, Department of Mathematics, University of Utah, August 1987.
 *
 */
template <typename MeshType, typename RealType = float, typename LO = int>
class LoopT : public SubdividerT<MeshType, RealType, LO> {

public:

  typedef RealType real_t;
  typedef MeshType mesh_t;
  typedef SubdividerT< mesh_t, real_t, LO > parent_t;

  typedef std::pair< real_t, real_t > weight_t;
  typedef std::vector< std::pair<real_t, real_t> > weights_t;

public:

  LoopT( ) : parent_t( ), _1over8( 1.0 / 8.0 ), _3over8( 3.0 / 8.0 ) {
    init_weights( );
  }

  LoopT(
      mesh_t& _m
      ) : parent_t( _m ), _1over8( 1.0 / 8.0 ), _3over8( 3.0 / 8.0 ) {
    init_weights( );
  }

  ~LoopT( ) {
  }


public:

  const char *name( ) const {
    return "Uniform Loop";
  }


  /// Pre-compute weights

  void init_weights(
      size_t _max_valence = 50
      ) {
    weights_.resize( _max_valence );
    std::generate( weights_.begin( ), weights_.end( ), compute_weight( ) );
  }

  //  bool subdivideAndCreateMatrix(
  //      mesh_t & _m,
  //      size_t _n,
  //      bem4i::SparseMatrix< int, RealType > * matrix,
  //      const bool _update_points = true
  //      ) {
  //
  //    ///TODO:Implement fixed positions   
  //
  //    typename mesh_t::FaceIter fit, f_end;
  //    typename mesh_t::EdgeIter eit, e_end;
  //    typename mesh_t::VertexIter vit;
  //
  //    int nNodes;
  //    int nEdges;
  //    int nnz;
  //    std::vector< int > * rowInd;
  //    std::vector< int > * colInd;
  //    std::vector< RealType > * values;
  //
  //    // Do _n subdivisions
  //    for ( size_t i = 0; i < _n; ++i ) {
  //      nNodes = _m.n_vertices( );
  //      nEdges = _m.n_edges( );
  //      nnz = 4 * nEdges + 7 * nNodes;
  //      rowInd = new std::vector< int >( );
  //      colInd = new std::vector< int >( );
  //      values = new std::vector< RealType >( );
  //      rowInd->reserve( nnz );
  //      colInd->reserve( nnz );
  //      values->reserve( nnz );
  //
  //      if ( _update_points ) {
  //        // compute new positions for old vertices
  //        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
  //          smooth( _m, *vit, *rowInd, *colInd, *values );
  //        }
  //      } else {
  //        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
  //          rowInd->push_back( vit->idx() );
  //          colInd->push_back( vit->idx() );
  //          values->push_back( 1.0 );
  //        }
  //      }
  //
  //      // Compute position for new vertices and store them in the edge property
  //      for ( eit = _m.edges_begin( ); eit != _m.edges_end( ); ++eit )
  //        compute_midpoint( _m, *eit, *rowInd, *colInd, *values, _update_points );
  //
  //      // Split each edge at midpoint and store precomputed positions (stored in
  //      // edge property ep_pos_) in the vertex property vp_pos_;
  //
  //      // Attention! Creating new edges, hence make sure the loop ends correctly.
  //      e_end = _m.edges_end( );
  //      for ( eit = _m.edges_begin( ); eit != e_end; ++eit )
  //        split_edge( _m, *eit );
  //
  //
  //      // Commit changes in topology and reconsitute consistency
  //
  //      // Attention! Creating new faces, hence make sure the loop ends correctly.
  //      f_end = _m.faces_end( );
  //      for ( fit = _m.faces_begin( ); fit != f_end; ++fit )
  //        split_face( _m, *fit );
  //
  //      if ( _update_points ) {
  //        // Commit changes in geometry
  //        for ( vit = _m.vertices_begin( );
  //            vit != _m.vertices_end( ); ++vit ) {
  //          _m.set_point( *vit, _m.property( vp_pos_, *vit ) );
  //        }
  //      }
  //
  //
  //#if defined(_DEBUG) || defined(DEBUG)
  //      // Now we have an consistent mesh!
  //      assert( OpenMesh::Utils::MeshCheckerT<mesh_t>( _m ).check( ) );
  //#endif
  //
  //      if ( i == 0 ) {
  //        matrix->setFromTriplets( nNodes + nEdges, nNodes, *rowInd, *colInd,
  //            *values );
  //      } else {
  //        bem4i::SparseMatrix< int, RealType > * thisLevel =
  //            new bem4i::SparseMatrix< int, RealType >( );
  //        thisLevel->setFromTriplets( nNodes + nEdges, nNodes, *rowInd, *colInd,
  //            *values );
  //
  //        bem4i::SparseMatrix< int, RealType > * prevLevels =
  //            new bem4i::SparseMatrix< int, RealType >( *matrix );
  //
  //        matrix->resize( thisLevel->getNRows( ), prevLevels->getNCols( ) );
  //        matrix->multiply( *thisLevel, *prevLevels );
  //
  //        delete thisLevel;
  //        delete prevLevels;
  //      }
  //      delete rowInd;
  //      delete colInd;
  //      delete values;
  //    }
  //
  //    return true;
  //  }

protected:

  bool prepare(
      mesh_t& _m
      ) {
    _m.add_property( vp_pos_ );
    _m.add_property( ep_pos_ );
    return true;
  }

  bool cleanup(
      mesh_t& _m
      ) {
    _m.remove_property( vp_pos_ );
    _m.remove_property( ep_pos_ );
    return true;
  }

  //  bool subdivide(
  //      mesh_t& _m,
  //      size_t _n,
  //      const bool _update_points = true
  //      ) {
  //
  //    ///TODO:Implement fixed positions   
  //
  //    typename mesh_t::FaceIter fit, f_end;
  //    typename mesh_t::EdgeIter eit, e_end;
  //    typename mesh_t::VertexIter vit;
  //
  //    // Do _n subdivisions
  //    for ( size_t i = 0; i < _n; ++i ) {
  //
  //      if ( _update_points ) {
  //        // compute new positions for old vertices
  //        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
  //          smooth( _m, *vit );
  //        }
  //      }
  //
  //      // Compute position for new vertices and store them in the edge property
  //      for ( eit = _m.edges_begin( ); eit != _m.edges_end( ); ++eit )
  //        compute_midpoint( _m, *eit );
  //
  //      // Split each edge at midpoint and store precomputed positions (stored in
  //      // edge property ep_pos_) in the vertex property vp_pos_;
  //
  //      // Attention! Creating new edges, hence make sure the loop ends correctly.
  //      e_end = _m.edges_end( );
  //      for ( eit = _m.edges_begin( ); eit != e_end; ++eit )
  //        split_edge( _m, *eit );
  //
  //
  //      // Commit changes in topology and reconsitute consistency
  //
  //      // Attention! Creating new faces, hence make sure the loop ends correctly.
  //      f_end = _m.faces_end( );
  //      for ( fit = _m.faces_begin( ); fit != f_end; ++fit )
  //        split_face( _m, *fit );
  //
  //      if ( _update_points ) {
  //        // Commit changes in geometry
  //        for ( vit = _m.vertices_begin( );
  //            vit != _m.vertices_end( ); ++vit ) {
  //          _m.set_point( *vit, _m.property( vp_pos_, *vit ) );
  //        }
  //      }
  //
  //
  //#if defined(_DEBUG) || defined(DEBUG)
  //      // Now we have an consistent mesh!
  //      assert( OpenMesh::Utils::MeshCheckerT<mesh_t>( _m ).check( ) );
  //#endif
  //    }
  //
  //    return true;
  //  }

  //  bool subdivide(
  //      mesh_t & _m,
  //      size_t _n,
  //      const bool _update_points = true
  //      ) {
  //
  //    typename mesh_t::FaceIter fit, f_end;
  //    typename mesh_t::EdgeIter eit, e_end;
  //    typename mesh_t::VertexIter vit;
  //
  //    // Do _n subdivisions
  //    for ( size_t i = 0; i < _n; ++i ) {
  //
  //      if ( _update_points ) {
  //        // compute new positions for old vertices
  //        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
  //          smooth( _m, *vit, nullptr, nullptr, nullptr );
  //        }
  //      } else {
  //        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
  //          rowInd->push_back( vit->idx( ) );
  //          colInd->push_back( vit->idx( ) );
  //          values->push_back( 1.0 );
  //        }
  //      }
  //
  //      // Compute position for new vertices and store them in the edge property
  //      for ( eit = _m.edges_begin( ); eit != _m.edges_end( ); ++eit )
  //        compute_midpoint( _m, *eit, nullptr, nullptr, nullptr, _update_points );
  //
  //      // Split each edge at midpoint and store precomputed positions (stored in
  //      // edge property ep_pos_) in the vertex property vp_pos_;
  //
  //      // Attention! Creating new edges, hence make sure the loop ends correctly.
  //      e_end = _m.edges_end( );
  //      for ( eit = _m.edges_begin( ); eit != e_end; ++eit )
  //        split_edge( _m, *eit );
  //
  //
  //      // Commit changes in topology and reconsitute consistency
  //
  //      // Attention! Creating new faces, hence make sure the loop ends correctly.
  //      f_end = _m.faces_end( );
  //      for ( fit = _m.faces_begin( ); fit != f_end; ++fit )
  //        split_face( _m, *fit );
  //
  //      if ( _update_points ) {
  //        // Commit changes in geometry
  //        for ( vit = _m.vertices_begin( );
  //            vit != _m.vertices_end( ); ++vit ) {
  //          _m.set_point( *vit, _m.property( vp_pos_, *vit ) );
  //        }
  //      }
  //
  //
  //#if defined(_DEBUG) || defined(DEBUG)
  //      // Now we have an consistent mesh!
  //      assert( OpenMesh::Utils::MeshCheckerT<mesh_t>( _m ).check( ) );
  //#endif
  //    }
  //
  //    return true;
  //  }

  bool subdivide(
      mesh_t & _m,
      size_t _n,
      std::vector< bem4i::SparseMatrix< LO, RealType > * > * matrices
      = nullptr,
      const bool _update_points = true
      ) {

    typename mesh_t::FaceIter fit, f_end;
    typename mesh_t::EdgeIter eit, e_end;
    typename mesh_t::VertexIter vit;

    LO nNodes = 0;
    LO nEdges = 0;
    LO nnz = 0;
    std::vector< LO > * rowInd = nullptr;
    std::vector< LO > * colInd = nullptr;
    std::vector< RealType > * values = nullptr;

    // Do _n subdivisions
    for ( size_t i = 0; i < _n; ++i ) {
      if ( matrices ) {
        nNodes = _m.n_vertices( );
        nEdges = _m.n_edges( );
        nnz = 4 * nEdges + 7 * nNodes;
        rowInd = new std::vector< LO >( );
        colInd = new std::vector< LO >( );
        values = new std::vector< RealType >( );
        rowInd->reserve( nnz );
        colInd->reserve( nnz );
        values->reserve( nnz );
      }

      if ( _update_points ) {
        // compute new positions for old vertices
        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
          smooth( _m, *vit, rowInd, colInd, values );
        }
      } else {
        for ( vit = _m.vertices_begin( ); vit != _m.vertices_end( ); ++vit ) {
          if ( matrices ) {
            rowInd->push_back( vit->idx( ) );
            colInd->push_back( vit->idx( ) );
            values->push_back( 1.0 );
          }
        }
      }

      // Compute position for new vertices and store them in the edge property
      for ( eit = _m.edges_begin( ); eit != _m.edges_end( ); ++eit )
        compute_midpoint( _m, *eit, rowInd, colInd, values, _update_points );

      // Split each edge at midpoint and store precomputed positions (stored in
      // edge property ep_pos_) in the vertex property vp_pos_;

      // Attention! Creating new edges, hence make sure the loop ends correctly.
      e_end = _m.edges_end( );
      for ( eit = _m.edges_begin( ); eit != e_end; ++eit )
        split_edge( _m, *eit );


      // Commit changes in topology and reconsitute consistency

      // Attention! Creating new faces, hence make sure the loop ends correctly.
      f_end = _m.faces_end( );
      for ( fit = _m.faces_begin( ); fit != f_end; ++fit )
        split_face( _m, *fit );

      if ( _update_points ) {
        // Commit changes in geometry
        for ( vit = _m.vertices_begin( );
            vit != _m.vertices_end( ); ++vit ) {
          _m.set_point( *vit, _m.property( vp_pos_, *vit ) );
        }
      }


#if defined(_DEBUG) || defined(DEBUG)
      // Now we have an consistent mesh!
      assert( OpenMesh::Utils::MeshCheckerT<mesh_t>( _m ).check( ) );
#endif

      if ( matrices ) {
        ( *matrices )[ i ]->setFromTriplets( nNodes + nEdges, nNodes, *rowInd,
            *colInd, *values );
        delete rowInd;
        delete colInd;
        delete values;
      }
    }

    return true;
  }

private:

  /// Helper functor to compute weights for Loop-subdivision
  /// \internal

  struct compute_weight {

    compute_weight( ) : valence( -1 ) {
    }

    weight_t operator( )(void) {
#if !defined(OM_CC_MIPS)
      using std::cos;
#endif
      //              1
      // alpha(n) = ---- * (40 - ( 3 + 2 cos( 2 Pi / n ) )� )
      //             64

      if ( ++valence ) {
        double inv_v = 1.0 / double( valence );
        double alpha;
        //#define WARREN
#ifdef WARREN
        if ( valence == 3 ) {
          alpha = 0.5625; // 9/16
        } else {
          alpha = 0.375; // 3/8
        }
#else
        double t = ( 3.0 + 2.0 * cos( 2.0 * M_PI * inv_v ) );
        alpha = ( 40.0 - t * t ) / 64.0;
#endif
        return weight_t( 1.0 - alpha, inv_v * alpha );
      }
      return weight_t( 0.0, 0.0 );
    }
    int valence;
  };

private: // topological modifiers

  void split_face(
      mesh_t& _m,
      const typename mesh_t::FaceHandle& _fh
      ) {
    typename mesh_t::HalfedgeHandle
    heh1( _m.halfedge_handle( _fh ) ),
        heh2( _m.next_halfedge_handle( _m.next_halfedge_handle( heh1 ) ) ),
        heh3( _m.next_halfedge_handle( _m.next_halfedge_handle( heh2 ) ) );

    // Cutting off every corner of the 6_gon
    corner_cutting( _m, heh1 );
    corner_cutting( _m, heh2 );
    corner_cutting( _m, heh3 );
  }

  void corner_cutting(
      mesh_t& _m,
      const typename mesh_t::HalfedgeHandle& _he
      ) {
    // Define Halfedge Handles
    typename mesh_t::HalfedgeHandle
    heh1( _he ),
        heh5( heh1 ),
        heh6( _m.next_halfedge_handle( heh1 ) );

    // Cycle around the polygon to find correct Halfedge
    for (; _m.next_halfedge_handle( _m.next_halfedge_handle( heh5 ) ) != heh1;
        heh5 = _m.next_halfedge_handle( heh5 ) ) {
    }

    typename mesh_t::VertexHandle
    vh1 = _m.to_vertex_handle( heh1 ),
        vh2 = _m.to_vertex_handle( heh5 );

    typename mesh_t::HalfedgeHandle
    heh2( _m.next_halfedge_handle( heh5 ) ),
        heh3( _m.new_edge( vh1, vh2 ) ),
        heh4( _m.opposite_halfedge_handle( heh3 ) );

    /* Intermediate result
     *
     *            *
     *         5 /|\
     *          /_  \
     *    vh2> *     *
     *        /|\3   |\
     *       /_  \|4   \
     *      *----\*----\*
     *          1 ^   6
     *            vh1 (adjust_outgoing halfedge!)
     */

    // Old and new Face
    typename mesh_t::FaceHandle fh_old( _m.face_handle( heh6 ) );
    typename mesh_t::FaceHandle fh_new( _m.new_face( ) );


    // Re-Set Handles around old Face
    _m.set_next_halfedge_handle( heh4, heh6 );
    _m.set_next_halfedge_handle( heh5, heh4 );

    _m.set_face_handle( heh4, fh_old );
    _m.set_face_handle( heh5, fh_old );
    _m.set_face_handle( heh6, fh_old );
    _m.set_halfedge_handle( fh_old, heh4 );

    // Re-Set Handles around new Face
    _m.set_next_halfedge_handle( heh1, heh3 );
    _m.set_next_halfedge_handle( heh3, heh2 );

    _m.set_face_handle( heh1, fh_new );
    _m.set_face_handle( heh2, fh_new );
    _m.set_face_handle( heh3, fh_new );

    _m.set_halfedge_handle( fh_new, heh1 );
  }

  void split_edge(
      mesh_t& _m,
      const typename mesh_t::EdgeHandle& _eh
      ) {
    bem4i::EdgeTag etag = _m.data( _eh ).tag( );

    typename mesh_t::HalfedgeHandle
    heh = _m.halfedge_handle( _eh, 0 ),
        opp_heh = _m.halfedge_handle( _eh, 1 );

    typename mesh_t::HalfedgeHandle new_heh, opp_new_heh, t_heh;
    typename mesh_t::VertexHandle vh;
    typename mesh_t::VertexHandle vh1( _m.to_vertex_handle( heh ) );
    typename mesh_t::Point midP( _m.point( _m.to_vertex_handle( heh ) ) );
    midP += _m.point( _m.to_vertex_handle( opp_heh ) );
    midP *= 0.5;

    // new vertex
    vh = _m.new_vertex( midP );

    // memorize position, will be set later
    _m.property( vp_pos_, vh ) = _m.property( ep_pos_, _eh );


    // Re-link mesh entities
    if ( _m.is_boundary( _eh ) ) {
      for ( t_heh = heh;
          _m.next_halfedge_handle( t_heh ) != opp_heh;
          t_heh =
          _m.opposite_halfedge_handle( _m.next_halfedge_handle( t_heh ) ) ) {
      }
    } else {
      for ( t_heh = _m.next_halfedge_handle( opp_heh );
          _m.next_halfedge_handle( t_heh ) != opp_heh;
          t_heh = _m.next_halfedge_handle( t_heh ) ) {
      }
    }

    new_heh = _m.new_edge( vh, vh1 );
    opp_new_heh = _m.opposite_halfedge_handle( new_heh );
    _m.set_vertex_handle( heh, vh );

    _m.set_next_halfedge_handle( t_heh, opp_new_heh );
    _m.set_next_halfedge_handle( new_heh, _m.next_halfedge_handle( heh ) );
    _m.set_next_halfedge_handle( heh, new_heh );
    _m.set_next_halfedge_handle( opp_new_heh, opp_heh );

    if ( _m.face_handle( opp_heh ).is_valid( ) ) {
      _m.set_face_handle( opp_new_heh, _m.face_handle( opp_heh ) );
      _m.set_halfedge_handle( _m.face_handle( opp_new_heh ), opp_new_heh );
    }

    _m.set_face_handle( new_heh, _m.face_handle( heh ) );
    _m.set_halfedge_handle( vh, new_heh );
    _m.set_halfedge_handle( _m.face_handle( heh ), heh );
    _m.set_halfedge_handle( vh1, opp_new_heh );

    // Never forget this, when playing with the topology
    _m.adjust_outgoing_halfedge( vh );
    _m.adjust_outgoing_halfedge( vh1 );

    // set tags
    _m.data( _m.edge_handle( new_heh ) ).tag( etag );
    if ( etag == bem4i::EDGE_CREASE ) {
      _m.data( vh ).tag( bem4i::VERTEX_CREASE );
    } else {
      _m.data( vh ).tag( bem4i::VERTEX_SMOOTH );
    }
  }

private: // geometry helper

  //  void compute_midpoint( mesh_t& _m, const typename mesh_t::EdgeHandle& _eh ) {
  //#define V( X ) vector_cast< typename mesh_t::Normal >( X )
  //    typename mesh_t::HalfedgeHandle heh, opp_heh;
  //
  //    heh = _m.halfedge_handle( _eh, 0 );
  //    opp_heh = _m.halfedge_handle( _eh, 1 );
  //
  //    typename mesh_t::Point
  //    pos( _m.point( _m.to_vertex_handle( heh ) ) );
  //
  //    pos += V( _m.point( _m.to_vertex_handle( opp_heh ) ) );
  //
  //    // boundary edge: just average vertex positions
  //    if ( _m.is_boundary( _eh ) ) {
  //      pos *= 0.5;
  //    } else // inner edge: add neighbouring Vertices to sum
  //    {
  //      pos *= real_t( 3.0 );
  //      pos += V( _m.point(
  //          _m.to_vertex_handle( _m.next_halfedge_handle( heh ) ) ) );
  //      pos += V( _m.point(
  //          _m.to_vertex_handle( _m.next_halfedge_handle( opp_heh ) ) ) );
  //      pos *= _1over8;
  //    }
  //    _m.property( ep_pos_, _eh ) = pos;
  //#undef V
  //  }

  void compute_midpoint(
      mesh_t& _m,
      const typename mesh_t::EdgeHandle& _eh,
      std::vector< LO > * rowInd,
      std::vector< LO > * colInd,
      std::vector< RealType > * values,
      bool _update_points
      ) {

    bem4i::EdgeTag etag = _m.data( _eh ).tag( );

    LO nNodes = _m.n_vertices( );

#define V( X ) vector_cast< typename mesh_t::Normal >( X )
    typename mesh_t::HalfedgeHandle heh, opp_heh;

    heh = _m.halfedge_handle( _eh, 0 );
    opp_heh = _m.halfedge_handle( _eh, 1 );

    typename mesh_t::Point
    pos( _m.point( _m.to_vertex_handle( heh ) ) );

    pos += V( _m.point( _m.to_vertex_handle( opp_heh ) ) );

    // crease edge: just average vertex positions
    if ( etag == bem4i::EDGE_CREASE || !_update_points ) {
      pos *= 0.5;

      if ( rowInd && colInd && values ) {
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle( heh ).idx( ) );
        values->push_back( 0.5 ); // 1/2
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle( opp_heh ).idx( ) );
        values->push_back( 0.5 ); // 1/2
      }
    } else // inner edge: add neighbouring Vertices to sum
    {
      pos *= real_t( 3.0 );
      pos += V( _m.point( _m.to_vertex_handle(
          _m.next_halfedge_handle( heh ) ) ) );
      pos += V( _m.point( _m.to_vertex_handle(
          _m.next_halfedge_handle( opp_heh ) ) ) );
      pos *= _1over8;

      if ( rowInd && colInd && values ) {
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle( heh ).idx( ) );
        values->push_back( 0.375 ); // 3/8
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle( opp_heh ).idx( ) );
        values->push_back( 0.375 ); // 3/8
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle(
            _m.next_halfedge_handle( heh ) ).idx( ) );
        values->push_back( 0.125 ); // 1/8
        rowInd->push_back( nNodes + _eh.idx( ) );
        colInd->push_back( _m.to_vertex_handle(
            _m.next_halfedge_handle( opp_heh ) ).idx( ) );
        values->push_back( 0.125 ); // 1/8
      }
    }
    _m.property( ep_pos_, _eh ) = pos;
#undef V
  }

  //  void smooth( mesh_t& _m, const typename mesh_t::VertexHandle& _vh ) {
  //    typename mesh_t::Point pos( 0.0, 0.0, 0.0 );
  //
  //    if ( _m.is_boundary( _vh ) ) // if boundary: Point 1-6-1
  //    {
  //      typename mesh_t::HalfedgeHandle heh, prev_heh;
  //      heh = _m.halfedge_handle( _vh );
  //
  //      if ( heh.is_valid( ) ) {
  //        assert( _m.is_boundary( _m.edge_handle( heh ) ) );
  //
  //        prev_heh = _m.prev_halfedge_handle( heh );
  //
  //        typename mesh_t::VertexHandle
  //        to_vh = _m.to_vertex_handle( heh ),
  //            from_vh = _m.from_vertex_handle( prev_heh );
  //
  //        // ( v_l + 6 v + v_r ) / 8
  //        pos = _m.point( _vh );
  //        pos *= real_t( 6.0 );
  //        pos += vector_cast< typename mesh_t::Normal > ( _m.point( to_vh ) );
  //        pos += vector_cast< typename mesh_t::Normal > ( _m.point( from_vh ) );
  //        pos *= _1over8;
  //
  //      } else
  //        return;
  //    } else // inner vertex: (1-a) * p + a/n * Sum q, q in one-ring of p
  //    {
  //      typedef typename mesh_t::Normal Vec;
  //      typename mesh_t::VertexVertexIter vvit;
  //      size_t valence( 0 );
  //
  //      // Calculate Valence and sum up neighbour points
  //      for ( vvit = _m.vv_iter( _vh ); vvit.is_valid( ); ++vvit ) {
  //        ++valence;
  //        pos += vector_cast< Vec >( _m.point( *vvit ) );
  //      }
  //      pos *= weights_[valence].second; // alpha(n)/n * Sum q, q in one-ring of p
  //      pos += weights_[valence].first
  //          * vector_cast<Vec>( _m.point( _vh ) ); // + (1-a)*p
  //    }
  //
  //    _m.property( vp_pos_, _vh ) = pos;
  //  }

  void smooth(
      mesh_t& _m,
      const typename mesh_t::VertexHandle& _vh,
      std::vector< LO > * rowInd,
      std::vector< LO > * colInd,
      std::vector< RealType > * values
      ) {

    typename mesh_t::Point pos( 0.0, 0.0, 0.0 );

    bem4i::VertexTag vtag = _m.data( _vh ).tag( );

    if ( vtag == bem4i::VERTEX_CORNER ) { // 1
      pos = _m.point( _vh );
      if ( rowInd && colInd && values ) {
        rowInd->push_back( _vh.idx( ) );
        colInd->push_back( _vh.idx( ) );
        values->push_back( 1.0 );
      }
    } else if ( vtag == bem4i::VERTEX_CREASE ) { // 1-6-1

      typename mesh_t::CVOHIter vh_it = _m.cvoh_iter( _vh );
      bem4i::EdgeTag etag;
      typename mesh_t::VertexHandle v1, v2;

      // find first EDGE_CREASE edge
      for (; vh_it.is_valid( ); ++vh_it ) {
        etag = _m.data( _m.edge_handle( *vh_it ) ).tag( );
        if ( etag == bem4i::EDGE_CREASE ) {
          v1 = _m.to_vertex_handle( *vh_it );
          ++vh_it;
          break;
        }
      }

      // find second EDGE_CREASE edge
      for (; vh_it.is_valid( ); ++vh_it ) {
        etag = _m.data( _m.edge_handle( *vh_it ) ).tag( );
        if ( etag == bem4i::EDGE_CREASE ) {
          v2 = _m.to_vertex_handle( *vh_it );
          break;
        }
      }

      // ( v_l + 6 v + v_r ) / 8
      pos = _m.point( _vh );
      pos *= real_t( 6.0 );
      pos += vector_cast< typename mesh_t::Normal > ( _m.point( v1 ) );
      pos += vector_cast< typename mesh_t::Normal > ( _m.point( v2 ) );
      pos *= _1over8;

      if ( rowInd && colInd && values ) {
        rowInd->push_back( _vh.idx( ) );
        colInd->push_back( _vh.idx( ) );
        values->push_back( 0.75 ); // 6/8
        rowInd->push_back( _vh.idx( ) );
        colInd->push_back( v1.idx( ) );
        values->push_back( 0.125 ); // 1/8
        rowInd->push_back( _vh.idx( ) );
        colInd->push_back( v2.idx( ) );
        values->push_back( 0.125 ); // 1/8
      }

    } else { //if ( vtag == bem4i::VERTEX_SMOOTH ) { // inner vertex: (1-a) * p + a/n * Sum q, q in one-ring of p
      typedef typename mesh_t::Normal Vec;
      typename mesh_t::VertexVertexIter vvit;
      size_t valence( 0 );

      // Calculate valence and sum up neighbour points
      for ( vvit = _m.vv_iter( _vh ); vvit.is_valid( ); ++vvit ) {
        ++valence;
        pos += vector_cast< Vec >( _m.point( *vvit ) );
        if ( rowInd && colInd ) {
          rowInd->push_back( _vh.idx( ) );
          colInd->push_back( vvit->idx( ) );
        }
      }
      if ( values ) {
        for ( size_t i = 0; i < valence; ++i ) {
          values->push_back( weights_[valence].second );
        }
      }
      pos *= weights_[valence].second; // alpha(n)/n * Sum q, q in one-ring of p
      pos += weights_[valence].first
          * vector_cast<Vec>( _m.point( _vh ) ); // + (1-a)*p

      if ( rowInd && colInd && values ) {
        rowInd->push_back( _vh.idx( ) );
        colInd->push_back( _vh.idx( ) );
        values->push_back( weights_[valence].first );
      }
    }

    _m.property( vp_pos_, _vh ) = pos;
  }

private: // data

  OpenMesh::VPropHandleT< typename mesh_t::Point > vp_pos_;
  OpenMesh::EPropHandleT< typename mesh_t::Point > ep_pos_;

  weights_t weights_;

  const real_t _1over8;
  const real_t _3over8;

};


//=============================================================================
} // END_NS_UNIFORM
} // END_NS_SUBDIVIDER
} // END_NS_OPENMESH
//=============================================================================
#endif // OPENMESH_SUBDIVIDER_UNIFORM_COMPOSITELOOPT_HH defined
//=============================================================================