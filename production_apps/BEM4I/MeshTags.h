/*!
 * @file    MeshTags.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    June 19, 2014
 * @brief   Tags used for subdivision meshes
 * 
 */

#ifndef MESHTAGS_H
#define	MESHTAGS_H

namespace bem4i {

enum VertexTag {

  VERTEX_SMOOTH, // standard subdivision rules apply
  VERTEX_CREASE, // vertex on crease
  VERTEX_CORNER, // fixed vertex
      
  NO_VERTEX_TAGS // Always last!
};

enum EdgeTag {

  EDGE_SMOOTH, // standard subdivision rules apply
  EDGE_CREASE, // creased edge
  
  NO_EDGE_TAGS // Always last!
};

}

#endif	/* MESHTAGS_H */

