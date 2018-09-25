/*!
 * @file    Tree.h
 * @author  Michal Merta 
 * @date    August 8, 2013
 * @brief   Header file for class Tree
 * 
 */

#ifndef TREE_H
#define	TREE_H

#include "Mesh.h"
#include <vector>
#include <iostream>
#include <limits.h>

namespace bem4i {

/*  
template<class C> class TreeMember;

template<class C> void DeleteFunc (C &c);
  
template<class C>
void TreeMemberDeleteFunc (TreeMember<C> * & treeMember)
{
  if (treeMember)
    if (treeMember->data)
      DeleteFunc<C>(treeMember->data);
};

 */

template<class C>
class TreeMember {
public:
  int level;
  C data;
  TreeMember<C> *parent;
  std::vector<TreeMember<C>*> sons;

  TreeMember( ) {
  };

  ~TreeMember( ) {
    delete data;
  }
};

//! class represents arbitrary N-ary tree

template<class C, class LO>
class Tree {
public:

  //! default constructor
  Tree( );

  //! copy constructor
  Tree( const Tree& orig );

  //! destructor
  virtual ~Tree( );

  //! creates a root entry with given data
  void createRoot( C data );

  //! sets pointer to the root

  inline void setPointerToRoot( ) {
    pointer = root;
  };

  //! sets pointer to the i-th son of current TreeMember

  inline void setPointerToSon( LO i ) {
    pointer = pointer->sons[i];
  }

  //! sets pointer to parent of current tree member

  inline void setPointerToParent( ) {
    pointer = pointer->parent;
  }

  //! returns current level in the tree

  inline LO getLevel( ) const {
    return pointer->level;
  }

  //! returns data of current TreeMember

  inline C get( ) const {
    return pointer->data;
  }

  //! returns number of sons of current TreeMember

  inline LO getNSons( ) const {
    return pointer->sons.size( );
  }

  //! adds son to current TreeMember

  inline void addSon( C data ) {
    TreeMember<C> *son = new TreeMember<C>;
    son->level = pointer->level + 1;
    son->data = data;
    son->parent = pointer;
    pointer->sons.push_back( son );
  }

  // method for "external" tree climbing

  inline TreeMember<C>* getRootNode( ) const {
    return root;
  }

  inline TreeMember<C>* getSonNode( TreeMember<C> * const node, LO i ) const {
    return node->sons[i];
  }

  inline TreeMember<C>* getParentNode( TreeMember<C> * const node ) const {
    return node->parent;
  }

  inline C get( TreeMember<C> const * node ) const {
    return node->data;
  }

  inline LO getNSons( TreeMember<C> const * node ) const {
    return node->sons.size( );
  }

  inline void addSon( TreeMember<C>* node, C data ) {
    TreeMember<C> *son = new TreeMember<C>;
    son->level = pointer->level + 1;
    son->data = data;
    son->parent = pointer;
    node->sons.push_back( son );
  }

  inline LO getNodeLevel(
      TreeMember<C> const * node
      ) const {
    return node->level;
  }

  //! returns minimal depth in the tree (of the shortest branch)  

  inline LO getMinDepth( TreeMember<C> const * node ) const {
    if( node == nullptr ) return INT_MAX;
    if( node->sons.size() == 0 ) return 0;
    return 1 + std::min(this->getMinDepth(node->sons[0]), 
        this->getMinDepth(node->sons[1]));
  }

  //! returns minimal depth in the tree (of the shortest branch)  

  inline LO getMaxDepth( TreeMember<C> const * node ) const {
    if( node == nullptr ) return INT_MAX;
    if( node->sons.size() == 0 ) return 0;
    return 1 + std::max(this->getMaxDepth(node->sons[0]), 
        this->getMaxDepth(node->sons[1]));
  }

  TreeMember<C>* search( LO elem ) {
    std::cout << "Only in specialized!" << std::endl;
  }

  

  void printLeaves(
      std::vector< LO > & cluster
      ) {
    std::cout << "Only in specialized!" << std::endl;
  }

private:

  //! root of the tree
  TreeMember<C>* root;

  //! current position of climber in tree
  TreeMember<C>* pointer;

  //! deletes the whole tree
  void deleteTree( );

  TreeMember<C>* search( LO elem, TreeMember< C > const * node ) {
    std::cout << "Only in specialized!" << std::endl;
    return nullptr;
  }

  void printLeaves(
      std::vector< std::vector< LO > > & data,
      TreeMember< C > const * node
      ) {
    std::cout << "Only in specialized!" << std::endl;
  }

};

}

// include .cpp file to overcome linking problems due to templates
#include "Tree.cpp"

#endif	/* TREE_H */
