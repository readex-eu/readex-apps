/*!
 * @file    Tree.cpp
 * @author  Michal Merta 
 * @date    August 8, 2013
 * 
 */

#ifdef TREE_H

namespace bem4i {

template<class C, class LO>
Tree<C, LO>::Tree( ) {
  root = nullptr;
  pointer = nullptr;
}

template<class C, class LO>
Tree<C, LO>::Tree( const Tree& orig ) {
}

template<class C, class LO>
Tree<C, LO>::~Tree( ) {
  setPointerToRoot( );

  // recursively delete all data in tree
  deleteTree( );
}

template<class C, class LO>
void Tree<C, LO>::createRoot( C data ) {
  root = new TreeMember<C>;
  root->level = 0;
  root->data = data;
  root->parent = nullptr;
  pointer = root;
}

template<class C, class LO>
void Tree<C, LO>::deleteTree( ) {
  for ( int i = 0; i < getNSons( ); i++ ) {
    setPointerToSon( i );
    deleteTree( );
  }
  TreeMember<C>* parent = pointer->parent;
  delete pointer;
  pointer = parent;
}

template<>
void Tree< BECluster< int, double > *, int >::printLeaves(
    std::vector< std::vector< int > > & data,
    TreeMember< BECluster< int, double > *> const * node
    ) {

  typedef int LO;

  LO nSons = node->sons.size( );

  if ( nSons != 0 ) {
    for ( LO i = 0; i < nSons; ++i ) {
      this->printLeaves( data, node->sons[ i ] );
    }
  } else {
    LO size = node->data->elems->size( );
    std::cout << size << ": ";
    std::vector< LO > v;
    data.push_back( v );
    for ( LO i = 0; i < size; ++i ) {
      std::cout << node->data->elems->at( i ) << " ";
      data.back( ).push_back( node->data->elems->at( i ) );
    }
    std::cout << std::endl;
  }
}

template<> 
TreeMember< BECluster< int, double > *> * 
Tree< BECluster< int, double > *, int >::search( int elem ) {
    std::cout << "Write your code bellow" << std::endl;
    // go throw all leaves and find the elem  


    return nullptr;
  }

template<>
void Tree< BECluster< int, double > *, int >::printLeaves(
    std::vector< int > & cluster
    ) {

  typedef int LO;

  std::vector< std::vector< LO > > data;
  this->printLeaves( data, this->root );

  for ( LO i = 0; i < data.size( ); ++i ) {
    for ( LO j = 0; j < data[ i ].size( ); ++j ) {
      cluster[ data[ i ][ j ] ] =  i;
    }
  }
}

template<>
void Tree< BECluster< long, double > *, long >::printLeaves(
    std::vector< std::vector< long > > & data,
    TreeMember< BECluster< long, double > *> const * node
    ) {

  typedef long LO;

  LO nSons = node->sons.size( );

  if ( nSons != 0 ) {
    for ( LO i = 0; i < nSons; ++i ) {
      this->printLeaves( data, node->sons[ i ] );
    }
  } else {
    LO size = node->data->elems->size( );
    //std::cout << size << " ";
    std::vector< LO > v;
    data.push_back( v );
    for ( LO i = 0; i < size; ++i ) {
      //std::cout << node->data->elems->at( i ) << " ";
      data.back( ).push_back( node->data->elems->at( i ) );
    }
    //std::cout << std::endl;
  }
}

template<>
void Tree< BECluster< long, double > *, long >::printLeaves(
    std::vector< long > & cluster
    ) {

  typedef long LO;

  std::vector< std::vector< LO > > data;
  this->printLeaves( data, this->root );

  for ( LO i = 0; i < data.size( ); ++i ) {
    for ( LO j = 0; j < data[ i ].size( ); ++j ) {
      cluster[ data[ i ][ j ] ] =  i;
    }
  }
}

}

#endif
