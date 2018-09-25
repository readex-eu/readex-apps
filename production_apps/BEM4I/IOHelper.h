/*!
 * @file    IOHelper.h
 * @author  Michal Merta 
 * @author  Jan Zapletal
 * @date    July 1, 2014
 * @brief   Helper for IO
 * 
 */

#ifndef IOHELPER_H
#define	IOHELPER_H

#include <vector>

namespace bem4i {

template< class LO, class SC >
class IOHelper {

public:

  static void print(
      const std::vector< SC > & v,
      std::ostream & stream = std::cout
      ) {
    
    for ( LO i = 0; i < v.size( ); ++i ) {
      stream << v[ i ] << " ";
    }
    stream << std::endl;
  }

};

}

#endif	/* IOHELPER_H */

