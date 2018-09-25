/*!
 * @file    Offloadable.h
 * @author  Jan Zapletal
 * @date    September 15, 2016
 * @brief   Header file for abstract class Offloadable
 * 
 */

#ifndef OFFLOADABLE_H
#define OFFLOADABLE_H

namespace bem4i {

/*! 
 * Abstract class representing arbitrary class offloadable to MIC
 * 
 */
class Offloadable {
public:

  //! transfer all relevant class data to MIC
  virtual void xferToMIC(
      int device = 0
      ) const = 0;

  //! transfer all relevant class data to host
  virtual void xferToHost(
      int device = 0
      ) const = 0;

  //! update data on MIC, includes xfer but no allocation!
  virtual void updateMIC(
      int device = 0
      ) const = 0;

  //! dealloc all relevant class data to MIC
  virtual void deleteMIC(
      int device = 0
      ) const = 0;

};

}

#endif /* OFFLOADABLE_H */
