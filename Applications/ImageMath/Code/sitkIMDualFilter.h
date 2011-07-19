#ifndef __sitkIMDualFilter_h
#define __sitkIMDualFilter_h

#include "sitkImage.h"

#include <string>

namespace sitkIM
{
class DualFilter
{
public:
  typedef DualFilter Self;

  /**
   * Default Constructor that takes no arguments and initializes
   * default parameters
   */
  DualFilter(){};

  /**
   * Default Destructor
   */
  ~DualFilter(){};

  // Print ourselves out
  virtual std::string ToString() const = 0;

  virtual itk::simple::Image* Execute ( itk::simple::Image*,
                                        itk::simple::Image* ) = 0;

private:

};


} // end namespace sitkIM
#endif
