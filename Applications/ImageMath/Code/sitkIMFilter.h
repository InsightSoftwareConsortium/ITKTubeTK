#ifndef __sitkIMFilter_h
#define __sitkIMFilter_h

#include "sitkImage.h"

#include <string>

namespace sitkIM
{
class Filter
{
public:
  typedef Filter Self;

  /**
   * Default Constructor that takes no arguments and initializes
   * default parameters
   */
  Filter(){};

  /**
   * Default Destructor
   */
  ~Filter(){};

  // Print ourselves out
  virtual std::string ToString() const = 0;

  virtual itk::simple::Image* Execute ( itk::simple::Image* ) = 0;

private:

};


} // end namespace sitkIM
#endif
