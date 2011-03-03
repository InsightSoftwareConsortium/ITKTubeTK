#ifndef __sitkIMBinaryProcessFilter_h
#define __sitkIMBinaryProcessFilter_h

#include "SimpleITK.h"

#include "sitkIMDualFilter.h"


namespace sitkIM
{

class BinaryProcessFilter : public DualFilter
{
public:
  typedef BinaryProcessFilter Self;

  /** Default Constructor */
  BinaryProcessFilter();

  /** Mode list typedef */
  typedef enum {MultMode} ModeType;

  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Mode */
  Self& SetMode(int m);
  int GetMode();

  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image1,
                              itk::simple::Image* image2,
                              int mode);
  itk::simple::Image* Execute(itk::simple::Image* image1,
                              itk::simple::Image* image2);

private:

  /** Process mode flag */
  int m_Mode;

};

} // end namespace sitkIM
#endif
