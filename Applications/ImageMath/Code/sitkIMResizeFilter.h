#ifndef __sitkIMResizeFilter_h
#define __sitkIMResizeFilter_h

#include "SimpleITK.h"
#include "sitkIMFilter.h"


namespace sitkIM
{

class ResizeFilter : public Filter
{
public:
  typedef ResizeFilter Self;

  /** Default Constructor */
  ResizeFilter();


  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Factor */
  Self& SetFactor(double f);
  double GetFactor();


  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image, double f);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Resize factor */
  double m_Factor;

};

} // end namespace sitkIM
#endif
