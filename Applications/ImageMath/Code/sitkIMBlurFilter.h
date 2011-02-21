#ifndef __sitkIMBlurFilter_h
#define __sitkIMBlurFilter_h

#include "SimpleITK.h"
#include "sitkIMFilter.h"


namespace sitkIM
{

class BlurFilter : public Filter
{
public:
  typedef BlurFilter Self;

  /** Default Constructor */
  BlurFilter();


  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Sigma */
  Self& SetSigma(float s);
  float GetSigma();


  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image, float sigma);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Gaussian filter sigma value */
  float m_Sigma;

};

} // end namespace sitkIM
#endif
