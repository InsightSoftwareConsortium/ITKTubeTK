#ifndef __sitkIMMaskingFilter_h
#define __sitkIMMaskingFilter_h

#include "SimpleITK.h"

#include "sitkIMDualFilter.h"


namespace sitkIM
{

class MaskingFilter : public DualFilter
{
public:
  typedef MaskingFilter Self;

  /** Default Constructor */
  MaskingFilter();


  /** Print outselves out */
  std::string ToString() const;


  /** Get/Set OutsideValue */
  Self& SetOutsideValue(float ov);
  float GetOutsideValue();

  /** Get/Set LowerThreshold */
  Self& SetLowerThreshold(float lt);
  float GetLowerThreshold();

  /** Get/Set UpperThreshold */
  Self& SetUpperThreshold(float ut);
  float GetUpperThreshold();


  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image1,
                              itk::simple::Image* image2,
                              float lowThresh, float upThresh, float valFalse);
  itk::simple::Image* Execute(itk::simple::Image* image1,
                              itk::simple::Image* image2);

private:

  /** Lower threshold value */
  float m_LowerThreshold;

  /** Upper threshold value */
  float m_UpperThreshold;

  /** Value applied to pisels outside the threshold range */
  float m_OutsideValue;

};

} // end namespace sitkIM
#endif
