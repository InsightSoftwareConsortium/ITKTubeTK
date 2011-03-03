#include "sitkIMMaskingFilter.h"
#include "sitkIMResampleFilter.h"
#include "sitkCastImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
MaskingFilter::MaskingFilter()
{
  m_LowerThreshold = 1;
  m_UpperThreshold = 255;
  m_OutsideValue = 0;
}


//
// ToString
//
std::string MaskingFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Masking" << std::endl
      << "  Lower Threshold: " << m_LowerThreshold << std::endl
      << "  Upper Threshold: " << m_UpperThreshold << std::endl
      << "  Outside Value: " << m_OutsideValue
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetLowerThreshold / SetLowerThreshold
//
float MaskingFilter::GetLowerThreshold()
{
  return this->m_LowerThreshold;
}
MaskingFilter::Self& MaskingFilter::SetLowerThreshold(float lt)
{
  this->m_LowerThreshold = lt;
  return *this;
}

//
// GetUpperThreshold / SetUpperThreshold
//
float MaskingFilter::GetUpperThreshold()
{
  return this->m_UpperThreshold;
}
MaskingFilter::Self& MaskingFilter::SetUpperThreshold( float ut )
{
  this->m_UpperThreshold = ut;
  return *this;
}

//
// GetOutsideValue / SetOutsideValue
//
float MaskingFilter::GetOutsideValue()
{
  return this->m_OutsideValue;
}
MaskingFilter::Self& MaskingFilter::SetOutsideValue(float ov)
{
  this->m_OutsideValue = ov;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* MaskingFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2,
                                         float lowThresh, float upThresh, float valFalse)
{
  this->SetLowerThreshold( lowThresh );
  this->SetUpperThreshold( upThresh );
  this->SetOutsideValue( valFalse );

  return this->Execute(image1, image2);
}

itk::simple::Image* MaskingFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2)
{
  // Cast the second image to the same type as the first
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( image1->GetPixelIDValue() );
  std::auto_ptr<itk::simple::Image> maskImage( castFilter.Execute( image2 ) );

  // Resample the second image to match the first
  sitkIM::ResampleFilter resampleFilter;
  maskImage.reset( resampleFilter.Execute(maskImage.get(), image1) );

  // Threshold the second image
  itk::simple::BinaryThresholdImageFilter thresholdFilter;
  thresholdFilter.SetLowerThreshold( this->GetLowerThreshold() );
  thresholdFilter.SetUpperThreshold( this->GetUpperThreshold() );
  thresholdFilter.SetOutsideValue( 0 );
  thresholdFilter.SetInsideValue( 1 );
  maskImage.reset( thresholdFilter.Execute( maskImage.get() ) );

  // Mask the first image
  itk::simple::MaskImageFilter maskFilter;
  maskFilter.SetOutsideValue( this->GetOutsideValue() );
  return maskFilter.Execute( image1, maskImage.get() );

}

} // end namespace sitkIM
