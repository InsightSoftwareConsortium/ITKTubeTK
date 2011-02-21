#include "sitkIMResizeFilter.h"
#include "sitkIMResampleFilter.h"
#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
ResizeFilter::ResizeFilter()
{
  m_Factor = 1;
}


//
// ToString
//
std::string ResizeFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Resize" << std::endl
      << "  Factor: " << m_Factor
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetFactor / SetFactor
//
double ResizeFilter::GetFactor()
{
  return this->m_Factor;
}
ResizeFilter::Self& ResizeFilter::SetFactor(double f)
{
  this->m_Factor = f;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* ResizeFilter::Execute(itk::simple::Image* image,
                                                double f)
{
  this->SetFactor(f);
  return this->Execute(image);
}

itk::simple::Image* ResizeFilter::Execute(itk::simple::Image* image)
{

  std::vector<unsigned int> size =  image->GetSize();
  std::vector<double> spacing =  image->GetSpacing();

  // Standard refactor
  if (this->m_Factor != 0)
    {
    for (unsigned int i = 0; i < image->GetDimension(); ++i)
      {
      size[i] = (unsigned int)((double)size[i] / this->m_Factor);
      spacing[i] *= this->m_Factor;
      }
    }
  // 0 => Make isotropic
  else
    {
    // Calculate mean spacing
    double meanSpacing = 0;
    for (unsigned int i = 0; i < image->GetDimension(); ++i)
      {
      meanSpacing += spacing[i];
      }
    meanSpacing /= image->GetDimension();

    // Set up size and spacing output
    for (unsigned int i = 0; i < image->GetDimension(); ++i)
      {
      double factor = meanSpacing/spacing[i];
      size[i] = (unsigned int)((double)size[i] / factor);
      spacing[i] = meanSpacing;
      }
    }

  // Create temporary image with the right specs
  std::auto_ptr<itk::simple::Image> targetImage;
  if (image->GetDimension() == 2)
    {
    targetImage.reset( new itk::simple::Image(size[0], size[1], itk::simple::sitkUInt8) );
    }
  else
    {
    targetImage.reset( new itk::simple::Image(size[0], size[1], size[2], itk::simple::sitkUInt8) );
    }
  targetImage->SetSpacing(spacing);
  targetImage->SetOrigin( image->GetOrigin() );

  // Resample the image
  sitkIM::ResampleFilter filter;
  itk::simple::Image* out = filter.Execute( image, targetImage.get() );

  // Return result
  return out;

}

} // end namespace sitkIM
