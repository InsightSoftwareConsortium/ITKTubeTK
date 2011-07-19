#include "sitkIMBlurFilter.h"
#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
BlurFilter::BlurFilter()
{
  m_Sigma = 1;
}


//
// ToString
//
std::string BlurFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Blur" << std::endl
      << "  Sigma: " << m_Sigma
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetSigma / SetSigma
//
float BlurFilter::GetSigma()
{
  return this->m_Sigma;
}
BlurFilter::Self& BlurFilter::SetSigma(float s)
{
  this->m_Sigma = s;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* BlurFilter::Execute(itk::simple::Image* image,
                                                float sigma)
{
  this->SetSigma(sigma);
  return this->Execute(image);
}

itk::simple::Image* BlurFilter::Execute(itk::simple::Image* image)
{
  itk::simple::RecursiveGaussianImageFilter filter;
  filter.SetSigma( this->GetSigma() );
  filter.SetOrder( itk::simple::RecursiveGaussianImageFilter::ZeroOrder );

  filter.SetDirection(0);
  std::auto_ptr<itk::simple::Image> out( filter.Execute(image) );
  for (unsigned int i = 1; i < image->GetDimension(); ++i)
    {
    filter.SetDirection( i );
    out.reset( filter.Execute(out.get()) );
    }

  // release the auto pointer and return
  return out.release();
}

} // end namespace sitkIM
