#include "sitkIMBinaryProcessFilter.h"
#include "sitkIMResampleFilter.h"
#include "sitkCastImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
BinaryProcessFilter::BinaryProcessFilter()
{
  m_Mode = MultMode;
}


//
// ToString
//
std::string BinaryProcessFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: BinaryProcess" << std::endl
      << "  Mode: " << m_Mode
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetMode / SetMode
//
int BinaryProcessFilter::GetMode()
{
  return this->m_Mode;
}
BinaryProcessFilter::Self& BinaryProcessFilter::SetMode(int m)
{
  this->m_Mode = m;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* BinaryProcessFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2,
                                         int mode)
{
  this->SetMode( mode );
  return this->Execute(image1, image2);
}

itk::simple::Image* BinaryProcessFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2)
{
  itk::simple::Image* out = image1;

  // Cast the second image to the same type as the first
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( image1->GetPixelIDValue() );
  std::auto_ptr<itk::simple::Image> image2rsmp( castFilter.Execute( image2 ) );

  // Resample the second image to match the first
  sitkIM::ResampleFilter resampleFilter;
  image2rsmp.reset( resampleFilter.Execute(image2rsmp.get(), image1) );

  // Mult
  if( this->GetMode() == MultMode )
    {
    itk::simple::MultiplyImageFilter filter;
    out = filter.Execute( image1, image2rsmp.get() );
    }
  else
    {
    sitkExceptionMacro("Invalid Mode -> Options: 0 = Mult");
    }


  // return the result
  return out;

}

} // end namespace sitkIM
