#include "sitkIMUnaryProcessFilter.h"
#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
UnaryProcessFilter::UnaryProcessFilter()
{
  m_Mode = 0;
}


//
// ToString
//
std::string UnaryProcessFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: UnaryProcess" << std::endl
      << "  Mode: " << m_Mode
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetMode / SetMode
//
int UnaryProcessFilter::GetMode()
{
  return this->m_Mode;
}
UnaryProcessFilter::Self& UnaryProcessFilter::SetMode(int m)
{
  this->m_Mode = m;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* UnaryProcessFilter::Execute(itk::simple::Image* image,
                                                        int mode)
{
  this->SetMode(mode);
  return this->Execute(image);
}

itk::simple::Image* UnaryProcessFilter::Execute(itk::simple::Image* image)
{
  std::auto_ptr<itk::simple::Image> out;

  // Abs
  if( this->GetMode() == AbsMode )
    {
    itk::simple::AbsImageFilter filter;
    out.reset( filter.Execute( image ) );
    }
  else
    {
    sitkExceptionMacro("Invalid Mode -> Options: 0 = Abs");
    }


  // return the result
  return out.release();
}

} // end namespace sitkIM
