#include "sitkIMAddFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
AddFilter::AddFilter()
{
  m_Weight1 = 1;
  m_Weight2 = 1;
}


//
// ToString
//
std::string AddFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Add" << std::endl
      << "  Weight1: " << m_Weight1
      << "  Weight2: " << m_Weight2
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetWeight1 / SetWeight1
//
float AddFilter::GetWeight1()
{
  return this->m_Weight1;
}
AddFilter::Self& AddFilter::SetWeight1(float w)
{
  this->m_Weight1 = w;
  return *this;
}

//
// GetWeight2 / SetWeight2
//
float AddFilter::GetWeight2()
{
  return this->m_Weight2;
}
AddFilter::Self& AddFilter::SetWeight2(float w)
{
  this->m_Weight2 = w;
  return *this;
}



//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* AddFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2,
                                         float weight1, float weight2)
{
  this->SetWeight1(weight1);
  this->SetWeight2(weight2);

  return this->Execute(image1, image2);
}

itk::simple::Image* AddFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2)
{
  // Apply the weights
  itk::simple::MultiplyByConstantImageFilter multFilter;
  multFilter.SetConstant( this->GetWeight1() );
  std::auto_ptr<itk::simple::Image> image1w( multFilter.Execute( image1 ) );
  multFilter.SetConstant( this->GetWeight2() );
  std::auto_ptr<itk::simple::Image> image2w( multFilter.Execute( image2 ) );

  // Add the images
  itk::simple::AddImageFilter addFilter;
  itk::simple::Image* out = addFilter.Execute( image1w.get(), image2w.get() );

  // Return
  return out;

}

} // end namespace sitkIM
