#include "sitkIMBlurOrderFilter.h"
#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
BlurOrderFilter::BlurOrderFilter()
{
  m_Sigma = 1;
  m_Order = 0;
  m_Direction = 0;
}


//
// ToString
//
std::string BlurOrderFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: BlurOrder" << std::endl
      << "  Sigma: " << m_Sigma << std::endl
      << "  Order: " << m_Order << std::endl
      << "  Direction: " << m_Direction
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetSigma / SetSigma
//
float BlurOrderFilter::GetSigma()
{
  return this->m_Sigma;
}
BlurOrderFilter::Self& BlurOrderFilter::SetSigma(float s)
{
  this->m_Sigma = s;
  return *this;
}

//
// GetOrder / SetOrder
//
int BlurOrderFilter::GetOrder()
{
  return this->m_Order;
}
BlurOrderFilter::Self& BlurOrderFilter::SetOrder(int o)
{
  this->m_Order = o;
  return *this;
}

//
// GetDirection / SetDirection
//
int BlurOrderFilter::GetDirection()
{
  return this->m_Direction;
}
BlurOrderFilter::Self& BlurOrderFilter::SetDirection(int d)
{
  this->m_Direction = d;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* BlurOrderFilter::Execute(itk::simple::Image* image,
                                                float sigma,
                                                int order,
                                                int direction)
{
  this->SetSigma(sigma);
  this->SetOrder(order);
  this->SetDirection(direction);
  return this->Execute(image);
}

itk::simple::Image* BlurOrderFilter::Execute(itk::simple::Image* image)
{
  itk::simple::RecursiveGaussianImageFilter filter;
  filter.SetSigma( this->GetSigma() );
  filter.SetDirection( this->GetDirection() );
  switch( this->GetOrder() )
    {
    case 0:
      filter.SetOrder( itk::simple::RecursiveGaussianImageFilter::ZeroOrder );
      break;
    case 1:
      filter.SetOrder( itk::simple::RecursiveGaussianImageFilter::FirstOrder );
      break;
    case 2:
      filter.SetOrder( itk::simple::RecursiveGaussianImageFilter::SecondOrder );
      break;
    default:
      sitkExceptionMacro("Invalid Order value -> Options: 0, 1, 2");
    }


  // run the filter and return
  return filter.Execute( image );
}

} // end namespace sitkIM
