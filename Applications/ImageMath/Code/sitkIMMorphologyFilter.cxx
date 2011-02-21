#include "sitkIMMorphologyFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
MorphologyFilter::MorphologyFilter()
{
  m_ForegroundValue = 1;
  m_BackgroundValue = 0;
  m_Radius = 1;
  m_Operation = Erode;
}


//
// ToString
//
std::string MorphologyFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Morphology"
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetOperation / SetOperation
//
MorphologyFilter::OperationType MorphologyFilter::GetOperation()
{
  return this->m_Operation;
}
MorphologyFilter::Self& MorphologyFilter::SetOperation(MorphologyFilter::OperationType op)
{
  this->m_Operation = op;
  return *this;
}

//
// GetForegroundValue / SetForegroundValue
//
float MorphologyFilter::GetForegroundValue()
{
  return this->m_ForegroundValue;
}
MorphologyFilter::Self& MorphologyFilter::SetForegroundValue(float fg)
{
  this->m_ForegroundValue = fg;
  return *this;
}

//
// GetBackgroundValue / SetBackgroundValue
//
float MorphologyFilter::GetBackgroundValue()
{
  return this->m_BackgroundValue;
}
MorphologyFilter::Self& MorphologyFilter::SetBackgroundValue(float bg)
{
  this->m_BackgroundValue = bg;
  return *this;
}

//
// GetRadius / SetRadius
//
float MorphologyFilter::GetRadius()
{
  return this->m_Radius;
}
MorphologyFilter::Self& MorphologyFilter::SetRadius(float r)
{
  this->m_Radius = r;
  return *this;
}



//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* MorphologyFilter::Execute(itk::simple::Image* image,
                                                OperationType op,
                                                float fg,
                                                float bg,
                                                float r)
{
  this->SetOperation(op);
  this->SetForegroundValue(fg);
  this->SetBackgroundValue(bg);
  this->SetRadius(r);

  return this->Execute(image);
}

itk::simple::Image* MorphologyFilter::Execute(itk::simple::Image* image)
{
  // Erode
  if (this->GetOperation() == Erode)
    {
    itk::simple::ErodeObjectMorphologyImageFilter filter;

    filter.SetKernelType( itk::simple::ErodeObjectMorphologyImageFilter::BinaryBallKernel );
    filter.SetObjectValue( this->GetForegroundValue() );
    filter.SetBackgroundValue( this->GetBackgroundValue() );
    filter.SetKernelRadius( 1 );

    std::auto_ptr<itk::simple::Image> out( filter.Execute(image) );
    for( int r = 1; r < this->GetRadius(); ++r)
      {
      out.reset( filter.Execute(out.get()) );
      }

    return out.release();
    }
  // Dilate
  else if (this->GetOperation() == Dilate)
    {
    itk::simple::DilateObjectMorphologyImageFilter filter;

    filter.SetKernelType( itk::simple::DilateObjectMorphologyImageFilter::BinaryBallKernel );
    filter.SetObjectValue( this->GetForegroundValue() );
    filter.SetKernelRadius( 1 );

    std::auto_ptr<itk::simple::Image> out( filter.Execute(image) );
    for( int r = 1; r < this->GetRadius(); ++r)
      {
      out.reset( filter.Execute(out.get()) );
      }

    return out.release();
    }
  // Unknown
  else
    {
    sitkExceptionMacro("Unknown Option Type");
    }
}

} // end namespace sitkIM
