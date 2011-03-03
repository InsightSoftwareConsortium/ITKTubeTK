#include "sitkIMAlgorithmReporter.h"
#include "sitkCastImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
AlgorithmReporter::AlgorithmReporter()
{
  m_LowerThreshold = 0;
  m_UpperThreshold = 255;
}


//
// ToString
//
std::string AlgorithmReporter::ToString() const
{
  std::stringstream out;
  out << "Reporter: Masking" << std::endl
      << "  Lower Threshold: " << m_LowerThreshold << std::endl
      << "  Upper Threshold: " << m_UpperThreshold << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetLowerThreshold / SetLowerThreshold
//
float AlgorithmReporter::GetLowerThreshold()
{
  return this->m_LowerThreshold;
}
AlgorithmReporter::Self& AlgorithmReporter::SetLowerThreshold(float lt)
{
  this->m_LowerThreshold = lt;
  return *this;
}

//
// GetUpperThreshold / SetUpperThreshold
//
float AlgorithmReporter::GetUpperThreshold()
{
  return this->m_UpperThreshold;
}
AlgorithmReporter::Self& AlgorithmReporter::SetUpperThreshold( float ut )
{
  this->m_UpperThreshold = ut;
  return *this;
}

//
// GetMode/ SetMode
//
AlgorithmReporter::ModeType AlgorithmReporter::GetMode()
{
  return this->m_Mode;
}
AlgorithmReporter::Self& AlgorithmReporter::SetMode( ModeType m )
{
  this->m_Mode = m;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
double AlgorithmReporter::Execute(itk::simple::Image* image1,
                                                       itk::simple::Image* image2,
                                                       float lowThresh, float upThresh,
                                                       ModeType mode)
{
  this->SetLowerThreshold( lowThresh );
  this->SetUpperThreshold( upThresh );
  this->SetMode( mode );

  return this->Execute(image1, image2);
}

double AlgorithmReporter::Execute(itk::simple::Image* inImage1,
                                  itk::simple::Image* inImage2)
{
  // Cast the image to a double image
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( itk::simple::sitkFloat64 );
  std::auto_ptr<itk::simple::Image> inImage1d( castFilter.Execute( inImage1 ) );
  std::auto_ptr<itk::simple::Image> inImage2d( castFilter.Execute( inImage2 ) );

  // Set up ITK types
  typedef itk::Image<double, 2> Image2DType;
  typedef itk::Image<double, 3> Image3DType;


  // Set up the computation variables  
  double sum = 0;
  double sumS = 0;
  unsigned int count = 0;
  double mean = 0;

  // 2D Case
  if( inImage1d->GetDimension() == 2 && inImage2d->GetDimension() == 2 )
    {
    // Get the ITK images
    Image2DType::Pointer image1 = dynamic_cast <Image2DType*>(
      inImage1d->GetImageBase().GetPointer() );
    Image2DType::Pointer image2 = dynamic_cast <Image2DType*>(
      inImage2d->GetImageBase().GetPointer() );
    if ( image1.IsNull() || image2.IsNull() )
      {
      sitkExceptionMacro( "Unexpected template dispatch error!" );
      }

    // Compute the statistics
    itk::ImageRegionIterator< Image2DType > it1( image1,
          image1->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< Image2DType > it2( image2,
          image2->GetLargestPossibleRegion() );
    it1.GoToBegin();
    it2.GoToBegin();
    while( !it1.IsAtEnd() && !it2.IsAtEnd() )
      {
      double maskV = it2.Get();
      if( maskV >= this->m_LowerThreshold && maskV <= this->m_UpperThreshold )
        {
        sum += it1.Get();
        sumS += it1.Get() * it1.Get();
        ++count;
        }
      ++it1;
      ++it2;
      }
    mean = sum/count;
    }

  // 3D Case
  else if( inImage1d->GetDimension() == 3 && inImage2d->GetDimension() == 3 )
    {
    // Get the ITK images
    Image3DType::Pointer image1 = dynamic_cast <Image3DType*>(
      inImage1d->GetImageBase().GetPointer() );
    Image3DType::Pointer image2 = dynamic_cast <Image3DType*>(
      inImage2d->GetImageBase().GetPointer() );
    if ( image1.IsNull() || image2.IsNull() )
      {
      sitkExceptionMacro( "Unexpected template dispatch error!" );
      }

    // Compute the statistics
    itk::ImageRegionIterator< Image3DType > it1( image1,
          image1->GetLargestPossibleRegion() );
    itk::ImageRegionIterator< Image3DType > it2( image2,
          image2->GetLargestPossibleRegion() );
    it1.GoToBegin();
    it2.GoToBegin();
    while( !it1.IsAtEnd() && !it2.IsAtEnd() )
      {
      double maskV = it2.Get();
      if( maskV >= this->m_LowerThreshold && maskV <= this->m_UpperThreshold )
        {
        sum += it1.Get();
        sumS += it1.Get() * it1.Get();
        ++count;
        }
      ++it1;
      ++it2;
      }
    mean = sum/count;
    }

  // Bad dimension
  else
    {
    sitkExceptionMacro("sitkIM::AlgorithmReporter -> Invalid Image Dimensions");
    }


  // Return the apropriate output
  if( this->m_Mode == MeanMode )
    {
    return mean;
    }
  else if( this->m_Mode == StdDevMode )
    {
    return (sumS - (sum*mean))/(count-1);
    }
  else
    {
    sitkExceptionMacro( "Unknown Mode" );
    }

}

} // end namespace sitkIM
