#include "sitkIMGaussianNoiseFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNormalVariateGenerator.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
GaussianNoiseFilter::GaussianNoiseFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_MinValue = 0;
  m_MaxValue = 255;
  m_NoiseMean = 0;
  m_NoiseStdDeviation = 10;
  m_Seed = NO_SEED;
}


//
// ToString
//
std::string GaussianNoiseFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: GaussianNoise" << std::endl
      << "  Min Value: " << m_MinValue << std::endl
      << "  Max Value: " << m_MaxValue << std::endl
      << "  Noise Mean: " << m_NoiseMean << std::endl
      << "  Noise Standard Deviation: " << m_NoiseStdDeviation << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetMinValue / SetMinValue
//
float GaussianNoiseFilter::GetMinValue()
{
  return this->m_MinValue;
}
GaussianNoiseFilter::Self& GaussianNoiseFilter::SetMinValue(float minv)
{
  this->m_MinValue = minv;
  return *this;
}

//
// GetMaxValue / SetMaxValue
//
float GaussianNoiseFilter::GetMaxValue()
{
  return this->m_MaxValue;
}
GaussianNoiseFilter::Self& GaussianNoiseFilter::SetMaxValue(float maxv)
{
  this->m_MaxValue = maxv;
  return *this;
}

//
// GetNoiseMean / SetNoiseMean
//
float GaussianNoiseFilter::GetNoiseMean()
{
  return this->m_NoiseMean;
}
GaussianNoiseFilter::Self& GaussianNoiseFilter::SetNoiseMean(float m)
{
  this->m_NoiseMean = m;
  return *this;
}

//
// GetNoiseStdDeviation / SetNoiseStdDeviation
//
float GaussianNoiseFilter::GetNoiseStdDeviation()
{
  return this->m_NoiseStdDeviation;
}
GaussianNoiseFilter::Self& GaussianNoiseFilter::SetNoiseStdDeviation(float s)
{
  this->m_NoiseStdDeviation = s;
  return *this;
}

//
// GetSeed / SetSeed
//
int GaussianNoiseFilter::GetSeed()
{
  return this->m_Seed;
}
GaussianNoiseFilter::Self& GaussianNoiseFilter::SetSeed(int s)
{
  this->m_Seed = s;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* GaussianNoiseFilter::Execute(itk::simple::Image* image,
                                                 float minv, float maxv,
                                                 float mean, float std, int seed)
{
  this->SetMinValue( minv );
  this->SetMaxValue( maxv );
  this->SetNoiseMean( mean );
  this->SetNoiseStdDeviation( std );
  this->SetSeed( seed );

  return this->Execute( image );
}

//
// Execute
//
itk::simple::Image* GaussianNoiseFilter::Execute(
                                          itk::simple::Image* image)
{
  // Prep to call ExecuteInternal
  itk::simple::PixelIDValueType type = image->GetPixelIDValue();
  unsigned int dimension = image->GetDimension();

  // Dispatch the proper member function
  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image );

}


//----------------------------------------------------------------------------
// ExecuteInternal
//
template <class TImageType>
itk::simple::Image* GaussianNoiseFilter::ExecuteInternal(
                                  itk::simple::Image* inImage )
  {
  // Set up ITK image types
  typedef TImageType     InputImageType;
  typedef InputImageType OutputImageType;

  // Cast the images and make sure they are the correct type
  typename InputImageType::Pointer image =
    dynamic_cast <InputImageType*> ( inImage->GetImageBase().GetPointer() );
  if( image.IsNull() )
    {
    sitkExceptionMacro( "Unexpected template dispatch error!" );
    }

  // Set up gaussian generator
  typedef itk::Statistics::NormalVariateGenerator gaussGenType;
  typename gaussGenType::Pointer gaussGen = gaussGenType::New();

  if (this->m_Seed != NO_SEED)
    {
    gaussGen->Initialize( this->m_Seed );
    }

  // Perform the guts of the filter
  itk::ImageRegionIterator< InputImageType > it2( image, image->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    if( tf >= this->m_MinValue && tf <= this->m_MaxValue )
      {
      tf += gaussGen->GetVariate() * m_NoiseStdDeviation + m_NoiseMean;
      it2.Set( ( typename InputImageType::PixelType )tf );
      }
    ++it2;
    }

  // Return the resulting SimpleITK image
  return new itk::simple::Image( image );
  }

} // end namespace sitkIM
