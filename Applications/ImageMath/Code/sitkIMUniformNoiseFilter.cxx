#include "sitkIMUniformNoiseFilter.h"
#include "sitkIMGaussianNoiseFilter.h"

#include "itkImageRegionIterator.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
UniformNoiseFilter::UniformNoiseFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_MinValue = 0;
  m_MaxValue = 255;
  m_NoiseMean = 0;
  m_NoiseRange = 10;
  m_Seed = GaussianNoiseFilter::NO_SEED;
}


//
// ToString
//
std::string UniformNoiseFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: UniformNoise" << std::endl
      << "  Min Value: " << m_MinValue << std::endl
      << "  Max Value: " << m_MaxValue << std::endl
      << "  Noise Mean: " << m_NoiseMean << std::endl
      << "  Noise Range: " << m_NoiseRange << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetMinValue / SetMinValue
//
float UniformNoiseFilter::GetMinValue()
{
  return this->m_MinValue;
}
UniformNoiseFilter::Self& UniformNoiseFilter::SetMinValue(float minv)
{
  this->m_MinValue = minv;
  return *this;
}

//
// GetMaxValue / SetMaxValue
//
float UniformNoiseFilter::GetMaxValue()
{
  return this->m_MaxValue;
}
UniformNoiseFilter::Self& UniformNoiseFilter::SetMaxValue(float maxv)
{
  this->m_MaxValue = maxv;
  return *this;
}

//
// GetNoiseMean / SetNoiseMean
//
float UniformNoiseFilter::GetNoiseMean()
{
  return this->m_NoiseMean;
}
UniformNoiseFilter::Self& UniformNoiseFilter::SetNoiseMean(float m)
{
  this->m_NoiseMean = m;
  return *this;
}

//
// GetNoiseRange / SetNoiseRange
//
float UniformNoiseFilter::GetNoiseRange()
{
  return this->m_NoiseRange;
}
UniformNoiseFilter::Self& UniformNoiseFilter::SetNoiseRange(float r)
{
  this->m_NoiseRange = r;
  return *this;
}

//
// GetSeed / SetSeed
//
int UniformNoiseFilter::GetSeed()
{
  return this->m_Seed;
}
UniformNoiseFilter::Self& UniformNoiseFilter::SetSeed(int s)
{
  this->m_Seed = s;
  return *this;
}

//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* UniformNoiseFilter::Execute(
                                          itk::simple::Image* image,
                                          float minv, float maxv,
                                          float mean, float range, int seed)
{
  this->SetMinValue( minv );
  this->SetMaxValue( maxv );
  this->SetNoiseMean( mean );
  this->SetNoiseRange( range );
  this->SetSeed( seed );

  return this->Execute( image );
}

//
// Execute
//
itk::simple::Image* UniformNoiseFilter::Execute(
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
itk::simple::Image* UniformNoiseFilter::ExecuteInternal(
                                                itk::simple::Image* inImage )
  {
  // Set up ITK image types
  typedef TImageType     InputImageType;
  typedef InputImageType OutputImageType;

  // Cast the images and make sure they are the correct type
  typename InputImageType::Pointer image =
    dynamic_cast <InputImageType*> ( inImage->GetImageBase().GetPointer() );
  if ( image.IsNull() )
    {
    sitkExceptionMacro( "Unexpected template dispatch error!" );
    }

  // Perform the guts of the filter
  if (this->m_Seed != GaussianNoiseFilter::NO_SEED)
    {
    srand( (unsigned int)this->m_Seed );
    }
  itk::ImageRegionIterator< InputImageType > it2( image, image->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    if( tf >= this->m_MinValue && tf <= this->m_MaxValue )
      {
      tf += ( 2.0*( rand()/( double )RAND_MAX )-1 ) * this->m_NoiseRange
            + this->m_NoiseMean;
      it2.Set( ( typename InputImageType::PixelType )tf );
      }
    ++it2;
    }

  // Return the resulting SimpleITK image
  return new itk::simple::Image( image );
  }

} // end namespace sitkIM
