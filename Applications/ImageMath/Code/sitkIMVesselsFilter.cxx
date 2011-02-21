#include "sitkIMVesselsFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

// tubetk includes
#include "itkTubeRidgeExtractor.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
VesselsFilter::VesselsFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_NumScales = 1;
  m_ScaleMin = 1;
  m_ScaleMax = 1;
}


//
// ToString
//
std::string VesselsFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Vessels" << std::endl
      << "  Scale Min: " << m_ScaleMin << std::endl
      << "  Scale Max: " << m_ScaleMax << std::endl
      << "  Num Scales: " << m_NumScales << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNumScales / SetNumScales
//
unsigned int VesselsFilter::GetNumScales()
{
  return this->m_NumScales;
}
VesselsFilter::Self& VesselsFilter::SetNumScales(unsigned int n)
{
  this->m_NumScales = n;
  return *this;
}

//
// GetScaleMin / SetScaleMin
//
double VesselsFilter::GetScaleMin()
{
  return this->m_ScaleMin;
}
VesselsFilter::Self& VesselsFilter::SetScaleMin(double min)
{
  this->m_ScaleMin = min;
  return *this;
}

//
// GetScaleMax / SetScaleMax
//
double VesselsFilter::GetScaleMax()
{
  return this->m_ScaleMax;
}
VesselsFilter::Self& VesselsFilter::SetScaleMax(double max)
{
  this->m_ScaleMax = max;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* VesselsFilter::Execute(itk::simple::Image* image,
                                                        double min, double max,
                                                        unsigned int n)
{
  this->SetNumScales( n );
  this->SetScaleMin( min );
  this->SetScaleMax( max );

  return this->Execute( image );
}

//
// Execute (Filter version)
//
itk::simple::Image* VesselsFilter::Execute( itk::simple::Image* image )
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
itk::simple::Image* VesselsFilter::ExecuteInternal(
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
  double logScaleStep = (vcl_log(this->m_ScaleMax) - vcl_log(this->m_ScaleMin))
    / (double)(this->m_NumScales-1);

  typedef itk::tube::RidgeExtractor< InputImageType > RidgeFuncType;
  typename RidgeFuncType::Pointer imFunc = RidgeFuncType::New();
  imFunc->SetInputImage( image );

  typename InputImageType::Pointer image2 = InputImageType::New();
  image2->SetRegions( image->GetLargestPossibleRegion() );
  image2->SetOrigin( image->GetOrigin() );
  image2->SetSpacing( image->GetSpacing() );
  image2->Allocate();

  itk::ImageRegionIteratorWithIndex< InputImageType > it1( image,
        image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< InputImageType > it2( image2,
        image2->GetLargestPossibleRegion() );

  double ridgeness = 0;
  double roundness = 0;
  double curvature = 0;
  double scale = this->m_ScaleMin;
  imFunc->SetScale( scale );
  std::cout << "   Processing scale " << scale << std::endl;
  it1.GoToBegin();
  it2.GoToBegin();
  typename RidgeFuncType::ContinuousIndexType cIndx;
  while( !it1.IsAtEnd() )
    {
    for( unsigned int d=0; d<InputImageType::ImageDimension; ++d )
      {
      cIndx[d] = it1.GetIndex()[d];
      }
    ridgeness = imFunc->Ridgeness( cIndx, roundness, curvature  );
    it2.Set( ( typename InputImageType::PixelType )ridgeness );
    ++it1;
    ++it2;
    }
  for( unsigned int i=1; i<this->m_NumScales; i++ )
    {
    scale = vcl_exp(vcl_log(this->m_ScaleMin) + i * logScaleStep);
    imFunc->SetScale( scale );
    std::cout << "   Processing scale " << scale << std::endl;
    it1.GoToBegin();
    it2.GoToBegin();
    while( !it1.IsAtEnd() )
      {
      for( unsigned int d=0; d<InputImageType::ImageDimension; ++d )
        {
        cIndx[d] = it1.GetIndex()[d];
        }
      ridgeness = imFunc->Ridgeness( cIndx, roundness, curvature  );
      if( ridgeness > it2.Get() )
        {
        it2.Set( ( typename InputImageType::PixelType )ridgeness );
        }
      ++it1;
      ++it2;
      }
    }
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    it1.Set( it2.Get() );
    ++it1;
    ++it2;
    }

  // return the result
  return new itk::simple::Image( image );

  }

} // end namespace sitkIM
