#include "sitkIMIntensityMultFilter.h"
#include "sitkIMResampleFilter.h"
#include "sitkCastImageFilter.h"

#include "itkImageRegionIterator.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
IntensityMultFilter::IntensityMultFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
}


//
// ToString
//
std::string IntensityMultFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: IntensityMultiplicativeCorrect"
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* IntensityMultFilter::Execute(
                              itk::simple::Image* image1,
                              itk::simple::Image* image2)
{

  // Cast the second image to the same type as the first
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( image1->GetPixelIDValue() );
  itk::simple::Image* image2rsmp = castFilter.Execute( image2 );

  // Resample the second image to match the first
  sitkIM::ResampleFilter resampleFilter;
  image2rsmp = resampleFilter.Execute(image2rsmp, image1);

  // Prep to call ExecuteInternal
  itk::simple::PixelIDValueType type = image1->GetPixelIDValue();
  unsigned int dimension = image1->GetDimension();

  // todo need better error handling and potential type conversion
  if( type != image2rsmp->GetPixelIDValue() ||
      dimension != image2rsmp->GetDimension() ||
      image1->GetWidth() != image2rsmp->GetWidth() ||
      image1->GetHeight() != image2rsmp->GetHeight() ||
      image1->GetDepth() != image2rsmp->GetDepth() )
    {
    sitkExceptionMacro( "Both image for IntensityMultFilter don't match type or dimension!" );
    }

  // Dispatch the proper member function
  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image1 , image2rsmp );

}


//----------------------------------------------------------------------------
// ExecuteInternal
//
template <class TImageType>
itk::simple::Image* IntensityMultFilter::ExecuteInternal(
                                itk::simple::Image* inImage1,
                                itk::simple::Image* inImage2 )
  {
  // Set up ITK image types
  typedef TImageType     InputImageType;
  typedef InputImageType OutputImageType;

  // Cast the images and make sure they are the correct type
  typename InputImageType::Pointer image1 =
    dynamic_cast <InputImageType*> ( inImage1->GetImageBase().GetPointer() );
  typename InputImageType::Pointer image2 =
    dynamic_cast <InputImageType*> ( inImage2->GetImageBase().GetPointer() );
  if ( image1.IsNull() || image2.IsNull() )
    {
    sitkExceptionMacro( "Unexpected template dispatch error!" );
    }

  // Perform the guts of the filter
  itk::ImageRegionIterator< InputImageType > it2( image2, image2->GetLargestPossibleRegion() );
  int count = 0;
  double mean = 0;
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    mean += tf;
    if( tf != 0 )
      {
      ++count;
      }
    ++it2;
    }
  mean /= count;
  itk::ImageRegionIterator< InputImageType > it3( image1, image1->GetLargestPossibleRegion() );
  it3.GoToBegin();
  it2.GoToBegin();
  while( !it3.IsAtEnd() )
    {
    double tf = it3.Get();
    double tf2 = it2.Get();
    if( tf2 != 0 )
      {
      double alpha = mean / tf2;
      tf = tf * alpha;
      it3.Set( ( typename InputImageType::PixelType )tf );
      }
    ++it3;
    ++it2;
    }

  // Return the resulting SimpleITK image
  return new itk::simple::Image( image1 );
  }

} // end namespace sitkIM
