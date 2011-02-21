#include "sitkIMFuseFilter.h"
#include "sitkIMResampleFilter.h"
#include "sitkCastImageFilter.h"

#include "itkImageRegionIterator.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
FuseFilter::FuseFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_Offset = 1;
}


//
// ToString
//
std::string FuseFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: FuseiplicativeCorrect" << std::endl
      << "  Offset: " << m_Offset
      << std::endl;

  return out.str();
}


//-----------------------------------------------------------------------------
// GetOffset / SetOffset
//
float FuseFilter::GetOffset()
{
  return this->m_Offset;
}
FuseFilter::Self& FuseFilter::SetOffset(float o)
{
  this->m_Offset = o;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* FuseFilter::Execute(itk::simple::Image* image1,
                                        itk::simple::Image* image2,
                                        float offset)
{
  this->SetOffset( offset );
  return this->Execute( image1, image2 );
}

//
// Execute
//
itk::simple::Image* FuseFilter::Execute(itk::simple::Image* image1,
                                                itk::simple::Image* image2)
{

  // Cast the second image to the same type as the first
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( image1->GetPixelIDValue() );
  std::auto_ptr<itk::simple::Image> image2rsmp( castFilter.Execute( image2 ) );

  // Resample the second image to match the first
  sitkIM::ResampleFilter resampleFilter;
  image2rsmp.reset( resampleFilter.Execute(image2rsmp.get(), image1) );

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
    sitkExceptionMacro("Both image for FuseFilter don't match type or dimension!");
    }

  // Dispatch the proper member function
  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image1 , image2rsmp.get() );

}


//----------------------------------------------------------------------------
// ExecuteInternal
//
template <class TImageType>
itk::simple::Image* FuseFilter::ExecuteInternal(itk::simple::Image* inImage1,
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
  itk::ImageRegionIterator< InputImageType > it1( image1,
        image1->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< InputImageType > it2( image2,
        image2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    double tf2 = it2.Get();
    if( tf2>tf1 )
      {
      double tf = this->m_Offset + tf2;
      it1.Set( ( typename InputImageType::PixelType )tf );
      }
    ++it1;
    ++it2;
    }

  // Return the resulting SimpleITK image
  return new itk::simple::Image( image1 );
  }

} // end namespace sitkIM
