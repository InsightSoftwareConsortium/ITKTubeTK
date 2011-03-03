#include "sitkIMResampleFilter.h"

#include "itkResampleImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor (currently does nothing)
//
ResampleFilter::ResampleFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();
}


//
// ToString
//
std::string ResampleFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Resample"
      << std::endl;

  return out.str();
}



//-----------------------------------------------------------------------------
// Execute
//
itk::simple::Image* ResampleFilter::Execute(itk::simple::Image* image1,
                                         itk::simple::Image* image2)
{

  // Cast the second image to the same type as the first
  itk::simple::CastImageFilter castFilter;
  castFilter.SetOutputPixelType( image1->GetPixelIDValue() );
  std::auto_ptr<itk::simple::Image> image2rsmp( castFilter.Execute( image2 ) );

  // Prep to call ExecuteInternal
  itk::simple::PixelIDValueType type = image1->GetPixelIDValue();
  unsigned int dimension = image1->GetDimension();

  // Dispatch the proper member function
  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image1 , image2rsmp.get() );

}



//----------------------------------------------------------------------------
// ExecuteInternal
//
template <class TImageType>
itk::simple::Image* ResampleFilter::ExecuteInternal(
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
  typename OutputImageType::Pointer output = image1;

  for( unsigned int i = 0; i < InputImageType::ImageDimension; i++ )
    {
    if( image1->GetLargestPossibleRegion().GetSize()[i]
          != image2->GetLargestPossibleRegion().GetSize()[i]
        || image1->GetLargestPossibleRegion().GetIndex()[i]
            != image2->GetLargestPossibleRegion().GetIndex()[i]
        || image1->GetSpacing()[i] != image2->GetSpacing()[i]
        || image1->GetOrigin()[i] != image2->GetOrigin()[i]  )
      {
      typedef typename itk::ResampleImageFilter< OutputImageType,
                OutputImageType> ResampleFilterType;
      typename ResampleFilterType::Pointer filter =
        ResampleFilterType::New();
      filter->SetInput( image1 );
      filter->SetSize( image2->GetLargestPossibleRegion().GetSize() );
      filter->SetOutputOrigin( image2->GetOrigin() );
      filter->SetOutputSpacing( image2->GetSpacing() );
      filter->SetDefaultPixelValue( 0 );
      filter->Update();
      output = filter->GetOutput();
      }
    }

  // Return the resulting SimpleITK image
  return new itk::simple::Image( output );
  }


} // end namespace sitkIM
