#include "sitkIMCorrectionFilter.h"
#include "sitkIMResampleFilter.h"

#include "itkImageRegionIterator.h"
#include "itkHistogramMatchingImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
CorrectionFilter::CorrectionFilter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_NumBins = 255;
  m_NumMatchPoints = 10;  // Might need to change this
}


//
// ToString
//
std::string CorrectionFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: Correction" << std::endl
      << "  Num Bins: " << m_NumBins << std::endl
      << "  Num Match Points: " << m_NumMatchPoints << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNumBins / SetNumBins
//
unsigned int CorrectionFilter::GetNumBins()
{
  return this->m_NumBins;
}
CorrectionFilter::Self& CorrectionFilter::SetNumBins(unsigned int n)
{
  this->m_NumBins = n;
  return *this;
}

//
// GetNumMatchPoints / SetNumMatchPoints
//
unsigned int CorrectionFilter::GetNumMatchPoints()
{
  return this->m_NumMatchPoints;
}
CorrectionFilter::Self& CorrectionFilter::SetNumMatchPoints(unsigned int n)
{
  this->m_NumMatchPoints = n;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* CorrectionFilter::Execute(itk::simple::Image* image1,
                                              itk::simple::Image* image2,
                                              unsigned int nbins, unsigned int nmatch)
{
  this->SetNumBins( nbins );
  this->SetNumMatchPoints( nmatch );

  return this->Execute( image1, image2 );
}

//
// Execute (Filter version)
//
itk::simple::Image* CorrectionFilter::Execute( itk::simple::Image* image1,
                                                       itk::simple::Image* image2 )
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
    sitkExceptionMacro("Both images for FuseFilter don't match type or dimension!");
    }

  // Dispatch the proper member function
  return this->m_MemberFactory->GetMemberFunction( type, dimension )( image1, image2rsmp.get() );

}


//----------------------------------------------------------------------------
// ExecuteInternal
//
template <class TImageType>
itk::simple::Image* CorrectionFilter::ExecuteInternal(
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
  typedef itk::HistogramMatchingImageFilter< InputImageType, InputImageType >
      HistogramMatchFilterType;
  typename HistogramMatchFilterType::Pointer matchFilter;
  matchFilter = HistogramMatchFilterType::New();
  matchFilter->SetReferenceImage( image2 );
  matchFilter->SetInput( image1 );
  matchFilter->SetNumberOfHistogramLevels( this->m_NumBins );
  matchFilter->SetNumberOfMatchPoints( this->m_NumMatchPoints );
  matchFilter->Update();
  image1 = matchFilter->GetOutput();

  // return the result
  return new itk::simple::Image( image1 );

  }

} // end namespace sitkIM
