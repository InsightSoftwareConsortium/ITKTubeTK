#include "sitkIMCorrectionSliceFilter.h"

#include "itkImageRegionIterator.h"
#include "itkHistogramMatchingImageFilter.h"

#include <sstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
CorrectionSliceFilter::CorrectionSliceFilter()
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
std::string CorrectionSliceFilter::ToString() const
{
  std::stringstream out;
  out << "Filter: CorrectionSlice" << std::endl
      << "  Num Bins: " << m_NumBins << std::endl
      << "  Num Match Points: " << m_NumMatchPoints << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNumBins / SetNumBins
//
unsigned int CorrectionSliceFilter::GetNumBins()
{
  return this->m_NumBins;
}
CorrectionSliceFilter::Self& CorrectionSliceFilter::SetNumBins(unsigned int n)
{
  this->m_NumBins = n;
  return *this;
}

//
// GetNumMatchPoints / SetNumMatchPoints
//
unsigned int CorrectionSliceFilter::GetNumMatchPoints()
{
  return this->m_NumMatchPoints;
}
CorrectionSliceFilter::Self& CorrectionSliceFilter::SetNumMatchPoints(unsigned int n)
{
  this->m_NumMatchPoints = n;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* CorrectionSliceFilter::Execute(itk::simple::Image* image,
                                                        unsigned int nbins, unsigned int nmatch)
{
  this->SetNumBins( nbins );
  this->SetNumMatchPoints( nmatch );

  return this->Execute( image );
}

//
// Execute (Filter version)
//
itk::simple::Image* CorrectionSliceFilter::Execute( itk::simple::Image* image )
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
itk::simple::Image* CorrectionSliceFilter::ExecuteInternal(
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

  unsigned int dimensionT = inImage->GetDimension();

  // Perform the guts of the filter
  typedef itk::Image<typename InputImageType::PixelType, 2> ImageType2D;
  typedef itk::HistogramMatchingImageFilter< ImageType2D, ImageType2D >
      HistogramMatchFilterType;
  typename HistogramMatchFilterType::Pointer matchFilter;
  typename ImageType2D::Pointer im2DRef = ImageType2D::New();
  typename ImageType2D::Pointer im2DIn = ImageType2D::New();
  typename ImageType2D::SizeType size2D;
  size2D[0] = image->GetLargestPossibleRegion().GetSize()[0];
  size2D[1] = image->GetLargestPossibleRegion().GetSize()[1];
  im2DRef->SetRegions( size2D );
  im2DRef->Allocate();
  im2DIn->SetRegions( size2D );
  im2DIn->Allocate();
  itk::ImageRegionIterator< InputImageType > it3D( image,
        image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< InputImageType > it3DSliceStart( image,
        image->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType2D > it2DRef( im2DRef,
        im2DRef->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType2D > it2DIn( im2DIn,
        im2DIn->GetLargestPossibleRegion() );
  unsigned int z, y, x;
  it3D.GoToBegin();
  unsigned int zMax = 1;
  if( dimensionT == 3 )
    {
    zMax = image->GetLargestPossibleRegion().GetSize()[2];
    }
  for( z=0; z<dimensionT && z<zMax; z++ )
    {
    it2DRef.GoToBegin();
    for( y=0; y<image->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<image->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it3D.Get() );
        ++it2DRef;
        ++it3D;
        }
      }
    }
  for(; z<zMax; z++ )
    {
    it2DIn.GoToBegin();
    it3DSliceStart = it3D;
    for( y=0; y<image->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<image->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DIn.Set( it3D.Get() );
        ++it2DIn;
        ++it3D;
        }
      }
    matchFilter = HistogramMatchFilterType::New();
    matchFilter->SetReferenceImage( im2DRef );
    matchFilter->SetInput( im2DIn );
    matchFilter->SetNumberOfHistogramLevels( this->m_NumBins );
    matchFilter->SetNumberOfMatchPoints( this->m_NumMatchPoints );
    matchFilter->Update();
    itk::ImageRegionIterator< ImageType2D > it2DOut(
          matchFilter->GetOutput(),
          im2DIn->GetLargestPossibleRegion() );
    it2DRef.GoToBegin();
    it2DOut.GoToBegin();
    it3D = it3DSliceStart;
    for( y=0; y<image->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<image->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it2DOut.Get() );
        it3D.Set( it2DOut.Get() );
        ++it2DRef;
        ++it2DOut;
        ++it3D;
        }
      }
    }

  // return the result
  return new itk::simple::Image( image );

  }

} // end namespace sitkIM
