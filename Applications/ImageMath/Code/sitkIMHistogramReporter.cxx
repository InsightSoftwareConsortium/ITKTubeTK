#include "sitkIMHistogramReporter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkArray.h"

#include <sstream>
#include <iostream>
#include <fstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
HistogramReporter::HistogramReporter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_NBins = 255;
  m_OutputFilename = "";
}


//
// ToString
//
std::string HistogramReporter::ToString() const
{
  std::stringstream out;
  out << "Reporter: Histogram" << std::endl
      << "  NBins: " << m_NBins << std::endl
      << "  Output Filename: " << m_OutputFilename << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNBins / SetNBins
//
unsigned int HistogramReporter::GetNBins()
{
  return this->m_NBins;
}
HistogramReporter::Self& HistogramReporter::SetNBins(unsigned int n)
{
  this->m_NBins = n;
  return *this;
}

//
// GetOutputFilename / SetOutputFilename
//
std::string HistogramReporter::GetOutputFilename()
{
  return this->m_OutputFilename;
}
HistogramReporter::Self& HistogramReporter::SetOutputFilename(std::string f)
{
  this->m_OutputFilename = f;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* HistogramReporter::Execute(itk::simple::Image* image,
                                                       unsigned int n, std::string f)
{
  this->SetNBins( n );
  this->SetOutputFilename( f );

  return this->Execute( image );
}

//
// Execute (Filter version)
//
itk::simple::Image* HistogramReporter::Execute( itk::simple::Image* image )
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
itk::simple::Image* HistogramReporter::ExecuteInternal(
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
  itk::ImageRegionIteratorWithIndex< InputImageType > it1( image,
        image->GetLargestPossibleRegion() );
  it1.GoToBegin();
  double binMin = it1.Get();
  double binMax = it1.Get();
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    if( tf < binMin )
      {
      binMin = tf;
      }
    else
      {
      if( tf > binMax )
        {
        binMax = tf;
        }
      }
    ++it1;
    }
  std::cout << "  binMin = " << binMin << std::endl;
  std::cout << "  binMax = " << binMax << std::endl;
  it1.GoToBegin();
  itk::Array<double> bin;
  bin.set_size( this->m_NBins );
  bin.Fill( 0 );
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    tf = ( tf-binMin )/( binMax-binMin ) * this->m_NBins;
    if( tf>this->m_NBins-1 )
      {
      tf = this->m_NBins-1;
      }
    else
      {
      if( tf<0 )
        {
        tf = 0;
        }
      }
    bin[( int )tf]++;
    ++it1;
    }
  std::ofstream writeStream;
  writeStream.open( this->m_OutputFilename.c_str(), std::ios::binary | std::ios::out );
  if( ! writeStream.rdbuf()->is_open() )
    {
    sitkExceptionMacro( "Cannot write to file : " << this->m_OutputFilename << std::endl);
    return NULL;
    }
  for( unsigned int i=0; i<this->m_NBins; i++ )
    {
    writeStream << ( i/( double )this->m_NBins )*( binMax-binMin )+binMin
                << " " << bin[i] << std::endl;
    }
  writeStream.close();

  // Just return the input image if successful
  return inImage;

  }

} // end namespace sitkIM
