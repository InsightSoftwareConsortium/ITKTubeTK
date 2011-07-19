#include "sitkIMHistogram2Reporter.h"

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
Histogram2Reporter::Histogram2Reporter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_NBins = 255;
  m_BinMin = 0;
  m_BinSize = 1;
  m_OutputFilename = "";
}


//
// ToString
//
std::string Histogram2Reporter::ToString() const
{
  std::stringstream out;
  out << "Reporter: Histogram2" << std::endl
      << "  NBins: " << m_NBins << std::endl
      << "  Bin Min: " << m_BinMin << std::endl
      << "  Bin Size: " << m_BinSize << std::endl
      << "  Output Filename: " << m_OutputFilename << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNBins / SetNBins
//
unsigned int Histogram2Reporter::GetNBins()
{
  return this->m_NBins;
}
Histogram2Reporter::Self& Histogram2Reporter::SetNBins(unsigned int n)
{
  this->m_NBins = n;
  return *this;
}

//
// GetBinMin / SetBinMin
//
double Histogram2Reporter::GetBinMin()
{
  return this->m_BinMin;
}
Histogram2Reporter::Self& Histogram2Reporter::SetBinMin(double min)
{
  this->m_BinMin = min;
  return *this;
}

//
// GetBinSize / SetBinSize
//
double Histogram2Reporter::GetBinSize()
{
  return this->m_BinSize;
}
Histogram2Reporter::Self& Histogram2Reporter::SetBinSize(double sz)
{
  this->m_BinSize = sz;
  return *this;
}

//
// GetOutputFilename / SetOutputFilename
//
std::string Histogram2Reporter::GetOutputFilename()
{
  return this->m_OutputFilename;
}
Histogram2Reporter::Self& Histogram2Reporter::SetOutputFilename(std::string f)
{
  this->m_OutputFilename = f;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* Histogram2Reporter::Execute(itk::simple::Image* image,
                                                        unsigned int n, double min, double sz,
                                                        std::string f)
{
  this->SetNBins( n );
  this->SetBinMin( min );
  this->SetBinSize( sz );
  this->SetOutputFilename( f );

  return this->Execute( image );
}

//
// Execute (Filter version)
//
itk::simple::Image* Histogram2Reporter::Execute( itk::simple::Image* image )
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
itk::simple::Image* Histogram2Reporter::ExecuteInternal(
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
  double binMax = this->m_BinMin + this->m_BinSize*this->m_NBins;

  itk::ImageRegionIteratorWithIndex< InputImageType > it1( image,
        image->GetLargestPossibleRegion() );
  it1.GoToBegin();
  itk::Array<double> bin;
  bin.set_size( this->m_NBins );
  bin.Fill( 0 );
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    tf = ( tf-this->m_BinMin )/( binMax-this->m_BinMin ) * this->m_NBins;
    if( tf<this->m_NBins && tf>0 )
      {
      bin[( int )tf]++;
      }
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
    writeStream << ( i/( double )this->m_NBins )*( binMax-this->m_BinMin )+this->m_BinMin
                << " " << bin[i] << std::endl;
    }
  writeStream.close();

  // Just return the input image if successful
  return inImage;

  }

} // end namespace sitkIM
