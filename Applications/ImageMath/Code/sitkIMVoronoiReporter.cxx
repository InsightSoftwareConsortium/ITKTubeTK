#include "sitkIMVoronoiReporter.h"

#include "itkImageRegionIteratorWithIndex.h"


// tubetk includes
#include "itkTubeCVTImageFilter.h"

#include <sstream>
#include <iostream>
#include <fstream>

namespace sitkIM
{

//-----------------------------------------------------------------------------
// Default Constructor
//
VoronoiReporter::VoronoiReporter()
{
  this->m_MemberFactory.reset(
    new itk::simple::detail::MemberFunctionFactory<MemberFunctionType>( this ) );
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 3 > ();
  this->m_MemberFactory->RegisterMemberFunctions< PixelIDTypeList, 2 > ();

  m_NumCentroids = 100;
  m_NumIterations = 100;
  m_NumSamples = 100;
  m_OutputFilename = "";
}

//
// ToString
//
std::string VoronoiReporter::ToString() const
{
  std::stringstream out;
  out << "Reporter: Voronoi" << std::endl
      << "  Num Centroids: " << m_NumCentroids << std::endl
      << "  Num Iterations: " << m_NumIterations << std::endl
      << "  Num Samples: " << m_NumSamples << std::endl
      << "  Output Filename: " << m_OutputFilename << std::endl;

  return out.str();
}

//-----------------------------------------------------------------------------
// GetNumCentroids / SetNumCentroids
//
unsigned int VoronoiReporter::GetNumCentroids()
{
  return this->m_NumCentroids;
}
VoronoiReporter::Self& VoronoiReporter::SetNumCentroids(unsigned int n)
{
  this->m_NumCentroids = n;
  return *this;
}

//
// GetNumIterations / SetNumIterations
//
unsigned int VoronoiReporter::GetNumIterations()
{
  return this->m_NumIterations;
}
VoronoiReporter::Self& VoronoiReporter::SetNumIterations(unsigned int n)
{
  this->m_NumIterations = n;
  return *this;
}

// GetNumSamples / SetNumSamples
//
unsigned int VoronoiReporter::GetNumSamples()
{
  return this->m_NumSamples;
}
VoronoiReporter::Self& VoronoiReporter::SetNumSamples(unsigned int n)
{
  this->m_NumSamples = n;
  return *this;
}

//
// GetOutputFilename / SetOutputFilename
//
std::string VoronoiReporter::GetOutputFilename()
{
  return this->m_OutputFilename;
}
VoronoiReporter::Self& VoronoiReporter::SetOutputFilename(std::string f)
{
  this->m_OutputFilename = f;
  return *this;
}


//-----------------------------------------------------------------------------
// Execute
itk::simple::Image* VoronoiReporter::Execute(itk::simple::Image* image,
                                                       unsigned int centroids,
                                                       unsigned int iters,
                                                       unsigned int samples,
                                                       std::string f)
{
  this->SetNumCentroids( centroids );
  this->SetNumIterations( iters );
  this->SetNumSamples( samples );
  this->SetOutputFilename( f );

  return this->Execute( image );
}

//
// Execute (Filter version)
//
itk::simple::Image* VoronoiReporter::Execute( itk::simple::Image* image )
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
itk::simple::Image* VoronoiReporter::ExecuteInternal(
                                                itk::simple::Image* inImage )
  {
  // Set up ITK image types
  typedef TImageType     InputImageType;

  // Cast the images and make sure they are the correct type
  typename InputImageType::Pointer image =
    dynamic_cast <InputImageType*> ( inImage->GetImageBase().GetPointer() );
  if ( image.IsNull() )
    {
    sitkExceptionMacro( "Unexpected template dispatch error!" );
    }

  // Perform the guts of the filter
  typedef itk::tube::CVTImageFilter<InputImageType, InputImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetNumberOfSamples( this->m_NumSamples );
  filter->SetNumberOfCentroids( this->m_NumCentroids );
  filter->SetNumberOfIterations( this->m_NumIterations );
  filter->SetNumberOfSamplesPerBatch( this->m_NumIterations );
  filter->Update();

  std::ofstream writeStream;
  writeStream.open( this->m_OutputFilename.c_str(),
    std::ios::binary | std::ios::out );
  if( ! writeStream.rdbuf()->is_open() )
    {
    sitkExceptionMacro( "Cannot write to file : " << this->m_OutputFilename );
    return NULL;
    }
  writeStream << this->m_NumCentroids << std::endl;
  for( unsigned int i=0; i<this->m_NumCentroids; i++ )
    {
    for( unsigned int j = 0; j<3; j++ )
      {
      writeStream << ( *( filter->GetCentroids() ) )[i][j];
      if( j<2 )
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  image = filter->GetOutput();
  typename InputImageType::SizeType size =
    image->GetLargestPossibleRegion().GetSize();

  this->m_OutputFilename = this->m_OutputFilename + ".mat";

  vnl_matrix<int> aMat( this->m_NumCentroids, this->m_NumCentroids );
  aMat.fill( 0 );

  itk::Index<InputImageType::ImageDimension> indx;
  itk::Index<InputImageType::ImageDimension> indx2;
  itk::Index<InputImageType::ImageDimension> indx3;
  indx.Fill( 0 );
  bool done = false;
  bool invalid = false;
  bool done2 = false;
  int c, n;
  while( !done )
    {
    c = ( int )( image->GetPixel( indx )-1 );
    indx2.Fill( 0 );
    indx2[0] = 1;
    invalid = false;
    done2 = false;
    while( !done2 )
      {
      invalid = false;
      for( unsigned int d=0; d<InputImageType::ImageDimension; d++ )
        {
        indx3[d] = indx[d] + indx2[d];
        if( indx3[d] >= ( int )size[d] )
          {
          invalid = true;
          break;
          }
        }
      if( !invalid )
        {
        n = ( int )( image->GetPixel( indx3 )-1 );
        if( c != n )
          {
          aMat( c, n ) = 1;
          aMat( n, c ) = 1;
          }
        }
      int i=0;
      indx2[i]++;
      while( !done2 && indx2[i]>=2 )
        {
        indx2[i] = 0;
        i++;
        if( i>2 )
          {
          done2 = true;
          }
        else
          {
          indx2[i]++;
          }
        }
      }
    int i = 0;
    indx[i]++;
    while( !done && indx[i]>=( int )size[i] )
      {
      indx[i] = 0;
      i++;
      if( i>2 )
        {
        done = true;
        }
      else
        {
        if( i == 2 )
          {
          std::cout << "Computing adjacency of slice : " << indx[2]
                    << std::endl;
          }
        indx[i]++;
        }
      }
    }

  writeStream.open( this->m_OutputFilename.c_str(),
    std::ios::binary | std::ios::out );
  if( ! writeStream.rdbuf()->is_open() )
    {
    sitkExceptionMacro( "Cannot write to file : " << this->m_OutputFilename );
    return NULL;
    }
  writeStream << this->m_NumCentroids << std::endl;
  for( unsigned int i=0; i<this->m_NumCentroids; i++ )
    {
    for( unsigned int j = 0; j<this->m_NumCentroids; j++ )
      {
      writeStream << aMat( i, j );
      if( j<this->m_NumCentroids-1 )
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  // Just return the input image if successful
  return inImage;

  }

} // end namespace sitkIM
