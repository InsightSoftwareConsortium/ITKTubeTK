/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "tubeMessage.h"

#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "ComputeRegionSignaturesCLP.h"

// STL includes

typedef enum
{
  SIGNATURE_CVT_CTR = 0,
  SIGNATURE_CVT_MIN = 1
} SignatureType;


template< class TPixel, unsigned int VImageDimension >
int DoIt( int argc, char * argv[] );

#include "tubeCLIHelperFunctions.h"


int main( int argc, char **argv )
{
  PARSE_ARGS;
  return tube::ParseArgsAndCallDoIt( argSegImageFileName, argc, argv );
}


/** Check for equality of VNL vectors.
 *
 * \param v 1st VNL vector
 * \param g 2nd VNL vector
 * \param tol Tolerance value
 * \return true in case of equality, false else
 *
 */
template< typename T >
bool CheckVNLVectorEquality( const vnl_vector< T > &v,
  const vnl_vector< T > &g, double tol = 1.0e-6 )
{
  if( v.size() != g.size() )
    {
    return false;
    }
  for( unsigned i=0; i<v.size(); ++i )
    {
    if( vnl_math_abs( g.get( i ) - v.get( i ) ) > tol )
      {
      return false;
      }
    }
  return true;
}


/** Check for equality of VNL matrices.
 *
 * \param V 1st VNL matrix
 * \param G 2nd VNL matrix
 * \param tol Tolerance value
 * \return true in case of equality, false otherwise
 *
 */
template< typename T >
bool CheckVNLMatrixEquality( const vnl_matrix< T > &V,
  const vnl_matrix< T > &G, double tol = 1.0e-6 )
{
  if( V.rows() != G.rows() || V.cols() != G.cols() )
    {
    return false;
    }
  for( unsigned int r=0; r<V.rows(); ++r )
    {
    for( unsigned int c=0; c<V.cols(); ++c )
      {
      if( vnl_math_abs( V.get( r, c ) - G.get( r, c ) ) > tol )
        {
        return false;
        }
      }
    }
  return true;
}


/** Check compatibility of imgages in terms of size, origin, spacing, direction.
 *
 *  \param imageA 1st input image of type ImageType
 *  \param imageB 2nd input image of type ImageType
 *  \returns true if image have equal spacing and size, false otherwise
 *
 */
template< typename ImageT >
bool CheckCompatibility( typename ImageT::Pointer imageA,
  typename ImageT::Pointer imageB )
{
  // Spacing tolerance is imageA's 1st coord. spacing * 1e-6
  const double spacingTol = imageA->GetSpacing()[0] * 1.0e-6;

  if( !CheckVNLVectorEquality( imageA->GetOrigin().GetVnlVector(),
    imageB->GetOrigin().GetVnlVector() ) )
    {
    tube::ErrorMessage( "Origin mismatch between input images" );
    return false;
    }
  if( !CheckVNLVectorEquality( imageA->GetSpacing().GetVnlVector(),
    imageB->GetSpacing().GetVnlVector(), spacingTol ) )
    {
    tube::ErrorMessage( "Spacing mismatch between input images" );
    return false;
    }
  if( !( imageA->GetLargestPossibleRegion().GetSize() ==
         imageB->GetLargestPossibleRegion().GetSize() ) )
    {
    tube::ErrorMessage( "Size mismatch between input images" );
    return false;
    }
  if( !CheckVNLMatrixEquality(
    imageA->GetDirection().GetVnlMatrix().as_ref(),
    imageB->GetDirection().GetVnlMatrix().as_ref() ) )
    {
    tube::ErrorMessage( "Directions mismatch between input images" );
    return false;
    }
  return true;
}


/** Create an empty image.
 *
 *  \param outImage Image that is about to be created ( has to exist )
 *  \param targetSize Desired size of the image
 *
 */
template< class ImageT >
void CreateEmptyImage( typename ImageT::Pointer & outImage,
  typename ImageT::SizeType targetSize )
{
  typename ImageT::IndexType start;
  start.Fill( 0 );
  typename ImageT::SizeType outImageSize = targetSize;

  typename ImageT::RegionType region( start, outImageSize );
  outImage->SetRegions( region );
  outImage->Allocate();
  outImage->FillBuffer( 0 );
}


/** Check if image contains only discrete values.
 *
 * \param fileName Image file name
 * \return 'true' if image, given by 'fileName' has a discrete-value type
 *
 */
bool IsDiscrete( const std::string & fileName )
{
  typedef itk::ImageIOBase::IOComponentType ScalarPixelType;

  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
    fileName.c_str(), itk::ImageIOFactory::ReadMode );

  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();
  const ScalarPixelType pixelType = imageIO->GetComponentType();

  if( pixelType == itk::ImageIOBase::FLOAT ||
      pixelType == itk::ImageIOBase::DOUBLE )
    {
    return false;
    }
  return true;
}


/** Map all S elements in a set to {0,...S-1}.
 *
 *  \param in STL set of values to be renumbered
 *  \param map STL map that is filled up
 */
template< typename TPixel >
void CreateRenumberingMap( const std::set< TPixel > & in,
  std::map< TPixel, unsigned int > & map )
{
  unsigned int cnt = 0;
  typename std::set< TPixel >::const_iterator it = in.begin();
  while( it != in.end() )
    {
    map[( *it )] = cnt++;
    ++it;
    }
}


template< class TPixel, unsigned int VImageDimension >
int DoIt( int argc, char **argv )
{
  PARSE_ARGS;

  // Commonly-used types
  typedef TPixel                                          InputPixelType;
  typedef itk::Image< InputPixelType, VImageDimension >   InputImageType;

  typedef InputPixelType                                  OutputPixelType;
  typedef itk::Image< OutputPixelType, VImageDimension >  OutputImageType;

  typedef itk::ImageFileReader< InputImageType >          InputReaderType;
  typedef itk::ImageFileWriter< OutputImageType >         OutputWriterType;

  typedef itk::ImageRegionIteratorWithIndex< InputImageType >
    ImageIteratorType;

  typedef std::vector< typename InputImageType::IndexType >
    InputIndexVectorType;

  typedef std::map< TPixel, InputIndexVectorType >    IdToIndexVectorType;


  // Read segmenation image and image of Central-Voronoi-Tesellation ( CVT )
  // cells
  typename InputImageType::Pointer  segImage, cvtImage;
  typename InputReaderType::Pointer segImageReader = InputReaderType::New();
  typename InputReaderType::Pointer cvtImageReader = InputReaderType::New();

  segImageReader->SetFileName( argSegImageFileName.c_str() );
  cvtImageReader->SetFileName( argCVTImageFileName.c_str() );

  try
    {
    segImageReader->Update();
    cvtImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  segImage = segImageReader->GetOutput();
  cvtImage = cvtImageReader->GetOutput();

  if( !IsDiscrete( argSegImageFileName ) ||
      !IsDiscrete( argCVTImageFileName ) )
    {
    tube::ErrorMessage( "Input images contain non-discrete values" );
    return EXIT_FAILURE;
    }

  // Read the CVT center file
  typename std::vector< typename InputImageType::IndexType > cvtCentersIdx;

  std::ifstream cvtCenterFile( argCVTCenterFileName.c_str() );
  if( cvtCenterFile.is_open() )
    {
    double c0 = 0.0;
    double c1 = 0.0;
    double c2 = 0.0;
    int nCenters;

    cvtCenterFile >> nCenters;
    cvtCentersIdx.resize( nCenters );

    tube::FmtInfoMessage( "Reading %d CVT centers ...",
                          nCenters );

    unsigned int centerCounter = 0;
    while( cvtCenterFile >> c0 >> c1 >> c2 )
      {
      itkAssertOrThrowMacro( VImageDimension == 3,
                             "No support for 2-D CVT center file" );

      itk::Point< float, VImageDimension > p;
      p[0] = c0;
      p[1] = c1;
      p[2] = c2;

      typename InputImageType::IndexType targetIndex;
      bool isInside = cvtImage->TransformPhysicalPointToIndex( p,
        targetIndex );
      if( !isInside )
        {
        tube::ErrorMessage( "CVT cell center outside image" );
        return EXIT_FAILURE;
        }

      cvtCentersIdx[centerCounter] = targetIndex;
      ++centerCounter;
      }
    cvtCenterFile.close();
    itkAssertOrThrowMacro( centerCounter = nCenters,
                           "#CVT centers mismatch" );
    }
  else
    {
    tube::FmtErrorMessage( "Could not read %s",
                           argCVTCenterFileName.c_str() );
    return EXIT_FAILURE;
    }

  // Ensure compatibility of CVT and segmenation image
  if( !CheckCompatibility< InputImageType >( segImage, cvtImage ) )
    {
    tube::FmtErrorMessage( "%s and %s are incompatible",
      argSegImageFileName.c_str(), argCVTImageFileName.c_str() );
    return EXIT_FAILURE;
    }


  // Determine the number of unqique segmentation regions
  typename itk::ImageRegionIteratorWithIndex< InputImageType > segImageIt(
    segImage, segImage->GetLargestPossibleRegion() );

  // Holds the IDs of the segmentation regions
  std::set< TPixel > segmentationRegions;

  segImageIt.GoToBegin();
  while( !segImageIt.IsAtEnd() )
    {
    segmentationRegions.insert( segImageIt.Get() );
    ++segImageIt;
    }

  tube::FmtInfoMessage( "Found %d distinct segmentation IDs",
    segmentationRegions.size() );

  // Map segmentation IDs to the ( artificial ) segmentation IDs {0, ..., S-1}
  std::map< TPixel, unsigned int > fwdSegMapper;
  std::map< unsigned int, TPixel > invSegMapper;

  unsigned int segmentationRegionCounter = 0;
  typename std::set< TPixel>::const_iterator regionIt =
    segmentationRegions.begin();

  while( regionIt != segmentationRegions.end() )
    {
    fwdSegMapper[ *regionIt ] = segmentationRegionCounter;
    invSegMapper[ segmentationRegionCounter++ ] = *regionIt;
    ++regionIt;
    }

  itkAssertOrThrowMacro( segmentationRegions.size() ==
    segmentationRegionCounter,
    "Size mismatch after seg. region renumbering" );

  // Create a map of ( artificial ) segmentation IDs -> Vector of voxel
  // indices
  IdToIndexVectorType segToIndexVector;

  segImageIt.GoToBegin();
  while( !segImageIt.IsAtEnd() )
    {
    const typename InputImageType::IndexType &index =
      segImageIt.GetIndex();
    segToIndexVector[ fwdSegMapper[ segImageIt.Get() ] ].push_back( index );
    ++segImageIt;
    }

  itkAssertOrThrowMacro( segToIndexVector.size() ==
    segmentationRegions.size(),
    "Size mismatch after region -> index vector mapping" );

  // Map CVT cell IDs to ( artificial ) CVT cell IDs {0, ..., C-1}
  ImageIteratorType cvtImageIt( cvtImage,
    cvtImage->GetLargestPossibleRegion() );

  std::set< TPixel > cvtCellSet;
  cvtImageIt.GoToBegin();
  while( !cvtImageIt.IsAtEnd() )
    {
    cvtCellSet.insert( cvtImageIt.Get() );
    ++cvtImageIt;
    }

  std::map< TPixel, unsigned int > fwdCVTMapper;
  CreateRenumberingMap< TPixel >( cvtCellSet, fwdCVTMapper );

  itkAssertOrThrowMacro( fwdCVTMapper.size() == cvtCellSet.size(),
    "Size mismatch after renumbering CVT cells" );

  // Create a map of ( artificial ) CVT cell IDs -> Vector of voxel indices
  IdToIndexVectorType cvtToIndexVector;

  cvtImageIt.GoToBegin();
  while( !cvtImageIt.IsAtEnd() )
    {
    TPixel id = cvtImageIt.Get();
    cvtToIndexVector[ fwdCVTMapper[id] ].push_back( cvtImageIt.GetIndex() );
    ++cvtImageIt;
    }

  TPixel largestKey = cvtToIndexVector.rbegin()->first;
  itkAssertOrThrowMacro( ( largestKey + 1 ) ==
    static_cast< TPixel >( cvtToIndexVector.size() ),
    "Size mismatch after CVT to index vector mapping" );

  // Build index vector for exclusion regions
  boost::dynamic_bitset<> excludeRegions( segmentationRegions.size() );
  for( unsigned int e = 0; e < argExcludeRegions.size(); ++e )
    {
    unsigned int mappedRegionID = fwdCVTMapper[ argExcludeRegions[e] ];
    itkAssertOrThrowMacro( mappedRegionID < segmentationRegions.size(),
      "Excluded region ID out of range" );

    excludeRegions[ mappedRegionID ] = 1;
    }

  // Create numberOfSegmentedRegions-dimensional distance signature vectors
  unsigned int numberOfCVTCells = cvtToIndexVector.size();
  unsigned int numberOfSegmentedRegions = segmentationRegions.size();

  vnl_matrix< TPixel > signatureMatrix ( numberOfCVTCells,
    numberOfSegmentedRegions );

  typename IdToIndexVectorType::iterator mapIt = segToIndexVector.begin();
  while( mapIt != segToIndexVector.end() )
    {
    TPixel id = ( *mapIt ).first;

    if( excludeRegions[id] )
      {
      tube::FmtInfoMessage( "Exclude segmentation region %d ( Orig: %d )",
        id, invSegMapper[id] );

      ++mapIt;
      continue;
      }

    std::stringstream mapFileName, prtFileName;
    mapFileName << boost::format( "%04i-map.mha" ) % invSegMapper[id];
    prtFileName << boost::format( "%04i-prt.mha" ) % invSegMapper[id];

    boost::filesystem::path dir( argOutputDirectory.c_str() );
    boost::filesystem::path mapFilePath ( mapFileName.str().c_str() );
    boost::filesystem::path prtFilePath ( prtFileName.str().c_str() );
    boost::filesystem::path fullMapFilePath = dir / mapFilePath;
    boost::filesystem::path fullPrtFilePath = dir / prtFilePath;

    typename OutputImageType::Pointer prtImage, mapImage;
    if( boost::filesystem::exists( fullPrtFilePath.string().c_str() ) &&
        boost::filesystem::exists( fullMapFilePath.string().c_str() ) )
      // Try to load existing data ...
      {
      tube::FmtInfoMessage( "Found %s and %s - Loading ...",
        fullPrtFilePath.string().c_str(),
        fullMapFilePath.string().c_str() );

      typename InputReaderType::Pointer singlePrtReader =
        InputReaderType::New();
      typename InputReaderType::Pointer singleMapReader =
        InputReaderType::New();
      singlePrtReader->SetFileName( fullPrtFilePath.string().c_str() );
      singleMapReader->SetFileName( fullMapFilePath.string().c_str() );
      try
        {
        singlePrtReader->Update();
        singleMapReader->Update();
        }
      catch( itk::ExceptionObject & ex )
        {
        tube::ErrorMessage( ex.what() );
        return EXIT_FAILURE;
        }
      mapImage = singleMapReader->GetOutput();
      prtImage = singlePrtReader->GetOutput();
      }
    else
      // Otherwise, recompute distance maps ...
      {
      tube::FmtInfoMessage(
        "Recomputing distance map for ( original ) seg. %d",
        invSegMapper[id] );

      typename OutputImageType::Pointer outImage = OutputImageType::New();

      CreateEmptyImage< OutputImageType >( outImage,
        segImage->GetLargestPossibleRegion().GetSize() );

      itkAssertOrThrowMacro( segImage->GetLargestPossibleRegion().GetSize()
        == outImage->GetLargestPossibleRegion().GetSize(),
        "Newly created output image differs in size" );

      const InputIndexVectorType &indices = ( *mapIt ).second;
      for( unsigned int v=0; v<indices.size(); ++v )
        {
        outImage->SetPixel( indices[v], 1 );
        }

      typedef typename itk::DanielssonDistanceMapImageFilter<
        InputImageType,
        OutputImageType > DistanceMapType;

      typename DistanceMapType::Pointer distanceMap =
        DistanceMapType::New();
      distanceMap->SetInputIsBinary( true );
      distanceMap->SetUseImageSpacing( true );
      distanceMap->SetInput( outImage );
      distanceMap->Update();

      mapImage = distanceMap->GetDistanceMap();
      prtImage = outImage;

      typename OutputWriterType::Pointer prtWriter =
        OutputWriterType::New();
      prtWriter->SetFileName( fullPrtFilePath.string().c_str() );
      prtWriter->SetInput( prtImage );
      prtWriter->SetUseCompression( true );

      typename OutputWriterType::Pointer mapWriter =
        OutputWriterType::New();
      mapWriter->SetFileName( fullMapFilePath.string().c_str() );
      mapWriter->SetInput( mapImage );
      mapWriter->SetUseCompression( true );

      try
        {
        mapWriter->Update();
        prtWriter->Update();
        }
      catch( itk::ExceptionObject & ex )
        {
        tube::ErrorMessage( ex.what() );
        return EXIT_FAILURE;
        }
      }

      switch( argSignatureSelector )
      {
      case SIGNATURE_CVT_CTR:
        {
        for( unsigned int c=0; c<numberOfCVTCells; ++c )
          {
          // Distance to c-th ( i.e., c-th CVT cell ) center
          TPixel dist = mapImage->GetPixel( cvtCentersIdx[c] );
          signatureMatrix.put( c, id, dist );
          }
        break;
        }
      case SIGNATURE_CVT_MIN:
        {
        namespace acc = boost::accumulators;
        typedef acc::stats< acc::tag::min >               StatsType;
        typedef acc::accumulator_set< TPixel, StatsType > AccumulatorType;

        typename IdToIndexVectorType::iterator cvtToIndexVectorIt =
          cvtToIndexVector.begin();

        while( cvtToIndexVectorIt != cvtToIndexVector.end() )
          {
          AccumulatorType cellAcc;

          // Artificial CVT cell ID
          TPixel cvtCell = ( *cvtToIndexVectorIt ).first;

          const std::vector< typename InputImageType::IndexType > &
            voxelIndices = ( *cvtToIndexVectorIt ).second;

          // Add distances to accumulator to keep track of the minimum
          for( unsigned int i=0; i<voxelIndices.size(); ++i )
            {
            TPixel dist = mapImage->GetPixel( voxelIndices[i] );
            cellAcc( dist );
            }

          // Get the minimum distance to the current segment
          TPixel minDist = acc::min( cellAcc );

          // We need to map the original CVT cell ID to the artificial one
          signatureMatrix.put( cvtCell, id, minDist );
          ++cvtToIndexVectorIt;
          }
        break;
        }
      }
      ++mapIt;
    }

  std::ofstream outSignatureFile( argDistanceSignatureFileName.c_str() );
  if( outSignatureFile.is_open() )
    {
    for( unsigned int r=0; r<signatureMatrix.rows(); ++r )
      {
      for( unsigned int c=0; c<signatureMatrix.columns()-1; ++c )
        {
        outSignatureFile << signatureMatrix.get( r, c ) << ",";
        }
      outSignatureFile << signatureMatrix.get( r,
        signatureMatrix.columns() - 1 ) << std::endl;
      }
    }
  else
    {
    tube::FmtErrorMessage( "Could not open %s for writing",
      argDistanceSignatureFileName.c_str() );
    return EXIT_FAILURE;
    }
  outSignatureFile.close();
  return EXIT_SUCCESS;
}
