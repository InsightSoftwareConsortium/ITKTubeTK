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

#include <boost/dynamic_bitset.hpp>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>

#include "TransferLabelsToRegionsCLP.h"

template< class TPixel, unsigned int VImageDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  return tube::ParseArgsAndCallDoIt( argInImageFileName, argc, argv );
}

/** Check if image contains discrete values.
 *
 * \param fileName Image file name
 * \return 'true' if image, given by 'fileName' has a discrete-value type
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


/** Check for equality of VNL vectors.
 *
 * \param v 1st VNL vector
 * \param g 2nd VNL vector
 * \param tol Tolerance value
 * \return true in case of equality, false else
 */
template< class T >
bool check_vnl_vector_equality( const vnl_vector<T> &v,
                                const vnl_vector<T> &g,
                                double tol = 1.0e-6 )
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
 * \return true in case of equality, false else
 */
template< class T >
bool check_vnl_matrix_equality( const vnl_matrix<T> &V,
                                const vnl_matrix<T> &G,
                                double tol = 1.0e-6 )
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


/** Create an empty image.
 *
 *  \param outImage Image that is about to be created ( has to exist )
 *  \param targetSize Desired size of the image
 */
template< class TImage >
void CreateEmptyImage(
  typename TImage::Pointer &outImage,
  typename TImage::SizeType targetSize )
{
  typename TImage::IndexType start;
  start.Fill( 0 );
  typename TImage::SizeType outImageSize = targetSize;

  typename TImage::RegionType region( start, outImageSize );
  outImage->SetRegions( region );
  outImage->Allocate();
  outImage->FillBuffer( 0 );
}


/** Evaluate if two images are compatible.
 *  Compatibility is defined in terms of spacing, directions, size and origin.
 *
 *  \param imageA 1st input image of type TImage
 *  \param imageB 2nd input image of type TImage
 *  \returns true if image have equal spacing and size, false else
 */
template< class TImage >
bool CheckCompatibility(
  typename TImage::Pointer imageA,
  typename TImage::Pointer imageB )
{
  // Spacing tolerance is imageA's 1st coord. spacing * 1e-6
  double spacingTol = imageA->GetSpacing()[0] * 1.0e-6;

  if( !check_vnl_vector_equality(
    imageA->GetOrigin().GetVnlVector(),
    imageB->GetOrigin().GetVnlVector() ) )
    {
    tube::ErrorMessage( "Origin mismatch between input images!" );
    return false;
    }
  if( !check_vnl_vector_equality(
    imageA->GetSpacing().GetVnlVector(),
    imageB->GetSpacing().GetVnlVector(), spacingTol ) )
    {
    tube::ErrorMessage( "Spacing mismatch between input images!" );
    return false;
    }
  if( !( imageA->GetLargestPossibleRegion().GetSize() ==
         imageB->GetLargestPossibleRegion().GetSize() ) )
    {
    tube::ErrorMessage( "Size mismatch between input images!" );
    return false;
    }
  if( !check_vnl_matrix_equality(
    imageA->GetDirection().GetVnlMatrix().as_ref(),
    imageB->GetDirection().GetVnlMatrix().as_ref() ) )
    {
    tube::ErrorMessage( "Directions mismatch between input images!" );
    return false;
    }
  return true;
}


template< class TPixel, unsigned int VImageDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Some typedefs - Our input image type equals the output image type, since
  // the images are essentially the same, except that the voxel values
  // represent different entities.
  typedef TPixel                                          InputPixelType;
  typedef itk::Image< InputPixelType, VImageDimension >   InputImageType;

  typedef InputPixelType                                  OutputPixelType;
  typedef itk::Image< OutputPixelType, VImageDimension >  OutputImageType;

  typedef itk::ImageFileReader< InputImageType >          InputReaderType;
  typedef itk::ImageFileWriter< OutputImageType >         OutputWriterType;

  typedef itk::ImageIOBase::IOComponentType               ScalarPixelType;

  // Check input images's type ( we assume discrete valued images )
  if( !IsDiscrete( argInImageFileName ) ||
      !IsDiscrete( argInLabelFileName ) )
    {
    tube::ErrorMessage( "Non-discrete value types in input files!" );
    return EXIT_FAILURE;
    }


  // First, we try to read in the input image( s ), i.e., an image of the
  // Voronoi tesselation and an image that contains discrete label values
  // for each voxel of the input image.
  typename InputImageType::Pointer inImage;
  typename InputReaderType::Pointer inImageReader =
    InputReaderType::New();

  inImageReader->SetFileName( argInImageFileName.c_str() );
  try
    {
    inImageReader->Update();
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  inImage = inImageReader->GetOutput();

  typename InputImageType::Pointer inLabelImage;
  typename InputReaderType::Pointer inLabelImageReader =
    InputReaderType::New();

  inLabelImageReader->SetFileName( argInLabelFileName.c_str() );
  try
    {
    inLabelImageReader->Update();
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  inLabelImage = inLabelImageReader->GetOutput();


  // Since we will use the voxel indices, extracted from the input image, to
  // index voxels in the label image, we have to ensure compatibility w.r.t.
  // spacing, direction, size and origin.
  if( !CheckCompatibility< InputImageType >( inImage, inLabelImage ) )
    {
    tube::ErrorMessage( "Exiting due to image incompatibility!" );
    return EXIT_FAILURE;
    }


  // Determine the set of unique CVT cell IDs.
  typename itk::ImageRegionIteratorWithIndex< InputImageType > inImageIt(
    inImage, inImage->GetLargestPossibleRegion() );

  std::set< TPixel > cvtCellIdentifiers;

  inImageIt.GoToBegin();
  while( !inImageIt.IsAtEnd() )
    {
    cvtCellIdentifiers.insert( inImageIt.Get() );
    ++inImageIt;
    }

  tube::FmtDebugMessage( "%d distinct CVT cells!",
    static_cast<int>( cvtCellIdentifiers.size() ) );

  TPixel maxCVTIndex = *std::max_element( cvtCellIdentifiers.begin(),
                                          cvtCellIdentifiers.end() );
  TPixel minCVTIndex = *std::min_element( cvtCellIdentifiers.begin(),
                                          cvtCellIdentifiers.end() );

  try
    {
    itkAssertOrThrowMacro( static_cast<int>( maxCVTIndex - minCVTIndex )
      == static_cast<int>( cvtCellIdentifiers.size()-1 ),
      "CVT cell numbering not continuous!" );
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }

  tube::FmtDebugMessage( "Maximum CVT identifier = %d!",
    static_cast<int>( maxCVTIndex ) );


  // Create a mapping from CVT cell ID to the vector of corresponding
  // voxel indices.
  typedef typename InputImageType::IndexType                    IndexType;
  typedef typename std::map< TPixel, std::vector< IndexType > > MapType;

  MapType cvtToIndex;

  inImageIt.GoToBegin();
  while( !inImageIt.IsAtEnd() )
    {
    IndexType index = inImageIt.GetIndex();
    cvtToIndex[inImageIt.Get()].push_back( index );
    ++inImageIt;
    }


  // Determine the set of unique labels in the label map image.
  std::set< TPixel > labelSet;
  typename itk::ImageRegionIteratorWithIndex< InputImageType > inLabelImageIt(
    inLabelImage, inLabelImage->GetLargestPossibleRegion() );

  inLabelImageIt.GoToBegin();
  while( !inLabelImageIt.IsAtEnd() )
    {
    labelSet.insert( inLabelImageIt.Get() );
    ++inLabelImageIt;
    }

  TPixel maxLabel = *std::max_element( labelSet.begin(), labelSet.end() );

  tube::FmtDebugMessage( "Identified %d unique labels ( max. = %d )!",
    labelSet.size(), static_cast<int>( maxLabel ) );


  // Map { label_0, ...., label_N } -> {0, ...,  N-1}. This is convenient,
  // since it allows us to easily build a histogram later on.
  unsigned int labelCounter = 0;
  std::map< TPixel, unsigned int > labelRemapFwd;
  std::map< unsigned int, TPixel > labelRemapInv;

  typename std::set< TPixel >::const_iterator labelSetIt = labelSet.begin();
  while( labelSetIt != labelSet.end() )
    {
    labelRemapFwd[*( labelSetIt )] = labelCounter;
    labelRemapInv[labelCounter] = *( labelSetIt );
    ++labelCounter;
    ++labelSetIt;
    }

  typename std::map< TPixel, unsigned int >::iterator u = labelRemapFwd.begin();
  while( u != labelRemapFwd.end() )
    {
    tube::FmtDebugMessage( "%.4d mapsto %.4d",
      ( *u ).first, ( *u ).second );
    ++u;
    }


  // Create an empty output image.
  typename OutputImageType::Pointer outImage = OutputImageType::New();
  CreateEmptyImage< OutputImageType >( outImage,
    inImage->GetLargestPossibleRegion().GetSize() );


  // Count newly created IDs for CVT cells, other than segmentation labels
  unsigned int omittedRegionCounter = maxLabel + 1;
  tube::FmtDebugMessage( "OmittedRegionCounter init to %d",
    omittedRegionCounter );


  // Relabeling algorithm:
  //
  //   1 ) Create histogram of anatomical labels occurring in CVT cell
  //   2 ) Determine the anatomical label occurring most frequently
  //   3 ) Set voxel values ( of CVT cell ) in output image to dominant label
  std::ofstream outMappingFile;
  outMappingFile.open( argOutMappingFileName.c_str() );

  if( outMappingFile.fail() )
    {
    tube::FmtErrorMessage( "Could not open %s for writing!",
      argOutMappingFileName.c_str() );
    return EXIT_FAILURE;
    }
  outMappingFile << cvtCellIdentifiers.size() << std::endl;


  // Vector, holding the frequencies of label occurrence
  std::vector< unsigned int > cellHist( labelSet.size() );
  boost::dynamic_bitset<> omitLabel( labelSet.size() );

  // Set OMIT flag for certain labels
  for( unsigned int o=0; o<argOmit.size(); ++o )
    {
    omitLabel[labelRemapFwd[argOmit[o]]] = 1;
    }

  typename MapType::iterator mapIt = cvtToIndex.begin();
  while( mapIt != cvtToIndex.end() )
    {
    TPixel cell = ( *mapIt ).first;

    const std::vector< IndexType > & indexVector =
      ( *mapIt ).second;

    tube::FmtDebugMessage( "Cell %d: Index vector has size = %d!",
      static_cast<int>( cell ),
      static_cast<int>( indexVector.size() ) );

    // Build histogram of anatomical labels in current CVT cell
    for( unsigned int v=0; v<indexVector.size(); ++v )
      {
      IndexType index = indexVector[v];
      TPixel label = inLabelImage->GetPixel( index );
      cellHist[labelRemapFwd[label]] += 1;
      }

    // Find the dominant label
    unsigned int dominantLabel = distance(
      cellHist.begin(),
      std::max_element(
        cellHist.begin(),
        cellHist.end() ) );

    // Map the artificial label to the original one
    TPixel originalLabel = labelRemapInv[dominantLabel];
    tube::FmtDebugMessage( "Dominant label in CVT cell %d = %d ( %d )",
      cell, dominantLabel, labelRemapInv[dominantLabel] );

    if( omitLabel[dominantLabel] )
      {
      // If we omit a region, give it a unique label - by that we
      // can ensure that omitted background ( for instance ) cells do
      // not get just one background label
      originalLabel = omittedRegionCounter++;
      }

    for( unsigned int v=0; v < indexVector.size(); ++v )
      {
      IndexType index = indexVector[v];
      outImage->SetPixel( index, originalLabel );
      }

    if( argOutputHistogram )
      {
      for( unsigned int i=0; i<cellHist.size(); ++i )
        {
        outMappingFile << cellHist[i] << " ";
        }
      }
    else
      {
      outMappingFile << originalLabel;
      }
    outMappingFile << std::endl;

    // Reset cell histogram
    std::fill( cellHist.begin(), cellHist.end(), 0 );
    ++mapIt;
    }

  outMappingFile.close();
  tube::FmtDebugMessage( "Introduced %d new labels!",
    omittedRegionCounter - maxLabel - 1 );

  // Write relabeled CVT image to file
  typename OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName( argOutImageFileName.c_str() );
  writer->SetInput( outImage );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
