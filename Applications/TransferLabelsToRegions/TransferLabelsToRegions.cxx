/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
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
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkMacro.h"

#include "itkResampleImageFilter.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkTranslationTransform.h"
#include "itkImageFileWriter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkOrientImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkImageDuplicator.h"


#include "TransferLabelsToRegionsCLP.h"

// STL includes
#include <algorithm>
#include <vector>
#include <set>
#include <map>

template< class TPixel, unsigned int VImageDimension >
int DoIt( int argc, char * argv[] );

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"


int main( int argc, char **argv )
{
  PARSE_ARGS;
  return tube::ParseArgsAndCallDoIt( inImageFileName, argc, argv );

}


/** Check if image contains discrete values.
 *
 * \param fileName Image file name
 * \return 'true' if image, given by 'fileName' has a discrete-value type
 *
 */
bool IsDiscrete( const std::string & fileName )
{
   typedef itk::ImageIOBase::IOComponentType ScalarTPixelype;

  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
    fileName.c_str(), itk::ImageIOFactory::ReadMode);

  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();
  const ScalarTPixelype pixelType = imageIO->GetComponentType();

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
 *
 */
template <typename T>
bool check_vnl_vector_equality( const vnl_vector<T> &v,
                                const vnl_vector<T> &g,
                                double tol = 1.0e-6 )
{
  if( v.size() != g.size() )
    {
    return false;
    }
  for( unsigned i=0; i<v.size(); ++i)
    {
    if( vnl_math_abs(g.get(i) - v.get(i)) > tol )
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
 *
 */
template <typename T>
bool check_vnl_matrix_equality( const vnl_matrix<T> &V,
                                const vnl_matrix<T> &G,
                                double tol = 1.0e-6 )
{
  if( V.rows() != G.rows() || V.cols() != G.cols() )
    {
    return false;
    }
  for( unsigned int r=0; r<V.rows(); ++r)
    {
    for( unsigned int c=0; c<V.cols(); ++c)
      {
      if( vnl_math_abs(V.get(r,c) - G.get(r,c)) > tol )
        {
        return false;
        }
      }
    }
  return true;
}


/** Create an empty image.
 *
 *  \param outImage Image that is about to be created (has to exist)
 *  \param targetSize Desired size of the image
 *
 */
template < class ImageType >
void CreateEmptyImage(
  typename ImageType::Pointer &outImage,
  typename ImageType::SizeType targetSize)
{
  typename ImageType::IndexType start; start.Fill(0);
  typename ImageType::SizeType outImageSize = targetSize;

  typename ImageType::RegionType region( start, outImageSize );
  outImage->SetRegions( region );
  outImage->Allocate();
  outImage->FillBuffer(0);
}


/** Evaluate if two images are compatible.
 *  Compatibility is defined in terms of spacing, directions, size and origin.
 *
 *  \param imageA 1st input image of type ImageType
 *  \param imageB 2nd input image of type ImageType
 *  \returns true if image have equal spacing and size, false else
 *
 */
template < typename ImageType >
bool CheckCompatibility(
  typename ImageType::Pointer imageA,
  typename ImageType::Pointer imageB )
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
int DoIt( int argc, char **argv )
{
  PARSE_ARGS;

  /*
   * Some typedefs - Our input image type equals the output image type, since
   * the images are essentially the same, except that the voxel values
   * represent different entities.
   */
  typedef TPixel                                          InputTPixelype;
  typedef itk::Image< InputTPixelype, VImageDimension >   InputImageType;

  typedef InputTPixelype                                  OutputTPixelype;
  typedef itk::Image< OutputTPixelype, VImageDimension >  OutputImageType;

  typedef itk::ImageFileReader< InputImageType >          InputReaderType;
  typedef itk::ImageFileWriter< OutputImageType >         OutputWriterType;

  typedef itk::ImageIOBase::IOComponentType               ScalarTPixelype;

  // Check input images's type (we assume discrete valued images)
  if( !IsDiscrete( inImageFileName ) ||
      !IsDiscrete( inLabelFileName ) )
    {
    tube::ErrorMessage( "Non-discrete value types in input files!" );
    return EXIT_FAILURE;
    }


  /*
   * First, we try to read in the input image(s), i.e., an image of the Voronoi
   * tesselation and an image that contains discrete label values for each
   * voxel of the input image.
   */
  typename InputImageType::Pointer inImage;
  typename InputReaderType::Pointer inImageReader =
    InputReaderType::New( );

  inImageReader->SetFileName( inImageFileName.c_str( ) );
  try
    {
    inImageReader->Update( );
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  inImage = inImageReader->GetOutput( );

  typename InputImageType::Pointer inLabelImage;
  typename InputReaderType::Pointer inLabelImageReader =
    InputReaderType::New( );

  inLabelImageReader->SetFileName( inLabelFileName.c_str() );
  try
    {
    inLabelImageReader->Update( );
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  inLabelImage = inLabelImageReader->GetOutput( );


  /*
   * Since we will use the voxel indices, extracted from the input image, to
   * index voxels in the label image, we have to ensure compatibility w.r.t.
   * spacing, direction and size.
   */
  if( !CheckCompatibility< InputImageType >( inImage, inLabelImage ) )
    {
    tube::ErrorMessage( "Exiting due to image incompatibility!" );
    return EXIT_FAILURE;
    }


  /*
   * Determine the set of unique CVT cell IDs.
   */
  typename itk::ImageRegionIteratorWithIndex< InputImageType > inImageIt(
    inImage, inImage->GetLargestPossibleRegion());

  std::set< TPixel > regions;

  inImageIt.GoToBegin();
  while( !inImageIt.IsAtEnd() )
    {
    regions.insert( inImageIt.Get() );
    ++inImageIt;
    }

  tube::FmtInfoMessage("%d distinct CVT cells!",
    static_cast<int>( regions.size() ) );

  TPixel maxCVTIndex = *std::max_element( regions.begin(), regions.end() );
  TPixel minCVTIndex = *std::min_element( regions.begin(), regions.end() );

  try
    {
    // Assert if CVT cells are numbered  1 ... C
    itkAssertOrThrowMacro( static_cast<int>( maxCVTIndex - minCVTIndex )
      == static_cast<int>(regions.size()-1),
      "CVT cell numbering not continuous!" );
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }

  tube::FmtInfoMessage("Maximum CVT identifier = %d!",
    static_cast<int>( maxCVTIndex ) );


  /*
   * Create a mapping from CVT cell ID to the vector of corresponding
   * voxel indices.
   */
  typedef typename InputImageType::IndexType indexType;
  typedef typename std::map< TPixel, std::vector< indexType > > mapType;

  mapType cvtToIndex;

  inImageIt.GoToBegin();
  while( !inImageIt.IsAtEnd() )
   {
   indexType index = inImageIt.GetIndex();
   cvtToIndex[inImageIt.Get()].push_back( index );
   ++inImageIt;
   }


  /*
   * Determine the set of unique labels in the label map image.
   */
  std::set< TPixel > labelSet;
  typename itk::ImageRegionIteratorWithIndex< InputImageType > inLabelImageIt(
    inLabelImage, inLabelImage->GetLargestPossibleRegion());

  inLabelImageIt.GoToBegin();
  while( !inLabelImageIt.IsAtEnd() )
    {
    labelSet.insert( inLabelImageIt.Get() );
    ++inLabelImageIt;
    }

  TPixel maxLabel = *std::max_element( labelSet.begin(), labelSet.end() );

  tube::FmtDebugMessage( "Identified %d unique labels (max. = %d)!",
    labelSet.size(), static_cast<int>( maxLabel ) );


  /*
   * Map { label_0, ...., label_N } -> {0, ...,  N-1}. This is convenient,
   * since it allows us to easily build a histogram later on.
   */
  unsigned int labelCounter = 0;
  std::map< TPixel, unsigned int > labelRemap;
  std::map< unsigned int, TPixel > labelRemapInverse;

  typename std::set< TPixel >::const_iterator labelSetIt = labelSet.begin();
  while( labelSetIt != labelSet.end() )
    {
    labelRemap[*(labelSetIt)] = labelCounter;
    labelRemapInverse[labelCounter] = *(labelSetIt);
    ++labelCounter;
    ++labelSetIt;
    }

  typename std::map< TPixel, unsigned int >::iterator u = labelRemap.begin();
  while( u != labelRemap.end() )
    {
    tube::FmtDebugMessage("%.4d mapsto %.4d",
      (*u).first, (*u).second);
    ++u;
    }


  /*
   * Create an empty output image.
   */
  typename OutputImageType::Pointer outImage = OutputImageType::New();
  CreateEmptyImage< OutputImageType >( outImage,
    inImage->GetLargestPossibleRegion().GetSize() );


  /*
   * Relabeling algorithm:
   *
   *   1) Create histogram of anatomical labels occurring in CVT cell
   *   2) Determine the anatomical label occurring most frequently
   *   3) Set voxel values (of CVT cell) in output image to dominant label
   */
  std::ofstream outMappingFile;
  outMappingFile.open( outMappingFileName.c_str() );

  if( outMappingFile.fail() )
    {
    tube::FmtErrorMessage( "Could not open %s for writing!",
      outMappingFileName.c_str() );
    return EXIT_FAILURE;
    }

  // Vector, holding the frequencies of label occurrence
  std::vector< unsigned int > cellHist( labelSet.size() );

  typename mapType::iterator mapIt = cvtToIndex.begin();
  while( mapIt != cvtToIndex.end() )
    {
    TPixel cell = (*mapIt).first;

    tube::FmtDebugMessage("Processing cell %d!",
      static_cast<int>(cell) );

    const std::vector< indexType > & indexVector =
      (*mapIt).second;

    tube::FmtDebugMessage("Index vector has size = %d!",
      static_cast<int>( indexVector.size() ) );

    // Build histogram of anatomical labels in CVT cell
    for(unsigned int v=0; v<indexVector.size(); ++v)
      {
      indexType index = indexVector[v];
      TPixel label = inLabelImage->GetPixel( index );
      cellHist[labelRemap[label]] += 1;
      }

    // Determine position of highest freq.
    unsigned int dominantLabel = distance(
      cellHist.begin(),
      std::max_element(
        cellHist.begin(),
        cellHist.end() ) );

    // Map label back to original one
    TPixel originalLabel = labelRemapInverse[dominantLabel];
    for(unsigned int v=0; v < indexVector.size(); ++v)
      {
      indexType index = indexVector[v];
      outImage->SetPixel( index, originalLabel );
      }

    outMappingFile << cell << " ";
    if( argOutputHistogram )
      {
      for( unsigned int i=0; i<cellHist.size(); ++i)
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
    std::fill( cellHist.begin(), cellHist.end(), 0);
    ++mapIt;
    }

  outMappingFile.close();

  typename OutputWriterType::Pointer writer = OutputWriterType::New( );
  writer->SetFileName( outImageFileName.c_str( ) );
  writer->SetInput( outImage );
  writer->SetUseCompression( true );
  try
    {
    writer->Update( );
    }
  catch( itk::ExceptionObject &ex )
    {
    tube::ErrorMessage( ex.what() );
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
