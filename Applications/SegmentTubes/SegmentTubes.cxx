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

#include "tubeSegmentTubes.h"
#include "itktubeTubeExtractorIO.h"

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "SegmentTubesCLP.h"

#include <sstream>

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and
//   forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "RidgeExtractor",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                     PixelType;
  typedef itk::Image< PixelType, VDimension >        ImageType;
  typedef itk::ImageFileReader< ImageType >          ReaderType;

  typedef tube::SegmentTubes< ImageType >           TubeOpType;
  typename TubeOpType::Pointer tubeOp = TubeOpType::New();

  typedef typename TubeOpType::TubeMaskImageType     MaskImageType;
  typedef itk::ImageFileReader< MaskImageType >      MaskReaderType;
  typedef itk::ImageFileWriter< MaskImageType >      MaskWriterType;

  typedef float                                      ScaleType;
  typedef itk::Image< ScaleType, VDimension >        ScaleImageType;
  typedef itk::ImageFileReader< ScaleImageType >     ScaleReaderType;

  typedef itk::ContinuousIndex< double, VDimension > IndexType;
  typedef std::vector< IndexType >                   IndexListType;

  typedef std::vector< ScaleType >                   ScaleListType;

  typedef itk::VesselTubeSpatialObject< VDimension > TubeType;

  typedef typename TubeType::TransformType           TransformType;

  typedef itk::SpatialObjectWriter< VDimension >     TubesWriterType;
  typedef itk::SpatialObjectReader< VDimension >     TubesReaderType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer inputImage = reader->GetOutput();
  tubeOp->SetInputImage( inputImage );

  if( !radiusInputVolume.empty() )
    {
    typename ReaderType::Pointer radiusReader = ReaderType::New();
    radiusReader->SetFileName( radiusInputVolume.c_str() );
    try
      {
      radiusReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading radius volume: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    tubeOp->SetRadiusInputImage( radiusReader->GetOutput() );
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  double scaleNorm = inputImage->GetSpacing()[0];
  if( scale / scaleNorm < 0.3 )
    {
    tube::ErrorMessage(
      "Error: Scale < 0.3 * voxel spacing is unsupported." );
    return EXIT_FAILURE;
    }
  tubeOp->SetRadius( scale / scaleNorm );

  IndexType seedIndex;
  IndexListType seedIndexList;
  seedIndexList.clear();

  ScaleListType seedRadiusList;
  seedRadiusList.clear();

  if( !seedI.empty() )
    {
    for( unsigned int seedINum = 0; seedINum < seedI.size(); ++seedINum )
      {
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        seedIndex[i] = seedI[seedINum][i];
        }
      seedIndexList.push_back( seedIndex );
      seedRadiusList.push_back( scale / scaleNorm );
      }
    }

  if( !seedP.empty() )
    {
    for( size_t seedNum = 0; seedNum < seedP.size(); ++seedNum )
      {
      typename ImageType::PointType point;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        point[i] = seedP[seedNum][i];
        }

      bool transformSuccess =
        inputImage->TransformPhysicalPointToContinuousIndex( point, seedIndex );
      if ( !transformSuccess )
        {
        std::cerr<<"Could not transform point #"
          << seedNum <<" to seed index." << std::endl;
        continue;
        }

      seedIndexList.push_back( seedIndex );
      seedRadiusList.push_back( scale / scaleNorm );
      }
    }

  if( !seedListFile.empty() )
    {
    std::ifstream readStream;
    readStream.open( seedListFile.c_str(), std::ios::binary |
      std::ios::in );
    std::string line;
    double seedScale;
    while( std::getline( readStream, line ) )
      {
      std::istringstream iss( line );
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        iss >> seedIndex[i];
        }
      iss >> seedScale;
      }
    seedIndexList.push_back( seedIndex );
    seedRadiusList.push_back( seedScale / scaleNorm );
    }

  if( !seedMask.empty() )
    {
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( seedMask.c_str() );
    try
      {
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading mask: Exception caught: "
                          + std::string(err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    typename MaskImageType::Pointer seedMaskImage = maskReader->GetOutput();

    if( !scaleMask.empty() )
      {
      typename ScaleReaderType::Pointer scaleReader =
        ScaleReaderType::New();
      scaleReader->SetFileName( scaleMask.c_str() );
      try
        {
        scaleReader->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        tube::ErrorMessage( "Reading scale: Exception caught: "
                            + std::string(err.GetDescription()) );
        timeCollector.Report();
        return EXIT_FAILURE;
        }
      typename ScaleImageType::Pointer scaleImage = scaleReader->GetOutput();

      itk::ImageRegionConstIteratorWithIndex< MaskImageType > iter(
        seedMaskImage, seedMaskImage->GetLargestPossibleRegion() );
      itk::ImageRegionConstIterator< ScaleImageType > iterS( scaleImage,
        scaleImage->GetLargestPossibleRegion() );
      int count = 0;
      while( !iter.IsAtEnd() )
        {
        if( iter.Get() )
          {
          if( ++count == seedMaskStride )
            {
            count = 0;
            seedIndexList.push_back( iter.GetIndex() );
            seedRadiusList.push_back( iterS.Get() / scaleNorm );
            }
          }
        ++iter;
        ++iterS;
        }
      }
    else
      {
      int count = 0;
      itk::ImageRegionConstIteratorWithIndex< MaskImageType > iter(
        seedMaskImage, seedMaskImage->GetLargestPossibleRegion() );
      while( !iter.IsAtEnd() )
        {
        if( iter.Get() )
          {
          if( ++count == seedMaskStride )
            {
            count = 0;
            seedIndexList.push_back( iter.GetIndex() );
            seedRadiusList.push_back( scale / scaleNorm );
            }
          }
        ++iter;
        }
      }
    }

  if( !existingVesselsMask.empty() )
    {
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( existingVesselsMask.c_str() );
    try
      {
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading vessels mask: Exception caught: "
                          + std::string(err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    typename MaskImageType::Pointer maskImage = maskReader->GetOutput();
    tubeOp->SetTubeMaskImage( maskImage );
    }

  if( !existingVessels.empty() )
    {
    typedef typename TubeType::ChildrenListType  ChildrenListType;
    typedef typename ChildrenListType::iterator  ChildrenIteratorType;

    typename TubesReaderType::Pointer tubeReader = TubesReaderType::New();

    try
      {
      tubeReader->SetFileName( existingVessels.c_str() );
      tubeReader->Update();
      }
    catch( ... )
      {
      std::cerr << "Error:: No readable Tubes found " << std::endl;
      }

    char tubeName[] = "Tube";
    ChildrenListType* tubeList = tubeReader->GetGroup()->GetChildren( 999999,
      tubeName );

    ChildrenIteratorType tubeIterator = tubeList->begin();
    while( tubeIterator != tubeList->end() )
      {
      tubeOp->AddTube( static_cast< TubeType * >(
          tubeIterator->GetPointer() ) );
      ++tubeIterator;
      }

    delete tubeList;
    }

  if( !parametersFile.empty() )
    {
    itk::tube::TubeExtractorIO< ImageType > teReader;
    teReader.SetTubeExtractor( tubeOp->GetTubeOp() );
    teReader.Read( parametersFile.c_str() );
    }

  tubeOp->SetDebug( false );
  tubeOp->GetRidgeOp()->SetDebug( false );
  tubeOp->GetRadiusOp()->SetDebug( false );

  if( border > 0 )
    {
    typename ImageType::IndexType minIndx = inputImage->
      GetLargestPossibleRegion().GetIndex();
    typename ImageType::SizeType size = inputImage->
      GetLargestPossibleRegion().GetSize();
    typename ImageType::IndexType maxIndx = minIndx + size;
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      minIndx[i] += border;
      maxIndx[i] -= border;
      }
    tubeOp->SetExtractBoundMin( minIndx );
    tubeOp->SetExtractBoundMax( maxIndx );
    }

  timeCollector.Start( "Ridge Extractor" );
  typename IndexListType::iterator seedIndexIter =
    seedIndexList.begin();
  ScaleListType::iterator seedRadiusIter =
    seedRadiusList.begin();
  unsigned int count = 1;
  bool foundOneTube = false;
  while( seedIndexIter != seedIndexList.end() )
    {
    tubeOp->SetRadius( *seedRadiusIter );

    std::cout << "Extracting from index point " << *seedIndexIter
      << " at radius " << *seedRadiusIter << std::endl;
    typename TubeType::Pointer xTube =
      tubeOp->ExtractTube( *seedIndexIter, count, true );
    if( !xTube.IsNull() )
      {
      tubeOp->AddTube( xTube );
      std::cout << "  Extracted " << xTube->GetPoints().size() << " points."
        << std::endl;
      foundOneTube = true;
      }
    else
      {
      std::stringstream ss;
      ss << "Error: Ridge not found for seed #" << count;
      tube::Message( ss.str() );
      }

    ++seedIndexIter;
    ++seedRadiusIter;
    ++count;
    }

  if ( !foundOneTube )
    {
    tube::ErrorMessage("No Ridge found at all");
    return EXIT_FAILURE;
    }

  // Update tubes transform
  typename TransformType::InputVectorType scaleVector;
  typename TransformType::OffsetType offsetVector;
  typename TransformType::MatrixType directionMatrix;
  typename ImageType::SpacingType spacing = inputImage->GetSpacing();
  typename ImageType::PointType origin = inputImage->GetOrigin();

  for (unsigned int i = 0; i < VDimension; ++i)
    {
    scaleVector[i] = spacing[i];
    offsetVector[i] = origin[i];
    }

  tubeOp->GetTubeGroup()->GetObjectToParentTransform()->SetScale(
    scaleVector );
  tubeOp->GetTubeGroup()->GetObjectToParentTransform()->SetOffset(
    offsetVector );
  tubeOp->GetTubeGroup()->GetObjectToParentTransform()->SetMatrix(
    inputImage->GetDirection() );
  tubeOp->GetTubeGroup()->ComputeObjectToWorldTransform();

  // Save Tubes
  typename TubesWriterType::Pointer soWriter = TubesWriterType::New();
  soWriter->SetFileName( outputTubeFile );
  soWriter->SetInput( tubeOp->GetTubeGroup().GetPointer() );
  soWriter->Update();

  // Save Tube Mask Image
  if( !outputTubeImage.empty() )
    {
    typename MaskWriterType::Pointer writer = MaskWriterType::New();
    writer->SetFileName( outputTubeImage );
    writer->SetInput( tubeOp->GetTubeMaskImage() );
    writer->SetUseCompression( true );
    writer->Update();
    }

  std::cout << "Ridge termination code counts:" << std::endl;
  for( unsigned int code = 0; code < tubeOp->GetRidgeOp()
    ->GetNumberOfFailureCodes(); ++code )
    {
    std::cout << "   " << tubeOp->GetRidgeOp()->GetFailureCodeName(
      typename TubeOpType::RidgeExtractorFilterType::FailureCodeEnum( code ) ) << " : "
      << tubeOp->GetRidgeOp()->GetFailureCodeCount(
      typename TubeOpType::RidgeExtractorFilterType::FailureCodeEnum(
      code ) ) << std::endl;
    }


  timeCollector.Stop("Ridge Extractor");

  progressReporter.Report( 1.0 );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
