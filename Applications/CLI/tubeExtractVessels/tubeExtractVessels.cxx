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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkTubeRidgeExtractor.h"
#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkVesselTubeSpatialObject.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "tubeExtractVesselsCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
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

  typedef float                                         PixelType;
  typedef itk::Image< PixelType,  dimensionT >          ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;

  typedef itk::SpatialObject< dimensionT >              SpatialObjectType;
  typedef typename SpatialObjectType::ChildrenListType  ObjectListType;
  typedef itk::GroupSpatialObject< dimensionT >         GroupType;
  typedef itk::VesselTubeSpatialObject< dimensionT >    TubeType;
  typedef typename TubeType::PointListType              PointListType;
  typedef typename TubeType::PointType                  PointType;
  typedef typename TubeType::TubePointType              TubePointType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename ImageType::Pointer inputImage = reader->GetOutput();

  if( scale < 0.3 )
    {
    tube::ErrorMessage( "Errror: Scale < 0.3 is unsupported." );
    return EXIT_FAILURE;
    }

  typedef itk::tube::RidgeExtractor< ImageType > RidgeOpType;
  typename RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  ridgeOp->SetInputImage( inputImage );
  ridgeOp->SetScale( scale );
  ridgeOp->SetExtent( 3.0 );
  ridgeOp->SetStepX( 0.75 );
  ridgeOp->SetDynamicScale( true );

  if( seedX.size() > 0 )
    {
    if( seedX.size() != dimensionT )
      {
      std::cout << "seedX = ";
      for( unsigned int i=0; i<seedX.size(); i++ )
        {
        std::cout << seedX[i] << " ";
        }
      std::cout << std::endl;
      std::cout << "dimensionT = " << dimensionT << std::endl;
      std::cout << "size = " << seedX.size() << std::endl;
      tube::ErrorMessage(
        "Errror: X vector must be specified to initiate an extraction." );
      return EXIT_FAILURE;
      }

    timeCollector.Start("Ridge Extractor");

    itk::ContinuousIndex< double, dimensionT > cIndx;
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      cIndx[i] = seedX[i];
      }

    typename TubeType::Pointer xTube = ridgeOp->Extract( cIndx, 1 );

    if( xTube.IsNull() )
      {
      tube::ErrorMessage( "Errror: Ridge not found. " );
      return EXIT_FAILURE;
      }
    }

  progressReporter.Report( 0.2 );
  timeCollector.Stop("Ridge Extractor");

  progressReporter.Report( 1.0 );
  progressReporter.End( );

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
