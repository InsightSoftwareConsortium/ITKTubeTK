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

#include "itktubeRidgeExtractor.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkGroupSpatialObject.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTubeSpatialObjectPoint.h>
#include <itkVesselTubeSpatialObject.h>

#include "SegmentTubesCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
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

  typedef float                                         PixelType;
  typedef itk::Image< PixelType, VDimension >           ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;

  typedef itk::SpatialObject< VDimension >              SpatialObjectType;
  typedef typename SpatialObjectType::ChildrenListType  ObjectListType;
  typedef itk::GroupSpatialObject< VDimension >         GroupType;
  typedef itk::VesselTubeSpatialObject< VDimension >    TubeType;
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
    if( seedX.size() != VDimension )
      {
      std::cout << "seedX = ";
      for( unsigned int i=0; i<seedX.size(); i++ )
        {
        std::cout << seedX[i] << " ";
        }
      std::cout << std::endl;
      std::cout << "VDimension = " << VDimension << std::endl;
      std::cout << "size = " << seedX.size() << std::endl;
      tube::ErrorMessage(
        "Errror: X vector must be specified to initiate an extraction." );
      return EXIT_FAILURE;
      }

    timeCollector.Start("Ridge Extractor");

    itk::ContinuousIndex< double, VDimension > cIndx;
    for( unsigned int i=0; i<VDimension; i++ )
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
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
