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

#include "tubeConvertTubesToDensityImage.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include <itkImageFileReader.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "ConvertTubesToDensityImageCLP.h"

int DoIt( int argc, char * argv[] );

/** Main work happens here */
template< unsigned int Dimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  /*Typedefs..*/
  typedef float                                         TPixel;
  typedef itk::Vector< TPixel, Dimension >              TangentPixelType;
  typedef itk::Image< TPixel, Dimension >               DensityImageType;
  typedef itk::Image< TPixel, Dimension >               RadiusImageType;
  typedef itk::Image< TangentPixelType, Dimension >     TangentImageType;
  typedef itk::Image< unsigned char, Dimension >        TemplateImageType;
  typedef itk::ImageFileReader< TemplateImageType >     TemplateImageReaderType;

  /** Max Intensity value */
  TPixel   max_densityIntensity = 2048;

  typedef tube::ConvertTubesToDensityImage<
  TPixel, Dimension > TubeToDensityImageBuilderType;

  typedef itk::SpatialObjectReader< Dimension >         TubesReaderType;


  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter( 
    "tubeDensityImageRadiusBuilder",
    CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  typename TubeToDensityImageBuilderType::Pointer
    builder = TubeToDensityImageBuilderType::New();

  builder->SetMaxDensityIntensity( max_densityIntensity ); // Const

  if( !inputTemplateImage.empty() )
    {
    std::cout << "Trying to use template image as constraints!" << std::endl;

    timeCollector.Start( "Loading template image" );

    typename TemplateImageReaderType::Pointer imTemplateReader;
    imTemplateReader = TemplateImageReaderType::New();
    imTemplateReader->SetFileName( inputTemplateImage.c_str() );
    imTemplateReader->Update();

    typename TemplateImageType::Pointer imT = imTemplateReader->GetOutput();
    typename TubeToDensityImageBuilderType::SizeType size;
    double spacing[Dimension];
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      size[i] = imT->GetLargestPossibleRegion().GetSize()[i];
      spacing[i] = imT->GetSpacing()[i];
      }
    builder->SetSize( size );
    builder->SetSpacing( spacing );

    timeCollector.Stop( "Loading template image" );
    }
  else
    {
    std::cout << "Trying to use user-specified constraints!" << std::endl;

    if( !outputSize.size() )
      {
      std::cerr << "Output size is missing!" << std::endl;
      return -1;
      }
    typename TubeToDensityImageBuilderType::SizeType sizeValue;
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      sizeValue[i] = outputSize[i];
      }
    builder->SetSize( sizeValue );

    if( !outputSpacing.size() )
      {
      std::cerr << "Output spacing is missing!" << std::endl;
      return -1;
      }
    double sp[Dimension];
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      sp[i] = outputSpacing[i];
      }
    builder->SetSpacing( sp );
    }

  builder->SetUseSquareDistance( useSquareDistance );
  typename TubesReaderType::Pointer reader = TubesReaderType::New();
  try
    {
    reader->SetFileName( inputTubeFile.c_str() );
    std::cout << "Reading Tube group... ";
    reader->Update();
    std::cout << "Done." << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Error:: No readable Tubes found " << std::endl;
    return EXIT_FAILURE;
    }
  builder->SetInputTubeGroup( reader->GetGroup() );

  progress = 0.1; // At about 10% done
  progressReporter.Report( progress );

  timeCollector.Start( "Update filter" );
  builder->Update();
  timeCollector.Stop( "Update filter" );

  progress = 0.8; // At about 80% done after filter
  progressReporter.Report( progress );

  timeCollector.Start( "Save data" );
  std::cout << "Writing image: " << outputDensityImage.c_str() << std::endl;
  typedef itk::ImageFileWriter< DensityImageType > WriterType_d;
  typename WriterType_d::Pointer  writer_d = WriterType_d::New();

  writer_d->SetFileName( outputDensityImage.c_str() );
  writer_d->SetInput( builder->GetDensityMapImage() );
  writer_d->SetUseCompression( true );
  writer_d->Update();

  std::cout << "Writing image: " << outputRadiusImage.c_str() << std::endl;
  typedef itk::ImageFileWriter< RadiusImageType > WriterType_r;
  typename WriterType_r::Pointer  writer_r = WriterType_r::New();

  writer_r->SetFileName( outputRadiusImage.c_str() );
  writer_r->SetInput( builder->GetRadiusMapImage() );
  writer_r->SetUseCompression( true );
  writer_r->Update();

  std::cout << "Writing image: " << outputTangentImage.c_str() << std::endl;
  typedef itk::ImageFileWriter< TangentImageType > WriterType_t;
  typename WriterType_t::Pointer  writer_t = WriterType_t::New();

  writer_t->SetFileName( outputTangentImage.c_str() );
  writer_t->SetInput( builder->GetTangentMapImage() );
  writer_t->SetUseCompression( true );
  writer_t->Update();

  timeCollector.Stop( "Save data" );

  progress = 1.0;
  progressReporter.Report( progress );
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

  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTubeFile.c_str() );

  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input Tube file has no spatial objects" );
    delete mScene;
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 2:
      {
      bool result = DoIt<2>( argc, argv );
      delete mScene;
      return result;
      break;
      }
    case 3:
      {
      bool result = DoIt<3>( argc, argv );
      delete mScene;
      return result;
      break;
      }
    default:
      {
      tubeErrorMacro( 
        << "Error: Only 2D and 3D data is currently supported." );
      delete mScene;
      return EXIT_FAILURE;
      break;
      }
    }
  return EXIT_FAILURE;
}
