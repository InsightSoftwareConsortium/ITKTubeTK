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

#include "itktubeTubeSpatialObjectToDensityImage.h"
#include "tubeCLIProgressReporter.h"

#include <itkImageFileReader.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include "ConvertTubesToDensityImageCLP.h"

enum { Dimension = 3 };

typedef unsigned int                                  DensityPixelType;
typedef float                                         RadiusPixelType;
typedef itk::Vector< float, Dimension >               TangentPixelType;

typedef itk::Image< unsigned char, Dimension >        TemplateImageType;
typedef itk::Image< DensityPixelType, Dimension >     DensityImageType;
typedef itk::Image< RadiusPixelType, Dimension >      RadiusImageType;
typedef itk::Image< TangentPixelType, Dimension >     TangentImageType;

typedef itk::ImageFileReader< TemplateImageType >     TemplateImageReaderType;

typedef itk::tube::TubeSpatialObjectToDensityImage<
  DensityImageType, RadiusImageType > TubeToDensityImageBuilderType;

typedef TubeToDensityImageBuilderType::TubeGroupType  TubesType;
typedef itk::SpatialObjectReader< Dimension >         TubesReaderType;

/** Max Intensity value */
DensityPixelType   max_densityIntensity = 2048;

/** Forward declarations */
TubesType::Pointer ReadTubes( const char * file );
void WriteImage( const char * file, DensityImageType::Pointer image );
void WriteImage( const char * file, RadiusImageType::Pointer image );
void WriteImage( const char * file, TangentImageType::Pointer image );
int DoIt( int argc, char * argv[] );


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}


/** Main work happens here */
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  double progress = 0.0;
  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter(
    "tubeDensityImageRadiusBuilder",
    CLPProcessInformation );

  progressReporter.Start();
  progressReporter.Report( progress );

  TubeToDensityImageBuilderType::Pointer
    builder = TubeToDensityImageBuilderType::New();

  builder->SetMaxDensityIntensity( max_densityIntensity ); // Const

  if( !inputTemplateImage.empty() )
    {
    std::cout << "Trying to use template image as constraints!" << std::endl;

    timeCollector.Start( "Loading template image" );

    TemplateImageReaderType::Pointer imTemplateReader;
    imTemplateReader = TemplateImageReaderType::New();
    imTemplateReader->SetFileName(inputTemplateImage.c_str());
    imTemplateReader->Update();

    TemplateImageType::Pointer imT = imTemplateReader->GetOutput();
    TubeToDensityImageBuilderType::SizeType size;
    double spacing[Dimension];
    for(int i = 0; i < Dimension; i++ )
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
    TubeToDensityImageBuilderType::SizeType sizeValue;
    for(int i = 0; i < Dimension; i++ )
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
    for(int i = 0; i < Dimension; i++ )
      {
      sp[i] = outputSpacing[i];
      }
    builder->SetSpacing( sp );
    }

  builder->UseSquareDistance( useSquareDistance );
  builder->SetTubes( ReadTubes( inputTubeFile.c_str() ) );

  progress = 0.1; // At about 10% done
  progressReporter.Report( progress );

  timeCollector.Start( "Update filter" );
  builder->Update();
  timeCollector.Stop( "Update filter" );

  progress = 0.8; // At about 80% done after filter
  progressReporter.Report( progress );

  timeCollector.Start( "Save data" );
  WriteImage( outputDensityImage.c_str(), builder->GetDensityMap() );
  WriteImage( outputRadiusImage.c_str(), builder->GetRadiusMap() );
  WriteImage( outputTangentImage.c_str(), builder->GetTangentMap() );
  timeCollector.Stop( "Save data" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}


TubesType::Pointer ReadTubes( const char * file )
{
  TubesReaderType::Pointer    reader = TubesReaderType::New();

  try
    {
    reader->SetFileName( file );
    std::cout << "Reading Tube group... ";
    reader->Update();
    std::cout << "Done." << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Error:: No readable Tubes found " << std::endl;
    return NULL;
    }

  return reader->GetGroup();
}


void WriteImage( const char * file, DensityImageType::Pointer image )
{
  std::cout << "Writing image: " << file << std::endl;
  typedef itk::ImageFileWriter<DensityImageType>    WriterType;
  WriterType::Pointer  writer = WriterType::New();

  writer->SetFileName( file );
  writer->SetInput( image );
  writer->SetUseCompression(true);
  writer->Update();
}


void WriteImage( const char * file, RadiusImageType::Pointer image )
{
  std::cout << "Writing image: " << file << std::endl;
  typedef itk::ImageFileWriter<RadiusImageType>    WriterType;
  WriterType::Pointer  writer = WriterType::New();

  writer->SetFileName( file );
  writer->SetInput( image );
  writer->SetUseCompression(true);
  writer->Update();
}


void WriteImage( const char * file, TangentImageType::Pointer image )
{
  std::cout << "Writing image: " << file <<std::endl;
  typedef itk::ImageFileWriter<TangentImageType>    WriterType;
  WriterType::Pointer  writer = WriterType::New();

  writer->SetFileName( file );
  writer->SetInput( image );
  writer->SetUseCompression(true);
  writer->Update();
}
