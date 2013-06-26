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

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeCropROI.h"
#include "tubeMessage.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "CropImageCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "Crop",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                PixelType;
  typedef itk::Image< PixelType,  TDimension >  ImageType;
  typedef itk::ImageFileReader< ImageType >     ReaderType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    tube::ErrorMessage( out.str() );
    timeCollector.Stop("Load data");
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  progressReporter.Report( 0.1 );

  if( size.size() > 0 || max.size() > 0 || min.size() > 0 )
    {
    if( size.size() > 0 && max.size() > 0 )
      {
      tube::ErrorMessage(
        "You must specify either --size or --max options.  Not both." );
      return EXIT_FAILURE;
      }

    if( center.size() > 0 && min.size() > 0 )
      {
      tube::ErrorMessage(
        "You must specify either --center or --min options.  Not both." );
      return EXIT_FAILURE;
      }

    timeCollector.Start("CropFilter");

    tube::CropROI< TPixel, TDimension > cropFilter;

    cropFilter.SetInput( reader->GetOutput() );
    if( min.size() > 0 )
      {
      typename ImageType::IndexType minI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        minI[i] = min[i];
        }
      cropFilter.SetMin( minI );
      }

    if( max.size() > 0 )
      {
      typename ImageType::IndexType maxI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        maxI[i] = max[i];
        }
      cropFilter.SetMax( maxI );
      }

    if( size.size() > 0 )
      {
      typename ImageType::SizeType sizeI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        sizeI[i] = size[i];
        }
      cropFilter.SetSize( sizeI );
      }

    if( center.size() > 0 )
      {
      typename ImageType::IndexType centerI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        centerI[i] = center[i];
        }
      cropFilter.SetCenter( centerI );
      }

    if( boundary.size() > 0 )
      {
      typename ImageType::IndexType boundaryI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        boundaryI[i] = boundary[i];
        }
      cropFilter.SetBoundary( boundaryI );
      }

    cropFilter.SetTimeCollector( &timeCollector );
    cropFilter.SetProgressReporter( &progressReporter, 0.1, 0.8 );

    try
      {
      cropFilter.Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::stringstream out;
      out << "Crop Filter: itk exception: ";
      out << e;
      tube::ErrorMessage( out.str() );
      timeCollector.Stop("CropFilter");
      throw( out.str() );
      }
    catch( const std::string & s )
      {
      std::cerr << "Error during crop filter: " << s << std::endl;
      timeCollector.Stop("CropFilter");
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr << "Error during crop filter" << std::endl;
      timeCollector.Stop("CropFilter");
      return EXIT_FAILURE;
      }
    timeCollector.Stop("CropFilter");

    typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

    timeCollector.Start("Save data");
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput( cropFilter.GetOutput() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Exception caught: " << err << std::endl;
      timeCollector.Stop("Save data");
      return EXIT_FAILURE;
      }
    timeCollector.Stop("Save data");
    }
  else if( split.size() == TDimension )
    {
    tube::CropROI< TPixel, TDimension > cropFilter;

    typename ImageType::Pointer inputImage = reader->GetOutput();
    typename ImageType::SizeType inputImageSize = inputImage->
                                          GetLargestPossibleRegion().
                                          GetSize();

    cropFilter.SetInput( inputImage );

    if( boundary.size() > 0 )
      {
      typename ImageType::IndexType boundaryI;
      for( unsigned int i=0; i<TDimension; i++ )
        {
        boundaryI[i] = boundary[i];
        }
      cropFilter.SetBoundary( boundaryI );
      }

    cropFilter.SetTimeCollector( &timeCollector );
    cropFilter.SetProgressReporter( &progressReporter, 0.1, 0.8 );

    typename ImageType::IndexType roiStep;
    typename ImageType::IndexType roiSize;
    for( unsigned int i=0; i<TDimension; i++ )
      {
      roiStep[i] = inputImageSize[i]/(split[i]+1);
      roiSize[i] = inputImageSize[i]/split[i];
      }
    typename ImageType::IndexType roiIndex;
    roiIndex.Fill( 0 );
    typename ImageType::IndexType roiMin;
    roiMin.Fill( 0 );
    typename ImageType::IndexType roiMax;
    roiMax.Fill( 0 );
    bool done = false;
    while( !done )
      {
      timeCollector.Start("CropFilter");

      for( unsigned int i=0; i<TDimension; i++ )
        {
        roiMin[i] = roiIndex[i] * roiSize[i];
        roiMax[i] = roiMin[i] + roiSize[i] - 1;
        if( roiIndex[i] == split[i]-1 )
          {
          roiMax[i] = inputImageSize[i]-1;
          }
        }

      cropFilter.SetMin( roiMin );
      cropFilter.SetMax( roiMax );
      try
        {
        cropFilter.Update();
        }
      catch( itk::ExceptionObject & e )
        {
        std::stringstream out;
        out << "Crop Filter: itk exception: ";
        out << e;
        tube::ErrorMessage( out.str() );
        timeCollector.Stop("CropFilter");
        throw( out.str() );
        }
      catch( const std::string & s )
        {
        std::stringstream out;
        out << "Error during crop filter: " << s << std::endl;
        tube::ErrorMessage( out.str() );
        timeCollector.Stop("CropFilter");
        return EXIT_FAILURE;
        }
      catch( ... )
        {
        std::cerr << "Error during crop filter" << std::endl;
        timeCollector.Stop("CropFilter");
        return EXIT_FAILURE;
        }

      timeCollector.Stop("CropFilter");

      timeCollector.Start("Save data");
      typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
      typename ImageWriterType::Pointer writer = ImageWriterType::New();

      std::stringstream out;
      out << outputVolume;
      out << "_";
      for( unsigned int i=0; i<TDimension; i++ )
        {
        out << roiIndex[i];
        }
      out << ".mha";
      writer->SetFileName( out.str() );

      writer->SetInput( cropFilter.GetOutput() );
      writer->SetUseCompression( true );
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Exception caught: " << err << std::endl;
        timeCollector.Stop("Save data");
        return EXIT_FAILURE;
        }
      timeCollector.Stop("Save data");

      unsigned int i=0;
      while( !done && ++roiIndex[i] >= split[i] )
        {
        roiIndex[i++] = 0;
        if( i >= TDimension )
          {
          done = true;
          }
        }
      }
    }

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
