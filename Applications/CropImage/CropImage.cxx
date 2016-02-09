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
#include "tubeCLIProgressReporter.h"

#include "itktubeCropImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkTimeProbesCollectorBase.h>

using tube::CLIProgressReporter;

#include "CropImageCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "Crop",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                               PixelType;
  typedef itk::Image< PixelType, VDimension >  ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;

  timeCollector.Start( "Load data" );
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
    timeCollector.Stop( "Load data" );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );
  progressReporter.Report( 0.1 );

  if( size.size() > 0 || max.size() > 0 || min.size() > 0 ||
    matchVolume.size() > 0 || matchMask.size() > 0 )
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

    timeCollector.Start( "CropFilter" );

    typedef itk::tube::CropImageFilter< ImageType, ImageType > CropFilterType;
    typename CropFilterType::Pointer cropFilter = CropFilterType::New();

    cropFilter->SetInput( reader->GetOutput() );

    if( matchVolume.size() > 0 )
      {
      typename ReaderType::Pointer matchReader = ReaderType::New();
      matchReader->SetFileName( matchVolume.c_str() );
      matchReader->UpdateOutputInformation();

      cropFilter->SetMatchVolume( matchReader->GetOutput() );
      }

    if( matchMask.size() > 0 )
      {
      timeCollector.Start( "Mask Bounding Box" );

      typename ReaderType::Pointer maskReader = ReaderType::New();
      maskReader->SetFileName( matchMask.c_str() );
      maskReader->Update();

      cropFilter->SetMatchMask( maskReader->GetOutput() );

      timeCollector.Stop( "Mask Bounding Box" );
      }

    if( min.size() > 0 )
      {
      typename ImageType::IndexType minI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        minI[i] = min[i];
        }
      cropFilter->SetMin( minI );
      }

    if( max.size() > 0 )
      {
      typename ImageType::IndexType maxI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        maxI[i] = max[i];
        }
      cropFilter->SetMax( maxI );
      }

    if( size.size() > 0 )
      {
      typename ImageType::SizeType sizeI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        sizeI[i] = size[i];
        }
      cropFilter->SetSize( sizeI );
      }

    if( center.size() > 0 )
      {
      typename ImageType::IndexType centerI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        centerI[i] = center[i];
        }
      cropFilter->SetCenter( centerI );
      }

    if( boundary.size() > 0 )
      {
      typename ImageType::IndexType boundaryI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        boundaryI[i] = boundary[i];
        }
      cropFilter->SetBoundary( boundaryI );
      }

    progressReporter.Report( 0.6 );

    try
      {
      cropFilter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::stringstream out;
      out << "Crop Filter: itk exception: ";
      out << e;
      tube::ErrorMessage( out.str() );
      timeCollector.Stop( "CropFilter" );
      throw( out.str() );
      }
    catch( const std::string & s )
      {
      std::cerr << "Error during crop filter: " << s << std::endl;
      timeCollector.Stop( "CropFilter" );
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr << "Error during crop filter" << std::endl;
      timeCollector.Stop( "CropFilter" );
      return EXIT_FAILURE;
      }

    timeCollector.Stop( "CropFilter" );
    progressReporter.Report( 0.8 );
    timeCollector.Start( "Save data" );

    typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput( cropFilter->GetOutput() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "Exception caught: " << err << std::endl;
      timeCollector.Stop( "Save data" );
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Save data" );
    }
  else if( split.size() == VDimension )
    {
    typename ImageType::IndexType splitI;
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      splitI[i] = split[i];
      }

    typename ImageType::Pointer inputImage = reader->GetOutput();

    typename ImageType::IndexType roiIndex;
    roiIndex.Fill( 0 );

    bool done = false;
    // Update filter until we get all split part of the inputImage
    while( !done )
      {
      //Create new instance of cropFilter each time we extract a portion of
      // the image to reset parameters
      typedef itk::tube::CropImageFilter<ImageType, ImageType> CropFilterType;
      typename CropFilterType::Pointer cropFilter = CropFilterType::New();
      cropFilter->SetInput( inputImage );

      if( boundary.size() > 0 )
        {
        typename ImageType::IndexType boundaryI;
        for( unsigned int i = 0; i < VDimension; i++ )
          {
          boundaryI[i] = boundary[i];
          }
        cropFilter->SetBoundary( boundaryI );
        }

      progressReporter.Report( 0.6 );

      timeCollector.Start( "CropFilter" );

      cropFilter->SetSplitInput( splitI, roiIndex );

      try
        {
        cropFilter->Update();
        }
      catch( itk::ExceptionObject & e )
        {
        std::stringstream out;
        out << "Crop Filter: itk exception: ";
        out << e;
        tube::ErrorMessage( out.str() );
        timeCollector.Stop( "CropFilter" );
        throw( out.str() );
        }
      catch( const std::string & s )
        {
        std::stringstream out;
        out << "Error during crop filter: " << s << std::endl;
        tube::ErrorMessage( out.str() );
        timeCollector.Stop( "CropFilter" );
        return EXIT_FAILURE;
        }
      catch( ... )
        {
        std::cerr << "Error during crop filter" << std::endl;
        timeCollector.Stop( "CropFilter" );
        return EXIT_FAILURE;
        }

      timeCollector.Stop( "CropFilter" );
      progressReporter.Report( 0.8 );
      timeCollector.Start( "Save data" );

      typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
      typename ImageWriterType::Pointer writer = ImageWriterType::New();

      std::stringstream out;
      out << outputVolume;
      out << "_";
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        out << roiIndex[i];
        }
      out << ".mha";
      writer->SetFileName( out.str() );

      writer->SetInput( cropFilter->GetOutput() );
      writer->SetUseCompression( true );
      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "Exception caught: " << err << std::endl;
        timeCollector.Stop( "Save data" );
        return EXIT_FAILURE;
        }
      timeCollector.Stop( "Save data" );

      //Update ROIIndex value to get the next split image
      unsigned int i=0;
      while( !done && ++roiIndex[i] >= split[i] )
        {
        roiIndex[i++] = 0;
        if( i >= VDimension )
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
