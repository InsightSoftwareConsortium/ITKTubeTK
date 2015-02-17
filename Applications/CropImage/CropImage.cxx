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

#include "tubeCropROI.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkImageRegionConstIterator.h>

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

    tube::CropROI< TPixel, VDimension > cropFilter;

    cropFilter.SetInput( reader->GetOutput() );
    if( matchVolume.size() > 0 )
      {
      typename ImageType::Pointer image = reader->GetOutput();

      typename ReaderType::Pointer matchReader = ReaderType::New();
      matchReader->SetFileName( matchVolume.c_str() );
      matchReader->UpdateOutputInformation();
      typename ImageType::ConstPointer match = matchReader->GetOutput();

      const typename ImageType::RegionType matchRegion =
        match->GetLargestPossibleRegion();
      const typename ImageType::IndexType matchIndex =
        matchRegion.GetIndex();
      const typename ImageType::SizeType matchSize =
        matchRegion.GetSize();
      const typename ImageType::PointType matchOrigin =
        match->GetOrigin();
      const typename ImageType::SpacingType matchSpacing =
        match->GetSpacing();

      const typename ImageType::RegionType imgRegion =
        image->GetLargestPossibleRegion();
      const typename ImageType::IndexType imgIndex =
        imgRegion.GetIndex();
      const typename ImageType::PointType imgOrigin =
        image->GetOrigin();
      const typename ImageType::SpacingType imgSpacing =
        image->GetSpacing();

      typename ImageType::IndexType minI;
      typename ImageType::SizeType sizeI;

      if( imgOrigin != matchOrigin || imgSpacing != matchSpacing )
        {
        for( unsigned int i = 0; i < VDimension; i++ )
          {
          minI[i] = vnl_math_rnd( ( ( matchOrigin[i]
            + matchIndex[i] * matchSpacing[i] ) - ( imgOrigin[i]
            + imgIndex[i] * imgSpacing[i] ) )
            / imgSpacing[i] );
          sizeI[i] = vnl_math_rnd( ( matchSize[i] * matchSpacing[i] )
            / imgSpacing[i] );
          }
        }
      else
        {
        for( unsigned int i = 0; i < VDimension; i++ )
          {
          minI[i] = matchIndex[i];
          sizeI[i] = matchSize[i];
          }
        }
      cropFilter.SetMin( minI );
      cropFilter.SetSize( sizeI );
      }

    if( matchMask.size() > 0 )
      {
      timeCollector.Start( "Mask Bounding Box" );

      typename ReaderType::Pointer maskReader = ReaderType::New();
      maskReader->SetFileName( matchMask.c_str() );
      maskReader->Update();
      typename ImageType::Pointer maskImage = maskReader->GetOutput();

      typename ImageType::IndexType minI;
      typename ImageType::IndexType maxI;

      itk::ImageRegionConstIterator< ImageType > it( maskImage,
        maskImage->GetLargestPossibleRegion() );
      while( !it.IsAtEnd() && it.Get() == 0 )
        {
        ++it;
        }
      minI = it.GetIndex();
      while( !it.IsAtEnd() && it.Get() != 0 )
        {
        ++it;
        }
      maxI = it.GetIndex();
      while( !it.IsAtEnd() )
        {
        while( !it.IsAtEnd() && it.Get() == 0 )
          {
          ++it;
          }
        if( !it.IsAtEnd() )
          {
          for( unsigned int i=0; i<VDimension; ++i )
            {
            if( it.GetIndex()[i] < minI[i] )
              {
              minI[i] = it.GetIndex()[i];
              }
            }
          }
        while( !it.IsAtEnd() && it.Get() != 0 )
          {
          ++it;
          }
        if( !it.IsAtEnd() )
          {
          for( unsigned int i=0; i<VDimension; ++i )
            {
            if( it.GetIndex()[i] > maxI[i] )
              {
              maxI[i] = it.GetIndex()[i];
              }
            }
          }
        }

      typename ImageType::SizeType sizeI;
      for( unsigned int i=0; i<VDimension; ++i )
        {
        sizeI[i] = maxI[i] - minI[i];
        }

      cropFilter.SetMin( minI );
      cropFilter.SetSize( sizeI );
      timeCollector.Stop( "Mask Bounding Box" );
      }

    if( min.size() > 0 )
      {
      typename ImageType::IndexType minI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        minI[i] = min[i];
        }
      cropFilter.SetMin( minI );
      }

    if( max.size() > 0 )
      {
      typename ImageType::IndexType maxI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        maxI[i] = max[i];
        }
      cropFilter.SetMax( maxI );
      }

    if( size.size() > 0 )
      {
      typename ImageType::SizeType sizeI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        sizeI[i] = size[i];
        }
      cropFilter.SetSize( sizeI );
      }

    if( center.size() > 0 )
      {
      typename ImageType::IndexType centerI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        centerI[i] = center[i];
        }
      cropFilter.SetCenter( centerI );
      }

    if( boundary.size() > 0 )
      {
      typename ImageType::IndexType boundaryI;
      for( unsigned int i = 0; i < VDimension; i++ )
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

    typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

    timeCollector.Start( "Save data" );
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
      timeCollector.Stop( "Save data" );
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Save data" );
    }
  else if( split.size() == VDimension )
    {
    tube::CropROI< TPixel, VDimension > cropFilter;

    typename ImageType::Pointer inputImage = reader->GetOutput();
    typename ImageType::SizeType inputImageSize = inputImage->
                                          GetLargestPossibleRegion().
                                          GetSize();

    cropFilter.SetInput( inputImage );

    if( boundary.size() > 0 )
      {
      typename ImageType::IndexType boundaryI;
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        boundaryI[i] = boundary[i];
        }
      cropFilter.SetBoundary( boundaryI );
      }

    cropFilter.SetTimeCollector( &timeCollector );
    cropFilter.SetProgressReporter( &progressReporter, 0.1, 0.8 );

    typename ImageType::IndexType roiStep;
    typename ImageType::IndexType roiSize;
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      roiStep[i] = inputImageSize[i] / ( split[i] + 1 );
      roiSize[i] = inputImageSize[i] / split[i];
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
      timeCollector.Start( "CropFilter" );

      for( unsigned int i = 0; i < VDimension; i++ )
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

      writer->SetInput( cropFilter.GetOutput() );
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
