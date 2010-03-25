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

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkImageRegionIteratorWithIndex.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "MergeCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "Merge",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef pixelT                                        PixelType;
  typedef itk::Image< PixelType,  dimensionT >          ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;
  
  timeCollector.Start("Load data");

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( inputVolume1.c_str() );
  try
    {
    reader1->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( inputVolume2.c_str() );
  try
    {
    reader2->Update();
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

  typename ImageType::Pointer curImage1 = reader1->GetOutput();
  typename ImageType::Pointer curImage2 = reader2->GetOutput();

  typename ImageType::Pointer outImage = ImageType::New();
  outImage->CopyInformation( curImage1 );
  typename ImageType::PointType origin = curImage1->GetOrigin();
  typename ImageType::SpacingType spacing = curImage1->GetSpacing();

  typename ImageType::PointType pointX;
  typename ImageType::IndexType indexX;

  typename ImageType::IndexType minX1;
  minX1 = curImage1->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType size1 = curImage1->
                                         GetLargestPossibleRegion().
                                         GetSize();
  typename ImageType::IndexType maxX1;
  maxX1 = minX1 + size1;
  
  typename ImageType::IndexType minX2;
  typename ImageType::IndexType minX2Org;
  minX2Org = curImage2->GetLargestPossibleRegion().GetIndex();
  curImage2->TransformIndexToPhysicalPoint( minX2Org, pointX );
  curImage1->TransformPhysicalPointToIndex( pointX, minX2 );
  typename ImageType::SizeType size2 = curImage2->
                                         GetLargestPossibleRegion().
                                         GetSize();
  typename ImageType::IndexType maxX2;
  typename ImageType::IndexType maxX2Org;
  maxX2Org = minX2Org + size2;
  curImage2->TransformIndexToPhysicalPoint( maxX2Org, pointX );
  curImage1->TransformPhysicalPointToIndex( pointX, maxX2 );

  typename ImageType::IndexType minXOut;
  minXOut = minX1;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    if( minX2[i] < minXOut[i] )
      {
      minXOut[i] = minX2[i];
      }
    if( maxX2[i] < minXOut[i] )
      {
      minXOut[i] = maxX2[i];
      }
    }
  typename ImageType::IndexType maxXOut;
  maxXOut = maxX1;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    if( minX2[i] > maxXOut[i] )
      {
      maxXOut[i] = minX2[i];
      }
    if( maxX2[i] > maxXOut[i] )
      {
      maxXOut[i] = maxX2[i];
      }
    }

  bool useBoundary = false;
  if( boundary.size() == dimensionT)
    {
    useBoundary = true;
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      minX2Org[i] += boundary[i];
      maxX2Org[i] -= boundary[i];
      }
    }

  typename ImageType::SizeType sizeOut;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    sizeOut[i] = maxXOut[i] - minXOut[i] + 1;
    }
  typename ImageType::RegionType regionOut = outImage->
                                               GetLargestPossibleRegion();
  regionOut.SetSize( sizeOut );
  regionOut.SetIndex( minXOut );
  outImage->SetRegions( regionOut );
  outImage->Allocate();

  itk::ImageRegionIteratorWithIndex< ImageType > iter( outImage,
    outImage->GetLargestPossibleRegion() );
  iter.GoToBegin();
  typename ImageType::IndexType indexX2;
  while( !iter.IsAtEnd() )
    {
    indexX = iter.GetIndex();
    outImage->TransformIndexToPhysicalPoint( indexX, pointX );
    double tf1 = background;
    bool inImage1 = false;
    if( curImage1->TransformPhysicalPointToIndex( pointX, indexX ) )
      {
      tf1 = curImage1->GetPixel( indexX );
      inImage1 = true;
      }
    bool inImage2 = false;
    if( curImage2->TransformPhysicalPointToIndex( pointX, indexX2 ) )
      {
      double tf2 = curImage2->GetPixel( indexX2 );
      bool useTF2 = true;
      if( useBoundary )
        {
        for( unsigned int i=0; i<dimensionT; i++ )
          {
          if( indexX2[i]<minX2Org[i] || indexX2[i]>maxX2Org[i] )
            {
            useTF2 = false;
            break;
            }
          }
        }
      if( useTF2 )
        {
        if( average && inImage1 )
          {
          tf1 = ( tf2 + tf1 ) / 2.0;
          }
        else
          {
          tf1 = tf2;
          }
        inImage2 = true;
        }
      }
    iter.Set( tf1 );
    ++iter;
    }

  timeCollector.Start("Save data");
  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( outImage );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
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
  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
