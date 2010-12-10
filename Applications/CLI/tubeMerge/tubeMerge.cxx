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
#include "tubeMergeCLP.h"

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

  timeCollector.Stop("Load data");

  typename ImageType::Pointer curImage1 = reader1->GetOutput();
  typename ImageType::Pointer outImage = ImageType::New();
  outImage->CopyInformation( curImage1 );
  typename ImageType::PointType origin = curImage1->GetOrigin();
  typename ImageType::SpacingType spacing = curImage1->GetSpacing();

  typename ImageType::PointType pointX;
  typename ImageType::IndexType indexX;

  typename ImageType::IndexType minX1Org;
  minX1Org = curImage1->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType size1 = curImage1->
                                         GetLargestPossibleRegion().
                                         GetSize();
  typename ImageType::IndexType maxX1Org;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    maxX1Org[i] = minX1Org[i] + size1[i] - 1;
    }

  bool useBoundary = false;
  if( boundary.size() == dimensionT)
    {
    useBoundary = true;
    }

  double progress = 0.1;
  progressReporter.Report( progress );
  typename ImageType::IndexType minXOut;
  typename ImageType::IndexType maxXOut;
  typename ImageType::SizeType  sizeOut;
  minXOut = minX1Org;
  maxXOut = maxX1Org;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    sizeOut[i] = maxXOut[i] - minXOut[i] + 1;
    }
  for( unsigned int imageNum=0; imageNum<inputVolume2.size(); imageNum++ )
    {
    timeCollector.Start("Load data");
    typename ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName( inputVolume2[imageNum].c_str() );
    try
      {
      reader2->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Stop("Load data");
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop("Load data");

    progress += 1.0/(double)inputVolume2.size() * 0.4;
    progressReporter.Report( progress );

    typename ImageType::Pointer curImage2 = reader2->GetOutput();

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
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      maxX2Org[i] = minX2Org[i] + size2[i] - 1;
      }
    curImage2->TransformIndexToPhysicalPoint( maxX2Org, pointX );
    curImage1->TransformPhysicalPointToIndex( pointX, maxX2 );

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

    for( unsigned int i=0; i<dimensionT; i++ )
      {
      sizeOut[i] = maxXOut[i] - minXOut[i] + 1;
      }
    }

  typename ImageType::RegionType regionOut = outImage->
                                               GetLargestPossibleRegion();
  regionOut.SetSize( sizeOut );
  regionOut.SetIndex( minXOut );
  outImage->SetRegions( regionOut );
  outImage->Allocate();
  outImage->FillBuffer( background );

  if( useBoundary )
    {
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      minX1Org[i] += boundary[i];
      maxX1Org[i] -= boundary[i];
      }
    }

  itk::ImageRegionIteratorWithIndex< ImageType > iter( curImage1,
    curImage1->GetLargestPossibleRegion() );
  while( !iter.IsAtEnd() )
    {
    indexX = iter.GetIndex();
    bool useTf1 = true;
    if( !firstImageSpecial && useBoundary )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        if( indexX[i]<minX1Org[i] || indexX[i]>maxX1Org[i] )
          {
          useTf1 = false;
          break;
          }
        }
      }
    if( useTf1 )
      {
      outImage->SetPixel( indexX, iter.Get() );
      }
    ++iter;
    }
  progress += 0.1;
  progressReporter.Report( progress );

  for( unsigned int imageNum=0; imageNum<inputVolume2.size(); imageNum++ )
    {
    timeCollector.Start("Load data");
    typename ReaderType::Pointer reader2 = ReaderType::New();
    reader2->SetFileName( inputVolume2[imageNum].c_str() );
    try
      {
      reader2->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Stop("Load data");
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    typename ImageType::Pointer curImage2 = reader2->GetOutput();
    timeCollector.Stop("Load data");

    typename ImageType::IndexType minX2Org;
    minX2Org = curImage2->GetLargestPossibleRegion().GetIndex();
    typename ImageType::SizeType size2 = curImage2->
                                           GetLargestPossibleRegion().
                                           GetSize();
    typename ImageType::IndexType maxX2Org;
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      maxX2Org[i] = minX2Org[i] + size2[i] - 1;
      }
    if( useBoundary )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        minX2Org[i] += boundary[i];
        maxX2Org[i] -= boundary[i];
        }
      }

    progress += 1.0/(double)inputVolume2.size() * 0.2;
    progressReporter.Report( progress );

    itk::ImageRegionIteratorWithIndex< ImageType > iter2( curImage2,
      curImage2->GetLargestPossibleRegion() );
    iter2.GoToBegin();
    typename ImageType::IndexType indexX2;
    while( !iter2.IsAtEnd() )
      {
      indexX2 = iter2.GetIndex();
      bool useTf2 = true;
      if( useBoundary )
        {
        for( unsigned int i=0; i<dimensionT; i++ )
          {
          if( indexX2[i]<minX2Org[i] || indexX2[i]>maxX2Org[i] )
            {
            useTf2 = false;
            break;
            }
          }
        }
      if( useTf2 )
        {
        double tf2 = iter2.Get();
        curImage2->TransformIndexToPhysicalPoint( indexX2, pointX );
        if( outImage->TransformPhysicalPointToIndex( pointX, indexX ) )
          {
          if( average )
            {
            double tf1 = outImage->GetPixel( indexX );
            if( tf1 != background )
              {
              tf1 = ( tf2 + tf1 ) / 2.0;
              }
            else
              {
              tf1 = tf2;
              }
            outImage->SetPixel( indexX, tf1 );
            }
          else
            {
            outImage->SetPixel( indexX, tf2 );
            }
          }
        }
      ++iter2;
      }
    progress += 1.0/(double)inputVolume2.size() * 0.19;
    progressReporter.Report( progress );
    }

  timeCollector.Start("Save data");
  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( outImage );
  writer->SetUseCompression( true );
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
