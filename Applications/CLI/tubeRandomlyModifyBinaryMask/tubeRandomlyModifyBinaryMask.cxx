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
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "tubeRandomlyModifyBinaryMaskCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "RandomlyModifyBinaryMask",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef pixelT                                        PixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >  ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;

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

  std::vector< int > regionSize;
  regionSize.resize( dimensionT );
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    regionSize[i] = maxEffectSize * 2;
    }

  int minChange = 1;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    minChange *= minEffectSize;
    }

  typename ImageType::Pointer curImage = reader->GetOutput();

  typename ImageType::IndexType index = curImage->
                                           GetLargestPossibleRegion().
                                           GetIndex();
  typename ImageType::SizeType size = curImage->
                                           GetLargestPossibleRegion().
                                           GetSize();
  typename ImageType::IndexType bufferSize;
  typename ImageType::IndexType regionCenterMin;
  typename ImageType::IndexType regionCenterMax;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    bufferSize[i] = maxEffectSize/2;
    regionCenterMin[i] = index[i] + regionSize[i]/2 + bufferSize[i] + 1;
    regionCenterMax[i] = index[i] + size[i] - 1 - regionSize[i]/2
                                  - bufferSize[i] - 1;
    }

  typename ImageType::SizeType tmpSize;
  typename ImageType::SizeType roiSize;
  typename ImageType::IndexType roiIndex;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    tmpSize[i] = regionSize[i]+bufferSize[i]*2+1;
    roiSize[i] = regionSize[i];
    roiIndex[i] = index[i] + bufferSize[i];
    }

  typename ImageType::RegionType tmpRegion;
  tmpRegion = curImage->GetLargestPossibleRegion();

  typename ImageType::Pointer changeImage = ImageType::New();
  changeImage->CopyInformation( curImage );
  changeImage->SetRegions( tmpRegion );
  changeImage->Allocate();
  changeImage->FillBuffer( 0 );

  tmpRegion.SetSize( tmpSize );
  typename ImageType::Pointer tmpImage = ImageType::New();
  tmpImage->CopyInformation( curImage );
  tmpImage->SetRegions( tmpRegion );
  tmpImage->Allocate();

  typename ImageType::IndexType tmpCenter;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    tmpCenter[i] = index[i] + tmpSize[i]/2;
    }

  typename ImageType::RegionType tmpROIRegion;
  tmpROIRegion = curImage->GetLargestPossibleRegion();
  tmpROIRegion.SetSize( roiSize );
  tmpROIRegion.SetIndex( roiIndex );

  typename ImageType::RegionType addSubRegion;
  addSubRegion = tmpImage->GetLargestPossibleRegion();
  typename ImageType::RegionType::SizeType addSubRegionSize;
  typename ImageType::IndexType addSubRegionIndex;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    addSubRegionSize[i] = maxEffectSize + 1;
    addSubRegionIndex[i] = roiIndex[i] + (roiSize[i]/2-maxEffectSize/2);
    }
  addSubRegion.SetSize( addSubRegionSize );
  addSubRegion.SetIndex( addSubRegionIndex );

  typename ImageType::IndexType roiCenter;
  typename ImageType::IndexType tmpMin;
  typename ImageType::IndexType tmpROIMin;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer
    randGen =
    itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  if( seed != -1 )
    {
    randGen->Initialize( seed );
    }
  else
    {
    randGen->Initialize();
    }

  for( int regionNum=0; regionNum<numberOfRegions; ++regionNum )
    {
    bool foundForegroundPixel = false;
    while( !foundForegroundPixel )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        roiCenter[i] = (int)( randGen->GetUniformVariate(
                                   regionCenterMin[i],
                                   regionCenterMax[i] ) );
        tmpMin[i] = roiCenter[i] - (tmpSize[i]-1)/2;
        tmpROIMin[i] = tmpMin[i] + bufferSize[i];
        }
      if( curImage->GetPixel( roiCenter ) == foreground )
        {
        foundForegroundPixel = true;
        }
      }

    typename ImageType::RegionType curRegion;
    curRegion = curImage->GetLargestPossibleRegion();
    curRegion.SetSize( tmpSize );
    curRegion.SetIndex( tmpMin );

    int mode = (int)( randGen->GetUniformVariate( 0, 4 ) );
    switch( mode )
      {
      case 0:
      case 1:
        {
        tmpImage->FillBuffer( 0 );
        typename ImageType::IndexType curTestIndex;
        curTestIndex = roiCenter;
        typename ImageType::IndexType tmpTestIndex;
        tmpTestIndex = tmpCenter;
        if( mode == 0 )
          {
          for( unsigned int i=0; i<dimensionT; i++ )
            {
            curTestIndex[i] = roiCenter[i] + maxEffectSize/2;
            if( curImage->GetPixel( curTestIndex ) == background )
              {
              tmpTestIndex[i] = tmpCenter[i] + maxEffectSize/2;
              tmpImage->SetPixel( tmpTestIndex, foreground );
              tmpTestIndex[i] = tmpCenter[i];
              }
            curTestIndex[i] = roiCenter[i] - maxEffectSize/2;
            if( curImage->GetPixel( curTestIndex ) == background )
              {
              tmpTestIndex[i] = tmpCenter[i] - maxEffectSize/2;
              tmpImage->SetPixel( tmpTestIndex, foreground );
              tmpTestIndex[i] = tmpCenter[i];
              }
            curTestIndex[i] = roiCenter[i];
            }
          }
        else
          {
          tmpImage->SetPixel( tmpTestIndex, foreground );
          }
        typedef itk::BinaryBallStructuringElement<PixelType, dimensionT >
          BallType;
        BallType ball;
        if( mode == 0 )
          {
          ball.SetRadius( maxEffectSize );
          }
        else
          {
          ball.SetRadius( (maxEffectSize+minEffectSize)/2 - 1 );
          }
        ball.CreateStructuringElement();
        typedef itk::BinaryDilateImageFilter< ImageType, ImageType,
          BallType > FilterType;
        typename FilterType::Pointer filter = FilterType::New();
        filter->SetDilateValue( foreground );
        filter->SetKernel( ball );
        filter->SetInput( tmpImage );
        filter->Update();
        typename ImageType::Pointer outImage = filter->GetOutput();
        itk::ImageRegionIterator< ImageType >
          curIter( curImage, curRegion );
        itk::ImageRegionIterator< ImageType >
          tmpIter( tmpImage, tmpRegion );
        itk::ImageRegionIterator< ImageType >
          outIter( outImage, tmpRegion );
        curIter.GoToBegin();
        tmpIter.GoToBegin();
        outIter.GoToBegin();
        if( mode == 0 )
          {
          std::cout << "Divot : " << roiCenter << std::endl;
          while( !curIter.IsAtEnd() )
            {
            if( outIter.Get() == foreground )
              {
              tmpIter.Set( background );
              }
            else
              {
              tmpIter.Set( curIter.Get() );
              }
            ++outIter;
            ++tmpIter;
            ++curIter;
            }
          }
        else
          {
          std::cout << "Bump : " << roiCenter << std::endl;
          while( !curIter.IsAtEnd() )
            {
            if( outIter.Get() == foreground )
              {
              tmpIter.Set( foreground );
              }
            else
              {
              tmpIter.Set( curIter.Get() );
              }
            ++outIter;
            ++tmpIter;
            ++curIter;
            }
          }
        break;
        }
      case 2:
        {
        std::cout << "Add : " << roiCenter << std::endl;
        itk::ImageRegionIterator< ImageType >
          curIter( curImage, curRegion );
        itk::ImageRegionIterator< ImageType >
          tmpIter( tmpImage, tmpRegion );
        curIter.GoToBegin();
        tmpIter.GoToBegin();
        while( !curIter.IsAtEnd() )
          {
          tmpIter.Set( curIter.Get() );
          ++tmpIter;
          ++curIter;
          }
        itk::ImageRegionIterator< ImageType > iter( tmpImage,
                                                    addSubRegion );
        iter.GoToBegin();
        while( !iter.IsAtEnd() )
          {
          iter.Set( foreground );
          ++iter;
          }
        break;
        }
      case 3:
        {
        std::cout << "Subtract : " << roiCenter << std::endl;
        itk::ImageRegionIterator< ImageType >
          curIter( curImage, curRegion );
        itk::ImageRegionIterator< ImageType >
          tmpIter( tmpImage, tmpRegion );
        curIter.GoToBegin();
        tmpIter.GoToBegin();
        while( !curIter.IsAtEnd() )
          {
          tmpIter.Set( curIter.Get() );
          ++tmpIter;
          ++curIter;
          }
        itk::ImageRegionIterator< ImageType > iter( tmpImage,
                                                    addSubRegion );
        iter.GoToBegin();
        while( !iter.IsAtEnd() )
          {
          iter.Set( background );
          ++iter;
          }
        break;
        }
      }

    typename ImageType::RegionType curROIRegion;
    curROIRegion = curImage->GetLargestPossibleRegion();
    curROIRegion.SetSize( roiSize );
    curROIRegion.SetIndex( tmpROIMin );

    itk::ImageRegionIterator< ImageType > changeROIIter( changeImage,
                                                      curROIRegion );
    itk::ImageRegionIterator< ImageType > curROIIter( curImage,
                                                      curROIRegion );
    itk::ImageRegionIterator< ImageType > tmpROIIter( tmpImage,
                                                      tmpROIRegion );
    changeROIIter.GoToBegin();
    curROIIter.GoToBegin();
    tmpROIIter.GoToBegin();
    int count = 0;
    while( !tmpROIIter.IsAtEnd() )
      {
      if( curROIIter.Get() != tmpROIIter.Get() )
        {
        ++count;
        }
      if( changeROIIter.Get() != 0 )
        {
        count = 0;
        break;
        }
      ++changeROIIter;
      ++curROIIter;
      ++tmpROIIter;
      }
    if( count > minChange )
      {
      changeROIIter.GoToBegin();
      curROIIter.GoToBegin();
      tmpROIIter.GoToBegin();
      while( !tmpROIIter.IsAtEnd() )
        {
        if( curROIIter.Get() != tmpROIIter.Get() ||
            ( curROIIter.Get() == foreground ) )
          {
          curROIIter.Set( tmpROIIter.Get() );
          if( mode == 0 )
            {
            changeROIIter.Set( 1 );
            }
          if( mode == 3 )
            {
            changeROIIter.Set( 2 );
            }
          if( mode == 1 )
            {
            changeROIIter.Set( 3 );
            }
          if( mode == 2 )
            {
            changeROIIter.Set( 4 );
            }
          }
        ++changeROIIter;
        ++curROIIter;
        ++tmpROIIter;
        }
      }
    else
      {
      std::cout << "Rejected: change insufficient: "
                << count << " < " << minChange << std::endl;
      --regionNum;
      }
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( curImage );
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

  typename ImageWriterType::Pointer changeWriter = ImageWriterType::New();
  std::stringstream out;
  out << outputVolume;
  out << "_change.mha";
  changeWriter->SetFileName( out.str() );
  changeWriter->SetInput( changeImage );
  changeWriter->SetUseCompression( true );
  try
    {
    changeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing change volume: Exception caught: "
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
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
