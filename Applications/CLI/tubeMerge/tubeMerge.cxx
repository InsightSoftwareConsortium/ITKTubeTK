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
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDanielssonDistanceMapImageFilter.h"

// Include Slicer3's registration method
#include "itkImageToImageRegistrationHelper.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "tubeMergeCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#define PARSE_ARGS_FLOAT_ONLY
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer's GUI
  tube::CLIProgressReporter    progressReporter( "Merge",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef pixelT                                        PixelType;
  typedef itk::Image< PixelType,  dimensionT >          ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;
  typedef itk::ImageFileWriter< ImageType  >            WriterType;

  typename WriterType::Pointer writer = WriterType::New();

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

  typename ImageType::PointType pointX;
  typename ImageType::IndexType indexX;

  typename ImageType::IndexType minX1Org;
  minX1Org = curImage1->GetLargestPossibleRegion().GetIndex();
  if( boundary.size() == dimensionT )
    {
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      minX1Org[i] -= boundary[i];
      }
    }
  typename ImageType::SizeType size1 = curImage1->
                                         GetLargestPossibleRegion().
                                         GetSize();
  typename ImageType::IndexType maxX1Org;
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    maxX1Org[i] = minX1Org[i] + size1[i] - 1;
    }
  if( boundary.size() == dimensionT )
    {
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      maxX1Org[i] += 2*boundary[i];
      }
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

    timeCollector.Start("Determine ROI");
    progress += 1.0/(double)inputVolume2.size() * 0.4;
    progressReporter.Report( progress );

    typename ImageType::Pointer curImage2 = reader2->GetOutput();

    typename ImageType::IndexType minX2;
    typename ImageType::IndexType minX2Org;
    minX2Org = curImage2->GetLargestPossibleRegion().GetIndex();
    if( boundary.size() == dimensionT )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        minX2Org[i] -= boundary[i];
        }
      }
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
    if( boundary.size() == dimensionT )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        maxX2Org[i] += 2*boundary[i];
        }
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
    timeCollector.Stop("Determine ROI");
    }

  timeCollector.Start("Allocate output image");
  typename ImageType::RegionType regionOut;
  regionOut.SetSize( sizeOut );
  regionOut.SetIndex( minXOut );

  typename ImageType::Pointer outImage = ImageType::New();
  outImage->CopyInformation( curImage1 );
  outImage->SetRegions( regionOut );
  outImage->Allocate();
  outImage->FillBuffer( background );

  typename ImageType::Pointer outImageMap = ImageType::New();
  outImageMap->CopyInformation( curImage1 );
  outImageMap->SetRegions( regionOut );
  outImageMap->Allocate();
  outImageMap->FillBuffer( 1 );

  itk::ImageRegionIteratorWithIndex< ImageType > iter( curImage1,
    curImage1->GetLargestPossibleRegion() );
  while( !iter.IsAtEnd() )
    {
    indexX = iter.GetIndex();
    double tf = iter.Get();
    if( !mask || tf != 0 )
      {
      outImage->SetPixel( indexX, tf );
      outImageMap->SetPixel( indexX, 0 );
      }
    ++iter;
    }
  progress += 0.1;
  progressReporter.Report( progress );
  timeCollector.Stop("Allocate output image");

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

    timeCollector.Start("Register images");
    typedef typename itk::ImageToImageRegistrationHelper< ImageType >
                                                              RegFilterType;
    typename RegFilterType::Pointer regOp = RegFilterType::New();
    regOp->SetFixedImage( curImage1 );
    regOp->SetMovingImage( curImage2 );
    regOp->SetSampleFromOverlap( true );
    regOp->SetEnableLoadedRegistration( false );
    regOp->SetEnableInitialRegistration( false );
    regOp->SetEnableRigidRegistration( true );
    regOp->SetRigidSamplingRatio( samplingRatio );
    regOp->SetRigidMaxIterations( iterations );
    regOp->SetEnableAffineRegistration( false );
    regOp->SetEnableBSplineRegistration( false );
    regOp->SetExpectedOffsetPixelMagnitude( expectedOffset );
    regOp->SetExpectedRotationMagnitude( expectedRotation );

    regOp->Initialize();
    regOp->Update();

    if( saveTransform.size() > 1 )
      {
      regOp->SaveTransform( saveTransform );
      }
    timeCollector.Stop("Register images");

    timeCollector.Start("Resample Image");
    regOp->SetFixedImage( outImage );
    typename ImageType::ConstPointer tmpImage = regOp->ResampleImage(
      RegFilterType::OptimizedRegistrationMethodType::
      NEAREST_NEIGHBOR_INTERPOLATION,
      NULL, NULL, NULL, background );
    typename ImageType::Pointer curImage2Reg = ImageType::New();
    curImage2Reg->CopyInformation( tmpImage );
    curImage2Reg->SetRegions( tmpImage->GetLargestPossibleRegion() );
    curImage2Reg->Allocate();
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterTmp(
      tmpImage, tmpImage->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iter2(
      curImage2Reg, curImage2Reg->GetLargestPossibleRegion() );
    while( !iter2.IsAtEnd() )
      {
      iter2.Set( iterTmp.Get() );
      ++iter2;
      ++iterTmp;
      }

    writer->SetFileName( "image2Reg.mha" );
    writer->SetInput( curImage2Reg );
    writer->SetUseCompression( true );
    writer->Update();

    typedef typename itk::DanielssonDistanceMapImageFilter< ImageType,
            ImageType>   MapFilterType;
    typename MapFilterType::Pointer mapFilter = MapFilterType::New();
    mapFilter->SetInput( outImageMap );
    mapFilter->SetInputIsBinary( true );
    mapFilter->SetUseImageSpacing( true );
    mapFilter->Update();
    typename ImageType::Pointer outImageDistMap =
      mapFilter->GetDistanceMap();
    writer->SetFileName( "outImageDistMap.mha" );
    writer->SetInput( outImageDistMap );
    writer->SetUseCompression( true );
    writer->Update();
    timeCollector.Stop("Resample Image");

    timeCollector.Start("Resample Map");
    typename ImageType::Pointer curImage2Tmp = ImageType::New();
    curImage2Tmp->CopyInformation( curImage2 );
    curImage2Tmp->SetRegions( curImage2->GetLargestPossibleRegion() );
    curImage2Tmp->Allocate();
    curImage2Tmp->FillBuffer( 1 );
    if( mask )
      {
      itk::ImageRegionIteratorWithIndex< ImageType > iter2(
        curImage2, curImage2->GetLargestPossibleRegion() );
      itk::ImageRegionIteratorWithIndex< ImageType > iter2Tmp(
        curImage2Tmp, curImage2Tmp->GetLargestPossibleRegion() );
      while( !iter2.IsAtEnd() )
        {
        double tf = iter2.Get();
        if( tf == 0 )
          {
          iter2Tmp.Set( 0 );
          }
        ++iter2;
        ++iter2Tmp;
        }
      }

    regOp->SetMovingImage( curImage2Tmp );
    typename ImageType::ConstPointer imageTmp2 = regOp->ResampleImage(
      RegFilterType::OptimizedRegistrationMethodType::
      NEAREST_NEIGHBOR_INTERPOLATION,
      NULL, NULL, NULL, background );
    typename ImageType::Pointer curImage2Map = ImageType::New();
    curImage2Map->CopyInformation( tmpImage );
    curImage2Map->SetRegions( tmpImage->GetLargestPossibleRegion() );
    curImage2Map->Allocate();
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterTmp2(
      imageTmp2, imageTmp2->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iter2Map(
      curImage2Map, curImage2Map->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterOutMap(
      outImageMap, outImageMap->GetLargestPossibleRegion() );
    while( !iter2Map.IsAtEnd() )
      {
      double tf2 = iterTmp2.Get();
      double tfOut = iterOutMap.Get();
      if( tf2 != 0 && tfOut == 0 )
        {
        iter2Map.Set( 0 );
        }
      else
        {
        if( tf2 != 0 )
          {
          iter2Map.Set( -2 );
          }
        else
          {
          iter2Map.Set( -3 );
          }
        }
      if( tf2 != 0 )
        {
        iterOutMap.Set( 0 );
        }
      ++iterTmp2;
      ++iterOutMap;
      ++iter2Map;
      }
    timeCollector.Stop("Resample Map");

    progress += 1.0/(double)inputVolume2.size() * 0.2;
    progressReporter.Report( progress );

    timeCollector.Start("Distance Map");
    typename MapFilterType::Pointer mapFilter2 = MapFilterType::New();
    mapFilter2->SetInput( curImage2Map );
    mapFilter2->SetInputIsBinary( true );
    mapFilter2->SetUseImageSpacing( true );
    mapFilter2->Update();
    typename ImageType::Pointer curImage2DistMap =
      mapFilter2->GetDistanceMap();
    writer->SetFileName( "curImage2DistMap.mha" );
    writer->SetInput( curImage2DistMap );
    writer->SetUseCompression( true );
    writer->Update();
    timeCollector.Stop("Distance Map");

    itk::ImageRegionConstIteratorWithIndex< ImageType > iter2Dist(
      curImage2DistMap, curImage2DistMap->GetLargestPossibleRegion() );
    double distMax = 0;
    while( !iter2Dist.IsAtEnd() )
      {
      double tf2D = iter2Dist.Get();
      if( tf2D > distMax )
        {
        distMax = tf2D;
        }
      ++iter2Dist;
      }

    timeCollector.Start("Blend");
    iter2.GoToBegin();
    iter2Map.GoToBegin();
    iter2Dist.GoToBegin();
    itk::ImageRegionIteratorWithIndex< ImageType > iterOut(
      outImage, outImage->GetLargestPossibleRegion() );
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterOutDist(
      outImageDistMap, outImageDistMap->GetLargestPossibleRegion() );
    while( !iter2.IsAtEnd() )
      {
      double tf2 = iter2.Get();
      double tf2M = iter2Map.Get();
      double tf2D = iter2Dist.Get();
      double tf1 = iterOut.Get();
      double tf1D = iterOutDist.Get();
      if( average || weighted )
        {
        if( tf2M != -3 )
          {
          if( tf2M != -2 )
            {
            if( weighted )
              {
              double ratio = tf2D/distMax;
              ratio *= 0.5;
              if( tf2D < tf1D )
                {
                tf1 = ( ratio * tf2 + (1-ratio) * tf1 );
                }
              else
                {
                tf1 = ( (1-ratio) * tf2 + ratio * tf1 );
                }
              }
            else
              {
              tf1 = ( tf1 + tf2 ) / 2;
              }
            }
          else
            {
            tf1 = tf2;
            }
          }
        }
      else if( tf2M != -3 )
        {
        if( tf2 > tf1 )
          {
          tf1 = tf2;
          }
        }
      iterOut.Set( tf1 );
      ++iter2;
      ++iter2Map;
      ++iter2Dist;
      ++iterOut;
      ++iterOutDist;
      }
    timeCollector.Stop("Blend");
    progress += 1.0/(double)inputVolume2.size() * 0.19;
    progressReporter.Report( progress );
    }

  timeCollector.Start("Save data");
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
