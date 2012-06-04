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
  outImageMap->FillBuffer( 0 );
  typename ImageType::IndexType center = curImage1->
    GetLargestPossibleRegion().GetIndex();
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    center[i] += size1[i] / 2;
    }
  outImageMap->SetPixel( center, 1 );
  typedef typename itk::DanielssonDistanceMapImageFilter< ImageType,
          ImageType>   MapFilterType;
  typename MapFilterType::Pointer mapFilter = MapFilterType::New();
  mapFilter->SetInput( outImageMap );
  mapFilter->SetInputIsBinary( true );
  mapFilter->SetUseImageSpacing( true );
  mapFilter->Update();
  typename ImageType::Pointer outImageDistMap =
    mapFilter->GetDistanceMap();

  itk::ImageRegionIteratorWithIndex< ImageType > iter( curImage1,
    curImage1->GetLargestPossibleRegion() );
  while( !iter.IsAtEnd() )
    {
    indexX = iter.GetIndex();
    double tf = iter.Get();
    if( !mask || tf != 0 )
      {
      outImage->SetPixel( indexX, tf );
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

    typename ImageType::SizeType size2 = curImage2->
                                           GetLargestPossibleRegion().
                                           GetSize();

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
    timeCollector.Stop("Resample Image");

    typename ImageType::Pointer curImage2Map = ImageType::New();
    curImage2Map->CopyInformation( outImage );
    curImage2Map->SetRegions( outImage->GetLargestPossibleRegion() );
    curImage2Map->Allocate();
    curImage2Map->FillBuffer( 0 );
    typename ImageType::IndexType center = curImage2->
      GetLargestPossibleRegion().GetIndex();
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      center[i] += size2[i] / 2;
      }
    typename ImageType::PointType pnt;
    curImage2->TransformIndexToPhysicalPoint( center, pnt );
    typename ImageType::PointType pnt2;
    pnt2 = regOp->GetCurrentMatrixTransform()->TransformPoint( pnt );
    curImage2Map->TransformPhysicalPointToIndex( pnt2, center );
    curImage2Map->SetPixel( center, 1 );
    timeCollector.Start("Distance Map 2");
    typename MapFilterType::Pointer mapFilter2 = MapFilterType::New();
    mapFilter2->SetInput( curImage2Map );
    mapFilter2->SetInputIsBinary( true );
    mapFilter2->SetUseImageSpacing( true );
    mapFilter2->Update();
    typename ImageType::Pointer curImage2DistMap =
      mapFilter2->GetDistanceMap();
    timeCollector.Stop("Distance Map 2");

    progress += 1.0/(double)inputVolume2.size() * 0.2;
    progressReporter.Report( progress );

    iter2.GoToBegin();
    itk::ImageRegionConstIteratorWithIndex< ImageType > iter2Dist(
      curImage2DistMap, curImage2DistMap->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterOut(
      outImage, outImage->GetLargestPossibleRegion() );
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterOutDist(
      outImageDistMap, outImageDistMap->GetLargestPossibleRegion() );
    bool first = true;
    double dist2Max = 0;
    double distOutMax = 0;
    double dist2Min = 0;
    double distOutMin = 0;
    while( !iter2Dist.IsAtEnd() )
      {
      double tf2 = iter2.Get();
      double tf2D = iter2Dist.Get();
      double tfOut = iterOut.Get();
      double tfOutD = iterOutDist.Get();
      if( tf2 != background && tfOut != background )
        {
        if( first )
          {
          dist2Min = tf2D;
          distOutMin = tfOutD;
          first = false;
          }
        if( tf2D > dist2Max )
          {
          dist2Max = tf2D;
          }
        if( tf2D < dist2Min )
          {
          dist2Min = tf2D;
          }
        if( tfOutD > distOutMax )
          {
          distOutMax = tfOutD;
          }
        if( tfOutD < distOutMin )
          {
          distOutMin = tfOutD;
          }
        }
      ++iter2;
      ++iter2Dist;
      ++iterOut;
      ++iterOutDist;
      }

    std::cout << "2: " << dist2Min << " - " << dist2Max << std::endl;
    std::cout << "Out: " << distOutMin << " - " << distOutMax << std::endl;

    timeCollector.Start("Blend");
    iter2.GoToBegin();
    iter2Dist.GoToBegin();
    iterOut.GoToBegin();
    iterOutDist.GoToBegin();
    while( !iter2.IsAtEnd() )
      {
      double tf2 = iter2.Get();
      double tf2D = iter2Dist.Get();
      double tfOut = iterOut.Get();
      double tfOutD = iterOutDist.Get();
      if( average || weighted )
        {
        if( tf2 != background )
          {
          if( tfOut != background )
            {
            if( weighted )
              {
              double d2 = 1 - (tf2D-dist2Min)/(dist2Max-dist2Min);
              d2 = d2 * 2.5 - 1;
              double dOut = 1 - (tfOutD-distOutMin)/(distOutMax-distOutMin);
              dOut = dOut * 2.5 - 1;
              if( d2 < 0 )
                {
                d2 = 0;
                dOut = 1;
                }
              if( dOut < 0 )
                {
                d2 = 1;
                dOut = 0;
                }
              if( d2 >= 1 )
                {
                d2 = 1;
                dOut = 0;
                }
              if( dOut >= 1 )
                {
                d2 = 0;
                dOut = 1;
                }
              double ratio2 = vcl_sqrt( d2 * ( 1 - dOut ) );
              double ratioOut = vcl_sqrt( dOut * ( 1 - d2 ) );
              double ratio = ( ratio2 + ( 1-ratioOut ) ) / 2;
              tfOut = ( ratio * tf2 + (1-ratio) * tfOut );
              }
            else
              {
              tfOut = ( tfOut + tf2 ) / 2;
              }
            }
          else
            {
            tfOut = tf2;
            }
          }
        }
      else if( tf2 != background )
        {
        if( tf2 > tfOut )
          {
          tfOut = tf2;
          }
        }
      iterOut.Set( tfOut );
      ++iter2;
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
