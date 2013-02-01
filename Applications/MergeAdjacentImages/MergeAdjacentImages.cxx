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
#include "itkSignedDanielssonDistanceMapImageFilter.h"

// Include Slicer3's registration method
#include "itkImageToImageRegistrationHelper.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "MergeAdjacentImagesCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#define PARSE_ARGS_FLOAT_ONLY 1
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
    timeCollector.Start("Load data");
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

  timeCollector.Start("Allocate output image");
  typename ImageType::RegionType regionOut;
  regionOut.SetSize( sizeOut );
  regionOut.SetIndex( minXOut );

  typename ImageType::Pointer outImage = ImageType::New();
  outImage->CopyInformation( curImage1 );
  outImage->SetRegions( regionOut );
  outImage->Allocate();
  outImage->FillBuffer( background );

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

    timeCollector.Start("Load data");
    reader2->SetFileName( inputVolume2.c_str() );
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
    curImage2 = reader2->GetOutput();
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
    curImage2Reg->FillBuffer( background );
    typename ImageType::Pointer outImageMap = ImageType::New();
    outImageMap->CopyInformation( tmpImage );
    outImageMap->SetRegions( tmpImage->GetLargestPossibleRegion() );
    outImageMap->Allocate();
    outImageMap->FillBuffer( 3 );
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterTmp(
      tmpImage, tmpImage->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iter2(
      curImage2Reg, curImage2Reg->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterOut(
      outImage, outImage->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterOutMap(
      outImageMap, outImageMap->GetLargestPossibleRegion() );
    while( !iter2.IsAtEnd() )
      {
      double iVal = iterTmp.Get();
      iter2.Set( iVal );
      bool image2Point = false;
      if( iVal != background && (!mask || iVal != 0) )
        {
        image2Point = true;
        }
      double oVal = iterOut.Get();
      bool imageOutPoint = false;
      if( oVal != background && (!mask || oVal != 0) )
        {
        imageOutPoint = true;
        }
      if( image2Point && imageOutPoint )
        {
        iterOutMap.Set( 0 );
        }
      else if( image2Point )
        {
        iterOutMap.Set( 2 );
        }
      else if( imageOutPoint )
        {
        iterOutMap.Set( 1 );
        }
      ++iterTmp;
      ++iter2;
      ++iterOut;
      ++iterOutMap;
      }
    timeCollector.Stop("Resample Image");

    timeCollector.Start("Out Distance Map");
    typedef typename itk::DanielssonDistanceMapImageFilter< ImageType,
      ImageType>   MapFilterType;
    typename MapFilterType::Pointer mapDistFilter = MapFilterType::New();
    mapDistFilter->SetInput( outImageMap );
    mapDistFilter->SetInputIsBinary( false );
    mapDistFilter->SetUseImageSpacing( false );
    mapDistFilter->Update();
    typename ImageType::Pointer outImageDistMap =
      mapDistFilter->GetDistanceMap();
    typename ImageType::Pointer outImageVoronoiMap =
      mapDistFilter->GetVoronoiMap();
    timeCollector.Stop("Out Distance Map");

    progress += 1.0/(double)inputVolume2.size() * 0.2;
    progressReporter.Report( progress );

    typename ImageType::Pointer vorImageMap = ImageType::New();
    vorImageMap->CopyInformation( outImage );
    vorImageMap->SetRegions( outImage->GetLargestPossibleRegion() );
    vorImageMap->Allocate();
    vorImageMap->FillBuffer( 1 );
    itk::ImageRegionConstIteratorWithIndex< ImageType > iterOutVor(
      outImageVoronoiMap, outImageVoronoiMap->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterVorMap(
      vorImageMap, vorImageMap->GetLargestPossibleRegion() );
    while( !iterOutVor.IsAtEnd() )
      {
      double tf = iterOutVor.Get();
      if( tf != 1 )
        {
        iterVorMap.Set( 0 );
        }
      ++iterOutVor;
      ++iterVorMap;
      }

    timeCollector.Start("Voronoi Distance Map");
    typedef typename itk::SignedDanielssonDistanceMapImageFilter< ImageType,
      ImageType>   SignedMapFilterType;
    typename SignedMapFilterType::Pointer mapVorFilter =
      SignedMapFilterType::New();
    mapVorFilter->SetInput( vorImageMap );
    mapVorFilter->SetUseImageSpacing( false );
    mapVorFilter->Update();
    typename ImageType::Pointer vorImageDistMap =
      mapVorFilter->GetDistanceMap();
    timeCollector.Stop("Voronoi Distance Map");

    iter2.GoToBegin();
    iterOut.GoToBegin();
    itk::ImageRegionIteratorWithIndex< ImageType > iterOutDistMap(
      outImageDistMap, outImageDistMap->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex< ImageType > iterVorDistMap(
      vorImageDistMap, vorImageDistMap->GetLargestPossibleRegion() );
    while( !iter2.IsAtEnd() )
      {
      double iVal = iter2.Get();
      bool image2Point = false;
      if( iVal != background && (!mask || iVal != 0) )
        {
        image2Point = true;
        }
      double oVal = iterOut.Get();
      bool imageOutPoint = false;
      if( oVal != background && (!mask || oVal != 0) )
        {
        imageOutPoint = true;
        }

      if( imageOutPoint )
        {
        if( image2Point )
          {
          double vDist = iterVorDistMap.Get();
          double oDist = iterOutDistMap.Get();

          if( vDist < 0 )
            {
            vDist = -vDist;
            double ratio = 0.5*vDist/(oDist+vDist) + 0.5;
            oVal = ratio * oVal + (1-ratio)*iVal;
            }
          else if( vDist > 0 )
            {
            double ratio = 0.5*vDist/(oDist+vDist) + 0.5;
            oVal = ratio * iVal + (1-ratio)*oVal;
            }
          else
            {
            oVal = 0.5 * iVal + 0.5 * oVal;
            }
          }
        }
      else if( image2Point )
        {
        oVal = iVal;
        }

      iterOut.Set( oVal );

      ++iter2;
      ++iterOut;
      ++iterOutDistMap;
      ++iterVorDistMap;
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
