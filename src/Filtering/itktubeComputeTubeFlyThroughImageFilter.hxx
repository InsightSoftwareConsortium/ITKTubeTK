/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeComputeTubeFlyThroughImageFilter_hxx
#define __itktubeComputeTubeFlyThroughImageFilter_hxx


#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkResampleImageFilter.h>

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TPixel, unsigned int Dimension >
ComputeTubeFlyThroughImageFilter< TPixel, Dimension >
::ComputeTubeFlyThroughImageFilter( void )
{
  m_OutputMask = NULL;
}

template< class TPixel, unsigned int Dimension >
void
ComputeTubeFlyThroughImageFilter< TPixel, Dimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  SuperClass::PrintSelf( os, indent );
  os << "TubeId: " << m_TubeId << std::endl;
}

template< class TPixel, unsigned int Dimension >
void
ComputeTubeFlyThroughImageFilter< TPixel, Dimension >
::GenerateData( void )
{
  itkDebugMacro( << "ComputeTubeFlyThroughImageFilter::Update() called." );

  // get input tubes
  typename TubeGroupType::ConstPointer inputTubeGroup = this->GetInput();

  // Find the user specified tybe
  itkDebugMacro( << "Finding user specified tube" );

  typedef typename TubeGroupType::ChildrenListType      TubeListType;

  TubeType * inputTube = NULL;
  bool blnTubeFound = false;
  char tubeName[] = "Tube";

  TubeListType * tubeList = inputTubeGroup->GetChildren(
    inputTubeGroup->GetMaximumDepth(), tubeName );

  typename TubeListType::const_iterator itTubes = tubeList->begin();

  while( itTubes != tubeList->end() )
    {
    // std::cout m_TubeId << " == " << ( *itTubes )->GetId() << std::endl;

    if( static_cast< unsigned long >( ( *itTubes )->GetId() ) == m_TubeId )
      {
      inputTube = dynamic_cast< TubeType * >( itTubes->GetPointer() );

      if( inputTube )
        {
        blnTubeFound = true;
        }

      break;
      }

    itTubes++;
    }

  delete tubeList;

  if( !blnTubeFound )
    {
    itkExceptionMacro( "Unable to find the specified tube" );
    }

  inputTube->Update();
  inputTube->RemoveDuplicatePointsInObjectSpace();
  inputTube->ComputeTangentsAndNormals();

  // Get list of tube points
  typedef typename TubeType::TubePointListType    TubePointListType;

  TubePointListType tubePointList = inputTube->GetPoints();

  if( tubePointList.size() <= 0 )
    {
    itkExceptionMacro( "The specified tube does not contain any points" );
    }

  itkDebugMacro( << "Num Tube Points = " << tubePointList.size() );

  // Determine maximum radius among all tube points
  typename TubePointListType::const_iterator itPts = tubePointList.begin();

  double maxTubeRadius = itPts->GetRadiusInWorldSpace();

  while( itPts != tubePointList.end() )
    {
    if( itPts->GetRadiusInWorldSpace() > maxTubeRadius )
      {
      maxTubeRadius = itPts->GetRadiusInWorldSpace();
      }
    itPts++;
    }

  itkDebugMacro( << "Max Radius = " << maxTubeRadius );

  // Determine the mean distance between consecutive tube points
  double meanTubePointDist = 0;

  for( unsigned int pid = 1; pid < tubePointList.size(); pid++ )
    {
    // compute distance between current and previous tube point
    double curDist = 0;
    typename TubeType::PointType p1 = tubePointList[pid-1].GetPositionInWorldSpace();
    typename TubeType::PointType p2 = tubePointList[pid].GetPositionInWorldSpace();

    for( unsigned int i = 0; i < Dimension; i++ )
      {
      curDist += ( p2[i] - p1[i] ) * ( p2[i] - p1[i] );
      }

    // update mean
    meanTubePointDist += curDist;
    }

  meanTubePointDist /= tubePointList.size();

  itkDebugMacro( << "Mean Tube Point Distance = " << meanTubePointDist );

  // Determine minimum spacing of the input image
  typename InputImageType::SpacingType inputSpacing =
    m_InputImage->GetSpacing();
  double minInputSpacing = inputSpacing[0];

  for( unsigned int i = 1; i < Dimension; i++ )
    {
    if( inputSpacing[i] < minInputSpacing )
      {
      minInputSpacing = inputSpacing[i];
      }
    }

  itkDebugMacro( << "Min Input Spacing = " << minInputSpacing );

  // allocate space for output fly through image
  typename OutputImageType::Pointer outputImage = this->GetOutput();

  // set spacing
  // For last dimension its set to mean consecutive point distance
  // For other dimensions its set to the minimum input spacing
  typename OutputImageType::SpacingType outputSpacing;

  for( unsigned int i = 0; i < Dimension-1; i++ )
    {
    outputSpacing[i] = minInputSpacing;
    }

  outputSpacing[Dimension-1] = meanTubePointDist;

  outputImage->SetSpacing( outputSpacing );

  // set start index1
  typename OutputImageType::IndexType startIndex;
  startIndex.Fill( 0 );

  // set size
  typename OutputImageType::SizeType size;
  for( unsigned int i = 0; i < Dimension-1; i++ )
    {
    size[i] = 2 * ( typename OutputImageType::SizeValueType )
      ( 0.5 + ( maxTubeRadius / outputSpacing[i] ) ) + 1;
    }
  size[Dimension-1] = tubePointList.size();

  // set regions
  typename OutputImageType::RegionType region;
  region.SetIndex( startIndex );
  region.SetSize( size );
  outputImage->SetRegions( region );

  // Allocate and initialize
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  // allocate space for output fly through mask
  m_OutputMask = OutputMaskType::New();
  m_OutputMask->SetRegions( region );
  m_OutputMask->CopyInformation( outputImage );
  m_OutputMask->Allocate();
  m_OutputMask->FillBuffer( 0 );

  // For each tube point, extract normal plane image
  // and fill into corresponding slice in the output image
  itkDebugMacro( "Generating fly through image" );

  typedef typename TubeType::TubePointType             TubePointType;
  typedef typename TubePointType::CovariantVectorType  TubeNormalType;
  typedef ImageRegionIteratorWithIndex<
    OutputImageType >                                  OutputImageIteratorType;
  typedef ImageRegionIterator< OutputMaskType >        OutputMaskIteratorType;

  typedef LinearInterpolateImageFunction< InputImageType, double >
    InterpolatorType;

  typedef MinimumMaximumImageFilter< InputImageType >
    MinMaxImageFilterType;

  typename MinMaxImageFilterType::Pointer minmaxFilter =
  MinMaxImageFilterType::New();
  minmaxFilter->SetInput( m_InputImage );
  minmaxFilter->Update();
  TPixel outsideVal = minmaxFilter->GetMinimum();

  typename InterpolatorType::Pointer pInterpolator = InterpolatorType::New();
  pInterpolator->SetInputImage( m_InputImage );

  unsigned int ptInd = 0;
  unsigned long tubePixelCount = 0;

  for( itPts = tubePointList.begin();
    itPts != tubePointList.end(); itPts++, ptInd++ )
    {
    // Get position, radius and frenet-serret basis of current tube point
    // in the world coordinate system
    typename TubeType::PointType curTubePosition =
      itPts->GetPositionInWorldSpace();

    TubeNormalType curTubeNormal1 = itPts->GetNormal1InWorldSpace();
    curTubeNormal1.Normalize();

    TubeNormalType curTubeNormal2 = itPts->GetNormal2InWorldSpace();
    curTubeNormal2.Normalize();

    double curTubeRadius = ( *itPts ).GetRadiusInWorldSpace();

    //std::cout << curTubeNormal1.GetNorm() << std::endl;
    //std::cout << curTubeNormal2.GetNorm() << std::endl;

    // Define slice region in the output image
    typename OutputImageType::RegionType sliceRegion;

    typename OutputImageType::IndexType sliceStartIndex;
    sliceStartIndex.Fill( 0 );
    sliceStartIndex[Dimension-1] = ptInd;

    typename OutputImageType::SizeType sliceSize;
    sliceSize = outputImage->GetLargestPossibleRegion().GetSize();
    sliceSize[Dimension-1] = 1;

    sliceRegion.SetIndex( sliceStartIndex );
    sliceRegion.SetSize( sliceSize );

    // Iterate through corresponding slice of output image and fill each pixel
    OutputImageIteratorType itOutSlice( outputImage, sliceRegion );
    OutputMaskIteratorType itMask( m_OutputMask, sliceRegion );

    for( itOutSlice.GoToBegin(), itMask.GoToBegin();
      !itOutSlice.IsAtEnd(); ++itOutSlice, ++itMask )
      {
      // get index of the current output pixel
      typename OutputImageType::IndexType curOutIndex = itOutSlice.GetIndex();

      // compute corresponding position in the input image
      //typename InterpolatorType::ContinuousIndexType curInputPoint;
      typename OutputImageType::PointType curInputPoint;

      double distToCenter = 0;

      for( unsigned int i = 0; i < Dimension; i++ )
        {
        curInputPoint[i] = curTubePosition[i];
        }

      if( Dimension == 2 )
        {
        double stepN1 = ( curOutIndex[0] - 0.5 * sliceSize[0] )
          * outputSpacing[0];

        for( unsigned int i = 0; i < Dimension; i++ )
          {
          curInputPoint[i] += stepN1 * curTubeNormal1[i];
          }

          distToCenter = stepN1;
        }
      else if( Dimension == 3 )
        {
        double stepN1 = ( curOutIndex[0] - 0.5 * sliceSize[0] )
          * outputSpacing[0];

        double stepN2 = ( curOutIndex[1] - 0.5 * sliceSize[1] )
          * outputSpacing[1];

        for( unsigned int i = 0; i < Dimension; i++ )
          {
          curInputPoint[i] += stepN1 * curTubeNormal1[i];
          curInputPoint[i] += stepN2 * curTubeNormal2[i];
          }

          distToCenter = std::sqrt( stepN1 * stepN1 + stepN2 * stepN2 );
        }

      // set pixel values in the output images
      if( pInterpolator->IsInsideBuffer( curInputPoint ) )
        {
        // set intensity value by getting it from input image using
        // interpolation
        itOutSlice.Set( pInterpolator->Evaluate( curInputPoint ) );

        // if point is within the tube set tube mask pixel to on
        if( distToCenter <= curTubeRadius )
          {
          itMask.Set( 1.0 );
          }

        tubePixelCount++;
        }
      else
        {
        itOutSlice.Set( outsideVal );
        }
      }
    }

  itkDebugMacro( << "ComputeTubeFlyThroughImageFilter::Update() finished." );
}

} // End namespace tube

} // End namespace itk

#endif
