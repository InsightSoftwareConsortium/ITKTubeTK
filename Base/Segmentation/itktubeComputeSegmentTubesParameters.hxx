/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __itktubeComputeSegmentTubesParameters_hxx
#define __itktubeComputeSegmentTubesParameters_hxx

#include <algorithm>

#include "itktubeMetaTubeExtractor.h"

#include "itktubeComputeSegmentTubesParameters.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_c_vector.h>

int m_SortColumn = 0;
struct ComputeSegmentTubesParametersSortFunctionType {
    bool operator()( const vnl_vector< double > & first,
      const vnl_vector< double > & second ) const
      {
      return( first[m_SortColumn] < second[m_SortColumn] );
      }
    };

namespace itk
{

namespace tube
{

/**
 * Constructor */
template< class TPixel, unsigned int VDimension >
ComputeSegmentTubesParameters< TPixel, VDimension >
::ComputeSegmentTubesParameters( void )
{
  m_InputImage = NULL;
  m_MaskInputImage = NULL;
  m_ScaleInputImage = NULL;
}

/** Destructor */
template< class TPixel, unsigned int VDimension >
ComputeSegmentTubesParameters< TPixel, VDimension >
::~ComputeSegmentTubesParameters( void )
{
}

/** Get Tube data list */
template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetTubeData( void )
{
  return m_TubeData;
}

/** Get Bkg data list */
template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetBkgData( void )
{
  return m_BkgData;
}

/** Get seed data list */
template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetSeedData( void )
{
  return m_SeedData;
}

/** Get Seed Index List */
template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetSeedDataIndexList( void )
{
  return m_SeedDataIndexList;
}

/**
 * Get Tube Index List */
template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetTubeDataIndexList( void )
{
  return m_TubeDataIndexList;
}

/**
 * Get Bkg Index List */
template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetBkgDataIndexList( void )
{
  return m_BkgDataIndexList;
}

/**
 * Update */
template< class TPixel, unsigned int VDimension >
void
ComputeSegmentTubesParameters< TPixel, VDimension >
::Update( void )
{
  typename InputImageType::RegionType region =
    m_InputImage->GetLargestPossibleRegion();

  itk::ImageRegionIteratorWithIndex< InputImageType > itI(
    m_InputImage, region );
  itk::ImageRegionIteratorWithIndex< MaskImageType > itM(
    m_MaskInputImage, region );
  itk::ImageRegionIteratorWithIndex< ScaleImageType > itS(
    m_ScaleInputImage, region );

  typename RidgeExtractorType::Pointer ridgeExtractor = RidgeExtractorType::New();
  ridgeExtractor->SetInputImage( m_InputImage );

  const double BIGD = 9999999999;

  double scale = 0;
  double scaleMin = BIGD;
  double scaleMax = 0;

  double dataMin = BIGD;
  double dataMax = 0;

  MetricVectorType instance( 5, 0 );
  MetricVectorType instanceMin( 5, BIGD );
  MetricVectorType instanceMax( 5, -BIGD );
  typename RidgeExtractorType::IndexType minIndx;
  typename RidgeExtractorType::IndexType maxIndx;
  for( unsigned int i = 0; i < VDimension; ++i )
    {
    minIndx[i] = region.GetIndex()[i] + 10;
    maxIndx[i] = region.GetIndex()[i] + region.GetSize()[i] - 10;
    }

  ridgeExtractor->SetMinRidgeness( 0.6 );
  ridgeExtractor->SetMinRidgenessStart( 0.6 );
  ridgeExtractor->SetMinRoundness( -1 );
  ridgeExtractor->SetMinRoundnessStart( -1 );
  ridgeExtractor->SetMinCurvature( -1 );
  ridgeExtractor->SetMinCurvatureStart( -1 );
  ridgeExtractor->SetMinLevelness( -1 );
  ridgeExtractor->SetMinLevelnessStart( -1 );

  while( !itM.IsAtEnd() )
    {
    if( itM.Get() == m_MaskBackGroundId
      || itM.Get() == m_MaskTubeId )
      {
      typename RidgeExtractorType::IndexType indx;
      typename RidgeExtractorType::ContinuousIndexType cIndx;
      indx = itM.GetIndex();
      bool outOfBounds = false;
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        cIndx[i] = indx[i];
        if( indx[i] < minIndx[i] || indx[i] > maxIndx[i] )
          {
          outOfBounds = true;
          }
        }

      if( !outOfBounds )
        {
        scale = itS.Get();
        ridgeExtractor->SetScale( scale );

        if( itM.Get() == m_MaskTubeId )
          {
          double intensity = 0;
          double ridgeness = 0;
          double roundness = 0;
          double curvature = 0;
          double levelness = 0;
          ridgeness = ridgeExtractor->Ridgeness( cIndx, intensity, roundness, curvature,
            levelness );
          instance[0] = intensity;
          instance[1] = ridgeness;
          instance[2] = roundness;
          instance[3] = curvature;
          instance[4] = levelness;
          m_SeedData.push_back( instance );
          m_SeedDataIndexList.push_back( cIndx );

          if ( ridgeExtractor->LocalRidge( cIndx ) == RidgeExtractorType::SUCCESS )
            {
            if ( scale < scaleMin )
              {
              scaleMin = scale;
              }
            else if ( scale > scaleMax )
              {
              scaleMax = scale;
              }
            ridgeness = ridgeExtractor->Ridgeness( cIndx, intensity, roundness,
              curvature, levelness );
            instance[0] = intensity;
            instance[1] = ridgeness;
            instance[2] = roundness;
            instance[3] = curvature;
            instance[4] = levelness;
            m_TubeData.push_back( instance );
            m_TubeDataIndexList.push_back( cIndx );
            }
          }
        else if( itM.Get() == m_MaskBackGroundId )
          {
          for( unsigned int i = 0; i < VDimension; ++i )
            {
            cIndx[i] = indx[i];
            }
          double intensity = 0;
          double ridgeness = 0;
          double roundness = 0;
          double curvature = 0;
          double levelness = 0;
          ridgeness = ridgeExtractor->Ridgeness( cIndx, intensity, roundness, curvature,
            levelness );
          instance[0] = intensity;
          instance[1] = ridgeness;
          instance[2] = roundness;
          instance[3] = curvature;
          instance[4] = levelness;
          m_BkgData.push_back( instance );
          m_BkgDataIndexList.push_back( cIndx );
          }
        }
      }
    double intensity = itI.Get();

    if( intensity < dataMin )
      {
      dataMin = intensity;
      }
    if( intensity > dataMax )
      {
      dataMax = intensity;
      }
    ++itI;
    ++itS;
    ++itM;
    }

  itk::tube::MetaTubeExtractor params;

  // Heuristics to identify common intensity ranges
  if( dataMin > -0.5 && dataMin < 0.5 &&
    dataMax > 0.5 && dataMax <= 1.5 )
    {
    // Synthetic: 0 to 1
    dataMin = 0;
    dataMax = 1;
    }
  else if( dataMin > -10 && dataMin < 50 &&
    dataMax > 200 && dataMax <= 300 )
    {
    // MRI or synthetic: 0 to 255
    dataMin = 0;
    dataMax = 255;
    }
  else if( dataMin > -20 && dataMin < 100 &&
    dataMax > 420 && dataMax <= 600 )
    {
    // MRI or synthetic: 0 to 512
    dataMin = 0;
    dataMax = 512;
    }
  else if( dataMin > -1100 && dataMin < -900 &&
    dataMax > 900 && dataMax <= 3100 )
    { // 900 = cancellous bone
    dataMin = -1000;  // Air in HU
    dataMax = 3000;    // Dense bone in HU
    }

  itk::tube::MetaTubeExtractor::VectorType tubeColor( 4, 0.0 );
  tubeColor[0] = 1.0;
  tubeColor[3] = 1.0;

  params.SetGeneralProperties( dataMin, dataMax, tubeColor );

  if( scaleMax == scaleMin )
    {
    scaleMax = 10 * scaleMin;
    }
  double scaleRange = scaleMax - scaleMin;
  double ridgeScale = scaleMin + 0.25 * scaleRange;
  double ridgeScaleKernelExtent = 2.5;

  bool ridgeDynamicScale = true;

  bool ridgeDynamicStepSize = false;

  double ridgeStepX = 0.1;

  double ridgeMaxTangentChange = 0.75;

  double ridgeMaxXChange = 3.0;

  double portion = 1.0 / 1000.0;
  int clippedMax = ( int )( m_TubeData.size() * portion );
  portion = 1.0 / 500.0;
  int clippedMaxStart = ( int )( m_TubeData.size() * portion );

  double ridgeMinRidgeness;
  double ridgeMinRidgenessStart;
  double ridgeMinRoundness;
  double ridgeMinRoundnessStart;
  double ridgeMinCurvature;
  double ridgeMinCurvatureStart;
  double ridgeMinLevelness;
  double ridgeMinLevelnessStart;
  for ( int index = 1; index < 5; index++ )
    {
    m_SortColumn = index;
    std::sort( m_TubeData.begin(), m_TubeData.end(),
      ComputeSegmentTubesParametersSortFunctionType() );
    if( index == 1 )
      {
      ridgeMinRidgeness = m_TubeData[clippedMax][index];
      ridgeMinRidgenessStart = m_TubeData[clippedMaxStart][index];
      }
    else if( index == 2 )
      {
      ridgeMinRoundness = m_TubeData[clippedMax][index];
      ridgeMinRoundnessStart = m_TubeData[clippedMaxStart][index];
      }
    else if( index == 3 )
      {
      ridgeMinCurvature = m_TubeData[clippedMax][index];
      ridgeMinCurvatureStart = m_TubeData[clippedMaxStart][index];
      }
    else if( index == 4 )
      {
      ridgeMinLevelness = m_TubeData[clippedMax][index];
      ridgeMinLevelnessStart = m_TubeData[clippedMaxStart][index];
      }
    }

  int ridgeMaxRecoveryAttempts = 3;

  params.SetRidgeProperties( ridgeScale, ridgeScaleKernelExtent,
    ridgeDynamicScale,
    ridgeDynamicStepSize,
    ridgeStepX,
    ridgeMaxTangentChange,
    ridgeMaxXChange,
    ridgeMinRidgeness, ridgeMinRidgenessStart,
    ridgeMinRoundness, ridgeMinRoundnessStart,
    ridgeMinCurvature, ridgeMinCurvatureStart,
    ridgeMinLevelness, ridgeMinLevelnessStart,
    ridgeMaxRecoveryAttempts );

  double radiusStart = ridgeScale / m_InputImage->GetSpacing()[0];
  double radiusMin = ( 0.2 * scaleMin ) / m_InputImage->GetSpacing()[0];
  if( radiusMin < 0.3 )
    {
    radiusMin = 0.3;
    }
  double radiusMax = ( scaleMax + 0.5 * scaleRange ) /
    m_InputImage->GetSpacing()[0];

  // Should be a function of curvature
  double radiusMinMedialness = ridgeMinCurvature / 100;
  double radiusMinMedialnessStart = ridgeMinCurvatureStart / 100;

  params.SetRadiusProperties( radiusStart,
    radiusMin, radiusMax,
    radiusMinMedialness, radiusMinMedialnessStart );

  params.Write( m_ParameterFile.c_str() );

}

/**
 * PrintSelf */
template< class TPixel, unsigned int VDimension >
void
ComputeSegmentTubesParameters< TPixel, VDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeComputeSegmentTubesParameters_hxx)
