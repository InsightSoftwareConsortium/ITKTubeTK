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

#ifndef __itktubeImageSegmentationMath_hxx
#define __itktubeImageSegmentationMath_hxx

#include "itktubeCVTImageFilter.h"
#include "itktubeRidgeExtractor.h"

#include <itkConnectedThresholdImageFilter.h>

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageSegmentationMath<VDimension>
::EnhanceVessels(
    typename ImageType::Pointer imIn,
    double scaleMin, double scaleMax, double numScales )
{
  double logScaleStep = (vcl_log(scaleMax) - vcl_log(scaleMin))
    / (numScales-1);

  typedef itk::tube::RidgeExtractor< ImageType > RidgeFuncType;
  typename RidgeFuncType::Pointer imFunc = RidgeFuncType::New();
  imFunc->SetInputImage( imIn );

  typename ImageType::Pointer imIn2 = ImageType::New();
  imIn2->SetRegions( imIn->GetLargestPossibleRegion() );
  imIn2->SetOrigin( imIn->GetOrigin() );
  imIn2->SetSpacing( imIn->GetSpacing() );
  imIn2->CopyInformation( imIn );
  imIn2->Allocate();

  itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );

  double intensity = 0;
  double ridgeness = 0;
  double roundness = 0;
  double curvature = 0;
  double linearity = 0;
  double scale = scaleMin;
  imFunc->SetScale( scale );
  std::cout << "   Processing scale " << scale << std::endl;
  it1.GoToBegin();
  it2.GoToBegin();
  typename RidgeFuncType::ContinuousIndexType cIndx;
  while( !it1.IsAtEnd() )
    {
    for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
      {
      cIndx[d] = it1.GetIndex()[d];
      }
    ridgeness = imFunc->Ridgeness( cIndx, intensity, roundness,
      curvature, linearity );
    it2.Set( ( PixelType )ridgeness );
    ++it1;
    ++it2;
    }
  for( unsigned int i=1; i<numScales; i++ )
    {
    scale = vcl_exp(vcl_log(scaleMin) + i * logScaleStep);
    imFunc->SetScale( scale );
    std::cout << "   Processing scale " << scale << std::endl;
    it1.GoToBegin();
    it2.GoToBegin();
    while( !it1.IsAtEnd() )
      {
      for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
        {
        cIndx[d] = it1.GetIndex()[d];
        }
      ridgeness = imFunc->Ridgeness( cIndx, intensity, roundness,
        curvature, linearity );
      if( ridgeness > it2.Get() )
        {
        it2.Set( ( PixelType )ridgeness );
        }
      ++it1;
      ++it2;
      }
    }
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    it1.Set( it2.Get() );
    ++it1;
    ++it2;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageSegmentationMath<VDimension>
::SegmentUsingConnectedThreshold(
    typename ImageType::Pointer & imIn,
    float threshLow, float threshHigh, float labelValue,
    float x, float y, float z )
{
  typedef itk::ConnectedThresholdImageFilter<ImageType, ImageType>
             FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typename ImageType::IndexType seed;
  seed[0] = ( long int )x;
  seed[1] = ( long int )y;

  if( VDimension == 3 )
    {
    seed[VDimension-1] = ( long int )z;
    }

  filter->SetInput( imIn );
  filter->SetLower( threshLow );
  filter->SetUpper( threshHigh );
  filter->AddSeed( seed );
  filter->SetReplaceValue( labelValue );
  filter->Update();

  imIn = filter->GetOutput();
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageSegmentationMath<VDimension>
::ComputeVoronoiTessellation(
    typename ImageType::Pointer & imIn,
    unsigned int numberOfCentroids,
    unsigned int numberOfIterations,
    unsigned int numberOfSamples,
    const std::string & centroidOutFilePath )
{
  std::string filename = centroidOutFilePath;

  typedef itk::tube::CVTImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( imIn );
  filter->SetNumberOfSamples( numberOfSamples );
  filter->SetNumberOfCentroids( numberOfCentroids );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetNumberOfSamplesPerBatch( numberOfIterations );
  filter->Update();

  std::ofstream writeStream;
  writeStream.open( filename.c_str(),
    std::ios::binary | std::ios::out );
  if( !writeStream.rdbuf()->is_open() )
    {
    std::cerr << "Cannot write to file : " << filename << std::endl;
    return false;
    }
  writeStream << numberOfCentroids << std::endl;
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    for( unsigned int j = 0; j<VDimension; j++ )
      {
      writeStream << ( *( filter->GetCentroids() ) )[i][j];
      if( j<2 )
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  imIn = filter->GetOutput();
  typename ImageType::SizeType size =
    imIn->GetLargestPossibleRegion().GetSize();

  filename = filename + ".mat";

  vnl_matrix<int> aMat( numberOfCentroids, numberOfCentroids );
  aMat.fill( 0 );

  itk::Index<VDimension> indx;
  itk::Index<VDimension> indx2;
  itk::Index<VDimension> indx3;
  indx.Fill( 0 );
  bool done = false;
  int n;
  while( !done )
    {
    int c = ( int )( imIn->GetPixel( indx )-1 );
    indx2.Fill( 0 );
    indx2[0] = 1;
    bool done2 = false;
    while( !done2 )
      {
      bool invalid = false;
      for( unsigned int d=0; d<VDimension; d++ )
        {
        indx3[d] = indx[d] + indx2[d];
        if( indx3[d] >= ( int )size[d] )
          {
          invalid = true;
          break;
          }
        }
      if( !invalid )
        {
        n = ( int )( imIn->GetPixel( indx3 )-1 );
        if( c != n )
          {
          aMat( c, n ) = 1;
          aMat( n, c ) = 1;
          }
        }
      int i=0;
      indx2[i]++;
      while( !done2 && indx2[i]>=2 )
        {
        indx2[i] = 0;
        i++;
        if( i > (int)(VDimension)-1 )
          {
          done2 = true;
          }
        else
          {
          indx2[i]++;
          }
        }
      }
    int i = 0;
    indx[i]++;
    while( !done && indx[i]>=( int )size[i] )
      {
      indx[i] = 0;
      i++;
      if( i>(int)(VDimension)-1 )
        {
        done = true;
        }
      else
        {
        if( i == (int)(VDimension)-1 )
          {
          std::cout << "Computing adjacency of slice : "
                    << indx[VDimension-1] << std::endl;
          }
        indx[i]++;
        }
      }
    }

  writeStream.open( filename.c_str(),
    std::ios::binary | std::ios::out );
  if( !writeStream.rdbuf()->is_open() )
    {
    std::cerr << "Cannot write to file : " << filename << std::endl;
    return 0;
    }
  writeStream << numberOfCentroids << std::endl;
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    for( unsigned int j = 0; j<numberOfCentroids; j++ )
      {
      writeStream << aMat( i, j );
      if( j<numberOfCentroids-1 )
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  return true;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeImageSegmentationMath_hxx)
