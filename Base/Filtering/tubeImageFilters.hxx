/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeImageFilters_hxx
#define __tubeImageFilters_hxx

#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkMirrorPadImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkNormalVariateGenerator.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkMedianImageFilter.h>

#include "itktubeCVTImageFilter.h"
#include "itktubeNJetImageFunction.h"

namespace tube
{

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>::
ApplyIntensityWindowing( typename ImageType::Pointer imIn,
  float valMin, float valMax, float outMin, float outMax )
{
  itk::ImageRegionIterator< ImageType > it2( imIn,
    imIn->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    tf = ( tf-valMin )/( valMax-valMin );
    if( tf<0 )
      {
      tf = 0;
      }
    if( tf>1 )
      {
      tf = 1;
      }
    tf = ( tf * ( outMax-outMin ) ) + outMin;
    it2.Set( ( PixelType )tf );
    ++it2;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>::
ApplyIntensityMultiplicativeWithBiasCorrection( 
  typename ImageType::Pointer imIn,
  const std::string & inMeanFieldFilePath )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( inMeanFieldFilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
      << std::endl;
    return false;
    }
  imIn2 = ResampleImage( imIn2, imIn );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
    imIn2->GetLargestPossibleRegion() );
  int count = 0;
  double mean = 0;
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    mean += tf;
    if( tf != 0 )
      {
      ++count;
      }
    ++it2;
    }
  mean /= count;
  itk::ImageRegionIterator< ImageType > it3( imIn,
    imIn->GetLargestPossibleRegion() );
  it3.GoToBegin();
  it2.GoToBegin();
  while( !it3.IsAtEnd() )
    {
    double tf = it3.Get();
    double tf2 = it2.Get();
    if( tf2 != 0 )
      {
      double alpha = mean / tf2;
      tf = tf * alpha;
      it3.Set( ( PixelType )tf );
      }
    ++it3;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
typename ImageFilters<VDimension>::ImageType::Pointer
ImageFilters<VDimension>
::ResampleImage( typename ImageType::Pointer a,
  typename ImageType::Pointer b )
{
  typename ImageType::Pointer output = a;

  bool doResample = false;
  for( unsigned int i = 0; i < VDimension; i++ )
    {
    if( a->GetLargestPossibleRegion().GetSize()[i]
      != b->GetLargestPossibleRegion().GetSize()[i]
      || a->GetLargestPossibleRegion().GetIndex()[i]
      != b->GetLargestPossibleRegion().GetIndex()[i]
      || a->GetSpacing()[i] != b->GetSpacing()[i]
      || a->GetOrigin()[i] != b->GetOrigin()[i]  )
      {
      doResample = true;
      break;
      }
    }

  if( doResample )
    {
    typedef typename itk::ResampleImageFilter< ImageType,
      ImageType> ResampleFilterType;
    typename ResampleFilterType::Pointer filter =
      ResampleFilterType::New();
    filter->SetInput( a );
    filter->SetUseReferenceImage( true );
    filter->SetReferenceImage( b );
    filter->Update();
    output = filter->GetOutput();
    }

  return output;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>::
AddUniformNoise( typename ImageType::Pointer imIn,
  float valMin, float valMax, float noiseMean, float noiseRange, int seed )
{
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator UniformGenType;
  typename UniformGenType::Pointer uniformGen = UniformGenType::New();
  std::srand( seed );
  uniformGen->Initialize( ( int )seed );

  itk::ImageRegionIterator< ImageType > it2( imIn,
        imIn->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    if( tf >= valMin && tf <= valMax )
      {
      tf += ( ( 2.0 * uniformGen->GetVariate() ) - 1 ) * noiseRange
        + noiseMean;
      it2.Set( ( PixelType )tf );
      }
    ++it2;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::AddGaussianNoise( typename ImageType::Pointer imIn, float valMin,
  float valMax, float noiseMean, float noiseStdDev, int seed )
{
  typedef itk::Statistics::NormalVariateGenerator GaussGenType;
  typename GaussGenType::Pointer gaussGen = GaussGenType::New();
  std::srand( seed );
  gaussGen->Initialize( ( int )seed );

  itk::ImageRegionIterator< ImageType > it2( imIn,
        imIn->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    if( tf >= valMin && tf <= valMax )
      {
      tf += gaussGen->GetVariate()*noiseStdDev+noiseMean;
      it2.Set( ( PixelType )tf );
      }
    ++it2;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::AddImages( typename ImageType::Pointer imIn,
  const std::string & imIn2FilePath, float weight1, float weight2 )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( imIn2FilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    double tf2 = it2.Get();
    double tf = weight1*tf1 + weight2*tf2;
    it1.Set( ( PixelType )tf );
    ++it1;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::MultiplyImages( 
    typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( imIn2FilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    it1.Set( it1.Get() * it2.Get() );
    ++it1;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::MirrorAndPadImage( 
    typename ImageType::Pointer & imIn,
    int numPadVoxels )
{
  typedef itk::MirrorPadImageFilter< ImageType, ImageType > PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput( imIn );
  typename PadFilterType::InputImageSizeType bounds;
  bounds.Fill( numPadVoxels );
  padFilter->SetPadBound( bounds );
  padFilter->Update();
  imIn = padFilter->GetOutput();
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
template< typename TPixel >
void
ImageFilters<VDimension>
::NormalizeImage( 
    typename ImageType::Pointer & imIn,
    int normType )
{
  if( normType == 0 )
    {
    typedef itk::NormalizeImageFilter< ImageType, ImageType >
                                                  NormFilterType;
    typename NormFilterType::Pointer normFilter = NormFilterType::New();
    normFilter->SetInput( imIn );
    normFilter->Update();
    imIn = normFilter->GetOutput();
    }
  else
    {
    unsigned int nBins = 50;

    itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
          imIn->GetLargestPossibleRegion() );
    it1.GoToBegin();
    double binMin = it1.Get();
    double binMax = it1.Get();
    while( !it1.IsAtEnd() )
      {
      double tf = it1.Get();
      if( tf < binMin )
        {
        binMin = tf;
        }
      else
        {
        if( tf > binMax )
          {
          binMax = tf;
          }
        }
      ++it1;
      }

    int loop = 0;
    double meanV = 0;
    double stdDevV = 1;
    itk::Array<double> bin;
    bin.set_size( nBins );
    while( loop++ < 4 )
      {
      std::cout << "binMin = " << binMin << " : binMax = " << binMax
        << std::endl;
      std::cout << "  Mean = " << meanV << " : StdDev = " << stdDevV
        << std::endl;

      if( ( binMax - binMin ) < nBins
        && 1 == static_cast< TPixel >( 1.1 ) )
        {
        std::cout << "Stretching represent int values" << std::endl;
        int binMid = ( binMax + binMin ) / 2.0;
        int binStep = ( nBins - 1 ) / 2;
        binMin = binMid - binStep;
        binMax = binMin + nBins;
        }

      it1.GoToBegin();
      bin.Fill( 0 );
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        tf = ( tf-binMin )/( binMax-binMin ) * ( nBins-1 );
        if( tf>=0 && tf<nBins )
          {
          bin[( int )tf]++;
          if( tf > 0 )
            {
            bin[( int )( tf - 1 )] += 0.5;
            }
          if( tf < nBins-1 )
            {
            bin[( int )( tf + 1 )] += 0.5;
            }
          }
        ++it1;
        }

      int maxBin = 0;
      double maxBinV = bin[0];
      for( unsigned int i=1; i<nBins; i++ )
        {
        if( bin[i] >= maxBinV )
          {
          maxBinV = bin[i];
          maxBin = i;
          }
        }
      double fwhm = maxBinV / 2;
      double binFWHMMin = maxBin;
      while( binFWHMMin>0 && bin[( int )binFWHMMin]>=fwhm )
        {
        --binFWHMMin;
        }
      std::cout << "  binfwhmmin = " << binFWHMMin
        << std::endl;
      binFWHMMin += ( fwhm - bin[( int )binFWHMMin] )
        / ( bin[( int )binFWHMMin+1] - bin[( int )binFWHMMin] );
      std::cout << "  tweak: binfwhmmin = " << binFWHMMin
        << std::endl;

      double binFWHMMax = maxBin;
      while( binFWHMMax<( int )nBins-1 && bin[( int )binFWHMMax]>=fwhm )
        {
        ++binFWHMMax;
        }
      std::cout << "  binfwhmmax = " << binFWHMMax
        << std::endl;
      binFWHMMax -= ( fwhm - bin[( int )binFWHMMax] )
        / ( bin[( int )binFWHMMax-1] - bin[( int )binFWHMMax] );
      std::cout << "  tweak: binfwhmmax = " << binFWHMMax
        << std::endl;
      if( binFWHMMax <= binFWHMMin )
        {
        binFWHMMin = maxBin-1;
        binFWHMMax = maxBin+1;
        }

      double minV = ( ( binFWHMMin + 0.5 ) / ( nBins - 1.0 ) )
        * ( binMax-binMin ) + binMin;
      double maxV = ( ( binFWHMMax + 0.5 ) / ( nBins - 1.0 ) )
        * ( binMax-binMin ) + binMin;

      meanV = ( maxV + minV ) / 2.0;

      // FWHM to StdDev relationship from
      //   http://mathworld.wolfram.com/GaussianFunction.html
      stdDevV = ( maxV - minV ) / 2.3548;

      binMin = meanV - 1.5 * stdDevV;
      binMax = meanV + 1.5 * stdDevV;
      }

    std::cout << "FINAL: binMin = " << binMin << " : binMax = "
      << binMax << std::endl;
    std::cout << "  Mean = " << meanV << " : StdDev = " << stdDevV
      << std::endl;

    it1.GoToBegin();
    if( normType == 1 )
      {
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        it1.Set( ( tf - meanV ) / stdDevV );
        ++it1;
        }
      }
    else
      {
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        it1.Set( tf - meanV );
        ++it1;
        }
      }
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::FuseImages( 
    typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath,
    float offset2 )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( imIn2FilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    double tf2 = it2.Get();
    if( tf2>tf1 )
      {
      double tf = offset2 + tf2;
      it1.Set( ( PixelType )tf );
      }
    ++it1;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::ThresholdImage( 
    typename ImageType::Pointer imIn,
    float threshLow, float threshHigh, float valTrue, float valFalse )
{
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  it1.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    if( tf1 >= threshLow && tf1 <= threshHigh )
      {
      it1.Set( ( PixelType )valTrue );
      }
    else
      {
      it1.Set( ( PixelType )valFalse );
      }
    ++it1;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::MedianImage( typename ImageType::Pointer & imIn, int filterSize )
{
  typedef itk::MedianImageFilter< ImageType, ImageType >
    FilterType;
  typename ImageType::Pointer imTemp;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( imIn );
  typename ImageType::SizeType radius;
  radius.Fill( filterSize );
  filter->SetRadius( radius );
  imTemp = filter->GetOutput();
  filter->Update();
  imIn = imTemp;
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
double
ImageFilters<VDimension>
::ComputeImageStdDevOrMeanWithinRangeUsingMask( 
    typename ImageType::Pointer imIn,
    const std::string & maskFilePath,
    float threshLow, float threshHigh, int mode, bool & success )
{
  success = true;
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( maskFilePath.c_str() );
  typename ImageType::Pointer imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    success = false;
    return 0.0;
    }

  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  double sum = 0;
  double sumS = 0;
  unsigned int count = 0;
  while( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    double maskV = it2.Get();
    if( maskV >= threshLow && maskV <= threshHigh )
      {
      sum += it1.Get();
      sumS += it1.Get() * it1.Get();
      ++count;
      }
    ++it1;
    ++it2;
    }
  double mean = sum/count;
  if( mode == 0 )
    {
    return mean;
    }
  else
    {
    double stdDev = ( sumS - ( sum*mean ) )/( count-1 );
    return stdDev;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::AbsoluteImage( typename ImageType::Pointer imIn )
{
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  it1.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    it1.Set( vnl_math_abs( it1.Get() ) );
    ++it1;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::MaskImageWithValueIfNotWithinSecondImageRange( 
    typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath,
    float threshLow, float threshHigh, float valFalse )
{

  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( imIn2FilePath.c_str() );

  typename ImageType::Pointer imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  imIn2 = ResampleImage( imIn2, imIn );
  itk::ImageRegionIterator< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imIn2,
        imIn2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf2 = it2.Get();
    if( tf2 >= threshLow && tf2 <= threshHigh )
      {
      //double tf1 = it1.Get();
      //it1.Set( ( PixelType )tf1 );
      }
    else
      {
      it1.Set( ( PixelType )valFalse );
      }
    ++it1;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::MorphImage( 
    typename ImageType::Pointer & imIn,
    int mode, float radius, float foregroundValue, float backgroundValue )
{
  typedef itk::BinaryBallStructuringElement<PixelType, VDimension>
    BallType;
  BallType ball;
  ball.SetRadius( radius );
  ball.CreateStructuringElement();

  typedef itk::ErodeObjectMorphologyImageFilter
               <ImageType, ImageType, BallType>       ErodeFilterType;
  typedef itk::DilateObjectMorphologyImageFilter
               <ImageType, ImageType, BallType>       DilateFilterType;
  switch( mode )
    {
    case 0:
      {
      typename ErodeFilterType::Pointer filter =
        ErodeFilterType::New();
      filter->SetBackgroundValue( backgroundValue );
      filter->SetKernel( ball );
      filter->SetObjectValue( foregroundValue );
      filter->SetInput( imIn );
      filter->Update();
      imIn = filter->GetOutput();
      break;
      }
    case 1:
      {
      typename DilateFilterType::Pointer filter =
        DilateFilterType::New();
      filter->SetKernel( ball );
      filter->SetObjectValue( foregroundValue );
      filter->SetInput( imIn );
      filter->Update();
      imIn = filter->GetOutput();
      break;
      }
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::OverwriteImage( 
    typename ImageType::Pointer imIn,
    const std::string & maskFilePath,
    float maskKeyVal, float imageKeyVal, float newImageVal )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( maskFilePath.c_str() );
  reader2->Update();
  typename ImageType::Pointer maskIm = reader2->GetOutput();

  itk::ImageRegionIterator< ImageType > itIm( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > itMask( maskIm,
        maskIm->GetLargestPossibleRegion() );
  while( !itIm.IsAtEnd() )
    {
    if( itMask.Get() == maskKeyVal )
      {
      if( itIm.Get() == imageKeyVal )
        {
        itIm.Set( newImageVal );
        }
      }
    ++itIm;
    ++itMask;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::BlurImage( 
    typename ImageType::Pointer & imIn,
    float sigma )
{
  typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
    filter;
  typename ImageType::Pointer imTemp;
  for( unsigned int i=0; i<VDimension; i++ )
    {
    filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
    filter->SetInput( imIn );
    //filter->SetNormalizeAcrossScale( true );
    filter->SetSigma( sigma );

    filter->SetOrder( 
             itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
    filter->SetDirection( i );

    imTemp = filter->GetOutput();
    filter->Update();
    imIn = imTemp;
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::BlurOrderImage( 
    typename ImageType::Pointer & imIn,
    float sigma, int order, int direction )
{
  typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
    filter;
  filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
  filter->SetInput( imIn );
  //filter->SetNormalizeAcrossScale( true );
  filter->SetSigma( sigma );
  filter->SetDirection( direction );
  switch( order )
    {
    case 0:
      filter->SetOrder( 
        itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      break;
    case 1:
      filter->SetOrder( 
        itk::RecursiveGaussianImageFilter<ImageType>::FirstOrder );
      break;
    case 2:
      filter->SetOrder( 
        itk::RecursiveGaussianImageFilter<ImageType>::SecondOrder );
      break;
    }
  typename ImageType::Pointer imTemp = filter->GetOutput();
  filter->Update();
  imIn = imTemp;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::ComputeImageHistogram( 
    typename ImageType::Pointer imIn,
    unsigned int nBins, const std::string & histOutputFilePath )
{
  itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  it1.GoToBegin();
  double binMin = it1.Get();
  double binMax = it1.Get();
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    if( tf < binMin )
      {
      binMin = tf;
      }
    else
      {
      if( tf > binMax )
        {
        binMax = tf;
        }
      }
    ++it1;
    }
  std::cout << "  binMin = " << binMin << std::endl;
  std::cout << "  binMax = " << binMax << std::endl;
  it1.GoToBegin();
  itk::Array<double> bin;
  bin.set_size( nBins );
  bin.Fill( 0 );
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    tf = ( tf-binMin )/( binMax-binMin ) * nBins;
    if( tf>nBins-1 )
      {
      tf = nBins-1;
      }
    else
      {
      if( tf<0 )
        {
        tf = 0;
        }
      }
    bin[( int )tf]++;
    ++it1;
    }
  std::ofstream writeStream;
  writeStream.open( histOutputFilePath.c_str(), std::ios::binary | std::ios::out );
  if( !writeStream.rdbuf()->is_open() )
    {
    std::cerr << "Cannot write to file : " << histOutputFilePath << std::endl;
    return false;
    }
  for( unsigned int i=0; i<nBins; i++ )
    {
    writeStream << ( i/( double )nBins )*( binMax-binMin )+binMin
                << " " << bin[i] << std::endl;
    }
  writeStream.close();
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::ComputeImageHistogram2( 
    typename ImageType::Pointer imIn,
    unsigned int nBins, double binMin, double binSize,
    const std::string & histOutputFilePath )
{
  double binMax = binMin + binSize*nBins;
  itk::ImageRegionIteratorWithIndex< ImageType > it1( imIn,
        imIn->GetLargestPossibleRegion() );
  it1.GoToBegin();
  itk::Array<double> bin;
  bin.set_size( nBins );
  bin.Fill( 0 );
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    tf = ( tf-binMin )/( binMax-binMin ) * nBins;
    if( tf<nBins && tf>0 )
      {
      bin[( int )tf]++;
      }
    ++it1;
    }
  std::ofstream writeStream;
  writeStream.open( histOutputFilePath.c_str(), std::ios::binary | std::ios::out );
  if( !writeStream.rdbuf()->is_open() )
    {
    std::cerr << "Cannot write to file : " << histOutputFilePath << std::endl;
    return false;
    }
  for( unsigned int i=0; i<nBins; i++ )
    {
    writeStream << ( i/( double )nBins )*( binMax-binMin )+binMin
                << " " << bin[i] << std::endl;
    }
  writeStream.close();
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::CorrectIntensitySliceBySliceUsingHistogramMatching( 
    typename ImageType::Pointer imIn,
    unsigned int numberOfBins, unsigned int numberOfMatchPoints )
{
  typedef itk::Image<PixelType, 2> ImageType2D;
  typedef itk::HistogramMatchingImageFilter< ImageType2D, ImageType2D >
      HistogramMatchFilterType;
  typename HistogramMatchFilterType::Pointer matchFilter;
  typename ImageType2D::Pointer im2DRef = ImageType2D::New();
  typename ImageType2D::Pointer im2DIn = ImageType2D::New();
  typename ImageType2D::SizeType size2D;
  size2D[0] = imIn->GetLargestPossibleRegion().GetSize()[0];
  size2D[1] = imIn->GetLargestPossibleRegion().GetSize()[1];
  im2DRef->SetRegions( size2D );
  im2DRef->Allocate();
  im2DIn->SetRegions( size2D );
  im2DIn->Allocate();
  itk::ImageRegionIterator< ImageType > it3D( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it3DSliceStart( imIn,
        imIn->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType2D > it2DRef( im2DRef,
        im2DRef->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType2D > it2DIn( im2DIn,
        im2DIn->GetLargestPossibleRegion() );
  unsigned int x;
  unsigned int y;
  unsigned int z;
  it3D.GoToBegin();
  unsigned int zMax = 1;
  if( VDimension == 3 )
    {
    zMax = imIn->GetLargestPossibleRegion().GetSize()[VDimension-1];
    }
  for( z=0; z<VDimension && z<zMax; z++ )
    {
    it2DRef.GoToBegin();
    for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it3D.Get() );
        ++it2DRef;
        ++it3D;
        }
      }
    }
  for( ; z<zMax; z++ )
    {
    it2DIn.GoToBegin();
    it3DSliceStart = it3D;
    for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DIn.Set( it3D.Get() );
        ++it2DIn;
        ++it3D;
        }
      }
    matchFilter = HistogramMatchFilterType::New();
    matchFilter->SetReferenceImage( im2DRef );
    matchFilter->SetInput( im2DIn );
    matchFilter->SetNumberOfHistogramLevels( numberOfBins );
    matchFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
    matchFilter->Update();
    itk::ImageRegionIterator< ImageType2D > it2DOut( 
          matchFilter->GetOutput(),
          im2DIn->GetLargestPossibleRegion() );
    it2DRef.GoToBegin();
    it2DOut.GoToBegin();
    it3D = it3DSliceStart;
    for( y=0; y<imIn->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<imIn->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it2DOut.Get() );
        it3D.Set( it2DOut.Get() );
        ++it2DRef;
        ++it2DOut;
        ++it3D;
        }
      }
    }
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::CorrectIntensityUsingHistogramMatching( 
    typename ImageType::Pointer & imIn,
    unsigned int numberOfBins, unsigned int numberOfMatchPoints,
    const std::string & referenceVolumeFilePath )
{
  typedef itk::HistogramMatchingImageFilter< ImageType, ImageType >
      HistogramMatchFilterType;
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( referenceVolumeFilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  typename HistogramMatchFilterType::Pointer matchFilter;
  matchFilter = HistogramMatchFilterType::New();
  matchFilter->SetReferenceImage( imIn2 );
  matchFilter->SetInput( imIn );
  matchFilter->SetNumberOfHistogramLevels( numberOfBins );
  matchFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
  matchFilter->Update();
  imIn = matchFilter->GetOutput();
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::Resize( 
    typename ImageType::Pointer & imIn,
    double factor )
{
  typename ImageType::Pointer imSub2 = ImageType::New();
  imSub2->CopyInformation( imIn );
  typename ImageType::SizeType size;
  typename ImageType::SpacingType spacing;
  if( factor != 0 )
    {
    for( unsigned int i=0; i<VDimension; i++ )
      {
      size[i] = ( long unsigned int )
                ( imIn->GetLargestPossibleRegion().GetSize()[i]
                  / factor );
      spacing[i] = imIn->GetSpacing()[i]*factor;
      }
    }
  else
    {
    for( unsigned int i=0; i<VDimension; i++ )
      {
      spacing[i] = imIn->GetSpacing()[i];
      }

    double meanSpacing = ( spacing[0] + spacing[1] ) / 2;
    if( VDimension == 3 )
      {
      meanSpacing = ( meanSpacing + spacing[VDimension-1] ) / 2;
      }
    factor = meanSpacing/spacing[0];
    size[0] = ( long unsigned int )
              ( imIn->GetLargestPossibleRegion().GetSize()[0]/factor );
    factor = meanSpacing/spacing[1];
    size[1] = ( long unsigned int )
              ( imIn->GetLargestPossibleRegion().GetSize()[1]/factor );
    spacing[0] = meanSpacing;
    spacing[1] = meanSpacing;
    if( VDimension == 3 )
      {
      factor = meanSpacing/spacing[VDimension-1];
      size[VDimension-1] = ( long unsigned int )
                ( imIn->GetLargestPossibleRegion().GetSize()[VDimension-1]
                  / factor );
      spacing[VDimension-1] = meanSpacing;
      }
    }
  imSub2->SetRegions( size );
  imSub2->SetSpacing( spacing );
  imSub2->Allocate();

  imIn = ResampleImage( imIn, imSub2 );
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::Resize( 
    typename ImageType::Pointer & imIn,
    const std::string & imIn2FilePath )
{
  typename VolumeReaderType::Pointer reader2 = VolumeReaderType::New();
  reader2->SetFileName( imIn2FilePath.c_str() );
  typename ImageType::Pointer imIn2;
  imIn2 = reader2->GetOutput();
  try
    {
    reader2->Update();
    }
  catch( ... )
    {
    std::cout << "Problems reading file format of inFile2."
              << std::endl;
    return false;
    }
  imIn = ResampleImage( imIn, imIn2 );
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageFilters<VDimension>
::ExtractSlice( 
    typename ImageType::Pointer & imIn,
    unsigned int dimension,
    unsigned int slice )
{
  typedef itk::ExtractImageFilter<ImageType, ImageType>
    ExtractSliceFilterType;

  typename ExtractSliceFilterType::Pointer filter =
    ExtractSliceFilterType::New();

  typename ImageType::SizeType size =
    imIn->GetLargestPossibleRegion().GetSize();

  typename ImageType::IndexType  extractIndex;
  typename ImageType::SizeType   extractSize;
  for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
    {
    extractIndex[d] = 0;
    extractSize[d] = size[d];
    }

  extractIndex[dimension] = slice;
  extractSize[dimension] = 1;

  typename ImageType::RegionType desiredRegion( extractIndex,
    extractSize );

  filter->SetInput( imIn );
  filter->SetExtractionRegion( desiredRegion );
  filter->UpdateLargestPossibleRegion();
  filter->Update();

  imIn = filter->GetOutput();
}

template< unsigned int VDimension >
void
ImageFilters<VDimension>
::EnhanceVessels( typename ImageType::Pointer imIn,
  double scaleMin, double scaleMax, double numScales )
{
  double logScaleStep = ( std::log( scaleMax ) - std::log( scaleMin ) )
    / ( numScales-1 );

  typedef itk::tube::NJetImageFunction< ImageType > ImageFunctionType;
  typename ImageFunctionType::Pointer imFunc = ImageFunctionType::New();
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

  double ridgeness = 0;
  double intensity = 0;
  double scale = scaleMin;
  std::cout << "   Processing scale " << scale << std::endl;
  it1.GoToBegin();
  it2.GoToBegin();
  typename ImageFunctionType::ContinuousIndexType cIndx;
  while( !it1.IsAtEnd() )
    {
    for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
      {
      cIndx[d] = it1.GetIndex()[d];
      }
    ridgeness = imFunc->RidgenessAtContinuousIndex( cIndx, scale );
    intensity = imFunc->GetMostRecentIntensity();
    double val = ridgeness * intensity;
    it2.Set( ( PixelType )val );
    ++it1;
    ++it2;
    }
  for( unsigned int i=1; i<numScales; i++ )
    {
    scale = std::exp( std::log( scaleMin ) + i * logScaleStep );
    std::cout << "   Processing scale " << scale << std::endl;
    it1.GoToBegin();
    it2.GoToBegin();
    while( !it1.IsAtEnd() )
      {
      for( unsigned int d=0; d<ImageType::ImageDimension; ++d )
        {
        cIndx[d] = it1.GetIndex()[d];
        }
      ridgeness = imFunc->RidgenessAtContinuousIndex( cIndx, scale );
      intensity = imFunc->GetMostRecentIntensity();
      double val = ridgeness * intensity;
      if( val > it2.Get() )
        {
        it2.Set( ( PixelType )val );
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

template< unsigned int VDimension >
void
ImageFilters<VDimension>
::SegmentUsingConnectedThreshold( typename ImageType::Pointer & imIn,
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

template< unsigned int VDimension >
bool
ImageFilters<VDimension>
::ComputeVoronoiTessellation( typename ImageType::Pointer & imIn,
  unsigned int numberOfCentroids, unsigned int numberOfIterations,
  unsigned int numberOfSamples, const std::string & centroidOutFilePath )
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
        if( i > ( int )( VDimension )-1 )
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
      if( i>( int )( VDimension )-1 )
        {
        done = true;
        }
      else
        {
        if( i == ( int )( VDimension )-1 )
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

#endif // End !defined( __tubeImageFilters_hxx )
