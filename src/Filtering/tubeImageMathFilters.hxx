/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __tubeImageMathFilters_hxx
#define __tubeImageMathFilters_hxx

#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
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

template< unsigned int VDimension >
ImageMathFilters<VDimension>::
ImageMathFilters()
{
  m_Input = nullptr;
}

template< unsigned int VDimension >
ImageMathFilters<VDimension>::
~ImageMathFilters()
{
  m_Input = nullptr;
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
ApplyIntensityWindowing( 
  float valMin, float valMax, float outMin, float outMax )
{
  itk::ImageRegionIterator< ImageType > it2( m_Input,
    m_Input->GetLargestPossibleRegion() );
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

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
ApplyIntensityMultiplicativeBiasCorrection(
  ImageType * inMeanFieldImage )
{
  itk::ImageRegionIterator< ImageType > it2( inMeanFieldImage,
    inMeanFieldImage->GetLargestPossibleRegion() );
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
  itk::ImageRegionIterator< ImageType > it3( m_Input,
    m_Input->GetLargestPossibleRegion() );
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
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
ResampleImage( ImageType * ref )
{
  bool doResample = false;
  for( unsigned int i = 0; i < VDimension; i++ )
    {
    if( m_Input->GetLargestPossibleRegion().GetSize()[i]
      != ref->GetLargestPossibleRegion().GetSize()[i]
      || m_Input->GetLargestPossibleRegion().GetIndex()[i]
      != ref->GetLargestPossibleRegion().GetIndex()[i]
      || m_Input->GetSpacing()[i] != ref->GetSpacing()[i]
      || m_Input->GetOrigin()[i] != ref->GetOrigin()[i] )
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
    typename ImageType::Pointer imTemp;
    filter->SetInput( m_Input );
    filter->SetUseReferenceImage( true );
    filter->SetReferenceImage( ref );
    imTemp = filter->GetOutput();
    filter->Update();
    m_Input = imTemp;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
AddUniformNoise( float valMin, float valMax, float noiseMin,
  float noiseMax, int seed )
{
  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator UniformGenType;
  typename UniformGenType::Pointer uniformGen = UniformGenType::New();
  std::srand( seed );
  uniformGen->Initialize( ( int )seed );

  itk::ImageRegionIterator< ImageType > it2( m_Input,
        m_Input->GetLargestPossibleRegion() );
  it2.GoToBegin();
  while( !it2.IsAtEnd() )
    {
    double tf = it2.Get();
    if( tf >= valMin && tf <= valMax )
      {
      tf += noiseMin + uniformGen->GetVariate() * (noiseMax-noiseMin);
      it2.Set( ( PixelType )tf );
      }
    ++it2;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
AddGaussianNoise( float valMin, float valMax, float noiseMean,
  float noiseStdDev, int seed )
{
  typedef itk::Statistics::NormalVariateGenerator GaussGenType;
  typename GaussGenType::Pointer gaussGen = GaussGenType::New();
  std::srand( seed );
  gaussGen->Initialize( ( int )seed );

  itk::ImageRegionIterator< ImageType > it2( m_Input,
        m_Input->GetLargestPossibleRegion() );
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

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
AddImages( ImageType * input2,
  float weight1, float weight2 )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( input2,
        input2->GetLargestPossibleRegion() );
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
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
MultiplyImages( ImageType * input2 )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( input2,
        input2->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    it1.Set( it1.Get() * it2.Get() );
    ++it1;
    ++it2;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
MirrorAndPadImage( int numPadVoxels )
{
  typedef itk::MirrorPadImageFilter< ImageType, ImageType > PadFilterType;
  typename PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput( m_Input );
  typename PadFilterType::InputImageSizeType bounds;
  bounds.Fill( numPadVoxels );
  padFilter->SetPadBound( bounds );
  padFilter->Update();
  m_Input = padFilter->GetOutput();
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
NormalizeImage( int normType )
{
  if( normType == 0 )
    {
    typedef itk::NormalizeImageFilter< ImageType, ImageType >
                                                  NormFilterType;
    typename NormFilterType::Pointer normFilter = NormFilterType::New();
    normFilter->SetInput( m_Input );
    normFilter->Update();
    m_Input = normFilter->GetOutput();
    }
  else
    {
    unsigned int nBins = 50;

    itk::ImageRegionIteratorWithIndex< ImageType > it1( m_Input,
          m_Input->GetLargestPossibleRegion() );
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
      //   https://mathworld.wolfram.com/GaussianFunction.html
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

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
FuseImages( ImageType * input2, float offset2 )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( input2,
        input2->GetLargestPossibleRegion() );
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
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
MedianImage( int filterSize )
{
  typedef itk::MedianImageFilter< ImageType, ImageType >
    FilterType;
  typename ImageType::Pointer imTemp;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( m_Input );
  typename ImageType::SizeType radius;
  radius.Fill( filterSize );
  filter->SetRadius( radius );
  imTemp = filter->GetOutput();
  filter->Update();
  m_Input = imTemp;
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>::
ThresholdImage( float threshLow, float threshHigh, float valTrue,
  float valFalse )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
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

//------------------------------------------------------------------------
template< unsigned int VDimension >
double
ImageMathFilters<VDimension>
::ComputeImageStatisticsWithinMaskRange(
    ImageType * mask, float maskThreshLow, float maskThreshHigh,
    int mode )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( mask,
        mask->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  double sum = 0;
  double sumS = 0;
  unsigned int count = 0;
  while( !it1.IsAtEnd() && !it2.IsAtEnd() )
    {
    double maskV = it2.Get();
    if( maskV >= maskThreshLow && maskV <= maskThreshHigh )
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

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::AbsoluteImage( void )
{
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  it1.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    it1.Set( std::fabs( it1.Get() ) );
    ++it1;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::ReplaceValuesOutsideMaskRange(
    ImageType * mask,
    float maskThreshLow, float maskThreshHigh, float valFalse )
{
  ImageMathFilters<VDimension> imf2 = ImageMathFilters<VDimension>();
  imf2.SetInput( mask );
  imf2.ResampleImage( m_Input );
  itk::ImageRegionIterator< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( imf2.GetOutput(),
        imf2.GetOutput()->GetLargestPossibleRegion() );
  it1.GoToBegin();
  it2.GoToBegin();
  while( !it1.IsAtEnd() )
    {
    double tf2 = it2.Get();
    if( ! ( tf2 >= maskThreshLow && tf2 <= maskThreshHigh ) )
      {
      it1.Set( ( PixelType )valFalse );
      }
    ++it1;
    ++it2;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::MorphImage( int mode, int radius, float foregroundValue,
  float backgroundValue )
{
  typedef itk::BinaryBallStructuringElement<PixelType, VDimension>
    BallType;
  BallType ball;
  ball.SetRadius( radius );
  ball.CreateStructuringElement();

  typedef itk::BinaryErodeImageFilter<ImageType, ImageType, BallType>
    ErodeFilterType;
  typedef itk::BinaryDilateImageFilter<ImageType, ImageType, BallType>
    DilateFilterType;
  switch( mode )
    {
    case 0:
      {
      typename ErodeFilterType::Pointer filter =
        ErodeFilterType::New();
      typename ImageType::Pointer imTemp;
      filter->SetBackgroundValue( backgroundValue );
      filter->SetKernel( ball );
      filter->SetErodeValue( foregroundValue );
      filter->SetInput( m_Input );
      imTemp = filter->GetOutput();
      filter->Update();
      m_Input = imTemp;
      break;
      }
    case 1:
      {
      typename DilateFilterType::Pointer filter =
        DilateFilterType::New();
      typename ImageType::Pointer imTemp;
      filter->SetKernel( ball );
      filter->SetDilateValue( foregroundValue );
      filter->SetInput( m_Input );
      imTemp = filter->GetOutput();
      filter->Update();
      m_Input = imTemp;
      break;
      }
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::ReplaceValueWithinMaskRange( ImageType * mask, float maskThreshLow,
  float maskThreshHigh, float imageVal, float newImageVal )
{
  itk::ImageRegionIterator< ImageType > itIm( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > itMask( mask,
        mask->GetLargestPossibleRegion() );
  while( !itIm.IsAtEnd() )
    {
    if( itMask.Get() >= maskThreshLow && itMask.Get() <= maskThreshHigh
      && itIm.Get() == imageVal )
      {
      itIm.Set( newImageVal );
      }
    ++itIm;
    ++itMask;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::BlurImage( float sigma )
{
  typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
    filter;
  typename ImageType::Pointer imTemp;
  for( unsigned int i=0; i<VDimension; i++ )
    {
    filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
    filter->SetInput( m_Input );
    filter->SetNormalizeAcrossScale( true );
    filter->SetSigma( sigma );

    filter->SetOrder(
             itk::GaussianOrderEnum::ZeroOrder );
    filter->SetDirection( i );

    imTemp = filter->GetOutput();
    filter->Update();
    m_Input = imTemp;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::BlurOrderImage( float sigma, int order, int direction )
{
  typename itk::RecursiveGaussianImageFilter< ImageType >::Pointer
    filter;
  filter = itk::RecursiveGaussianImageFilter< ImageType >::New();
  filter->SetInput( m_Input );
  filter->SetNormalizeAcrossScale( true );
  filter->SetSigma( sigma );
  filter->SetDirection( direction );
  switch( order )
    {
    case 0:
      filter->SetOrder( itk::GaussianOrderEnum::ZeroOrder );
      break;
    case 1:
      filter->SetOrder( itk::GaussianOrderEnum::FirstOrder );
      break;
    case 2:
      filter->SetOrder( itk::GaussianOrderEnum::SecondOrder );
      break;
    }
  typename ImageType::Pointer imTemp = filter->GetOutput();
  filter->Update();
  m_Input = imTemp;
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::CopyImageInformation( ImageType * sourceImage )
{
  m_Input->SetOrigin( sourceImage->GetOrigin() );
  m_Input->SetSpacing( sourceImage->GetSpacing() );
  m_Input->SetDirection( sourceImage->GetDirection() );
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
std::vector<double>
ImageMathFilters<VDimension>
::ComputeImageHistogram( unsigned int nBins, float & binMin, float & binSize )
{
  itk::ImageRegionIteratorWithIndex< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  it1.GoToBegin();
  std::vector<double> bin( nBins, 0 );

  double binMax = binMin + binSize*nBins;
  if( binMin == 0 && binSize == 0 )
    {
    binMin = it1.Get();
    binMax = it1.Get();
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
    }
  binSize = (binMax - binMin) / nBins;
  std::cout << "  binMin = " << binMin << std::endl;
  std::cout << "  binMax = " << binMax << std::endl;
  std::cout << "  binSize = " << binSize << std::endl;
  while( !it1.IsAtEnd() )
    {
    double tf = it1.Get();
    tf = ( tf-binMin )/( binMax-binMin ) * nBins;
    if( (int)tf<(int)nBins && (int)tf>=0 )
      {
      bin[( int )tf]++;
      }
    ++it1;
    }
  return bin;
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::CorrectIntensitySliceBySliceUsingHistogramMatching(
    unsigned int numberOfBins, unsigned int numberOfMatchPoints )
{
  typedef itk::Image<PixelType, 2> ImageType2D;
  typedef itk::HistogramMatchingImageFilter< ImageType2D, ImageType2D >
      HistogramMatchFilterType;
  typename HistogramMatchFilterType::Pointer matchFilter;
  typename ImageType2D::Pointer im2DRef = ImageType2D::New();
  typename ImageType2D::Pointer im2DIn = ImageType2D::New();
  typename ImageType2D::SizeType size2D;
  size2D[0] = m_Input->GetLargestPossibleRegion().GetSize()[0];
  size2D[1] = m_Input->GetLargestPossibleRegion().GetSize()[1];
  im2DRef->SetRegions( size2D );
  im2DRef->Allocate();
  im2DIn->SetRegions( size2D );
  im2DIn->Allocate();
  itk::ImageRegionIterator< ImageType > it3D( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it3DSliceStart( m_Input,
        m_Input->GetLargestPossibleRegion() );
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
    zMax = m_Input->GetLargestPossibleRegion().GetSize()[VDimension-1];
    }
  for( z=0; z<VDimension && z<zMax; z++ )
    {
    it2DRef.GoToBegin();
    for( y=0; y<m_Input->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<m_Input->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it3D.Get() );
        ++it2DRef;
        ++it3D;
        }
      }
    }
  while( z<zMax )
    {
    it2DIn.GoToBegin();
    it3DSliceStart = it3D;
    for( y=0; y<m_Input->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<m_Input->GetLargestPossibleRegion().GetSize()[0]; x++ )
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
    for( y=0; y<m_Input->GetLargestPossibleRegion().GetSize()[1]; y++ )
      {
      for( x=0; x<m_Input->GetLargestPossibleRegion().GetSize()[0]; x++ )
        {
        it2DRef.Set( it2DOut.Get() );
        it3D.Set( it2DOut.Get() );
        ++it2DRef;
        ++it2DOut;
        ++it3D;
        }
      }
    ++z;
    }
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::CorrectIntensityUsingHistogramMatching( unsigned int numberOfBins,
  unsigned int numberOfMatchPoints, ImageType * ref )
{
  typedef itk::HistogramMatchingImageFilter< ImageType, ImageType >
      HistogramMatchFilterType;
  typename HistogramMatchFilterType::Pointer matchFilter;
  matchFilter = HistogramMatchFilterType::New();
  matchFilter->SetReferenceImage( ref );
  matchFilter->SetInput( m_Input );
  matchFilter->SetNumberOfHistogramLevels( numberOfBins );
  matchFilter->SetNumberOfMatchPoints( numberOfMatchPoints );
  matchFilter->Update();
  m_Input = matchFilter->GetOutput();
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::Resize( double factor )
{
  typename ImageType::Pointer imSub2 = ImageType::New();
  imSub2->CopyInformation( m_Input );
  typename ImageType::SizeType size;
  typename ImageType::SpacingType spacing;
  if( factor != 0 )
    {
    for( unsigned int i=0; i<VDimension; i++ )
      {
      size[i] = ( long unsigned int )
                ( m_Input->GetLargestPossibleRegion().GetSize()[i]
                  / factor );
      spacing[i] = m_Input->GetSpacing()[i]*factor;
      }
    }
  else
    {
    for( unsigned int i=0; i<VDimension; i++ )
      {
      spacing[i] = m_Input->GetSpacing()[i];
      }

    double meanSpacing = ( spacing[0] + spacing[1] ) / 2;
    if( VDimension == 3 )
      {
      meanSpacing = ( meanSpacing + spacing[VDimension-1] ) / 2;
      }
    factor = meanSpacing/spacing[0];
    size[0] = ( long unsigned int )
              ( m_Input->GetLargestPossibleRegion().GetSize()[0]/factor );
    factor = meanSpacing/spacing[1];
    size[1] = ( long unsigned int )
              ( m_Input->GetLargestPossibleRegion().GetSize()[1]/factor );
    spacing[0] = meanSpacing;
    spacing[1] = meanSpacing;
    if( VDimension == 3 )
      {
      factor = meanSpacing/spacing[VDimension-1];
      size[VDimension-1] = ( long unsigned int )
                ( m_Input->GetLargestPossibleRegion().GetSize()[VDimension-1]
                  / factor );
      spacing[VDimension-1] = meanSpacing;
      }
    }
  imSub2->SetRegions( size );
  imSub2->SetSpacing( spacing );
  imSub2->Allocate();

  this->ResampleImage( imSub2 );
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::Resize( ImageType * ref )
{
  this->ResampleImage( ref );
}

//------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::ExtractSlice( unsigned int dimension, unsigned int slice )
{
  typedef itk::ExtractImageFilter<ImageType, ImageType>
    ExtractSliceFilterType;

  typename ExtractSliceFilterType::Pointer filter =
    ExtractSliceFilterType::New();

  typename ImageType::SizeType size =
    m_Input->GetLargestPossibleRegion().GetSize();

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

  filter->SetInput( m_Input );
  filter->SetExtractionRegion( desiredRegion );
  filter->UpdateLargestPossibleRegion();
  filter->Update();

  m_Input = filter->GetOutput();
}

template< unsigned int VDimension >
void
ImageMathFilters<VDimension>
::EnhanceVessels( double scaleMin, double scaleMax, double numScales )
{
  double logScaleStep = ( std::log( scaleMax ) - std::log( scaleMin ) )
    / ( numScales-1 );

  typedef itk::tube::NJetImageFunction< ImageType > ImageFunctionType;
  typename ImageFunctionType::Pointer imFunc = ImageFunctionType::New();
  imFunc->SetInputImage( m_Input );

  typename ImageType::Pointer input2 = ImageType::New();
  input2->SetRegions( m_Input->GetLargestPossibleRegion() );
  input2->SetOrigin( m_Input->GetOrigin() );
  input2->SetSpacing( m_Input->GetSpacing() );
  input2->CopyInformation( m_Input );
  input2->Allocate();

  itk::ImageRegionIteratorWithIndex< ImageType > it1( m_Input,
        m_Input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > it2( input2,
        input2->GetLargestPossibleRegion() );

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
ImageMathFilters<VDimension>
::SegmentUsingConnectedThreshold( float threshLow, float threshHigh,
  float labelValue, float x, float y, float z )
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

  filter->SetInput( m_Input );
  filter->SetLower( threshLow );
  filter->SetUpper( threshHigh );
  filter->AddSeed( seed );
  filter->SetReplaceValue( labelValue );
  filter->Update();

  m_Input = filter->GetOutput();
}

template< unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ImageMathFilters<VDimension>
::ComputeVoronoiTessellation( unsigned int numberOfCentroids,
  unsigned int numberOfIterations, unsigned int numberOfSamples )
{
  typedef itk::tube::CVTImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( m_Input );
  filter->SetNumberOfSamples( numberOfSamples );
  filter->SetNumberOfCentroids( numberOfCentroids );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetNumberOfSamplesPerBatch( numberOfIterations );
  filter->Update();

  std::vector< itk::ContinuousIndex< double, VDimension > > centroids;
  centroids = filter->GetCentroids();

  m_VoronoiTessellationAdjacencyMatrix = filter->GetAdjacencyMatrix();

  return centroids;
}

} // End namespace tube

#endif // End !defined( __tubeImageMathFilters_hxx )
