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

#ifndef __itktubeImageMath_hxx
#define __itktubeImageMath_hxx

#include <itkBinaryBallStructuringElement.h>
#include <itkCastImageFilter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkMetaImageIO.h>
#include <itkMirrorPadImageFilter.h>
#include <itkNormalizeImageFilter.h>
#include <itkNormalVariateGenerator.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkResampleImageFilter.h>

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMath<VDimension>
::ApplyIntensityWindowing(
    typename ImageType::Pointer imIn,
    float valMin,
    float valMax,
    float outMin,
    float outMax )
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
ImageMath<VDimension>
::ApplyIntensityMultiplicativeWithBiasCorrection(
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
typename ImageMath<VDimension>::ImageType::Pointer
ImageMath<VDimension>
::ResampleImage(
  typename ImageType::Pointer a,
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
    filter->SetUseReferenceImage(true);
    filter->SetReferenceImage(b);
    filter->Update();
    output = filter->GetOutput();
    }

  return output;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMath<VDimension>
::AddUniformNoise(
    typename ImageType::Pointer imIn,
    float valMin,
    float valMax,
    float noiseMean,
    float noiseRange,
    int seed )
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
ImageMath<VDimension>
::AddGaussianNoise(
    typename ImageType::Pointer imIn,
    float valMin,
    float valMax,
    float noiseMean,
    float noiseStdDev,
    int seed )
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
ImageMath<VDimension>
::AddImages(
    typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath,
    float weight1,
    float weight2 )
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
ImageMath<VDimension>
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
  while( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    double tf2 = it2.Get();
    double tf = tf1*tf2;
    it1.Set( ( PixelType )tf );
    ++it1;
    ++it2;
    }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMath<VDimension>
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
ImageMath<VDimension>
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

      if( (binMax - binMin) < nBins
        && 1 == static_cast< TPixel >( 1.1 ) )
        {
        std::cout << "Stretching represent int values" << std::endl;
        int binMid = (binMax + binMin ) / 2.0;
        int binStep = ( nBins - 1 ) / 2;
        binMin = binMid - binStep;
        binMax = binMin + nBins;
        }

      it1.GoToBegin();
      bin.Fill( 0 );
      while( !it1.IsAtEnd() )
        {
        double tf = it1.Get();
        tf = ( tf-binMin )/( binMax-binMin ) * (nBins-1);
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
      while( binFWHMMin>0 && bin[(int)binFWHMMin]>=fwhm )
        {
        --binFWHMMin;
        }
      std::cout << "  binfwhmmin = " << binFWHMMin
        << std::endl;
      binFWHMMin += ( fwhm - bin[(int)binFWHMMin] )
        / ( bin[(int)binFWHMMin+1] - bin[(int)binFWHMMin] );
      std::cout << "  tweak: binfwhmmin = " << binFWHMMin
        << std::endl;

      double binFWHMMax = maxBin;
      while( binFWHMMax<(int)nBins-1 && bin[(int)binFWHMMax]>=fwhm )
        {
        ++binFWHMMax;
        }
      std::cout << "  binfwhmmax = " << binFWHMMax
        << std::endl;
      binFWHMMax -= ( fwhm - bin[(int)binFWHMMax] )
        / ( bin[(int)binFWHMMax-1] - bin[(int)binFWHMMax] );
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
ImageMath<VDimension>
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
    ++it2; }
  return true;
}

//----------------------------------------------------------------------------
template< unsigned int VDimension >
void
ImageMath<VDimension>
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
double
ImageMath<VDimension>
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
    double stdDev = (sumS - (sum*mean))/(count-1);
    return stdDev;
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeImageMath_hxx)
