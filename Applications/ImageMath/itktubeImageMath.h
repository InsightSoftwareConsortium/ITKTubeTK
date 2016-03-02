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

#ifndef __itktubeImageMath_h
#define __itktubeImageMath_h

#include <itkImageFileReader.h>

namespace itk
{

namespace tube
{
template< unsigned int VDimension >
class ImageMath
{
public:
  typedef float                                    PixelType;
  typedef itk::Image< PixelType, VDimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >        VolumeReaderType;

  /** Intensity window inVal range to outValRange. */
  static void ApplyIntensityWindowing(
      typename ImageType::Pointer imIn,
      float inValMin, float inValMax, float outMin, float outMax );

  /** Intensity multiplicative correct using inMeanField. */
  static bool ApplyIntensityMultiplicativeWithBiasCorrection(
      typename ImageType::Pointer imIn,
      const std::string & inMeanFieldFilePath );

  /** Resamples image a to b if they are different, returns resampled_a */
  static typename ImageType::Pointer ResampleImage(
      typename ImageType::Pointer a, typename ImageType::Pointer b );

  /** Adds uniform noise to all pixels within inVal range */
  static void AddUniformNoise(
      typename ImageType::Pointer imIn,
      float valMin, float valMax,
      float noiseMean, float noiseRange, int seed );

  /** Adds Gaussian noise to all pixels within inVal range */
  static void AddGaussianNoise(
      typename ImageType::Pointer imIn,
      float valMin, float valMax,
      float noiseMean, float noiseStdDev, int seed );

  /** I( x ) = weight1*I( x ) + weight2*inFile2( x ) */
  static bool AddImages(
      typename ImageType::Pointer imIn, const std::string & imIn2FilePath,
      float weight1, float weight2 );

  /** I( x ) = I( x ) * inFile2( x ) */
  static bool MultiplyImages(
      typename ImageType::Pointer imIn, const std::string & imIn2FilePath );

  /** Use mirroring to pad an image. */
  static void MirrorAndPadImage(
      typename ImageType::Pointer & imIn, int numPadVoxels );

  /** Normalize: type0 = data's mean/std; 1 = FWHM estimate; 2 = FWHM mean (shift) only */
  template< typename TPixel >
  static void NormalizeImage(
      typename ImageType::Pointer & imIn, int normType );

  /** Fuse two images by max, applying offset to second image. */
  static bool FuseImages(
      typename ImageType::Pointer imIn, const std::string & imIn2FilePath,
      float offset2 );

  /** If I(x) in [tLow,tHigh] then I(x)=vTrue else I(x)=vFalse. */
  static void ThresholdImage(
      typename ImageType::Pointer imIn,
      float threshLow, float threshHigh, float valTrue, float valFalse );

  /** Return image value within masked region (mode: 0=mean, 1=stdDev) */
  static double ComputeImageStdDevOrMeanWithinRangeUsingMask(
      typename ImageType::Pointer imIn,
      const std::string & maskFilePath,
      float threshLow, float threshHigh, int mode, bool & success );

private:
  ImageMath();
  ~ImageMath();

}; // End class ImageMath

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeImageMath.hxx"
#endif

#endif // End !defined(__itktubeImageMath_h)
