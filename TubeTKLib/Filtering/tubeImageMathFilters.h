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

#ifndef __tubeImageMathFilters_h
#define __tubeImageMathFilters_h

#include <itkImageFileReader.h>

namespace tube
{
template< unsigned int VDimension >
class ImageMathFilters
{
public:
  typedef float                                    PixelType;
  typedef itk::Image< PixelType, VDimension >      ImageType;
  typedef itk::ImageFileReader< ImageType >        VolumeReaderType;

  /** Intensity window inVal range to outValRange. */
  static void ApplyIntensityWindowing( typename ImageType::Pointer imIn,
    float inValMin, float inValMax, float outMin, float outMax );

  /** Intensity multiplicative correct using inMeanField. */
  static bool ApplyIntensityMultiplicativeWithBiasCorrection(
    typename ImageType::Pointer imIn,
    const std::string & inMeanFieldFilePath );

  /** Resamples image a to b if they are different, returns resampled_a */
  static typename ImageType::Pointer ResampleImage(
    typename ImageType::Pointer a, typename ImageType::Pointer b );

  /** Adds uniform noise to all pixels within inVal range */
  static void AddUniformNoise( typename ImageType::Pointer imIn,
    float valMin, float valMax, float noiseMean, float noiseRange,
    int seed );

  /** Adds Gaussian noise to all pixels within inVal range */
  static void AddGaussianNoise( typename ImageType::Pointer imIn,
    float valMin, float valMax, float noiseMean, float noiseStdDev,
    int seed );

  /** I( x ) = weight1*I( x ) + weight2*inFile2( x ) */
  static bool AddImages( typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath, float weight1, float weight2 );

  /** I( x ) = I( x ) * inFile2( x ) */
  static bool MultiplyImages( typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath );

  /** Use mirroring to pad an image. */
  static void MirrorAndPadImage( typename ImageType::Pointer & imIn,
    int numPadVoxels );

  /** Normalize: = mean/std; 1 = FWHM ; 2 = FWHM mean ( shift ) only */
  template< typename TPixel >
  static void NormalizeImage( typename ImageType::Pointer & imIn,
    int normType );

  /** Fuse two images by max, applying offset to second image. */
  static bool FuseImages( typename ImageType::Pointer imIn,
    const std::string & imIn2FilePath, float offset2 );

  /** Apply median filter to the image. */
  static bool MedianImage( typename ImageType::Pointer & imIn,
    int size );

  /** If I( x ) in [tLow,tHigh] then I( x )=vTrue else I( x )=vFalse. */
  static void ThresholdImage( typename ImageType::Pointer imIn,
    float threshLow, float threshHigh, float valTrue, float valFalse );

  /** Return image value within masked region ( mode: 0=mean, 1=stdDev ) */
  static double ComputeImageStdDevOrMeanWithinRangeUsingMask(
    typename ImageType::Pointer imIn, const std::string & maskFilePath,
    float threshLow, float threshHigh, int mode, bool & success );

  /** Update image applying 'abs' unary operation. */
  static void AbsoluteImage( typename ImageType::Pointer imIn );

  /** If inFile( x ) in [tLow, tHigh] then I( x )=I( x ) else I( x )=vFalse */
  static bool MaskImageWithValueIfNotWithinSecondImageRange(
    typename ImageType::Pointer imIn, const std::string & imIn2FilePath,
    float threshLow, float threshHigh, float valFalse );

  /** Mathematical morphology using a sphere. Mode: 0=erode, 1=dilate. */
  static void MorphImage( typename ImageType::Pointer & imIn, int mode,
    float radius, float foregroundValue, float backgroundValue );

  /** Replace values within the image, with a mask. */
  static void OverwriteImage( typename ImageType::Pointer imIn,
    const std::string & maskFilePath, float maskKeyVal, float imageKeyVal,
    float newImageVal );

  /** Gaussian blur the image using the given sigma */
  static void BlurImage( typename ImageType::Pointer & imIn, float sigma );

  /** Gaussian blur the image using the given sigma, order and direction. */
  static void BlurOrderImage( typename ImageType::Pointer & imIn,
    float sigma, int order, int direction );

  /** Write the image's histogram to the designated file */
  static bool ComputeImageHistogram( typename ImageType::Pointer imIn,
    unsigned int nBins, const std::string & histOutputFilePath );

  /** Write the image's histogram to the designated file. */
  static bool ComputeImageHistogram2( typename ImageType::Pointer imIn,
    unsigned int nBins, double binMin, double binSize,
    const std::string & histOutputFilePath );

  /** Correct intensity slice-by-slice using HistogramMatchingFilter. */
  static void CorrectIntensitySliceBySliceUsingHistogramMatching(
    typename ImageType::Pointer imIn, unsigned int numberOfBins,
    unsigned int numberOfMatchPoints );

  /** Match intensity to another volume using HistogramMatchingFilter. */
  static bool CorrectIntensityUsingHistogramMatching(
    typename ImageType::Pointer & imIn, unsigned int numberOfBins,
    unsigned int numberOfMatchPoints,
    const std::string & referenceVolumeFilePath );

  /** Resample to reduce by a factor ( factor==0 means make isotropic ) */
  static void Resize( typename ImageType::Pointer & imIn,
    double factor );

  /** Resample to match inFile2. */
  static bool Resize( typename ImageType::Pointer & imIn,
    const std::string & imIn2FilePath );

  /** Extract a single slice from the image. */
  static void ExtractSlice(
    typename ImageType::Pointer & imIn, unsigned int dimension,
    unsigned int slice );

  /** Compute ridgness/vesselness for specified scales. */
  static void EnhanceVessels( typename ImageType::Pointer imIn,
    double scaleMin, double scaleMax, double numScales );

  /** Segment using ( inclusive ) threshold connected components. */
  static void SegmentUsingConnectedThreshold(
    typename ImageType::Pointer & imIn,
    float threshLow, float threshHigh, float labelValue,
    float x, float y, float z );

  /** Run centroid voronoi tessellation on the image. */
  static bool ComputeVoronoiTessellation(
    typename ImageType::Pointer & imIn, unsigned int numberOfCentroids,
    unsigned int numberOfIterations, unsigned int numberOfSamples,
    const std::string & centroidOutFilePath );

private:
  ImageMathFilters();
  ~ImageMathFilters();

}; // End class ImageMathFilters

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeImageMathFilters.hxx"
#endif

#endif // End !defined( __tubeImageMathFilters_h )
