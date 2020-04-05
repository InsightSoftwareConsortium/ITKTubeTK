/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.
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
#include <itkContinuousIndex.h>

#include <itkImageDuplicator.h>

namespace tube
{

template< unsigned int VDimension >
class ImageMathFilters 
{
public:
  typedef float                                    PixelType;
  typedef itk::Image< PixelType, VDimension >      ImageType;

  ImageMathFilters();
  ~ImageMathFilters();

  void SetInput( ImageType * input )
    {
    typename ImageType::Pointer tmpImage = ImageType::New();
    tmpImage = input;
    typename itk::ImageDuplicator< ImageType >::Pointer dupFilter =
      itk::ImageDuplicator< ImageType >::New();
    dupFilter->SetInputImage( tmpImage );
    dupFilter->Update();
    m_Input = dupFilter->GetOutput();
    }

  ImageType * GetInput( void )
    {
    return m_Input;
    }

  ImageType * GetOutput( void )
    {
    return m_Input;
    }

  /** Intensity window inVal range to outValRange. */
  void ApplyIntensityWindowing( float inValMin, float inValMax,
    float outMin, float outMax );

  /** Intensity multiplicative correct using inMeanField. */
  void ApplyIntensityMultiplicativeBiasCorrection(
    ImageType * inMeanFieldImage );

  /** Resamples image a to b if they are different, returns resampled_a */
  void ResampleImage( ImageType * ref );

  /** Adds uniform noise to all pixels within inVal range */
  void AddUniformNoise( float valMin, float valMax, float noiseMean,
    float noiseRange, int seed );

  /** Adds Gaussian noise to all pixels within inVal range */
  void AddGaussianNoise( float valMin, float valMax, float noiseMean,
    float noiseStdDev, int seed );

  /** I( x ) = weight1*I( x ) + weight2*inFile2( x ) */
  void AddImages( ImageType * input2,
    float weight1, float weight2 );

  /** I( x ) = I( x ) * inFile2( x ) */
  void MultiplyImages( ImageType * input2 );

  /** Use mirroring to pad an image. */
  void MirrorAndPadImage( int numPadVoxels );

  /** Normalize: = mean/std; 1 = FWHM ; 2 = FWHM mean ( shift ) only */
  void NormalizeImage( int normType );

  /** Fuse two images by max, applying offset to second image. */
  void FuseImages( ImageType * input2, float offset2 );

  /** Apply median filter to the image. */
  void MedianImage( int size );

  /** If I( x ) in [tLow,tHigh] then I( x )=vTrue else I( x )=vFalse. */
  void ThresholdImage( float threshLow, float threshHigh, float valTrue,
    float valFalse );

  /** Return image value within masked region ( mode: 0=mean, 1=stdDev ) */
  double ComputeImageStatisticsWithinMaskRange( ImageType * mask,
    float threshLow, float threshHigh, int mode );

  /** Update image applying 'abs' unary operation. */
  void AbsoluteImage( void );

  /** If inFile( x ) in [tLow, tHigh] then I( x )=I( x ) else I( x )=vFalse */
  void ReplaceValuesOutsideMaskRange( ImageType * mask, float maskThreshLow,
    float maskThreshHigh, float valFalse );

  /** Mathematical morphology using a sphere. Mode: 0=erode, 1=dilate. */
  void MorphImage( int mode, int radius, float foregroundValue,
    float backgroundValue );

  /** Replace values within the image, with a mask. */
  void ReplaceValueWithinMaskRange( ImageType * mask, float maskThreshLow,
    float maskThreshHigh, float imageVal, float newImageVal );

  /** Gaussian blur the image using the given sigma */
  void BlurImage( float sigma );

  /** Gaussian blur the image using the given sigma, order and direction. */
  void BlurOrderImage( float sigma, int order, int direction );

  /** Write the image's histogram to the designated file. */
  std::vector<double> ComputeImageHistogram( unsigned int nBins,
    float & binMin, float & binSize );

  /** Correct intensity slice-by-slice using HistogramMatchingFilter. */
  void CorrectIntensitySliceBySliceUsingHistogramMatching(
    unsigned int numberOfBins, unsigned int numberOfMatchPoints );

  /** Match intensity to another volume using HistogramMatchingFilter. */
  void CorrectIntensityUsingHistogramMatching( unsigned int numberOfBins,
    unsigned int numberOfMatchPoints, ImageType * ref );

  /** Resample to reduce by a factor ( factor==0 means make isotropic ) */
  void Resize( double factor );

  /** Resample to match inFile2. */
  void Resize( ImageType * ref );

  /** Extract a single slice from the image. */
  void ExtractSlice( unsigned int dimension, unsigned int slice );

  /** Compute ridgness/vesselness for specified scales. */
  void EnhanceVessels( double scaleMin, double scaleMax,
    double numScales );

  /** Segment using ( inclusive ) threshold connected components. */
  void SegmentUsingConnectedThreshold( float threshLow,
    float threshHigh, float labelValue, float x, float y, float z );

  /** Run centroid voronoi tessellation on the image. */
  std::vector< itk::ContinuousIndex< double, VDimension > >
    ComputeVoronoiTessellation( unsigned int numberOfCentroids,
      unsigned int numberOfIterations, unsigned int numberOfSamples );

  itk::VariableSizeMatrix< double > GetVoronoiTessellationAdjacencyMatrix()
    { return m_VoronoiTessellationAdjacencyMatrix; };

private:

  typename ImageType::Pointer       m_Input;

  itk::VariableSizeMatrix< double > m_VoronoiTessellationAdjacencyMatrix;

}; // End class ImageMathFilters

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeImageMathFilters.hxx"
#endif

#endif // End !defined( __tubeImageMathFilters_h )
