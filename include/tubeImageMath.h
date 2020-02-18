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
#ifndef __tubeImageMath_h
#define __tubeImageMath_h

// ITK includes
#include "itkProcessObject.h"
#include "itkCastImageFilter.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "tubeImageMathFilters.h"

namespace tube
{
/** \class ImageMath
 *  \brief Common image processing functions.
 *
 *  \ingroup TubeTK
 */

template< class TInputImage, class TOutputImage=TInputImage >
class ImageMath:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageMath                                  Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef tube::ImageMathFilters< TInputImage::ImageDimension >  FilterType;

  typedef TInputImage                                InputImageType;
  typedef typename FilterType::ImageType             ImageType;
  typedef TOutputImage                               OutputImageType;

  typedef itk::Image< unsigned char, TInputImage::ImageDimension> ImageTypeUChar;
  typedef itk::Image< short, TInputImage::ImageDimension>         ImageTypeShort;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ImageMath, ProcessObject );

  /** Set/Get input image.
   *  This is an in-place filter, so Input will change
   *  with each function call */
  void SetInput( InputImageType * input )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( input );
    castFilter->Update();
    typename ImageType::Pointer output = castFilter->GetOutput();
    output->Register();
    m_Filter.SetInput( output );
    this->Modified(); }

  /** Get current result */

  ImageType * GetInput( void )
  { return m_Filter.GetInput(); }

  /** Get current result */
  OutputImageType * GetOutput( void )
  { typedef itk::CastImageFilter< ImageType,OutputImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( m_Filter.GetOutput() );
    castFilter->Update();
    typename OutputImageType::Pointer output = castFilter->GetOutput();
    output->Register();
    return output; }

  /** Get current result */
  ImageType * GetOutputFloat( void )
  { return m_Filter.GetOutput(); }

  /** Get current result */
  ImageTypeUChar * GetOutputUChar( void )
  { typedef itk::CastImageFilter< ImageType, ImageTypeUChar > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( m_Filter.GetOutput() );
    castFilter->Update();
    typename ImageTypeUChar::Pointer output = castFilter->GetOutput();
    output->Register();
    return output; }

  /** Get current result */
  ImageTypeShort * GetOutputShort( void )
  { typedef itk::CastImageFilter< ImageType, ImageTypeShort > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( m_Filter.GetOutput() );
    castFilter->Update();
    typename ImageTypeShort::Pointer output = castFilter->GetOutput();
    output->Register();
    return output; }

  void IntensityWindow( float inValMin, float inValMax,
    float outMin, float outMax )
  { m_Filter.ApplyIntensityWindowing( inValMin, inValMax, outMin, outMax );
    this->Modified(); };

  void IntensityMultiplicativeBiasCorrection(
    ImageType * inMeanFieldImage )
  { m_Filter.ApplyIntensityMultiplicativeBiasCorrection( inMeanFieldImage );
    this->Modified(); };

  void Resample( ImageType * referenceImage )
  { m_Filter.ResampleImage( referenceImage ); this->Modified(); };

  void AddUniformNoise( float valMin, float valMax, float noiseMin,
    float noiseMax, int seed )
  { m_Filter.AddUniformNoise( valMin, valMax, noiseMin, noiseMax, seed );
    this->Modified(); };

  void AddGaussianNoise( float valMin, float valMax, float noiseMean,
    float noiseRange, int seed )
  { m_Filter.AddGaussianNoise( valMin, valMax, noiseMean, noiseRange, seed );
    this->Modified(); };

  void AddImages( InputImageType * input2, float weight1, float weight2 )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( input2 );
    castFilter->Update();
    m_Filter.AddImages( castFilter->GetOutput(), weight1, weight2 );
    this->Modified(); };

  void MultiplyImages( InputImageType * input2 )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( input2 );
    castFilter->Update();
    m_Filter.MultiplyImages( castFilter->GetOutput() );
    this->Modified(); };

  void PadUsingMirroring( int numPadVoxels )
  { m_Filter.MirrorAndPadImage( numPadVoxels ); this->Modified(); };

  void NormalizeMeanStdDev()
  { m_Filter.NormalizeImage( 0 ); this->Modified(); };

  void NormalizeFWHM()
  { m_Filter.NormalizeImage( 1 ); this->Modified(); };

  void NormalizeMeanShift()
  { m_Filter.NormalizeImage( 2 ); this->Modified(); };

  void FuseUsingMax( InputImageType * input2, float offset2 )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( input2 );
    castFilter->Update();
    m_Filter.FuseImages( castFilter->GetOutput(), offset2 );
    this->Modified(); };

  void MedianFilter( int size )
  { m_Filter.MedianImage( size ); this->Modified(); };

  void Threshold( float threshLow, float threshHigh, float valTrue,
    float valFalse )
  { m_Filter.ThresholdImage( threshLow, threshHigh, valTrue, valFalse);
  this->Modified(); };

  double MeanWithinMaskRange( InputImageType * mask, float maskThreshLow,
    float maskThreshHigh )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( mask );
    castFilter->Update();
    return m_Filter.ComputeImageStatisticsWithinMaskRange( castFilter->GetOutput(),
      maskThreshLow, maskThreshHigh, 0 ); };

  double StdDevWithinMaskRange( InputImageType * mask, float maskThreshLow,
    float maskThreshHigh )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( mask );
    castFilter->Update();
    return m_Filter.ComputeImageStatisticsWithinMaskRange( castFilter->GetOutput(),
      maskThreshLow, maskThreshHigh, 1 ); };

  void AbsoluteValue()
  { m_Filter.AbsoluteImage(); this->Modified(); };

  void ReplaceValuesOutsideMaskRange( InputImageType * mask, float maskThreshLow,
    float maskThreshHigh, float valFalse )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( mask );
    castFilter->Update();
    m_Filter.ReplaceValuesOutsideMaskRange( castFilter->GetOutput(),
      maskThreshLow, maskThreshHigh, valFalse ); this->Modified(); };

  void Erode( int radius, float fgVal, float bkgVal )
  { m_Filter.MorphImage( 0, radius, fgVal, bkgVal ); this->Modified(); };

  void Dilate( int radius, float fgVal, float bkgVal )
  { m_Filter.MorphImage( 1, radius, fgVal, bkgVal ); this->Modified(); };

  void ReplaceValueWithinMaskRange( InputImageType * mask, float maskThreshLow,
    float maskThreshHigh, float imageVal, float newImageVal )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( mask );
    castFilter->Update();
    m_Filter.ReplaceValueWithinMaskRange( castFilter->GetOutput(),
      maskThreshLow, maskThreshHigh, imageVal, newImageVal ); this->Modified(); };

  void Blur( float sigma )
  { m_Filter.BlurImage( sigma ); this->Modified(); };

  void BlurOrder( float sigma, int order, int direction )
  { m_Filter.BlurOrderImage( sigma, order, direction ); this->Modified(); };

  std::vector<double> Histogram( unsigned int nBins )
  {
  m_HistogramBinMin = 0;
  m_HistogramBinSize = 0;
  this->Modified();
  return m_Filter.ComputeImageHistogram( nBins, m_HistogramBinMin,
    m_HistogramBinSize );
  };

  float HistogramBinMin( void )
  { return m_HistogramBinMin; }

  float HistogramBinSize( void )
  { return m_HistogramBinSize; }

  std::vector<double> Histogram( unsigned int nBins, float binMin,
    float binSize )
  {
  m_HistogramBinMin = binMin;
  m_HistogramBinSize = binSize;
  this->Modified();
  return m_Filter.ComputeImageHistogram( nBins, m_HistogramBinMin,
    m_HistogramBinSize );
  }

  void IntensityCorrectionBySlice( unsigned int nBins,
    unsigned int nMatchPoints )
  { m_Filter.CorrectIntensitySliceBySliceUsingHistogramMatching(
  nBins, nMatchPoints ); this->Modified(); };

  void IntensityCorrection( unsigned int nBins, unsigned int nMatchPoints,
    InputImageType * referenceImage )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( referenceImage );
    castFilter->Update();
    m_Filter.CorrectIntensityUsingHistogramMatching( nBins, nMatchPoints,
      castFilter->GetOutput() ); this->Modified(); };

  void Resize( double factor )
  { m_Filter.Resize( factor ); this->Modified(); };

  void Resize( InputImageType * referenceImage )
  { typedef itk::CastImageFilter< InputImageType, ImageType > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( referenceImage );
    castFilter->Update();
    m_Filter.Resize( castFilter->GetOutput() ); this->Modified(); };

  void ExtractSlice( unsigned int dimension, unsigned int slice )
  { m_Filter.ExtractSlice( dimension, slice ); this->Modified(); };

  void EnhanceVessels( double scaleMin, double scaleMax, int numScales )
  { m_Filter.EnhanceVessels( scaleMin, scaleMax, numScales );
  this->Modified(); };

  void ConnectedComponents( float threshLow, float threshHigh, float labelVal,
  float x, float y, float z )
  { m_Filter.SegmentUsingConnectedThreshold( threshLow, threshHigh, labelVal,
  x, y, z ); this->Modified(); };

  std::vector< itk::ContinuousIndex< double, TInputImage::ImageDimension > >
  VoronoiTessellation( unsigned int nCentroids, unsigned int nIters,
    unsigned int nSamples )
  { this->Modified();
  return m_Filter.ComputeVoronoiTessellation( nCentroids, nIters, nSamples ); };

  itk::VariableSizeMatrix<double> GetVoronoiTessellationAdjacencyMatrix( void )
  { return m_Filter.GetVoronoiTessellationAdjacencyMatrix(); };

protected:
  ImageMath( void );
  ~ImageMath() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** ImageMath parameters **/
  ImageMath( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  FilterType m_Filter;

  float m_HistogramBinMin;
  float m_HistogramBinSize;
};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeImageMath.hxx"
#endif

#endif // End !defined( __tubeImageMath_h )
