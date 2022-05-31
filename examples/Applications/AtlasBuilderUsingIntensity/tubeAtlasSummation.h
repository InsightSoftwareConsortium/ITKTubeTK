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
#ifndef __tubeAtlasSummation_h
#define __tubeAtlasSummation_h

#include "itktubeCompleteImageResampleFilter.h"
#include "itktubeMeanAndSigmaImageBuilder.h"
#include "itktubeMinimizeImageSizeFilter.h"
#include "itktubeRobustMeanAndSigmaImageBuilder.h"

#include <itkAffineTransform.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkResampleImageFilter.h>

namespace tube
{

class AtlasSummation
{
public:
  /**
   * DEFAULT_PIXEL_FILL is the fill value considered not part of the
   * image when image is resampled. This tells the processor not to count
   * that value
   */
  enum { Dimension = 3, DEFAULT_PIXEL_FILL = 0 };

  typedef float                                      InputPixelType;
  typedef float                                      CountPixelType;
  typedef float                                      OutputPixelType;
  typedef float                                      MeanPixelType;
  typedef float                                      VariancePixelType;

  typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  typedef itk::Image< MeanPixelType, Dimension >     MeanImageType;
  typedef itk::Image< VariancePixelType, Dimension > VarianceImageType;
  typedef itk::Image< CountPixelType, Dimension >    CountImageType;
  typedef itk::AffineTransform<>                     TransformType;

  typedef InputImageType::Pointer                    InputImagePointer;
  typedef InputImageType::ConstPointer               InputImageConstPointer;
  typedef MeanImageType::Pointer                     MeanImagePointer;
  typedef VarianceImageType::Pointer                 VarianceImagePointer;
  typedef TransformType::Pointer                     TransformPointer;

  typedef InputImageType::SpacingType                SpacingType;
  typedef InputImageType::SizeType                   SizeType;
  typedef InputImageType::PointType                  PointType;

private:

  /** Pixel and Image Type for processing transition images */
  typedef double                                     ProcessPixelType;
  typedef itk::Image< ProcessPixelType, Dimension >  ProcessImageType;
  typedef ProcessImageType::Pointer                  ProcessImagePointer;

  typedef itk::ImageRegionConstIterator< InputImageType >
    InputConstIteratorType;

  typedef itk::ImageRegionIterator< InputImageType >    InputIteratorType;
  typedef itk::ImageRegionIterator< ProcessImageType >  ProcessIteratorType;
  typedef itk::ImageRegionIterator< MeanImageType >     MeanIteratorType;
  typedef itk::ImageRegionIterator< VarianceImageType > VarianceIteratorType;
  typedef itk::ImageRegionIterator< CountImageType >    CountIteratorType;

  typedef std::vector<InputImagePointer>                MedianImageListType;

  typedef itk::tube::MeanAndSigmaImageBuilder< InputImageType,
    MeanImageType, VarianceImageType >               RobustMeanBuilderType;

public:

  /** CTOR, DTOR */
  AtlasSummation( void );
  ~AtlasSummation( void );

  /** Initiate image for new atlas images to be inputed */
  void Clear( void ) { m_IsProcessing = false; }

  /** Add image with or without a transform */
  void AddImage( InputImageType::Pointer );

  /** Add image WITH transform -- Receives Moving -> Fixed Image
   * Transform */
  void AddImage( InputImageType::Pointer, TransformType::Pointer );

  /** Build Mean and variance image & end AddImage() addition abilities */
  void Finalize( void );

  /**
   * Return final Summation products-Mean ( or median ) & Variance
   * ( or standard deviation & image count for # of valid images
   */
  MeanImageType * GetMeanImage( void ) const
    { return m_MeanBuilder->GetOutputMeanImage(); }

  VarianceImageType * GetVarianceImage( void ) const
    { return m_MeanBuilder->GetOutputSigmaImage(); }

  CountImageType * GetValidCountImage( void ) const
    { return m_MeanBuilder->GetValidCountImage(); }

  /**
   * OPTIONAL PARAMETERS
   *
   * Set output size and spacing values to have images resampled as well.
   *
   * Note: Parameters need to be set before adding any images!
   */
  void SetOutputSpacing( const SpacingType& outputSpacing )
    {
    m_OutputSpacing = outputSpacing;
    m_OutputSpacingSet = true;
    }

  void SetOutputSize( const SizeType& outputSize )
    {
    m_OutputSize = outputSize;
    m_OutputSizeSet = true;
    }

  void SetOutputOrigin( const PointType& outputOrigin )
    {
    m_OutputOrigin = outputOrigin;
    m_OutputOriginSet = true;
    }

  const SpacingType GetOutputSpacing( void ) const
    { return m_OutputSpacing; }

  const SizeType GetOutputSize( void ) const
    { return m_OutputSize; }

  const PointType GetOutputOrigin( void ) const
    { return m_OutputOrigin; }

  /**
   * Do we want to use the variance ( S^2 ) or standard deviation ( S )
   * Default is to use standard deviation
   */
  bool GetUseStdDeviation( void ) const
    { return m_UseStdDeviation; }

  void SetUseStdDeviation( bool useStdDeviation )
    { m_UseStdDeviation = useStdDeviation; }

  /**
   * Set the minimum number of contributing images to a pixel
   * to consider that pixel valid for the mean and variance images,
   * default is 1
   */
  void SetImageCountThreshold( unsigned int imageCountThreshold )
    { m_ImageCountThreshold = imageCountThreshold; }

  unsigned int GetImageCountThreshold( void ) const
    { return m_ImageCountThreshold; }

  /**
   * Use the median as a location estimate instead of the mean;
   * Note: Needs to be called BEFORE adding the first image!!!
   */
  void UseMedian( unsigned int numOfImages )
    { m_NumOfImages = numOfImages; }

  bool UseMedian( void ) const
    { return ( m_NumOfImages > 0 ); }

  /**
   * Adjust all the resampled images origins and size ( if not already
   * defined ) so that no elements are cut off due to transforming &
   * resampling the image off of the screen.
   *
   * Note: This can require significant more processing time, but insures
   * that the first image entered does not dictate the final product image.
   * Otherwise the output mean is assumed to be identical to the first
   * image as far as size and origin is at ( 0,0,0 ).
   */
  void AdjustResampledImageSize( bool adjustResampledImageSize )
    { m_AdjustResampledImageSize = adjustResampledImageSize; }

  void AdjustResampledImageOrigin( bool adjustResampledImageOrigin )
    { m_AdjustResampledImageOrigin = adjustResampledImageOrigin; }


private:
  /**
   * Builds appropriate bounding box for size and origin so that
   * none of original image is cut off
   */
  void GetProperRegion( InputImageType::Pointer, ProcessImageType::Pointer,
    InputImageType::RegionType& );

  /** Returns the clipped images */
  InputImagePointer GetClippedImage( InputImagePointer image,
    TransformType::Pointer t );

  /**
   * Update the set output parameters ( but not spacing ) to include the area
   * given by inputed parameters and return true if the output variables
   * changed.
   */
  bool UpdateOutputProperties( SizeType inputSize, PointType inputOrigin,
    SpacingType inputSpacing );

  /**
   * Update the input size parameter to match the output size
   * and return true if the output size variable must be changed.
   */
  bool UpdateOutputSizeParameter( SizeType& inputSize );

  /** Resample the given image with transform & parameters */
  InputImagePointer TransformInputImage( InputImagePointer image,
    TransformPointer trans, SizeType size, SpacingType spacing,
    PointType origin );

  void Start( InputImageType::Pointer );
  void SumImage( InputImageType::Pointer );
  void WriteImage( MeanImageType::Pointer, const std::string & );
  void WriteImage( ProcessImagePointer, const std::string & );

  /** Median specific functions */
  MedianImageListType&  GetInputImageList( void )
    { return m_MedianList; }
  void SetupImageList( InputImagePointer example );
  void UpdateImageImageList( void );
  void AddMedianImage( InputImagePointer image );

  RobustMeanBuilderType::Pointer m_MeanBuilder;

  /** Median Image calculation variables */
  MedianImageListType            m_MedianList;
  InputPixelType                 m_MedianDefaultPixelValue;
  unsigned int                   m_NumOfImages;

  /** Output size, spacing & origin values */
  SizeType                       m_OutputSize;
  SpacingType                    m_OutputSpacing;
  PointType                      m_OutputOrigin;

  bool                           m_OutputSizeSet;
  bool                           m_OutputSpacingSet;
  bool                           m_OutputOriginSet;

  /** State markers and Setting values */
  unsigned int                   m_ImageNumber;
  unsigned int                   m_ImageCountThreshold;

  bool                           m_UseStdDeviation; // Defaults to TRUE
  bool                           m_IsProcessing;   // Indicates processing
  bool                           m_AdjustResampledImageSize;
  bool                           m_AdjustResampledImageOrigin;

  int                            m_Count;

}; // End class AtlasSummation

} // End namespace tube

#endif // End !defined( __tubeAtlasSummation_h )
