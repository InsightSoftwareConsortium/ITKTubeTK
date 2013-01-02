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


#ifndef __AtlasSummation_h
#define __AtlasSummation_h

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include "itkCompleteImageResampleFilter.h"
#include "itkMinimizeImageSizeFilter.h"
#include "itkRobustMeanAndSigmaImageBuilder.h"
#include "itkMeanAndSigmaImageBuilder.h"

#include <itkImageFileWriter.h>

namespace tube
{

class AtlasSummation
{
  public:

    /*
     * DEFAULT_PIXEL_FILL is the fill value considered not part of the
     * image when image is resampled. This tells the processor not to count
     * that value
     */
    enum { TDimensions = 3, DEFAULT_PIXEL_FILL = 0 };

    typedef float                                          InputPixelType;
    typedef float                                          CountPixelType;
    typedef float                                          OutputPixelType;
    typedef float                                          MeanPixelType;
    typedef float                                          VariancePixelType;

    typedef itk::Image< InputPixelType, TDimensions >      InputImageType;
    typedef itk::Image< MeanPixelType, TDimensions >       MeanImageType;
    typedef itk::Image< VariancePixelType, TDimensions >   VarianceImageType;
    typedef itk::Image< CountPixelType, TDimensions >      CountImageType;
    typedef itk::AffineTransform<>                         TransformType;

    typedef InputImageType::Pointer                        InputImagePointer;
    typedef InputImageType::ConstPointer                   InputImageConstPointer;
    typedef MeanImageType::Pointer                         MeanImagePointer;
    typedef VarianceImageType::Pointer                     VarianceImagePointer;
    typedef TransformType::Pointer                         TransformPointer;

    typedef InputImageType::SpacingType                    SpacingType;
    typedef InputImageType::SizeType                       SizeType;
    typedef InputImageType::PointType                      PointType;


  private:

    /**  Pixel and Image Type for processing transition images */
    typedef double                                        ProcessPixelType;
    typedef itk::Image< ProcessPixelType, TDimensions >   ProcessImageType;
    typedef ProcessImageType::Pointer                     ProcessImagePointer;

    typedef itk::ImageRegionConstIterator<InputImageType> InputConstIteratorType;
    typedef itk::ImageRegionIterator< InputImageType >    InputIteratorType;
    typedef itk::ImageRegionIterator< ProcessImageType >  ProcessIteratorType;
    typedef itk::ImageRegionIterator< MeanImageType >     MeanIteratorType;
    typedef itk::ImageRegionIterator< VarianceImageType > VarianceIteratorType;
    typedef itk::ImageRegionIterator< CountImageType >    CountIteratorType;

    typedef std::vector<InputImagePointer>                MedianImageListType;

    typedef itk::tube::MeanAndSigmaImageBuilder<
      InputImageType, MeanImageType, VarianceImageType >  RobustMeanBuilderType;


  public:

    /** CTOR, DTOR */
    AtlasSummation();
    ~AtlasSummation();

    /** Initiate image for new atlas images to be inputted */
    void Clear(void) { m_isProcessing = false; }

    /** Add image with or without a transform */
    void AddImage( InputImageType::Pointer );

    /** Add image WITH transform -- Receives Moving -> Fixed Image Transform */
    void AddImage( InputImageType::Pointer, TransformType::Pointer );

    /** Build Mean and variance image & end AddImage() addition abilities */
    void Finalize();

    /*
     * Return final Summation products-Mean (or median) & Variance
     * (or standard deviation & image count for # of valid images
     */
    MeanImageType * GetMeanImage() const
      { return m_meanBuilder->GetOutputMeanImage(); }
    VarianceImageType * GetVarianceImage() const
      { return m_meanBuilder->GetOutputSigmaImage(); }
    CountImageType * GetValidCountImage() const
      { return m_meanBuilder->GetValidCountImage(); }

    /*
     * OPTIONAL PARAMETERS
     *
     * Set output size and spacing values to have images resampled as well.
     *
     * Note: Parameters need to be set before adding any images!
     */
    void SetOutputSpacing( const SpacingType& s )
    {
      m_OutputSpacing = s;
      m_OutSpacingSet = true;
    }

    void SetOutputSize( const SizeType& s )
    {
      m_OutputSize = s;
      m_OutSizeSet = true;
    }

    void SetOutputOrigin( const PointType& o )
    {
      m_OutputOrigin = o;
      m_OutOriginSet = true;
    }

    const SpacingType GetOutputSpacing() const
      { return m_OutputSpacing; }
    const SizeType GetOutputSize() const
      { return m_OutputSize; }
    const PointType GetOutputOrigin() const
      { return m_OutputOrigin; }

    /*
     * Do we want to use the variance (S^2) or std. deviation (S)
     * Default is to use std. deviation
     */
    bool GetUseStdDeviation() const
      { return m_isStdDeviation; }
    void SetUseStdDeviation( bool b )
      { m_isStdDeviation = b; }

    /*
     * Set the minimum number of contributing images to a pixel
     * to consider that pixel valid for the mean and variance images,
     * default is 1
     */
    void SetImageCountThreshold( unsigned int t )
      { m_ImageCountThreshold = t; }
    unsigned int GetImageCountThreshold() const
      { return m_ImageCountThreshold; }

    /*
     * Use the median as a location estimate instead of the mean;
     * Note: Needs to be called BEFORE adding the first image!!!
     */
    void UseMedian( unsigned int numOfImages )
      { m_numOfImages = numOfImages; }
    bool UseMedian() const
      { return (m_numOfImages > 0); }

    /*
     * Adjust all the resampled images origins and size (if not already defined)
     * so that no elements are cut off due to transforming & resampling the image
     * off of the screen.
     *
     * Note: This can require significant more processing time, but insures that
     * the first image entered does not dictate the final product image. Otherwise
     * the output mean is assumed to be identical to the first image as far as size
     * and origin is at (0,0,0).
     */
    void AdjustResampledImageSize( bool b )
      { m_AdjustResampledImageSize = b; }
    void AdjustResampledImageOrigin( bool b )
      { m_AdjustResampledImageOrigin = b; }


  private:

    /*
     * Builds appropriate bounding box for size and origin so that
     * none of original image is cut off
     */
    void GetProperRegion( InputImageType::Pointer,
                          ProcessImageType::Pointer,
                          InputImageType::RegionType& );

    /** Returns the clipped images */
    InputImagePointer GetClippedImage( InputImagePointer image,
                                       TransformType::Pointer t );

    /*
     * Update the set output parameters (but not spacing) to include the area
     * given by inputted parameters and return true if the output variables changed.
     */
    bool UpdateOutputProperties( SizeType inputSize,
                                 PointType inputOrigin,
                                 SpacingType inputSpacing );

    /*
     * Update the input size parameter to match the output size
     * and return true if the output size variable must be changed.
     */
    bool UpdateOutputSizeParameter( SizeType& inputSize );

    /** Resample the given image with transform & parameters */
    InputImagePointer TransformInputImage( InputImagePointer image,
                                           TransformPointer trans,
                                           SizeType size,
                                           SpacingType spacing,
                                           PointType origin );

    void Start( InputImageType::Pointer );
    void SumImage( InputImageType::Pointer );
    void WriteImage( MeanImageType::Pointer, const char *);
    void WriteImage( ProcessImagePointer, const char *);

    /** Median specific functions */
    MedianImageListType&  GetInputImageList() { return m_medianList; }
    void SetupImageList( InputImagePointer example );
    void UpdateImageImageList();
    void AddMedianImage( InputImagePointer image );

    RobustMeanBuilderType::Pointer   m_meanBuilder;

    /** Median Image calculation variables */
    MedianImageListType       m_medianList;
    InputPixelType            TMedianDefaultPixelValue;
    unsigned int              m_numOfImages;

    /** Output size, spacing & origin values */
    SizeType                  m_OutputSize;
    SpacingType               m_OutputSpacing;
    PointType                 m_OutputOrigin;

    bool                      m_OutSizeSet;
    bool                      m_OutSpacingSet;
    bool                      m_OutOriginSet;

    /* State markers and Setting values */
    unsigned int              m_image_number;
    unsigned int              m_ImageCountThreshold;

    bool                      m_isStdDeviation; // Defaults to TRUE
    bool                      m_isProcessing;   // Indicates processing
    bool                      m_AdjustResampledImageSize;
    bool                      m_AdjustResampledImageOrigin;

    int count;

};


} // End namespace tube
#endif
