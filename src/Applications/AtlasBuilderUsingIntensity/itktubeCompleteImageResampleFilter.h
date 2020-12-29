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

#ifndef __itktubeCompleteImageResampleFilter_h
#define __itktubeCompleteImageResampleFilter_h

#include <itkFixedArray.h>
#include <itkImageFunction.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkSize.h>
#include <itkTransform.h>

namespace itk
{

namespace tube
{

/** \class CompleteImageResampleFilter
 * \brief Resample an image using ResampleImageFilter and return
 * image with modified size and spacing to cover the entire original image
 *
 * Due to output region variability, this is not effect for pipelining
 *
 * Also expects transforms in the identical manner as ResampleImageFilter
 *  ( i.e., Fixed->Moving )
 */
template< class TInputImage, class TOutputImage, class TNonSingularTransform,
          class TInterpolatorPrecisionType = double >
class CompleteImageResampleFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  /** Standard class typedefs. */
  typedef CompleteImageResampleFilter               Self;
  typedef ProcessObject                             Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  typedef TInputImage                               InputImageType;
  typedef TOutputImage                              OutputImageType;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename InputImageType::RegionType       InputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( CompleteImageResampleFilter, ImageToImageFilter );

  /** Number of dimensions. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  /** ResampleImageFilter typedef */
  typedef ResampleImageFilter< InputImageType,
                               OutputImageType,
                               TInterpolatorPrecisionType >
                                                  ResampleImageFilterType;

  /** Transform typedef
   *
   * Note: The transform needs to have an inverse to determine
   * the proper outcome dimensions
   */
  typedef TNonSingularTransform                   TransformType;
  typedef typename TransformType::ConstPointer    TransformPointerType;

  /** Interpolator typedef. */
  typedef InterpolateImageFunction< InputImageType,
                                    TInterpolatorPrecisionType>
                                                  InterpolatorType;
  typedef typename InterpolatorType::Pointer      InterpolatorPointerType;

  /** Image size typedef. */
  typedef Size< itkGetStaticConstMacro( ImageDimension ) >
                                                  SizeType;

  /** Image index typedef. */
  typedef typename TOutputImage::IndexType        IndexType;

  /** Image point typedef. */
  typedef typename InterpolatorType::PointType    PointType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType        PixelType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType       OutputImageRegionType;

  /** Image spacing typedef */
  typedef typename TOutputImage::SpacingType      SpacingType;
  typedef typename TOutputImage::PointType        OriginPointType;

  /*
   * Set the coordinate transformation.
   * Set the coordinate transform to use for resampling.  Note that this
   * must be in index coordinates and is the output-to-input transform,
   * NOT the input-to-output transform that you might naively expect.
   * The default is itk::AffineTransform<TInterpolatorPrecisionType,
   * ImageDimension>.
   */
  itkSetConstObjectMacro( Transform, TransformType );

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( Transform, TransformType );

  /**
   * Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType,
   * TInterpolatorPrecisionType>. Some other options are
   * itk::NearestNeighborInterpolateImageFunction
   * ( useful for binary masks and other images with a small number of
   * possible pixel values ), and itk::BSplineInterpolateImageFunction
   * ( which provides a higher order of interpolation ).  */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro( DefaultPixelValue, PixelType );

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetMacro( DefaultPixelValue, PixelType );

  /** Set the output image spacing. */
  itkSetMacro( OutputSpacing, SpacingType );
  virtual void SetOutputSpacing( const double values[ImageDimension] );

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  itkGetConstObjectMacro( Input, InputImageType );
  itkGetModifiableObjectMacro( Output, OutputImageType );

protected:

  /** This is where the work is done! */
  virtual void GenerateData( void ) override;

  /** Determine the output bounding box for the 3D case */
  void FindOutput3DParameters(
    InputImageConstPointer image,
    TransformPointerType transform,
    SizeType &outputSize,
    PointType& origin ) const;

  CompleteImageResampleFilter( void );
  ~CompleteImageResampleFilter( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const override;

private:

  CompleteImageResampleFilter( const Self & ); //purposely not implemented
  void operator=( const Self & ); //purposely not implemented

  InputImageConstPointer      m_Input;
  OutputImagePointer          m_Output;
  TransformPointerType        m_Transform;         // Coordinate transform
  InterpolatorPointerType     m_Interpolator;      // Interpolation function
  PixelType                   m_DefaultPixelValue; // Default pixel value if
                                                   // the point is outside of
                                                   // the image
  SpacingType                 m_OutputSpacing;     // Output image spacing

}; // End class CompleteImageResampleFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeCompleteImageResampleFilter.hxx"
#endif

#endif // End !defined( __itktubeCompleteImageResampleFilter_h )
