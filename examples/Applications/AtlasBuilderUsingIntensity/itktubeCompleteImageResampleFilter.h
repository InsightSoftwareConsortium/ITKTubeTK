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

#ifndef itktubeCompleteImageResampleFilter_h
#define itktubeCompleteImageResampleFilter_h

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
template <class TInputImage, class TOutputImage, class TNonSingularTransform, class TInterpolatorPrecisionType = double>
class CompleteImageResampleFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = CompleteImageResampleFilter;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(CompleteImageResampleFilter);

  /** Number of dimensions. */
  static constexpr unsigned int ImageDimension = TOutputImage::ImageDimension;

  /** ResampleImageFilter type alias */
  using ResampleImageFilterType = ResampleImageFilter<InputImageType, OutputImageType, TInterpolatorPrecisionType>;

  /** Transform type alias
   *
   * Note: The transform needs to have an inverse to determine
   * the proper outcome dimensions
   */
  using TransformType = TNonSingularTransform;
  using TransformPointerType = typename TransformType::ConstPointer;

  /** Interpolator type alias. */
  using InterpolatorType = InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>;
  using InterpolatorPointerType = typename InterpolatorType::Pointer;

  /** Image size type alias. */
  using SizeType = Size<Self::ImageDimension>;

  /** Image index type alias. */
  using IndexType = typename TOutputImage::IndexType;

  /** Image point type alias. */
  using PointType = typename InterpolatorType::PointType;

  /** Image pixel value type alias. */
  using PixelType = typename TOutputImage::PixelType;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Image spacing type alias */
  using SpacingType = typename TOutputImage::SpacingType;
  using OriginPointType = typename TOutputImage::PointType;

  /*
   * Set the coordinate transformation.
   * Set the coordinate transform to use for resampling.  Note that this
   * must be in index coordinates and is the output-to-input transform,
   * NOT the input-to-output transform that you might naively expect.
   * The default is itk::AffineTransform<TInterpolatorPrecisionType,
   * ImageDimension>.
   */
  itkSetConstObjectMacro(Transform, TransformType);

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro(Transform, TransformType);

  /**
   * Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType,
   * TInterpolatorPrecisionType>. Some other options are
   * itk::NearestNeighborInterpolateImageFunction
   * ( useful for binary masks and other images with a small number of
   * possible pixel values ), and itk::BSplineInterpolateImageFunction
   * ( which provides a higher order of interpolation ).  */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

  /** Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro(DefaultPixelValue, PixelType);

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetMacro(DefaultPixelValue, PixelType);

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void
  SetOutputSpacing(const double values[ImageDimension]);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  itkGetConstObjectMacro(Input, InputImageType);
  itkGetModifiableObjectMacro(Output, OutputImageType);

protected:
  /** This is where the work is done! */
  virtual void
  GenerateData(void) override;

  /** Determine the output bounding box for the 3D case */
  void
  FindOutput3DParameters(InputImageConstPointer image,
                         TransformPointerType   transform,
                         SizeType &             outputSize,
                         PointType &            origin) const;

  CompleteImageResampleFilter(void);
  ~CompleteImageResampleFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  CompleteImageResampleFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  InputImageConstPointer  m_Input;
  OutputImagePointer      m_Output;
  TransformPointerType    m_Transform;         // Coordinate transform
  InterpolatorPointerType m_Interpolator;      // Interpolation function
  PixelType               m_DefaultPixelValue; // Default pixel value if
                                               // the point is outside of
                                               // the image
  SpacingType m_OutputSpacing;                 // Output image spacing

}; // End class CompleteImageResampleFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeCompleteImageResampleFilter.hxx"
#endif

#endif // End !defined( __itktubeCompleteImageResampleFilter_h )
