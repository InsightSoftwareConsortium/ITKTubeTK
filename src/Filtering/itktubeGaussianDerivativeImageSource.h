/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itktubeGaussianDerivativeImageSource_h
#define itktubeGaussianDerivativeImageSource_h

#include <itkImageFunction.h>
#include <itkMatrix.h>
#include <itkParametricImageSource.h>
#include <itkSize.h>
#include <itkVector.h>

#include <vnl/vnl_c_vector.h>
#include <vnl/vnl_vector.h>

namespace itk
{

namespace tube
{

/** \class GaussianDerivativeImageSource
 * \brief Generate an n-dimensional image of a Gaussian.
 *
 * GaussianImageSource generates an image of a Gaussian.
 *
 * The output image may be of any dimension.
 *
 * \ingroup DataSources
 * \ingroup ITKImageSources
 */

template <typename TOutputImage>
class GaussianDerivativeImageSource : public ParametricImageSource<TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = GaussianDerivativeImageSource;
  using Superclass = ParametricImageSource<TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Typedef for the output image type. */
  using OutputImageType = TOutputImage;

  /** Typedef for the output image PixelType. */
  using OutputImagePixelType = typename TOutputImage::PixelType;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Spacing type alias support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
  using SpacingType = typename TOutputImage::SpacingType;

  /** Origin type alias support.  The origin is the geometric coordinates
   * of the index ( 0,0 ). */
  using PointType = typename TOutputImage::PointType;

  /** Direction type alias support.  The direction is the direction
   * cosines of the image. */
  using DirectionType = typename TOutputImage::DirectionType;

  /** Dimensionality of the output image */
  static constexpr unsigned int NDimensions = TOutputImage::ImageDimension;

  using SigmasType = Vector<double, TOutputImage::ImageDimension>;
  using OrdersType = Vector<int, TOutputImage::ImageDimension>;

  /** Size type matches that used for images */
  using IndexType = typename TOutputImage::IndexType;
  using SizeType = typename TOutputImage::SizeType;
  using SizeValueType = typename TOutputImage::SizeValueType;

  /** Types for parameters. */
  using ParametersValueType = typename Superclass::ParametersValueType;
  using ParametersType = typename Superclass::ParametersType;

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(GaussianDerivativeImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(Index, IndexType);
  itkGetConstReferenceMacro(Index, IndexType);

  itkSetMacro(Sigmas, SigmasType);
  itkGetConstReferenceMacro(Sigmas, SigmasType);

  itkSetMacro(Mean, PointType);
  itkGetConstReferenceMacro(Mean, PointType);

  itkSetMacro(Orders, OrdersType);
  itkGetConstReferenceMacro(Orders, OrdersType);

  /** Set/get the parameters for this source. When this source is
   * templated over an N-dimensional output image type, the first N
   * values in the parameter array are the Sigma parameters in each
   * dimension, the next N values are the Mean parameters in each
   * dimension, and the last value is the Scale. */
  virtual void
  SetParameters(const ParametersType & parameters) override;
  virtual ParametersType
  GetParameters() const override;

  /** Get the number of parameters for this image source. When this
   * source is templated over an N-dimensional output image type, the
   * number of parameters is 2*N+1. */
  virtual unsigned int
  GetNumberOfParameters() const override;

protected:
  GaussianDerivativeImageSource();

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  virtual void
  GenerateOutputInformation() override;

  virtual void
  GenerateData() override;

private:
  GaussianDerivativeImageSource(const GaussianDerivativeImageSource &);
  // purposely not implemented
  void
  operator=(const GaussianDerivativeImageSource &);
  // purposely not implemented

  IndexType m_Index;

  /** The standard deviation in each direction. */
  SigmasType m_Sigmas;

  /** The mean in each direction. */
  PointType m_Mean;

  OrdersType m_Orders;
}; // End class GaussianDerivativeImageSource

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeGaussianDerivativeImageSource.hxx"
#endif

#endif
