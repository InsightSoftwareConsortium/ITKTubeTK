/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#ifndef itktubeAnisotropicHybridDiffusionImageFilter_h
#define itktubeAnisotropicHybridDiffusionImageFilter_h

#include "itktubeAnisotropicDiffusionTensorFunction.h"
#include "itktubeAnisotropicDiffusionTensorImageFilter.h"
#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"
#include "itktubeSymmetricEigenVectorAnalysisImageFilter.h"

#include <itkDiffusionTensor3D.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>

namespace itk
{

namespace tube
{
/** \class AnisotropicHybridDiffusionImageFilter
 *  This class is an implementation of anisotropic hybrid diffusion with
 *  continuous switch:
 *   Mendrik et al., Noise reduction in computed tomography scans using 3-D
 *   anisotropic hybrid diffusion with continuous switch. IEEE Transactions on
 *   Medical Imaging 28( 10 ), pp. 1585-1594, 2009.
 *
 * \warning Does not handle image directions.  Re-orient images to axial
 * ( direction cosines = identity matrix ) before using this function,
 *
 * \sa AnisotropicDiffusionTensorImageFilter
 * \sa AnisotropicCoherenceEnhancingDiffusionImageFilter
 *
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */


template <class TInputImage, class TOutputImage>
class AnisotropicHybridDiffusionImageFilter : public AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias */
  using Self = AnisotropicHybridDiffusionImageFilter;

  using Superclass = AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>;

  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;


  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ) */
  itkOverrideGetNameOfClassMacro(AnisotropicHybridDiffusionImageFilter);

  /** Convenient type alias */
  using InputImageType = typename Superclass::InputImageType;
  using OutputImageType = typename Superclass::OutputImageType;
  using PixelType = typename Superclass::PixelType;

  typedef typename Superclass::DiffusionTensorImageType DiffusionTensorImageType;

  // Structure tensor type
  using StructureTensorFilterType = StructureTensorRecursiveGaussianImageFilter<InputImageType>;

  /** Dimensionality of input and output data is assumed to be the same.
   * It is inherited from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  using MatrixType = Matrix<double, ImageDimension, ImageDimension>;

  // Define image of matrix pixel type
  using OutputMatrixImageType = Image<MatrixType, ImageDimension>;

  // Define the symmetric tensor pixel type
  using TensorPixelType = SymmetricSecondRankTensor<double, ImageDimension>;
  using TensorImageType = Image<TensorPixelType, ImageDimension>;

  // Define the type for storing the eigenvalue
  using EigenValueArrayType = FixedArray<double, ImageDimension>;

  // Declare the types of the output images
  using EigenAnalysisOutputImageType = Image<EigenValueArrayType, ImageDimension>;

  /** The container type for the update buffer. */
  using UpdateBufferType = OutputImageType;

  /** Define diffusion image nbd type */
  typedef typename Superclass::DiffusionTensorNeighborhoodType DiffusionTensorNeighborhoodType;
  /** Set the contrast parameter for EED */
  itkSetMacro(ContrastParameterLambdaEED, double);

  /** Set the contrast parameter for CED */
  itkSetMacro(ContrastParameterLambdaCED, double);

  /** Set the contrast parameter for Hybrid */
  itkSetMacro(ContrastParameterLambdaHybrid, double);

  /** Set threshold parameter C */
  itkSetMacro(ThresholdParameterC, double);

  /** Set the sigma value for structure tensor computation */
  itkSetMacro(Sigma, double);
  itkSetMacro(SigmaOuter, double);

  /** Set the alpha value for structure tensor computation */
  itkSetMacro(Alpha, double);

  /** Get methods */
  itkGetMacro(ContrastParameterLambdaEED, double);
  itkGetMacro(ContrastParameterLambdaCED, double);
  itkGetMacro(ContrastParameterLambdaHybrid, double);
  itkGetMacro(ThresholdParameterC, double);
  itkGetMacro(Sigma, double);
  itkGetMacro(SigmaOuter, double);

protected:
  AnisotropicHybridDiffusionImageFilter(void);
  ~AnisotropicHybridDiffusionImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Update diffusion tensor image */
  void virtual UpdateDiffusionTensorImage(void) override;

private:
  // purposely not implemented
  AnisotropicHybridDiffusionImageFilter(const Self &);
  void
  operator=(const Self &); // purposely not implemented

  double m_ContrastParameterLambdaEED;
  double m_ContrastParameterLambdaCED;
  double m_ContrastParameterLambdaHybrid;
  double m_ThresholdParameterC;
  double m_Sigma;
  double m_SigmaOuter;
  double m_Alpha;

}; // End class AnisotropicHybridDiffusionImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeAnisotropicHybridDiffusionImageFilter.hxx"
#endif

#endif // End !defined( __itktubeAnisotropicHybridDiffusionImageFilter_h )
