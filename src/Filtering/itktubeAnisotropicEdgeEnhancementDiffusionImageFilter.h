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

#ifndef itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h
#  define itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h

#  include "itktubeAnisotropicDiffusionTensorFunction.h"
#  include "itktubeAnisotropicDiffusionTensorImageFilter.h"
#  include "itktubeStructureTensorRecursiveGaussianImageFilter.h"
#  include "itktubeSymmetricEigenVectorAnalysisImageFilter.h"

#  include <itkDiffusionTensor3D.h>
#  include <itkHessianRecursiveGaussianImageFilter.h>
#  include <itkSymmetricEigenAnalysisImageFilter.h>

namespace itk
{

namespace tube
{
/** \class AnisotropicEdgeEnhancementDiffusionImageFilter
 *  This class is an implementation of Edge-enhancing diffusion:
 *   Mendrik et al., Noise reduction in computed tomography scans using 3-D
 *   anisotropic hybrid diffusion with continuous switch. IEEE Transactions on
 *   Medical Imaging 28( 10 ), pp. 1585-1594, 2009.
 *
 * \warning Does not handle image directions.  Re-orient images to axial
 * ( direction cosines = identity matrix ) before using this function.
 *
 * \sa AnisotropicDiffusionTensorImageFilter
 * \sa AnisotropicCoherenceEnhancingDiffusionImageFilter
 * \sa AnisotropicHybridDiffusionImageFilter
 *
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */


template <class TInputImage, class TOutputImage>
class AnisotropicEdgeEnhancementDiffusionImageFilter
  : public AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias */
  using Self = AnisotropicEdgeEnhancementDiffusionImageFilter;

  using Superclass = AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>;

  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;


  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ) */
  itkOverrideGetNameOfClassMacro(AnisotropicEdgeEnhancementDiffusionImageFilter);

  /** Convenient type alias */
  using InputImageType = typename Superclass::InputImageType;
  using OutputImageType = typename Superclass::OutputImageType;
  using PixelType = typename Superclass::PixelType;

  typedef typename Superclass::DiffusionTensorImageType DiffusionTensorImageType;

  // Structure tensor type
  using StructureTensorFilterType = StructureTensorRecursiveGaussianImageFilter<InputImageType>;

  /** Dimensionality of input and output data is assumed to be the same.
   * It is inherited from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  using MatrixType = itk::Matrix<double, ImageDimension, ImageDimension>;

  // Define image of matrix pixel type
  using OutputMatrixImageType = itk::Image<MatrixType, ImageDimension>;

  // Define the symmetric tensor pixel type
  using TensorPixelType = itk::SymmetricSecondRankTensor<double, ImageDimension>;
  using TensorImageType = itk::Image<TensorPixelType, ImageDimension>;

  // Define the type for storing the eigenvalue
  using EigenValueArrayType = itk::FixedArray<double, ImageDimension>;

  // Declare the types of the output images
  using EigenAnalysisOutputImageType = itk::Image<EigenValueArrayType, ImageDimension>;

  /** The container type for the update buffer. */
  using UpdateBufferType = OutputImageType;

  /** Define diffusion image nbd type */
  typedef typename Superclass::DiffusionTensorNeighborhoodType DiffusionTensorNeighborhoodType;

  /** Set the contrast parameter */
  itkSetMacro(ContrastParameterLambdaE, double);

  /** Set threshold parameter C */
  itkSetMacro(ThresholdParameterC, double);

  /** Set the sigma and outer sigma value for structure tensor computation */
  itkSetMacro(Sigma, double);
  itkSetMacro(SigmaOuter, double);

  /**Get methods */
  itkGetMacro(ContrastParameterLambdaE, double);
  itkGetMacro(ThresholdParameterC, double);
  itkGetMacro(Sigma, double);
  itkGetMacro(SigmaOuter, double);


protected:
  AnisotropicEdgeEnhancementDiffusionImageFilter(void);
  ~AnisotropicEdgeEnhancementDiffusionImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Update diffusion tensor image */
  void
  UpdateDiffusionTensorImage(void) override;

private:
  // purposely not implemented
  AnisotropicEdgeEnhancementDiffusionImageFilter(const Self &);
  void
  operator=(const Self &); // purposely not implemented

  double m_ContrastParameterLambdaE;
  double m_ThresholdParameterC;
  double m_Sigma;
  double m_SigmaOuter;

}; // End class AnisotropicEdgeEnhancementDiffusionImageFilter

} // End namespace tube

} // End namespace itk

#  ifndef ITK_MANUAL_INSTANTIATION
#    include "itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.hxx"
#  endif

#endif
// End !defined( __itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h )
