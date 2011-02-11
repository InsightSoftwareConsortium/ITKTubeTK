/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkAnisotropicDiffusiveRegistrationFilter_h
#define __itkAnisotropicDiffusiveRegistrationFilter_h

#include "itkDiffusiveRegistrationFilter.h"

namespace itk
{

/** \class itkAnisotropicDiffusiveRegistrationFilter
 * \brief Algorithm for registration of images depicting sliding organs, using
 * an anisotropic diffusive regularization term.
 *
 * Traditional deformable image registration imposes a uniform
 * smoothness constraint on the deformation field. This
 * is not appropriate when registering images visualizing organs
 * that slide relative to each other, and therefore leads to registration
 * inaccuracies.
 *
 * This algorithm includes a deformation field regularization term that is
 * based on anisotropic diffusion and accommodates the deformation field
 * discontinuities that are expected when considering sliding motion.
 *
 * The update term is composed of two parts: an intensity distance term and a
 * regularization term.  The intensity distance term uses the sum of square
 * difference metric, so this registration algorithm is appropriate for
 * monomodal image registration term only.  The regularization term uses a
 * specified border between the organs (stored as a vtkPolyData *) and enforces
 * coupling between the organs while allowing the motion field to exhibit
 * sliding motion at the organ interface.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa itkAnisotropicDiffusiveRegistrationFunction
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT AnisotropicDiffusiveRegistrationFilter
  : public DiffusiveRegistrationFilter< TFixedImage,
                                        TMovingImage,
                                        TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveRegistrationFilter            Self;
  typedef DiffusiveRegistrationFilter< TFixedImage,
                                       TMovingImage,
                                       TDeformationField >  Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, DiffusiveRegistrationFilter);

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::MovingImagePointer       MovingImagePointer;
  typedef typename Superclass::DeformationFieldType     DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer  DeformationFieldPointer;
  typedef typename Superclass::TimeStepType             TimeStepType;
  typedef typename Superclass::FiniteDifferenceFunctionType
      FiniteDifferenceFunctionType;

  /** Typedefs used in multithreading */
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::UpdateBufferType         UpdateBufferType;
  typedef typename Superclass::ThreadRegionType         ThreadRegionType;

  /** Output image and update buffer types */
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  typedef typename Superclass::NeighborhoodType         NeighborhoodType;
  typedef typename Superclass::UpdateBufferRegionType   UpdateBufferRegionType;

  /** The registration function type */
  typedef typename Superclass::RegistrationFunctionType
      RegistrationFunctionType;
  typedef typename Superclass::RegistrationFunctionPointer
      RegistrationFunctionPointer;
  typedef typename Superclass::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename Superclass::SpacingType              SpacingType;

  /** Deformation field types. */
  typedef typename Superclass::DeformationVectorType    DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename Superclass::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename Superclass::DeformationVectorComponentImagePointer
      DeformationVectorComponentImagePointer;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodArrayType
      DeformationVectorComponentNeighborhoodArrayType;
  typedef typename Superclass::ThreadDeformationVectorComponentImageRegionType
      ThreadDeformationVectorComponentImageRegionType;

  /** Normal vector types */
  typedef typename Superclass::NormalVectorType         NormalVectorType;
  typedef typename Superclass::NormalVectorImageType    NormalVectorImageType;
  typedef typename Superclass::NormalVectorImagePointer
      NormalVectorImagePointer;
  typedef typename Superclass::NormalVectorNeighborhoodType
      NormalVectorNeighborhoodType;
  typedef typename Superclass::NormalVectorImageRegionType
      NormalVectorImageRegionType;
  typedef typename Superclass::ThreadNormalVectorImageRegionType
      ThreadNormalVectorImageRegionType;

  /** Diffusion tensor types */
  typedef typename Superclass::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename Superclass::DiffusionTensorImagePointer
      DiffusionTensorImagePointer;
  typedef typename Superclass::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;
  typedef typename Superclass::ThreadDiffusionTensorImageRegionType
      ThreadDiffusionTensorImageRegionType;

  /** The derivative matrix types */
  typedef typename Superclass::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename Superclass::TensorDerivativeImagePointer
      TensorDerivativeImagePointer;
  typedef typename Superclass::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;
  typedef typename Superclass::ThreadTensorDerivativeImageRegionType
      ThreadTensorDerivativeImageRegionType;

  /** Types for weighting between the anisotropic and diffusive (Gaussian)
    * regularization */
  typedef typename Superclass::WeightType               WeightType;
  typedef typename Superclass::WeightImageType          WeightImageType;
  typedef typename Superclass::WeightImagePointer       WeightImagePointer;
  typedef typename Superclass::WeightImageRegionType    WeightImageRegionType;

  /** Organ boundary surface types */
  typedef typename Superclass::BorderSurfacePointer     BorderSurfacePointer;

  /** Types for vector component extractor */
  typedef typename Superclass::VectorIndexSelectionFilterType
      VectorIndexSelectionFilterType;

protected:
  AnisotropicDiffusiveRegistrationFilter();
  virtual ~AnisotropicDiffusiveRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveRegistrationFilter.txx"
#endif

#endif
