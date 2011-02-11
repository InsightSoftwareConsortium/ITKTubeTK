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
#ifndef __itkAnisotropicDiffusiveRegistrationFunction_h
#define __itkAnisotropicDiffusiveRegistrationFunction_h

#include "itkDiffusiveRegistrationFunction.h"

namespace itk
{

/** \class itkAnisotropicDiffusiveRegistrationFunction
 * \brief Implements the update function for the
 * itkAnisotropicDiffusiveRegistrationFilter
 *
 * Uses the MeanSquareRegistrationFunction as the similarity metric between
 * images.  Takes advantage of the implementation in
 * itkAnisotropicDiffusionTensorFunction to calculate the update term for the
 * regularization.
 *
 * This class uses an anisotropic diffusive regularizer for registration of
 * images depicting sliding organs.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkAnisotropicDiffusiveRegistrationFilter
 * \sa itkDiffusiveRegistrationFunction
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */

template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT AnisotropicDiffusiveRegistrationFunction
  : public DiffusiveRegistrationFunction<TFixedImage,
                                         TMovingImage,
                                         TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveRegistrationFunction            Self;
  typedef DiffusiveRegistrationFunction<TFixedImage,
                                        TMovingImage,
                                        TDeformationField>    Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, DiffusiveRegistrationFunction);

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::MovingImagePointer       MovingImagePointer;
  typedef typename Superclass::DeformationFieldType     DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer
      DeformationFieldTypePointer;
  typedef typename Superclass::TimeStepType             TimeStepType;
  typedef typename Superclass::NeighborhoodType         NeighborhoodType;
  typedef typename Superclass::PixelType                PixelType;
  typedef typename Superclass::FloatOffsetType          FloatOffsetType;

  /** Deformation field types */
  typedef typename Superclass::DeformationVectorType    DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename Superclass::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename
      Superclass::DeformationVectorComponentImageBoundaryConditionType
      DeformationVectorComponentImageBoundaryConditionType;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodArrayType
      DeformationVectorComponentNeighborhoodArrayType;

  /** Typedefs for the intensity-based distance function */
  typedef typename Superclass::IntensityDistanceFunctionType
      IntensityDistanceFunctionType;
  typedef typename Superclass::IntensityDistanceFunctionPointer
      IntensityDistanceFunctionPointer;

  /** Typedefs for the regularization function */
  typedef typename Superclass::RegularizationFunctionType
      RegularizationFunctionType;
  typedef typename Superclass::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename Superclass::SpacingType              SpacingType;

  /** Normal vector types */
  typedef double
      NormalVectorComponentType;
  typedef itk::Vector< NormalVectorComponentType, ImageDimension >
      NormalVectorType;
  typedef itk::Image< NormalVectorType, ImageDimension >
      NormalVectorImageType;
  typedef ZeroFluxNeumannBoundaryCondition< NormalVectorImageType >
      NormalVectorImageBoundaryConditionType;
  typedef ConstNeighborhoodIterator
      < NormalVectorImageType, NormalVectorImageBoundaryConditionType >
      NormalVectorNeighborhoodType;

  /** Typedefs for the diffusion tensor image */
  typedef typename Superclass::DiffusionTensorType      DiffusionTensorType;
  typedef typename Superclass::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename Superclass::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;

  /** Typedefs for the matrices of derivatives */
  typedef typename Superclass::TensorDerivativeType     TensorDerivativeType;
  typedef typename Superclass::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename Superclass::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;

  /** Typedef for the global data type for this class of equations */
  typedef typename Superclass::GlobalDataStruct         GlobalDataStruct;

  /** Set/get whether to use the anisotropic diffusive regularization.  If
   *  false, the weighting term w=0 and Gaussian regularization is used.
   *  Default: true */
  void SetUseAnisotropicRegularization( bool diffuse )
    { m_UseAnisotropicRegularization = diffuse; }
  bool GetUseAnisotropicRegularization() const
    { return m_UseAnisotropicRegularization; }

  /** Compute the update value. */
  virtual PixelType ComputeUpdate(
      const NeighborhoodType &neighborhood,
      const DiffusionTensorNeighborhoodType
          &tangentialTensorNeighborhood,
      const TensorDerivativeImageRegionType
          &tangentialTensorDerivativeRegion,
      const DeformationVectorComponentNeighborhoodArrayType
          &tangentialDeformationComponentNeighborhoods,
      const NormalVectorNeighborhoodType
              &normalVectorNeighborhood,
      const DiffusionTensorNeighborhoodType
          &normalTensorNeighborhood,
      const TensorDerivativeImageRegionType
          &normalTensorDerivativeRegion,
      const DeformationVectorComponentNeighborhoodArrayType
          &normalDeformationComponentNeighborhoods,
      const SpacingType &spacing,
      void *globalData,
      const FloatOffsetType& = FloatOffsetType(0.0) );

protected:
  AnisotropicDiffusiveRegistrationFunction();
  virtual ~AnisotropicDiffusiveRegistrationFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  // We don't want this function to be made public any longer, but we are
  // allowed to use it internally
  PixelType ComputeUpdate(
      const NeighborhoodType &neighborhood,
      const DiffusionTensorNeighborhoodType
          &tensorNeighborhood,
      const TensorDerivativeImageRegionType
          &tensorDerivativeRegion,
      const DeformationVectorComponentNeighborhoodArrayType
          &deformationComponentNeighborhoods,
      const SpacingType &spacing,
      void *globalData,
      const FloatOffsetType& = FloatOffsetType(0.0) )
    { return Superclass::ComputeUpdate( neighborhood,
                                        tensorNeighborhood,
                                        tensorDerivativeRegion,
                                        deformationComponentNeighborhoods,
                                        spacing,
                                        globalData ); }

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFunction(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /* Whether or not to use the anisotropic regularization */
  bool                                  m_UseAnisotropicRegularization;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveRegistrationFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveRegistrationFunction.txx"
#endif

#endif
