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
#ifndef __itkAnisotropicDiffusiveSparseRegistrationFilter_h
#define __itkAnisotropicDiffusiveSparseRegistrationFilter_h

#include "itkAnisotropicDiffusiveRegistrationFilter.h"

namespace itk
{

/** \class itkAnisotropicDiffusiveSparseRegistrationFilter
 * \brief Registration filter for registrations using anisotropic diffusive
 * regularizers, for example for sliding organ registration.
 *
 * Traditional deformable image registration imposes a uniform
 * smoothness constraint on the deformation field. This
 * is not appropriate when registering images visualizing organs
 * that slide relative to each other, and therefore leads to registration
 * inaccuracies.
 *
 * This filter includes a regularization term based on anisotropic diffusion
 * that accomodates deformation field discontinuities that are expected when
 * considering sliding motion.
 *
 * TODO
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa itkAnisotropicDiffusiveRegistrationFilter
 * \sa itkAnisotropicDiffusiveRegistrationFunction
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT AnisotropicDiffusiveSparseRegistrationFilter
  : public AnisotropicDiffusiveRegistrationFilter< TFixedImage,
                                                   TMovingImage,
                                                   TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveSparseRegistrationFilter      Self;
  typedef AnisotropicDiffusiveRegistrationFilter<
      TFixedImage, TMovingImage, TDeformationField >        Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, AnisotropicDiffusiveRegistrationFilter);

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
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;

  /** Deformation field types. */
  typedef typename Superclass::DeformationVectorType    DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;

  /** Diffusion tensor image types */
  typedef typename Superclass::DiffusionTensorType      DiffusionTensorType;
  typedef typename Superclass::DiffusionTensorImageType
      DiffusionTensorImageType;

  /** Scalar derivative image types */
  typedef typename Superclass::ScalarDerivativeImageType
      ScalarDerivativeImageType;
  typedef typename Superclass::ScalarDerivativeImageArrayType
      ScalarDerivativeImageArrayType;

  /** Tensor derivative matrix image types */
  typedef typename Superclass::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename Superclass::TensorDerivativeImageArrayType
      TensorDerivativeImageArrayType;

  /** Typedefs for the multiplication vectors */
  typedef typename Superclass::DeformationVectorImageArrayType
      DeformationVectorImageArrayType;
  typedef typename Superclass::DeformationVectorImageRegionType
      DeformationVectorImageRegionType;

  /** Normal vector types */
  typedef typename Superclass::NormalVectorComponentType
      NormalVectorComponentType;
  typedef typename Superclass::NormalVectorType         NormalVectorType;
  typedef typename Superclass::NormalVectorImageType    NormalVectorImageType;
  typedef typename Superclass::NormalVectorImagePointer
      NormalVectorImagePointer;
  typedef typename Superclass::NormalVectorImageRegionType
      NormalVectorImageRegionType;

  /** Types for weighting between the anisotropic and diffusive (Gaussian)
    * regularization */
  typedef typename Superclass::WeightType               WeightType;
  typedef typename Superclass::WeightImageType          WeightImageType;
  typedef typename Superclass::WeightImagePointer       WeightImagePointer;
  typedef typename Superclass::WeigthImageRegionType    WeightImageRegionType;

  /** The number of div(Tensor \grad u)v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms() const
    { return 4; }

//  /** Get the normal components of the deformation field. */
//  virtual const DeformationFieldType * GetNormalDeformationComponentImage()
//      const
//    {
//    return this->GetDeformationComponentImage( NORMAL );
//    }

protected:
  AnisotropicDiffusiveSparseRegistrationFilter();
  virtual ~AnisotropicDiffusiveSparseRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Handy for array indexing. */
  enum DivTerm { SMOOTH_TANGENTIAL,
                 SMOOTH_NORMAL,
                 PROP_TANGENTIAL,
                 PROP_NORMAL };

  /** Allocate the deformation component images and their derivative images.
   *  (which may be updated throughout the registration). Reimplement in derived
   *  classes. */
  virtual void InitializeDeformationComponentAndDerivativeImages();

//  /** Allocate and populate the diffusion tensor images.
//   *  Reimplement in derived classes. */
//  virtual void ComputeDiffusionTensorImages();

//  /** Allocate and populate the images of multiplication vectors that the
//   *  div(T \grad(u)) values are multiplied by.  Allocate and populate all or
//   *  some of the multiplication vector images in derived classes.  Otherwise,
//   *  default to e_l, where e_l is the lth canonical unit vector. */
//  virtual void ComputeMultiplicationVectorImages();

//  /** Updates the deformation vector component images on each iteration. */
//  virtual void UpdateDeformationComponentImages();

//  /** If needed, allocates and computes the normal vector and weight images. */
//  virtual void SetupNormalVectorAndWeightImages();

//  /** Compute the normals for the border surface. */
//  void ComputeBorderSurfaceNormals();

//  /** Computes the normal vector image and weighting factors w given the
//   *  surface border polydata. */
//  virtual void ComputeNormalVectorAndWeightImages( bool computeNormals,
//                                                   bool computeWeights );

//  /** Computes the weighting factor w from the distance to the border.  The
//   *  weight should be 1 near the border and 0 away from the border. */
//  virtual WeightType ComputeWeightFromDistance( const WeightType distance )
//      const;

private:
  // Purposely not implemented
  AnisotropicDiffusiveSparseRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveSparseRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveSparseRegistrationFilter.txx"
#endif

#endif
