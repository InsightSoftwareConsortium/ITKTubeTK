/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 ( Tue, 10 Jul 2007 ) $
  Version:   $Revision: 0 $

  Copyright ( c ) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or https://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef itkRigidImageToImageRegistrationMethod_h
#define itkRigidImageToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkRigid2DTransform.h"

#include "itkOptimizedImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
class RigidImageToImageRegistrationMethod : public OptimizedImageToImageRegistrationMethod<TImage>
{

public:
  using Self = RigidImageToImageRegistrationMethod;
  using Superclass = OptimizedImageToImageRegistrationMethod<TImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(RigidImageToImageRegistrationMethod);

  itkNewMacro(Self);

  //
  // Typedefs from Superclass
  //

  static constexpr unsigned int ImageDimension = TImage::ImageDimension;

  // Overrides the superclass' TransformType typedef
  // We must use MatrixOffsetTransformBase since no itk rigid transform is
  //   templated over ImageDimension.
  using RigidTransformType =
    MatrixOffsetTransformBase<double, Self::ImageDimension, Self::ImageDimension>;
  using TransformType = RigidTransformType;

  //
  //  Custom Typedefs
  //
  using Rigid2DTransformType = Rigid2DTransform<double>;
  using Rigid3DTransformType = VersorRigid3DTransform<double>;

  using AffineTransformType = AffineTransform<double, Self::ImageDimension>;

  using AffineTransformPointer = typename AffineTransformType::Pointer;

  //
  //  Superclass Methods
  //
  void
  GenerateData(void) override;

  //
  // Custom Methods
  //

  /**
   * The function performs the casting.  This function should only appear
   *   once in the class hierarchy.  It is provided so that member
   *   functions that exist only in specific transforms ( e.g., SetIdentity )
   *   can be called without the caller having to do the casting. */
  TransformType *
  GetTypedTransform(void);

  const TransformType *
  GetTypedTransform(void) const;

  /**
   * This function creates a new affine transforms that implements the
   * current registration transform.   Provided to help with transform
   * composition. The transform is initialized with the current results
   * available in the GetTypedTransform() method. The returned transform is
   * not a member variable, and therefore, must be received into a
   * SmartPointer to prevent it from being destroyed by depletion of its
   * reference counting. */
  AffineTransformPointer
  GetAffineTransform(void) const;

  /** Initialize the transform parameters from an AffineTransform This method
   * is intended as an alternative to calling SetInitialTransformParameters()
   * and SetInitialTransformFixedParameters(). These later methods require
   * you to have a rigid transform at hand, and this is not always the case,
   * specially when a transform initializer is being used. The method below
   * facilitates to use the AffineTransform returned by the
   * InitialImageToImageRegistrationMethod to directly initialize this rigid
   * registration method. The received Affine transform will be approximated
   * to its closest rigid transform by using Polar decomposition. */
  void
  SetInitialTransformParametersFromAffineTransform(const AffineTransformType * affine);

protected:
  RigidImageToImageRegistrationMethod(void);
  virtual ~RigidImageToImageRegistrationMethod(void);

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // Purposely not implemented
  RigidImageToImageRegistrationMethod(const Self &);
  // Purposely not implemented
  void
  operator=(const Self &);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkRigidImageToImageRegistrationMethod.hxx"
#endif

#endif // __ImageToImageRegistrationMethod_h
