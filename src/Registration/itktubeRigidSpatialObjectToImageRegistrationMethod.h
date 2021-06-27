/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeRigidSpatialObjectToImageRegistrationMethod_h
#define __itktubeRigidSpatialObjectToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkRigid2DTransform.h"

#include "itktubeOptimizedSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
class RigidSpatialObjectToImageRegistrationMethod
  : public OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
{

public:

  typedef RigidSpatialObjectToImageRegistrationMethod             Self;
  typedef OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkTypeMacro( RigidSpatialObjectToImageRegistrationMethod,
                OptimizedSpatialObjectToImageRegistrationMethod );

  itkNewMacro( Self );

  //
  // Typedefs from Superclass
  //

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  // Overrides the superclass' TransformType typedef
  // We must use MatrixOffsetTransformBase since no itk rigid transform is
  //   templated over ImageDimension.
  typedef MatrixOffsetTransformBase<double,
                                    itkGetStaticConstMacro( ImageDimension ),
                                    itkGetStaticConstMacro( ImageDimension )>
  RigidTransformType;
  typedef RigidTransformType TransformType;

  //
  //  Custom Typedefs
  //
  typedef Rigid2DTransform<double>       Rigid2DTransformType;
  typedef VersorRigid3DTransform<double> Rigid3DTransformType;

  typedef AffineTransform<double,
                          itkGetStaticConstMacro( ImageDimension )>
  AffineTransformType;

  typedef typename AffineTransformType::Pointer AffineTransformPointer;

  //
  //  Superclass Methods
  //
  void GenerateData( void ) override;

  //
  // Custom Methods
  //

  /**
   * The function performs the casting.  This function should only appear
   *   once in the class hierarchy.  It is provided so that member
   *   functions that exist only in specific transforms ( e.g., SetIdentity )
   *   can be called without the caller having to do the casting. */
  TransformType * GetTypedTransform( void );

  const TransformType * GetTypedTransform( void ) const;

  /**
   * This function creates a new affine transforms that implements the
   * current registration transform.   Provided to help with transform
   * composition. The transform is initialized with the current results
   * available in the GetTypedTransform() method. The returned transform is
   * not a member variable, and therefore, must be received into a
   * SmartPointer to prevent it from being destroyed by depletion of its
   * reference counting. */
  AffineTransformPointer GetAffineTransform( void ) const;

  /** Initialize the transform parameters from an AffineTransform This method
   * is intended as an alternative to calling SetInitialTransformParameters()
   * and SetInitialTransformFixedParameters(). These later methods require
   * you to have a rigid transform at hand, and this is not always the case,
   * specially when a transform initializer is being used. The method below
   * facilitates to use the AffineTransform returned by the
   * InitialSpatialObjectToImageRegistrationMethod to directly initialize this rigid
   * registration method. The received Affine transform will be approximated
   * to its closest rigid transform by using Polar decomposition. */
  void SetInitialTransformParametersFromAffineTransform(
    const AffineTransformType * affine );

protected:

  RigidSpatialObjectToImageRegistrationMethod( void );
  virtual ~RigidSpatialObjectToImageRegistrationMethod( void );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  RigidSpatialObjectToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

};

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRigidSpatialObjectToImageRegistrationMethod.hxx"
#endif

#endif // __SpatialObjectToImageRegistrationMethod_h
