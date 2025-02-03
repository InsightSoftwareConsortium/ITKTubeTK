/*=========================================================================

Library:   TubeTK

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

#ifndef __itkAffineImageToImageRegistrationMethod_h
#define __itkAffineImageToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkAffineTransform.h"

#include "itkOptimizedImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
class AffineImageToImageRegistrationMethod
  : public OptimizedImageToImageRegistrationMethod<TImage>
{

public:

  typedef AffineImageToImageRegistrationMethod            Self;
  typedef OptimizedImageToImageRegistrationMethod<TImage> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkTypeMacro( AffineImageToImageRegistrationMethod,
                OptimizedImageToImageRegistrationMethod );

  itkNewMacro( Self );

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  //
  // Typedefs from Superclass
  //

  // Overrides the superclass' TransformType typedef
  typedef AffineTransform<double, Self:: ImageDimension >
                                                AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;
  typedef AffineTransformType                   TransformType;

  //
  // Superclass Methods
  //

  void GenerateData( void ) override;

  //
  // Custom Methods
  //

  /**
   * The function performs the casting.  This function should only appear
   *   once in the class hierarchy.  It is provided so that member
   *   functions that exist only in specific transforms
   *   ( e.g., SetIdentity )
   *   can be called without the caller having to do the casting.
   */
  TransformType * GetTypedTransform( void );

  const TransformType * GetTypedTransform( void ) const;

  /**
   * This function creates a new affine transforms that implements the
   * current registration transform.   Provided to help with transform
   * composition. The transform is initialized with the current results
   * available in the GetTypedTransform() method. The returned transform is
   * not a member variable, and therefore, must be received into a
   * SmartPointer to prevent it from being destroyed by depletion of its
   * reference counting.
   */
  AffineTransformPointer GetAffineTransform( void ) const;

  /** Initialize the transform parameters from an AffineTransform.
   * This method is intended as an alternative to calling
   * SetInitialTransformParameters() and
   * SetInitialTransformFixedParameters(). The method below facilitates to
   * use the AffineTransform returned by the
   * InitialImageToImageRegistrationMethod
   * to directly initialize this rigid registration method.
   */
  void SetInitialTransformParametersFromAffineTransform(
    const AffineTransformType * affine );

protected:

  AffineImageToImageRegistrationMethod( void );
  virtual ~AffineImageToImageRegistrationMethod( void );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  AffineImageToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAffineImageToImageRegistrationMethod.hxx"
#endif

#endif // __ImageToImageRegistrationMethod_h
