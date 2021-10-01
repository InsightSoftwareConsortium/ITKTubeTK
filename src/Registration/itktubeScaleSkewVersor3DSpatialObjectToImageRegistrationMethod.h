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

#ifndef __itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod_h
#define __itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkComposeScaleSkewVersor3DTransform.h"

#include "itktubeOptimizedSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

namespace tube
{

template <class TImage>
class ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod
  : public OptimizedSpatialObjectToImageRegistrationMethod<
    3, Image< typename TImage::PixelType, 3 > >
{

public:

  typedef ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod Self;
  typedef OptimizedSpatialObjectToImageRegistrationMethod<
    3, Image< typename TImage::PixelType, 3 > >           Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkTypeMacro( ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod,
                OptimizedSpatialObjectToImageRegistrationMethod );

  itkNewMacro( Self );

  itkStaticConstMacro( ImageDimension, unsigned int, 3 );
  itkStaticConstMacro( ObjectDimension, unsigned int, 3 );

  //
  // Typedefs from Superclass
  // 
  // Overrides the superclass' TransformType typedef
  typedef ::itk::ComposeScaleSkewVersor3DTransform< double >
            ScaleSkewVersor3DTransformType;
  typedef typename ScaleSkewVersor3DTransformType::Pointer
            ScaleSkewVersor3DTransformPointer;
  typedef ScaleSkewVersor3DTransformType
            TransformType;

  typedef AffineTransform<double, 3>
            AffineTransformType;
  typedef typename AffineTransformType::Pointer
            AffineTransformPointer;

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
   * InitialSpatialObjectToImageRegistrationMethod
   * to directly initialize this registration method.
   */
  void SetInitialTransformParametersFromAffineTransform(
    const AffineTransformType * transform );

protected:

  ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod( void );
  virtual ~ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod( void );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

};

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod.hxx"
#endif

#endif // __SpatialObjectToImageRegistrationMethod_h
