/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 ( Tue, 10 Jul 2007 ) $
  Version:   $Revision: 0 $

  Copyright ( c ) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScaleSkewAngle2DImageToImageRegistrationMethod_h
#define __itkScaleSkewAngle2DImageToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkTubeScaleSkewAngle2DTransform.h"

#include "itkOptimizedImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
class ScaleSkewAngle2DImageToImageRegistrationMethod
  : public OptimizedImageToImageRegistrationMethod<
    Image< typename TImage::PixelType, 2 > >
{

public:

  typedef ScaleSkewAngle2DImageToImageRegistrationMethod  Self;
  typedef OptimizedImageToImageRegistrationMethod< 
    Image< typename TImage::PixelType, 2 > >              Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkTypeMacro( ScaleSkewAngle2DImageToImageRegistrationMethod,
                OptimizedImageToImageRegistrationMethod );

  itkNewMacro( Self );

  itkStaticConstMacro( ImageDimension, unsigned int, 2 );

  //
  // Typedefs from Superclass
  //

  // Overrides the superclass' TransformType typedef
  typedef ::itk::TubeScaleSkewAngle2DTransform< double >
            ScaleSkewAngle2DTransformType;
  typedef typename ScaleSkewAngle2DTransformType::Pointer
            ScaleSkewAngle2DTransformPointer;
  typedef ScaleSkewAngle2DTransformType
            TransformType;

  typedef AffineTransform<double, 2>
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
   * InitialImageToImageRegistrationMethod
   * to directly initialize this registration method.
   */
  void SetInitialTransformParametersFromAffineTransform(
    const AffineTransformType * transform );

protected:

  ScaleSkewAngle2DImageToImageRegistrationMethod( void );
  virtual ~ScaleSkewAngle2DImageToImageRegistrationMethod( void );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  ScaleSkewAngle2DImageToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScaleSkewAngle2DImageToImageRegistrationMethod.txx"
#endif

#endif // __ImageToImageRegistrationMethod_h
