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

#ifndef __itkAnisotropicSimilarityLandmarkBasedTransformInitializer_h
#define __itkAnisotropicSimilarityLandmarkBasedTransformInitializer_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkAnisotropicSimilarity3DTransform.h"
#include "itkRigid2DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include <vector>
#include <iostream>

namespace itk
{

/** \brief AnisotropicSimilarityLandmarkBasedTransformInitializer
 * The class computes the transform that aligns the fixed and moving images
 * given a set of landmarks. The class is templated over the Transform type.
 *    The transform computed gives the best fit transform that maps the fixed
 * and moving images in a least squares sense. The indices are taken to
 * correspond, so point 1 in the first set will get mapped close to point
 * 1 in the second set, etc. An equal number of fixed and moving landmarks
 * need to be specified using SetFixedLandmarks() SetMovingLandmarks().
 * Any number of landmarks may be specified.
 * Call InitializeTransform() to initialize the transform.
 *
 * Currently, the  following transforms are supported by the class:
 *    AnisotropicSimilarity3DTransform
 *
 * The class is based in part on Hybrid/vtkLandmarkTransform originally
 * implemented in python by David G. Gobbi.
 *
 * The solution is based on
 * Berthold K. P. Horn ( 1987 ), "Closed-form solution of absolute orientation
 * using unit quaternions,"
 * http://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
 *
 *
 * \ingroup Transforms
 */
template <class TTransform,
          class TFixedImage,
          class TMovingImage>
class AnisotropicSimilarityLandmarkBasedTransformInitializer :
  public Object
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicSimilarityLandmarkBasedTransformInitializer Self;
  typedef Object                                                 Superclass;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( AnisotropicSimilarityLandmarkBasedTransformInitializer,
    Object );

  /** Type of the transform to initialize */
  typedef TTransform                      TransformType;
  typedef typename TransformType::Pointer TransformPointer;

  /** Dimension of parameters. */
  itkStaticConstMacro( InputSpaceDimension, unsigned int,
                      TransformType::InputSpaceDimension );
  itkStaticConstMacro( OutputSpaceDimension, unsigned int,
                      TransformType::OutputSpaceDimension );

  /** Set the transform to be initialized */
  itkSetObjectMacro( Transform, TransformType );

  /** Image Types to use in the initialization of the transform */
  typedef   TFixedImage  FixedImageType;
  typedef   TMovingImage MovingImageType;

  typedef   typename FixedImageType::ConstPointer  FixedImagePointer;
  typedef   typename MovingImageType::ConstPointer MovingImagePointer;

  /** \deprecated
   * Set the fixed image.
   * The method really doesn't do anything. The goal of this class is to compute
   * the optimal transform, for the templated TransformType between the fixed
   * and moving image grid, given a set of landmarks. Nothing is done with the
   * images themselves. The method will therefore be deprecated and removed */
  void SetFixedImage( const FixedImageType * image )
    {
    this->m_FixedImage = image;
    itkLegacyBodyMacro( SetFixedImage, 2.2 );
    }

  /** \deprecated
   * Set the moving image.
   * The method really doesn't do anything. The goal of this class is to compute
   * the optimal transform, for the templated TransformType between the fixed
   * and moving image grid, given a set of landmarks. Nothing is done with the
   * images themselves. The method will therefore be deprecated and removed. */
  void SetMovingImage( const MovingImageType * image )
    {
    this->m_MovingImage = image;
    itkLegacyBodyMacro( SetMovingImage, 2.2 );
    }

  /** Determine the image dimension. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    FixedImageType::ImageDimension );

  /** Convenience typedefs */
  typedef typename TransformType::InputPointType
    InputPointType;
  typedef typename TransformType::OutputVectorType
    OutputVectorType;
  typedef Point<double, itkGetStaticConstMacro( ImageDimension )>
    LandmarkPointType;
  typedef std::vector<LandmarkPointType>
    LandmarkPointContainer;
  typedef typename  LandmarkPointContainer::const_iterator
    PointsContainerConstIterator;
  typedef typename TransformType::ParametersType
    ParametersType;
  typedef typename ParametersType::ValueType
    ParameterValueType;

  /** Set the Fixed landmark point containers */
  void SetFixedLandmarks( const LandmarkPointContainer & fixedLandmarks )
    {
    this->m_FixedLandmarks = fixedLandmarks;
    }

  /** Set the Moving landmark point containers */
  void SetMovingLandmarks( const LandmarkPointContainer & movingLandmarks )
    {
    this->m_MovingLandmarks = movingLandmarks;
    }

  /**  Supported Transform typedefs */
  typedef AnisotropicSimilarity3DTransform<ParameterValueType>
  AnisotropicSimilarity3DTransformType;
  typedef VersorRigid3DTransform<ParameterValueType> VersorRigid3DTransformType;
  typedef Rigid2DTransform<ParameterValueType> Rigid2DTransformType;

  /** Initialize the transform from the landmarks */
  virtual void InitializeTransform();

protected:
  AnisotropicSimilarityLandmarkBasedTransformInitializer();
  ~AnisotropicSimilarityLandmarkBasedTransformInitializer()
    {
    }

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  // Supported Transform types
  typedef enum
    {
    AnisotropicSimilarity3Dtransform = 1,
    VersorRigid3Dtransform,
    Rigid2Dtransform,
    Else
    } InputTransformType;
private:
  // purposely not implemented
  AnisotropicSimilarityLandmarkBasedTransformInitializer( const Self & );

  // purposely not implemented
  void operator=( const Self & );

  FixedImagePointer  m_FixedImage;
  MovingImagePointer m_MovingImage;

  LandmarkPointContainer m_FixedLandmarks;
  LandmarkPointContainer m_MovingLandmarks;

  TransformPointer m_Transform;

}; // class LandmarkBasedTransformInitializer

}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnisotropicSimilarityLandmarkBasedTransformInitializer.hxx"
#endif

#endif /* __itkAnisotropicSimilarityLandmarkBasedTransformInitializer_h */
