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

#ifndef __itkBSplineImageToImageRegistrationMethod_h
#define __itkBSplineImageToImageRegistrationMethod_h

#include "itkImage.h"
#include "itkBSplineTransform.h"

#include "itkOptimizedImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
class BSplineImageToImageRegistrationMethod
  : public OptimizedImageToImageRegistrationMethod<TImage>
{

public:

  typedef BSplineImageToImageRegistrationMethod           Self;
  typedef OptimizedImageToImageRegistrationMethod<TImage> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkTypeMacro( BSplineImageToImageRegistrationMethod,
                OptimizedImageToImageRegistrationMethod );

  itkNewMacro( Self );

  //
  // Typedefs from Superclass
  //
  typedef TImage ImageType;
  static constexpr unsigned int ImageDimension = TImage::ImageDimension ;

  // Overrides the superclass' TransformType typedef
  typedef BSplineTransform<double, Self:: ImageDimension , 3>
      BSplineTransformType;

  typedef typename BSplineTransformType::Pointer BSplineTransformPointer;

  typedef BSplineTransformType TransformType;

  typedef typename BSplineTransformType::ParametersType ParametersType;

  //
  // Methods from Superclass
  //

  virtual void GenerateData( void ) override;

  //
  // Custom Methods
  //

  /**
   * The function performs the casting.  This function should only appear
   * once in the class hierarchy.  It is provided so that member
   * functions that exist only in specific transforms ( e.g., SetIdentity )
   * can be called without the caller having to do the casting. */
  virtual TransformType * GetTypedTransform( void );

  virtual const TransformType * GetTypedTransform( void ) const;

  itkSetMacro( ExpectedDeformationMagnitude, double );
  itkGetConstMacro( ExpectedDeformationMagnitude, double );

  itkSetClampMacro( NumberOfControlPoints, unsigned int, 3, 2000 );
  itkGetConstMacro( NumberOfControlPoints, unsigned int );

  itkSetClampMacro( NumberOfLevels, unsigned int, 1, 5 );
  itkGetConstMacro( NumberOfLevels, unsigned int );

  BSplineTransformPointer GetBSplineTransform( void ) const;

  void ComputeGridRegion( int numberOfControlPoints,
    typename TransformType::MeshSizeType & regionSize,
    typename TransformType::PhysicalDimensionsType & regionPhysicalDimensions,
    typename TransformType::OriginType & regionOrigin,
    typename TransformType::DirectionType & regionDirection );

  void ResampleControlGrid( int newNumberOfControlPoints,
    ParametersType & newParameters );

  itkSetMacro( GradientOptimizeOnly, bool );
  itkGetMacro( GradientOptimizeOnly, bool );

protected:

  BSplineImageToImageRegistrationMethod( void );
  virtual ~BSplineImageToImageRegistrationMethod( void );

  typedef InterpolateImageFunction<TImage, double> InterpolatorType;
  typedef ImageToImageMetric<TImage, TImage>       MetricType;

  virtual void Optimize( MetricType * metric, InterpolatorType * interpolator )
    override;

  virtual void GradientOptimize( MetricType * metric,
                                 InterpolatorType * interpolator );

  virtual void MultiResolutionOptimize( MetricType * metric,
                                        InterpolatorType * interpolator );

  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  BSplineImageToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  double m_ExpectedDeformationMagnitude;

  unsigned int m_NumberOfControlPoints;

  unsigned int m_NumberOfLevels;

  bool m_GradientOptimizeOnly;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBSplineImageToImageRegistrationMethod.hxx"
#endif

#endif // __ImageToImageRegistrationMethod_h
