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

#ifndef __itktubeOptimizedSpatialObjectToImageRegistrationMethod_h
#define __itktubeOptimizedSpatialObjectToImageRegistrationMethod_h

#include "itkImage.h"

#include "itktubeSpatialObjectToImageRegistrationMethod.h"
#include "itktubePointBasedSpatialObjectToImageMetric.h"

namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
class OptimizedSpatialObjectToImageRegistrationMethod
  : public SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
{

public:

  typedef OptimizedSpatialObjectToImageRegistrationMethod
                                                  Self;
  typedef SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
                                                  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  itkTypeMacro( OptimizedSpatialObjectToImageRegistrationMethod,
                SpatialObjectToImageRegistrationMethod );

  itkNewMacro( Self );

  //
  // Typedefs from Superclass
  //
  typedef TImage                         ImageType;
  typedef SpatialObject<ObjectDimension> SpatialObjectType;

  typedef typename ImageType::PixelType PixelType;

  typedef typename Superclass::TransformType TransformType;

  typedef typename TransformType::ParametersType TransformParametersType;

  typedef typename TransformType::ParametersType TransformParametersScalesType;

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  //
  // Custom Typedefs
  //
  enum TransformMethodEnumType { RIGID_TRANSFORM,
                                 AFFINE_TRANSFORM,
                                 PERSPECITVE_RIGID_TRANSFORM,
                                 PERSPECITVE_AFFINE_TRANSFORM };

  enum MetricMethodEnumType { IMAGE_INTENSITY_METRIC };

  //
  // Methods from Superclass
  //
  virtual void GenerateData( void ) override;

  //
  // Custom Methods
  //
  itkSetMacro( InitialTransformParameters, TransformParametersType );
  itkGetConstMacro( InitialTransformParameters, TransformParametersType );

  itkSetMacro( InitialTransformFixedParameters, TransformParametersType );
  itkGetConstMacro( InitialTransformFixedParameters, TransformParametersType );

  itkSetMacro( LastTransformParameters, TransformParametersType );
  itkGetConstMacro( LastTransformParameters, TransformParametersType );

  itkSetMacro( TransformParametersScales, TransformParametersScalesType );
  itkGetConstMacro( TransformParametersScales, TransformParametersScalesType );

  itkSetMacro( MaxIterations, unsigned int );
  itkGetConstMacro( MaxIterations, unsigned int );

  itkSetMacro( UseEvolutionaryOptimization, bool );
  itkGetConstMacro( UseEvolutionaryOptimization, bool );

  itkSetMacro( NumberOfSamples, unsigned int );
  itkGetConstMacro( NumberOfSamples, unsigned int );

  itkSetMacro( SamplingRatio, double );
  itkGetConstMacro( SamplingRatio, double );

  itkSetMacro( TargetError, double );
  itkGetConstMacro( TargetError, double );

  itkSetMacro( RandomNumberSeed, int );
  itkGetConstMacro( RandomNumberSeed, int );

  itkGetConstMacro( TransformMethodEnum, TransformMethodEnumType );

  itkSetMacro( MetricMethodEnum, MetricMethodEnumType );
  itkGetConstMacro( MetricMethodEnum, MetricMethodEnumType );

  itkGetMacro( FinalMetricValue, double );

protected:

  OptimizedSpatialObjectToImageRegistrationMethod( void );
  virtual ~OptimizedSpatialObjectToImageRegistrationMethod( void );

  itkSetMacro( FinalMetricValue, double );

  itkSetMacro( TransformMethodEnum, TransformMethodEnumType );

  typedef SpatialObjectToImageMetric<ObjectDimension, TImage>       MetricType;

  virtual void Optimize( MetricType * metric );

  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  OptimizedSpatialObjectToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  TransformParametersType m_InitialTransformParameters;
  TransformParametersType m_InitialTransformFixedParameters;

  TransformParametersType m_LastTransformParameters;

  TransformParametersScalesType m_TransformParametersScales;

  unsigned int m_MaxIterations;

  bool m_UseEvolutionaryOptimization;

  unsigned int m_NumberOfSamples;
  double       m_SamplingRatio;

  double m_TargetError;

  int m_RandomNumberSeed;

  TransformMethodEnumType m_TransformMethodEnum;

  MetricMethodEnumType m_MetricMethodEnum;

  double m_FinalMetricValue;
};

} // tube

} // itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeOptimizedSpatialObjectToImageRegistrationMethod.hxx"
#endif

#endif
