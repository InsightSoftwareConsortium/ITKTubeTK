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

#ifndef __itkOptimizedImageToImageRegistrationMethod_h
#define __itkOptimizedImageToImageRegistrationMethod_h

#include "itkImage.h"

#include "itkImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
class OptimizedImageToImageRegistrationMethod
  : public ImageToImageRegistrationMethod<TImage>
{

public:

  typedef OptimizedImageToImageRegistrationMethod Self;
  typedef ImageToImageRegistrationMethod<TImage>  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  itkTypeMacro( OptimizedImageToImageRegistrationMethod,
                ImageToImageRegistrationMethod );

  itkNewMacro( Self );

  //
  // Typedefs from Superclass
  //
  typedef TImage ImageType;

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
                                 BSPLINE_TRANSFORM };

  enum MetricMethodEnumType { MATTES_MI_METRIC,
                              NORMALIZED_CORRELATION_METRIC,
                              MEAN_SQUARED_ERROR_METRIC };

  enum InterpolationMethodEnumType { NEAREST_NEIGHBOR_INTERPOLATION,
                                     LINEAR_INTERPOLATION,
                                     BSPLINE_INTERPOLATION,
                                     SINC_INTERPOLATION };

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

  itkSetMacro( SampleFromOverlap, bool );
  itkGetConstMacro( SampleFromOverlap, bool );

  itkSetMacro( MinimizeMemory, bool );
  itkGetConstMacro( MinimizeMemory, bool );

  itkSetMacro( MaxIterations, unsigned int );
  itkGetConstMacro( MaxIterations, unsigned int );

  itkSetMacro( UseEvolutionaryOptimization, bool );
  itkGetConstMacro( UseEvolutionaryOptimization, bool );

  itkSetMacro( NumberOfSamples, unsigned int );
  itkGetConstMacro( NumberOfSamples, unsigned int );

  itkSetMacro( UseFixedImageSamplesIntensityThreshold, bool );
  itkGetConstMacro( UseFixedImageSamplesIntensityThreshold, bool );
  void SetFixedImageSamplesIntensityThreshold( PixelType val );

  itkGetConstMacro( FixedImageSamplesIntensityThreshold, PixelType );

  itkSetMacro( TargetError, double );
  itkGetConstMacro( TargetError, double );

  itkSetMacro( RandomNumberSeed, int );
  itkGetConstMacro( RandomNumberSeed, int );

  itkGetConstMacro( TransformMethodEnum, TransformMethodEnumType );

  itkSetMacro( MetricMethodEnum, MetricMethodEnumType );
  itkGetConstMacro( MetricMethodEnum, MetricMethodEnumType );

  itkSetMacro( InterpolationMethodEnum, InterpolationMethodEnumType );
  itkGetConstMacro( InterpolationMethodEnum, InterpolationMethodEnumType );

  itkGetMacro( FinalMetricValue, double );
protected:

  OptimizedImageToImageRegistrationMethod( void );
  virtual ~OptimizedImageToImageRegistrationMethod( void );

  itkSetMacro( FinalMetricValue, double );

  itkSetMacro( TransformMethodEnum, TransformMethodEnumType );

  typedef InterpolateImageFunction<TImage, double> InterpolatorType;
  typedef ImageToImageMetric<TImage, TImage>       MetricType;

  virtual void Optimize( MetricType * metric, InterpolatorType * interpolator );

  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Purposely not implemented
  OptimizedImageToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  TransformParametersType m_InitialTransformParameters;
  TransformParametersType m_InitialTransformFixedParameters;

  TransformParametersType m_LastTransformParameters;

  TransformParametersScalesType m_TransformParametersScales;

  bool m_SampleFromOverlap;

  bool m_MinimizeMemory;

  unsigned int m_MaxIterations;

  bool m_UseEvolutionaryOptimization;

  unsigned int m_NumberOfSamples;

  bool      m_UseFixedImageSamplesIntensityThreshold;
  PixelType m_FixedImageSamplesIntensityThreshold;

  double m_TargetError;

  int m_RandomNumberSeed;

  TransformMethodEnumType m_TransformMethodEnum;

  MetricMethodEnumType m_MetricMethodEnum;

  InterpolationMethodEnumType m_InterpolationMethodEnum;

  double m_FinalMetricValue;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOptimizedImageToImageRegistrationMethod.hxx"
#endif

#endif
