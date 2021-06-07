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

#ifndef __itktubeTubesToImageRegistrationMethod_h
#define __itktubeTubesToImageRegistrationMethod_h

#include "itktubeTubesToImageMetric.h"

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkSpatialObjectToImageRegistrationMethod.h>
#include <itkTransform.h>
#include <itkTubeSpatialObject.h>
#include <itkVectorContainer.h>

namespace itk
{

namespace tube
{

/** \class TubesToImageRegistrationMethod
 * \brief Register a hierarchy of tubes to an image.
 *
 * This class provides basic registration of a hierarchy of tubes to an image.
 *
 * \tparam TFixedImage Type of the image to register against.
 * \tparam TMovingSpatialObject Type of the moving spatial.  This could be a
 * TubeSpatialObject or perhaps a hierarchy of tubes contained in a
 * GroupSpatialObject.
 * \tparam TMovingTube Type of the tubes. Should be some type of
 * TubeSpatialObject.
 *
 * The basic elements of a registration method are:
 *   - Metric to compare the image and the tubes.
 *   - Transformation used to register the image against the tubes.
 *
 *  \ingroup AffineImageRegistration
 */

template< class TFixedImage, class TMovingSpatialObject >
class TubesToImageRegistration
  : public ImageToSpatialObjectRegistrationMethod< TFixedImage,
                                                   TMovingSpatialObject >
{
public:
  typedef TubesToImageRegistration                  Self;
  typedef SpatialObjectToImageRegistrationMethod< TFixedImage,
    TMovingSpatialObject >
                                                        Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingSpatialObjectType  MovingSpatialObjectType;
  typedef TMovingSpatialObject                          MovingTubeType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubesToImageRegistrationMethod,
                SpatialObjectToImageRegistrationMethod );

  typedef typename Superclass::MetricType  MetricType;

  typedef TubesToImageMetric< MovingSpatialObjectType, FixedImageType >
    DefaultMetricType;

  typedef typename DefaultMetricType::TransformParametersType
    ParametersType;
  typedef typename DefaultMetricType::TransformType
    TransformType;

  /**  Dimension of the images.  */
  enum { ImageDimension = FixedImageType::ImageDimension,
    ParametersDimension = TransformType::ParametersDimension };

  typedef typename Superclass::InterpolatorType   InterpolatorType;
  typedef typename Superclass::OptimizerType      OptimizerType;

  /** Method that initiates the registration. */
  void StartRegistration( void );

  /** Initialize the registration */
  void Initialize( void );

protected:
  TubesToImageRegistrationMethod( void );
  virtual ~TubesToImageRegistrationMethod( void ) {}

private:
  TubesToImageRegistrationMethod( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  bool                                     m_IsInitialized;

  FeatureWeightsType m_FeatureWeights;
}; // End class TubesToImageRegistrationMethod

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubesToImageRegistrationMethod.hxx"
#endif

#endif // End !defined( __itktubeTubesToImageRegistrationMethod_h )
