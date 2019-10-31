/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeImageToTubeRigidRegistration_h
#define __itktubeImageToTubeRigidRegistration_h

#include "itktubeImageToTubeRigidMetric.h"
#include "itktubeTubeExponentialResolutionWeightFunction.h"

#include <itkEuler3DTransform.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageToSpatialObjectRegistrationMethod.h>
#include <itkTransform.h>
#include <itkTubeSpatialObject.h>
#include <itkVectorContainer.h>

namespace itk
{

namespace tube
{

/** \class ImageToTubeRigidRegistration
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
 * \todo Below is out of date and needs to be updated.
 *
 *   - Optimization method used to search for the best transformation.
 *
 * Registration is not limited to Images, and for this reason
 * this class is templated over the type of the FixedImage object,
 * the TMovingTube object and the transformation. This types are obtained
 * from the Metric type, to reduce the number of redundant
 * template parameters
 *
 *  \ingroup AffineImageRegistration
 */

template< class TFixedImage, class TMovingSpatialObject, class TMovingTube >
class ImageToTubeRigidRegistration
  : public ImageToSpatialObjectRegistrationMethod< TFixedImage,
                                                   TMovingSpatialObject >
{
public:
  typedef ImageToTubeRigidRegistration                  Self;
  typedef ImageToSpatialObjectRegistrationMethod< TFixedImage,
    TMovingSpatialObject >
                                                        Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingSpatialObjectType  MovingSpatialObjectType;
  typedef TMovingTube                                   MovingTubeType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ImageToTubeRigidRegistration,
                ImageToSpatialObjectRegistrationMethod );

  typedef typename Superclass::MetricType  MetricType;

  typedef ImageToTubeRigidMetric< FixedImageType,
    MovingSpatialObjectType,
    MovingTubeType >
    DefaultMetricType;

  typedef typename DefaultMetricType::TransformParametersType
    ParametersType;
  typedef typename DefaultMetricType::TransformType
    TransformType;
  typedef ParametersType
    FeatureWeightsType;

  /**  Dimension of the images.  */
  enum { ImageDimension = FixedImageType::ImageDimension,
    ParametersDimension = TransformType::ParametersDimension };

  typedef typename Superclass::InterpolatorType   InterpolatorType;
  typedef typename Superclass::OptimizerType      OptimizerType;

  /** Method that initiates the registration. */
  void StartRegistration( void );

  /** Initialize the registration */
  void Initialize( void );

  /** Set/Get the scalar weights associated with every point in the tube.
   * The index of the point weights should correspond to "standard tube tree
   * interation". */
  void SetFeatureWeights( FeatureWeightsType & featureWeights );
  itkGetConstReferenceMacro( FeatureWeights, FeatureWeightsType )

protected:
  ImageToTubeRigidRegistration( void );
  virtual ~ImageToTubeRigidRegistration( void ) {}

private:
  ImageToTubeRigidRegistration( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  bool                                     m_IsInitialized;

  FeatureWeightsType m_FeatureWeights;
}; // End class ImageToTubeRigidRegistration

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeImageToTubeRigidRegistration.hxx"
#endif

#endif // End !defined( __itktubeImageToTubeRigidRegistration_h )
