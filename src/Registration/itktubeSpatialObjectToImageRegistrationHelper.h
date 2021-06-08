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

#ifndef __itkSpatialObjectToImageRegistrationHelper_h
#define __itkSpatialObjectToImageRegistrationHelper_h

#include "itkImage.h"
#include "itkCommand.h"

#include "itkSpatialObjectToImageRegistrationMethod.h"
#include "itkInitialSpatialObjectToImageRegistrationMethod.h"
#include "itkOptimizedSpatialObjectToImageRegistrationMethod.h"
#include "itkRigidSpatialObjectToImageRegistrationMethod.h"
#include "itkAffineSpatialObjectToImageRegistrationMethod.h"
#include "itkScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.h"
#include "itkScaleSkewVersor3DSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

template <class TSpatialObject, class TImage>
class SpatialObjectToImageRegistrationHelper : public Object
{

public:

  typedef SpatialObjectToImageRegistrationHelper Self;
  typedef Object                                 Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  itkTypeMacro( SpatialObjectToImageRegistrationHelper, Object );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef TImage ImageType;

  typedef typename TImage::PixelType PixelType;

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  itkStaticConstMacro( ObjectDimension, unsigned int,
                       TSpatialObject::ObjectDimension );

  //
  // Available Registration Methods
  //
  typedef SpatialObjectToImageRegistrationMethod<TImage>
  RegistrationMethodType;

  typedef InitialSpatialObjectToImageRegistrationMethod<TImage>
  InitialRegistrationMethodType;

  typedef OptimizedSpatialObjectToImageRegistrationMethod<TImage>
  OptimizedRegistrationMethodType;

  typedef RigidSpatialObjectToImageRegistrationMethod<TImage>
  RigidRegistrationMethodType;

  typedef AffineSpatialObjectToImageRegistrationMethod<TImage>
  AffineRegistrationMethodType;

  typedef ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TImage>
  Affine2DRegistrationMethodType;

  typedef ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>
  Affine3DRegistrationMethodType;

  //
  // Typedefs for the parameters of the registration methods
  //
  typedef typename RegistrationMethodType::GroupType
  GroupType;

  typedef typename RegistrationMethodType::ImageMaskObjectType
  ImageMaskObjectType;

  typedef typename RegistrationMethodType::SpatialObjectMaskObjectType
  SpatialObjectMaskObjectType;

  typedef typename OptimizedRegistrationMethodType::MetricMethodEnumType
  MetricMethodEnumType;

  typedef typename OptimizedRegistrationMethodType::InterpolationMethodEnumType
  InterpolationMethodEnumType;

  enum InitialMethodEnumType { INIT_WITH_NONE,
                               INIT_WITH_CURRENT_RESULTS,
                               INIT_WITH_IMAGE_CENTERS,
                               INIT_WITH_CENTERS_OF_MASS,
                               INIT_WITH_SECOND_MOMENTS,
                               INIT_WITH_LANDMARKS };

  enum RegistrationStageEnumType { PRE_STAGE,
                                   LOAD_STAGE,
                                   INIT_STAGE,
                                   RIGID_STAGE,
                                   AFFINE_STAGE };

  enum RegistrationMethodEnumType { NONE,
                                    INITIAL,
                                    RIGID,
                                    AFFINE,
                                    PIPELINE_RIGID,
                                    PIPELINE_AFFINE };

  typedef typename InitialRegistrationMethodType::TransformType
  InitialTransformType;

  typedef std::vector<std::vector<float>  > LandmarkVectorType;

  typedef typename TImage::PointType PointType;

  typedef typename RigidRegistrationMethodType::TransformType
  RigidTransformType;

  typedef typename AffineRegistrationMethodType::TransformType
  AffineTransformType;

  typedef typename Affine2DRegistrationMethodType::TransformType
  Affine2DTransformType;

  typedef typename Affine3DRegistrationMethodType::TransformType
  Affine3DTransformType;

  typedef AffineTransformType MatrixTransformType;

  //
  // Custom Methods
  //

  // **************
  // **************
  //  Specify the fixed and moving images
  // **************
  // **************
  void LoadFixedImage( const std::string & filename );

  itkSetConstObjectMacro( FixedImage, TImage );
  itkGetConstObjectMacro( FixedImage, TImage );

  void LoadMovingSpatialObject( const std::string & filename );

  void SetMovingSpatialObject( const TSpatialObject * MovingSpatialObject );

  itkSetConstObjectMacro( MovingGroupSpatialObject, GroupType );
  itkGetConstObjectMacro( MovingGroupSpatialObject, GroupType );

  // **************
  //  Generic file-save function
  // **************
  void SaveSpatialObject( const std::string & filename,
    const SpatialObject * spatialObject );

  itkSetMacro( RandomNumberSeed, unsigned int );
  itkGetMacro( RandomNumberSeed, unsigned int );

  // **************
  // **************
  //  Specify how the fixed image should be sampled when computing the
  //  metric and what ROI of the moving image is valid
  // **************
  // **************
  itkSetMacro( UseFixedImageMaskObject, bool );
  itkGetConstMacro( UseFixedImageMaskObject, bool );
  itkBooleanMacro( UseFixedImageMaskObject );

  void SetFixedImageMaskObject( const ImageMaskObjectType * mask );

  itkGetConstObjectMacro( FixedImageMaskObject, ImageMaskObjectType );

  itkSetMacro( UseMovingSpatialObjectMaskObject, bool );
  itkGetConstMacro( UseMovingSpatialObjectMaskObject, bool );
  itkBooleanMacro( UseMovingSpatialObjectMaskObject );

  void SetMovingSpatialObjectMaskObject(
    const SpatialObjectMaskObjectType * mask );

  itkGetConstObjectMacro( MovingSpatialObjectMaskObject,
    SpatialObjectMaskObjectType );

  itkSetMacro( UseRegionOfInterest, bool );
  itkGetMacro( UseRegionOfInterest, bool );
  itkSetMacro( RegionOfInterestPoint1, PointType );
  itkGetMacro( RegionOfInterestPoint1, PointType );
  itkSetMacro( RegionOfInterestPoint2, PointType );
  itkGetMacro( RegionOfInterestPoint2, PointType );
  void SetRegionOfInterest( const PointType & point1,
    const PointType & point2 );

  void SetRegionOfInterest( const std::vector<float> & points );

  // **************
  // **************
  //  Update
  // **************
  // **************
  void Initialize( void );

  /** This class provides an Update() method to fit the appearance of a
   * ProcessObject API, but it is not a ProcessObject.  */
  void Update( void );

  // **************
  // **************
  //  Resample
  // **************
  // **************
  const GroupType * ResampleGroupSpatialObject(
    const GroupType * movingGroupSpatialObject = NULL,
    const MatrixTransformType * matrixTransform = NULL,
    double portion = 1.0 );

  // Returns the moving image resampled into the space of the fixed image
  typename GroupType::ConstPointer  GetFinalMovingGroupSpatialObject( void );

  // **************
  // **************
  // Process Control
  // **************
  // **************

  // **************
  // Control which steps of the registration pipeline are applied
  // **************
  itkSetMacro( EnableLoadedRegistration, bool );
  itkGetConstMacro( EnableLoadedRegistration, bool );
  itkBooleanMacro( EnableLoadedRegistration );

  itkSetMacro( EnableInitialRegistration, bool );
  itkGetConstMacro( EnableInitialRegistration, bool );
  itkBooleanMacro( EnableInitialRegistration );

  itkSetMacro( EnableRigidRegistration, bool );
  itkGetConstMacro( EnableRigidRegistration, bool );
  itkBooleanMacro( EnableRigidRegistration );

  itkSetMacro( EnableAffineRegistration, bool );
  itkGetConstMacro( EnableAffineRegistration, bool );
  itkBooleanMacro( EnableAffineRegistration );

  void SetRegistration( RegistrationMethodEnumType reg );
  void SetInterpolation( InterpolationMethodEnumType interp );
  void SetMetric( MetricMethodEnumType metric );

  // **************
  // Specify the optimizer
  // **************
  itkSetMacro( UseEvolutionaryOptimization, bool );
  itkGetMacro( UseEvolutionaryOptimization, bool );
  // **************
  // Specify the expected magnitudes within the transform.  Used to
  //   guide the operating space of the optimizers
  // **************
  itkSetMacro( ExpectedOffsetMagnitude, double );
  itkGetConstMacro( ExpectedOffsetMagnitude, double );

  itkSetMacro( ExpectedRotationMagnitude, double );
  itkGetConstMacro( ExpectedRotationMagnitude, double );

  itkSetMacro( ExpectedScaleMagnitude, double );
  itkGetConstMacro( ExpectedScaleMagnitude, double );

  itkSetMacro( ExpectedSkewMagnitude, double );
  itkGetConstMacro( ExpectedSkewMagnitude, double );

  itkSetMacro( ExpectedDeformationMagnitude, double );
  itkGetConstMacro( ExpectedDeformationMagnitude, double );

  // **************
  //  Return the current product of the registration pipeline
  // **************
  itkGetConstObjectMacro( CurrentMatrixTransform, MatrixTransformType );

  // The image used for registration is updated at certain points in the
  //   registration pipeline for speed and transform composition.
  // Specifically, the SpatialObject is resmpled using the loaded transforms
  // prior to running the initial registration method and the spatialobject
  // is resampled after the affine registration / prior to running bspline
  // registration. The result of these resamplings is available as the
  // CurrentMovingSpatialObject.
  itkGetConstObjectMacro( CurrentMovingSpatialObject, TSpatialObject );
  itkGetConstObjectMacro( LoadedTransformResampledSpatialObject, TSpatialObject );
  itkGetConstObjectMacro( MatrixTransformResampledSpatialObject, TSpatialObject );

  // **************
  //  Not implemented at this time :(
  // **************
  void LoadParameters( const std::string & filename );

  void SaveParameters( const std::string & filename );

  // **************
  //  Final metric value after the pipeline has completed
  // **************
  itkGetMacro( FinalMetricValue, double );

  // **************
  //  Determine if progress messages should be sent to cout
  // **************
  itkSetMacro( ReportProgress, bool );
  itkGetMacro( ReportProgress, bool );
  itkBooleanMacro( ReportProgress );

  itkSetMacro( MinimizeMemory, bool );
  itkGetMacro( MinimizeMemory, bool );
  itkBooleanMacro( MinimizeMemory );

  //
  // Loaded transforms parameters
  //
  void LoadTransform( const std::string & filename, bool invert=false );

  void SaveTransform( const std::string & filename );

  void SaveDisplacementField( const std::string &filename );

  void SetLoadedMatrixTransform( const MatrixTransformType & tfm,
    bool invert=false );

  itkGetConstObjectMacro( LoadedMatrixTransform, MatrixTransformType );

  //
  // Initial Parameters
  //
  itkSetMacro( InitialMethodEnum, InitialMethodEnumType );
  itkGetConstMacro( InitialMethodEnum, InitialMethodEnumType );

  void SetFixedLandmarks( const LandmarkVectorType & fixedLandmarks );
  void SetMovingLandmarks( const LandmarkVectorType & movingLandmarks );

  //
  // Rigid Parameters
  //
  itkSetMacro( RigidSamplingRatio, double );
  itkGetConstMacro( RigidSamplingRatio, double );

  itkSetMacro( RigidTargetError, double );
  itkGetConstMacro( RigidTargetError, double );

  itkSetMacro( RigidMaxIterations, unsigned int );
  itkGetConstMacro( RigidMaxIterations, unsigned int );

  itkSetMacro( RigidMetricMethodEnum, MetricMethodEnumType );
  itkGetConstMacro( RigidMetricMethodEnum, MetricMethodEnumType );

  itkSetMacro( RigidInterpolationMethodEnum, InterpolationMethodEnumType );
  itkGetConstMacro( RigidInterpolationMethodEnum,
    InterpolationMethodEnumType );

  itkGetConstObjectMacro( RigidTransform, RigidTransformType );
  itkGetMacro( RigidMetricValue, double );

  //
  // Affine Parameters
  //
  itkSetMacro( AffineSamplingRatio, double );
  itkGetConstMacro( AffineSamplingRatio, double );

  itkSetMacro( AffineTargetError, double );
  itkGetConstMacro( AffineTargetError, double );

  itkSetMacro( AffineMaxIterations, unsigned int );
  itkGetConstMacro( AffineMaxIterations, unsigned int );

  itkSetMacro( AffineMetricMethodEnum, MetricMethodEnumType );
  itkGetConstMacro( AffineMetricMethodEnum, MetricMethodEnumType );

  itkSetMacro( AffineInterpolationMethodEnum, InterpolationMethodEnumType );
  itkGetConstMacro( AffineInterpolationMethodEnum,
    InterpolationMethodEnumType );

  itkGetConstObjectMacro( AffineTransform, AffineTransformType );
  itkGetMacro( AffineMetricValue, double );

protected:

  SpatialObjectToImageRegistrationHelper( void );
  virtual ~SpatialObjectToImageRegistrationHelper( void );

  void PrintSelfHelper( std::ostream & os, Indent indent,
    const std::string & basename, MetricMethodEnumType metric,
    InterpolationMethodEnumType interpolation ) const;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  template< int tmpImageDimension >
  void AffineRegND() 
   { typename Image< double, tmpImageDimension >::Pointer t; AffineRegND( t ); }

  void AffineRegND( Image< double, 2 > * t );

  void AffineRegND( Image< double, 3 > * t );

  typedef typename InitialRegistrationMethodType::LandmarkPointType
  LandmarkPointType;
  typedef typename InitialRegistrationMethodType::LandmarkPointContainer
  LandmarkPointContainer;

  // Purposely not implemented
  SpatialObjectToImageRegistrationHelper( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  //  Data
  typename TImage::ConstPointer              m_FixedImage;
  typename TSpatialObject::ConstPointer      m_MovingSpatialObject;

  bool                                       m_UseFixedImageMaskObject;
  typename ImageMaskObjectType::ConstPointer m_FixedImageMaskObject;

  bool                                       m_UseMovingSpatialObjectMaskObject;
  typename SpatialObjectMaskObjectType::ConstPointer
                                             m_MovingSpatialObjectMaskObject;

  bool      m_UseRegionOfInterest;
  PointType m_RegionOfInterestPoint1;
  PointType m_RegionOfInterestPoint2;

  unsigned int m_RandomNumberSeed;

  //  Process
  bool m_EnableLoadedRegistration;
  bool m_EnableInitialRegistration;
  bool m_EnableRigidRegistration;
  bool m_EnableAffineRegistration;

  double m_ExpectedOffsetMagnitude;
  double m_ExpectedRotationMagnitude;
  double m_ExpectedScaleMagnitude;
  double m_ExpectedSkewMagnitude;
  double m_ExpectedDeformationMagnitude;

  bool                      m_CompletedInitialization;
  RegistrationStageEnumType m_CompletedStage;
  bool                      m_CompletedResampling;

  typename TSpatialObject::ConstPointer   m_CurrentMovingSpatialObject;
  typename MatrixTransformType::Pointer   m_CurrentMatrixTransform;

  typename TSpatialObject::ConstPointer m_LoadedTransformResampledSpatialObject;
  typename TSpatialObject::ConstPointer m_MatrixTransformResampledSpatialObject;

  double m_FinalMetricValue;

  bool m_ReportProgress;

  bool m_MinimizeMemory;

  //  Optimizer
  bool m_UseEvolutionaryOptimization;

  //  Loaded Tansform
  typename MatrixTransformType::Pointer   m_LoadedMatrixTransform;

  //  Initial Parameters
  InitialMethodEnumType                   m_InitialMethodEnum;
  typename InitialTransformType::Pointer  m_InitialTransform;

  LandmarkPointContainer m_FixedLandmarks;
  LandmarkPointContainer m_MovingLandmarks;

  //  Rigid Parameters
  double       m_RigidSamplingRatio;
  double       m_RigidTargetError;
  unsigned int m_RigidMaxIterations;

  typename RigidTransformType::Pointer    m_RigidTransform;
  MetricMethodEnumType                    m_RigidMetricMethodEnum;
  InterpolationMethodEnumType             m_RigidInterpolationMethodEnum;

  double m_RigidMetricValue;

  //  Affine Parameters
  double       m_AffineSamplingRatio;
  double       m_AffineTargetError;
  unsigned int m_AffineMaxIterations;

  typename AffineTransformType::Pointer   m_AffineTransform;
  MetricMethodEnumType                    m_AffineMetricMethodEnum;
  InterpolationMethodEnumType             m_AffineInterpolationMethodEnum;

  double m_AffineMetricValue;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatialObjectToImageRegistrationHelper.hxx"
#endif

#endif
