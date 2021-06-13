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

#ifndef __itkSpatialObjectToImageRegistrationHelper_txx
#define __itkSpatialObjectToImageRegistrationHelper_txx

#include "itkSpatialObjectToImageRegistrationHelper.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkResampleImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkVector.h"
#include "itkAffineTransform.h"

namespace itk
{

template <class TSpatialObject, class TImage>
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SpatialObjectToImageRegistrationHelper()
{
  // Data
  m_FixedImage = NULL;
  m_MovingSpatialObject = NULL;

  // Masks
  m_UseFixedImageMaskObject = false;
  m_FixedImageMaskObject = NULL;
  m_UseMovingSpatialObjectMaskObject = false;
  m_MovingSpatialObjectMaskObject = NULL;

  m_UseFixedImageRegionOfInterest = false;
  m_FixedImageRegionOfInterestPoint1.Fill(0);
  m_FixedImageRegionOfInterestPoint2.Fill(0);

  m_RandomNumberSeed = 0;

  // Process
  m_EnableLoadedRegistration = true;
  m_EnableInitialRegistration = true;
  m_EnableRigidRegistration = true;
  m_EnableAffineRegistration = false;

  // Expected transform magnitude
  m_ExpectedOffsetMagnitude = 5;
  m_ExpectedRotationMagnitude = 0.01;
  m_ExpectedScaleMagnitude = 0.01;
  m_ExpectedSkewMagnitude = 0.0001;

  // Current state of the registration pipeline
  m_CompletedInitialization = false;
  m_CompletedStage = PRE_STAGE;
  m_CompletedResampling = false;

  m_CurrentMovingSpatialObject = NULL;
  m_CurrentMatrixTransform = NULL;

  m_LoadedTransformResampledImage = NULL;
  m_MatrixTransformResampledImage = NULL;

  // Results
  m_FinalMetricValue = 0.0;

  // Progress
  m_ReportProgress = false;

  // Optimizer
  m_UseEvolutionaryOptimization = true ;

  // Loaded
  m_LoadedMatrixTransform = NULL;

  // Initial
  m_InitialMethodEnum = INIT_WITH_NONE;
  m_InitialTransform = NULL;

  // Rigid
  m_RigidSamplingRatio = 0.05;
  m_RigidTargetError = 0.0001;
  m_RigidMaxIterations = 1000;
  m_RigidTransform = NULL;
  m_RigidMetricMethodEnum =
    OptimizedRegistrationMethodType::TUBE_TO_IMAGE_METRIC;
  m_RigidInterpolationMethodEnum =
    OptimizedRegistrationMethodType::LINEAR_INTERPOLATION;
  m_RigidMetricValue = 0.0;

  // Affine
  m_AffineSamplingRatio = 0.1;
  m_AffineTargetError = 0.0001;
  m_AffineMaxIterations = 500;
  m_AffineTransform = NULL;
  m_AffineMetricMethodEnum =
    OptimizedRegistrationMethodType::TUBE_TO_IMAGE_METRIC;
  m_AffineInterpolationMethodEnum =
    OptimizedRegistrationMethodType::LINEAR_INTERPOLATION;
  m_AffineMetricValue = 0.0;

}

template <class TSpatialObject, class TImage>
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::~SpatialObjectToImageRegistrationHelper()
{
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetMovingSpatialObject( const SpatialObjectType * movingSpatialObject )
{
  if( movingSpatialObject->GetTypeName() != "GroupSpatialObject" )
    {
    this->m_MovingGroupSpatialObject = GroupSpatialObjectType::New();
    this->m_MovingGroupSpatialObject->AddChild(movingSpatialObject);
    }
  else
    {
    this->m_MovingGroupSpatialObject = movingSpatialObject;
    }

  this->ProcessObject::SetNthInput(1, const_cast<GroupSpatialObjectType *>(
      m_MovingGroupSpatialObject ) );

  this->Modified();

  m_CompletedStage = PRE_STAGE;

  m_CompletedInitialization = false;
  m_CompletedResampling = false;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetMovingGroupSpatialObject(
  const GroupSpatialObjectType * movingGroupSpatialObject )
{
  this->m_MovingGroupSpatialObject = movingGroupSpatialObject;

  this->ProcessObject::SetNthInput(1, const_cast<GroupSpatialObjectType *>(
      m_MovingGroupSpatialObject ) );

  this->Modified();

  m_CompletedStage = PRE_STAGE;

  m_CompletedInitialization = false;
  m_CompletedResampling = false;
}


template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetFixedImageMaskObject( const MaskObjectType * maskObject )
{
  if( this->m_FixedImageMaskObject.GetPointer() != maskObject )
    {
    this->m_FixedImageMaskObject = maskObject;

    this->Modified();

    if( maskObject != NULL )
      {
      m_UseFixedImageMaskObject = true;
      }
    else
      {
      m_UseFixedImageMaskObject = false;
      }
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetMovingSpatialObjectMaskObject( const MaskObjectType * maskObject )
{
  if( this->m_MovingSpatialObjectMaskObject.GetPointer() != maskObject )
    {
    this->m_MovingSpatialObjectMaskObject = maskObject;

    this->Modified();

    if( maskObject != NULL )
      {
      m_UseMovingSpatialObjectMaskObject = true;
      }
    else
      {
      m_UseMovingSpatialObjectMaskObject = false;
      }
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetRegistration( RegistrationMethodEnumType reg )
{
  switch(reg)
    {
    default:
    case NONE:
      {
      this->SetEnableInitialRegistration(false);
      this->SetEnableRigidRegistration(false);
      this->SetEnableAffineRegistration(false);
      break;
      }
    case INITIAL:
      {
      this->SetEnableInitialRegistration(true);
      this->SetEnableRigidRegistration(false);
      this->SetEnableAffineRegistration(false);
      break;
      }
    case RIGID:
      {
      this->SetEnableInitialRegistration(false);
      this->SetEnableRigidRegistration(true);
      this->SetEnableAffineRegistration(false);
      break;
      }
    case AFFINE:
      {
      this->SetEnableInitialRegistration(false);
      this->SetEnableRigidRegistration(false);
      this->SetEnableAffineRegistration(true);
      break;
      }
    case PIPELINE_RIGID:
      {
      this->SetEnableInitialRegistration(true);
      this->SetEnableRigidRegistration(true);
      this->SetEnableAffineRegistration(false);
      break;
      }
    case PIPELINE_AFFINE:
      {
      this->SetEnableInitialRegistration(true);
      this->SetEnableRigidRegistration(true);
      this->SetEnableAffineRegistration(true);
      break;
      }
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetInterpolation( InterpolationMethodEnumType interp )
{
  this->SetRigidInterpolationMethodEnum( interp );
  this->SetAffineInterpolationMethodEnum( interp );
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetMetric( MetricMethodEnumType metric )
{
  this->SetRigidMetricMethodEnum( metric );
  this->SetAffineMetricMethodEnum( metric );
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::Initialize( void )
{
  // m_LoadedTransform = 0;  Not Initialized - since it is a user parameter
  m_InitialTransform = 0;
  m_RigidTransform = 0;
  m_AffineTransform = 0;

  m_CompletedStage = PRE_STAGE;

  m_CompletedInitialization = true;
  m_CompletedResampling = false;

  m_CurrentMatrixTransform = 0;

  m_FinalMetricValue = 0;
  m_RigidMetricValue = 0;
  m_AffineMetricValue = 0;

  if( m_InitialMethodEnum == INIT_WITH_CURRENT_RESULTS )
    {
    m_CurrentMovingGroupSpatialObject = GetFinalMovingGroupSpatialObject();
    }
  else
    {
    m_CurrentMovingGroupSpatialObject = m_MovingSpatialObject;
    }

  // Eventually these should only be reset if necessary - that is, if the
  //   only difference is enable BSpline registration, it shouldn't be
  //   necessary to re-run the entire registration pipeline
  m_LoadedTransformResampledImage = 0;
  m_MatrixTransformResampledImage = 0;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::AffineRegND( Image< double, 2 > * itkNotUsed( tmpImage ) )
{
  if( this->GetReportProgress() )
    {
    std::cout << "*** AFFINE REGISTRATION ***" << std::endl;
    }

  std::cout << "Start affine" << std::endl;
  unsigned long fixedImageNumPixels = m_FixedImage->GetLargestPossibleRegion()
    .GetNumberOfPixels();

  typename Affine2DRegistrationMethodType::Pointer regAff
    = Affine2DRegistrationMethodType::New();
  regAff->SetRandomNumberSeed( m_RandomNumberSeed );
  regAff->SetReportProgress( m_ReportProgress );
  regAff->SetMovingImage( m_CurrentMovingImage );
  regAff->SetFixedImage( m_FixedImage );
  regAff->SetNumberOfSamples( (unsigned int)(m_AffineSamplingRatio
    * fixedImageNumPixels) );
  if( m_UseFixedImageRegionOfInterest )
    {
    regAff->SetFixedImageRegionOfInterest( m_FixedImageRegionOfInterestPoint1,
      m_FixedImageRegionOfInterestPoint2 );
    }
  regAff->SetSampleFromOverlap( m_SampleFromOverlap );
  regAff->SetMaxIterations( m_AffineMaxIterations );
  regAff->SetTargetError( m_AffineTargetError );
  if( m_EnableRigidRegistration || !m_UseEvolutionaryOptimization )
    {
    regAff->SetUseEvolutionaryOptimization( false );
    }
  regAff->SetTargetError( m_AffineTargetError );
  if( m_UseFixedImageMaskObject )
    {
    if( m_FixedImageMaskObject.IsNotNull() )
      {
      regAff->SetFixedImageMaskObject( m_FixedImageMaskObject );
      }
    }
  if( m_UseMovingSpatialObjectMaskObject )
    {
    if( m_MovingSpatialObjectMaskObject.IsNotNull() )
      {
      regAff->SetMovingSpatialObjectMaskObject( m_MovingSpatialObjectMaskObject );
      }
    }
  if( m_SampleIntensityPortion > 0 )
    {
    typedef MinimumMaximumImageCalculator<ImageType> MinMaxCalcType;
    typename MinMaxCalcType::Pointer calc = MinMaxCalcType::New();
    calc->SetImage( m_FixedImage );
    calc->Compute();
    PixelType fixedImageMax = calc->GetMaximum();
    PixelType fixedImageMin = calc->GetMinimum();

    regAff->SetFixedImageSamplesIntensityThreshold(
      static_cast<PixelType>( ( m_SampleIntensityPortion
        * (fixedImageMax - fixedImageMin) ) + fixedImageMin ) );
    }
  regAff->SetMetricMethodEnum( m_AffineMetricMethodEnum );
  regAff->SetInterpolationMethodEnum( m_AffineInterpolationMethodEnum );
  typename AffineTransformType::ParametersType scales;
  scales.set_size( 7 );
  unsigned int scaleNum = 0;
  scales[scaleNum++] = 1.0 / ( m_ExpectedRotationMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedOffsetMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedOffsetMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedScaleMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedScaleMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedSkewMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedSkewMagnitude );
  regAff->SetTransformParametersScales( scales );

  if( m_CurrentMatrixTransform.IsNotNull() )
    {
    regAff->GetTypedTransform()->SetCenter(
      m_CurrentMatrixTransform->GetCenter() );
    //regAff->GetTypedTransform()->SetMatrix(
      //m_CurrentMatrixTransform->GetMatrix() );
    regAff->GetTypedTransform()->SetOffset(
      m_CurrentMatrixTransform->GetOffset() );
    regAff->SetInitialTransformParameters(
      regAff->GetTypedTransform()->GetParameters() );
    regAff->SetInitialTransformFixedParameters(
      regAff->GetTypedTransform()->GetFixedParameters() );
    }

  regAff->Update();

  m_AffineTransform = regAff->GetAffineTransform();
  m_CurrentMatrixTransform = m_AffineTransform;

  m_FinalMetricValue = regAff->GetFinalMetricValue();
  m_AffineMetricValue = m_FinalMetricValue;

  m_CompletedStage = AFFINE_STAGE;
  m_CompletedResampling = false;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::AffineRegND( Image< double, 3 > * itkNotUsed( tmpImage ) )
{
  if( this->GetReportProgress() )
    {
    std::cout << "*** AFFINE REGISTRATION ***" << std::endl;
    }

  unsigned long fixedImageNumPixels = m_FixedImage->GetLargestPossibleRegion()
    .GetNumberOfPixels();

  typename Affine3DRegistrationMethodType::Pointer regAff =
    Affine3DRegistrationMethodType::New();
  regAff->SetRandomNumberSeed( m_RandomNumberSeed );
  regAff->SetReportProgress( m_ReportProgress );
  regAff->SetMovingImage( m_CurrentMovingImage );
  regAff->SetFixedImage( m_FixedImage );
  regAff->SetNumberOfSamples( (unsigned int)(m_AffineSamplingRatio
    * fixedImageNumPixels) );
  if( m_UseFixedImageRegionOfInterest )
    {
    regAff->SetFixedImageRegionOfInterest( m_FixedImageRegionOfInterestPoint1,
      m_FixedImageRegionOfInterestPoint2 );
    }
  regAff->SetSampleFromOverlap( m_SampleFromOverlap );
  regAff->SetMaxIterations( m_AffineMaxIterations );
  regAff->SetTargetError( m_AffineTargetError );
  if( m_EnableRigidRegistration || !m_UseEvolutionaryOptimization )
    {
    regAff->SetUseEvolutionaryOptimization( false );
    }
  regAff->SetTargetError( m_AffineTargetError );
  if( m_UseFixedImageMaskObject )
    {
    if( m_FixedImageMaskObject.IsNotNull() )
      {
      regAff->SetFixedImageMaskObject( m_FixedImageMaskObject );
      }
    }
  if( m_UseMovingSpatialObjectMaskObject )
    {
    if( m_MovingSpatialObjectMaskObject.IsNotNull() )
      {
      regAff->SetMovingSpatialObjectMaskObject( m_MovingSpatialObjectMaskObject );
      }
    }
  if( m_SampleIntensityPortion > 0 )
    {
    typedef MinimumMaximumImageCalculator<ImageType> MinMaxCalcType;
    typename MinMaxCalcType::Pointer calc = MinMaxCalcType::New();
    calc->SetImage( m_FixedImage );
    calc->Compute();
    PixelType fixedImageMax = calc->GetMaximum();
    PixelType fixedImageMin = calc->GetMinimum();

    regAff->SetFixedImageSamplesIntensityThreshold(
      static_cast<PixelType>( ( m_SampleIntensityPortion
        * (fixedImageMax - fixedImageMin) ) + fixedImageMin ) );
    }
  regAff->SetMetricMethodEnum( m_AffineMetricMethodEnum );
  regAff->SetInterpolationMethodEnum( m_AffineInterpolationMethodEnum );
  typename AffineTransformType::ParametersType scales;

  scales.set_size( 12 );
  unsigned int scaleNum = 0;
  scales[scaleNum++] = 1.0 / ( m_ExpectedRotationMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedRotationMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedRotationMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedOffsetMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedOffsetMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedOffsetMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedScaleMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedScaleMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedScaleMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedSkewMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedSkewMagnitude );
  scales[scaleNum++] = 1.0 / ( m_ExpectedSkewMagnitude );
  regAff->SetTransformParametersScales( scales );

  if( m_CurrentMatrixTransform.IsNotNull() )
    {
    regAff->GetTypedTransform()->SetCenter(
      m_CurrentMatrixTransform->GetCenter() );
    regAff->GetTypedTransform()->SetMatrix(
      m_CurrentMatrixTransform->GetMatrix() );
    regAff->GetTypedTransform()->SetOffset(
      m_CurrentMatrixTransform->GetOffset() );
    regAff->SetInitialTransformParameters(
      regAff->GetTypedTransform()->GetParameters() );
    regAff->SetInitialTransformFixedParameters(
      regAff->GetTypedTransform()->GetFixedParameters() );
    }

  regAff->Update();

  m_AffineTransform = regAff->GetAffineTransform();
  m_CurrentMatrixTransform = m_AffineTransform;

  m_FinalMetricValue = regAff->GetFinalMetricValue();
  m_AffineMetricValue = m_FinalMetricValue;

  m_CompletedStage = AFFINE_STAGE;
  m_CompletedResampling = false;
}

/** This class provides an Update() method to fit the appearance of a
 * ProcessObject API, but it is not a ProcessObject.  */
template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::Update( void )
{
  if( !(this->m_CompletedInitialization) )
    {
    this->Initialize();
    }

  if( m_EnableLoadedRegistration
      && m_LoadedMatrixTransform.IsNotNull() )
    {
    if( m_LoadedTransformResampledImage.IsNotNull() )
      {
      m_CurrentMovingImage = m_LoadedTransformResampledImage;
      if( this->GetReportProgress() )
        {
        std::cout << "*** Using existing loaded transform ***" << std::endl;
        }
      }

    m_MatrixTransformResampledImage = 0;

    m_CompletedStage = LOAD_STAGE;
    m_CompletedResampling = true;

    m_CurrentMatrixTransform = 0;
    }

  if( this->GetReportProgress() )
    {
    std::cout << "*** INITIAL REGISTRATION ***" << std::endl;
    }

  typename InitialRegistrationMethodType::Pointer regInit =
    InitialRegistrationMethodType::New();
  regInit->SetReportProgress( m_ReportProgress );
  regInit->SetMovingImage( m_CurrentMovingImage );
  regInit->SetFixedImage( m_FixedImage );
  if( m_UseFixedImageMaskObject )
    {
    if( m_FixedImageMaskObject.IsNotNull() )
      {
      regInit->SetFixedImageMaskObject( m_FixedImageMaskObject );
      }
    }
  if( m_UseMovingSpatialObjectMaskObject )
    {
    if( m_MovingSpatialObjectMaskObject.IsNotNull() )
      {
      regInit->SetMovingSpatialObjectMaskObject( m_MovingSpatialObjectMaskObject );
      }
    }
  if( m_EnableInitialRegistration )
    {
    switch( m_InitialMethodEnum )
      {
      case INIT_WITH_NONE:
        regInit->SetComputeCenterOfRotationOnly( true );
        break;
      case INIT_WITH_IMAGE_CENTERS:
        regInit->SetNumberOfMoments( 0 );
        break;
      case INIT_WITH_CENTERS_OF_MASS:
        regInit->SetNumberOfMoments( 1 );
        break;
      case INIT_WITH_SECOND_MOMENTS:
        regInit->SetNumberOfMoments( 2 );
        break;
      case INIT_WITH_LANDMARKS:
        regInit->SetUseLandmarks( true );
        regInit->SetFixedLandmarks( m_FixedLandmarks );
        regInit->SetMovingLandmarks( m_MovingLandmarks );
        break;
      case INIT_WITH_CURRENT_RESULTS:
      default:
        break;
      }
    regInit->Update();
    m_InitialTransform = regInit->GetAffineTransform();
    }
  else
    {
    if( m_EnableLoadedRegistration
      && m_LoadedMatrixTransform.IsNotNull() )
      {
      m_InitialTransform = m_LoadedMatrixTransform;
      }
    else
      {
      regInit->SetComputeCenterOfRotationOnly( true );
      regInit->Update();
      m_InitialTransform = regInit->GetAffineTransform();
      }
    }

  m_CurrentMatrixTransform = m_InitialTransform;

  m_CompletedStage = INIT_STAGE;
  m_CompletedResampling = false;

  typename TImage::SizeType fixedImageSize;
  fixedImageSize = m_FixedImage->GetLargestPossibleRegion().GetSize();
  unsigned long fixedImageNumPixels = m_FixedImage->GetLargestPossibleRegion()
    .GetNumberOfPixels();

  if( m_EnableRigidRegistration )
    {
    if( this->GetReportProgress() )
      {
      std::cout << "*** RIGID REGISTRATION ***" << std::endl;
      }

    typename RigidRegistrationMethodType::Pointer regRigid;
    regRigid = RigidRegistrationMethodType::New();
    regRigid->SetRandomNumberSeed( m_RandomNumberSeed );
    if( !m_UseEvolutionaryOptimization )
      {
      regRigid->SetUseEvolutionaryOptimization( false );
      }
    regRigid->SetReportProgress( m_ReportProgress );
    regRigid->SetMovingImage( m_CurrentMovingImage );
    regRigid->SetFixedImage( m_FixedImage );
    regRigid->SetNumberOfSamples( (unsigned int)( m_RigidSamplingRatio
                                                  * fixedImageNumPixels ) );
    regRigid->SetSampleFromOverlap( m_SampleFromOverlap );
    regRigid->SetMaxIterations( m_RigidMaxIterations );
    regRigid->SetTargetError( m_RigidTargetError );
    if( m_UseFixedImageMaskObject )
      {
      if( m_FixedImageMaskObject.IsNotNull() )
        {
        regRigid->SetFixedImageMaskObject( m_FixedImageMaskObject );
        }
      }
    if( m_UseMovingSpatialObjectMaskObject )
      {
      if( m_MovingSpatialObjectMaskObject.IsNotNull() )
        {
        regRigid->SetMovingSpatialObjectMaskObject( m_MovingSpatialObjectMaskObject );
        }
      }
    if( m_SampleIntensityPortion > 0 )
      {
      typedef MinimumMaximumImageCalculator<ImageType> MinMaxCalcType;
      typename MinMaxCalcType::Pointer calc = MinMaxCalcType::New();
      calc->SetImage( m_FixedImage );
      calc->Compute();
      PixelType fixedImageMax = calc->GetMaximum();
      PixelType fixedImageMin = calc->GetMinimum();

      regRigid->SetFixedImageSamplesIntensityThreshold(
        static_cast<PixelType>( ( m_SampleIntensityPortion *
        (fixedImageMax - fixedImageMin) ) + fixedImageMin ) );
      }
    if( m_UseFixedImageRegionOfInterest )
      {
      regRigid->SetFixedImageRegionOfInterest( m_FixedImageRegionOfInterestPoint1,
        m_FixedImageRegionOfInterestPoint2 );
      }
    regRigid->SetSampleFromOverlap( m_SampleFromOverlap );
    regRigid->SetMetricMethodEnum( m_RigidMetricMethodEnum );
    regRigid->SetInterpolationMethodEnum( m_RigidInterpolationMethodEnum );
    typename RigidTransformType::ParametersType scales;
    if( ImageDimension == 2 )
      {
      scales.set_size( 3 );
      scales[0] = 1.0 / m_ExpectedRotationMagnitude;
      scales[1] = 1.0 / m_ExpectedOffsetMagnitude;
      scales[2] = 1.0 / m_ExpectedOffsetMagnitude;
      }
    else if( ImageDimension == 3 )
      {
      scales.set_size( 6 );
      scales[0] = 1.0 / m_ExpectedRotationMagnitude;
      scales[1] = 1.0 / m_ExpectedRotationMagnitude;
      scales[2] = 1.0 / m_ExpectedRotationMagnitude;
      scales[3] = 1.0 / m_ExpectedOffsetMagnitude;
      scales[4] = 1.0 / m_ExpectedOffsetMagnitude;
      scales[5] = 1.0 / m_ExpectedOffsetMagnitude;
      }
    else
      {
      std::cerr
        << "ERROR: Only 2 and 3 dimensional images are supported."
        << std::endl;
      }
    regRigid->SetTransformParametersScales( scales );

    if( m_CurrentMatrixTransform.IsNotNull() )
      {
      regRigid->GetTypedTransform()->SetCenter(
        m_CurrentMatrixTransform->GetCenter() );
      regRigid->GetTypedTransform()->SetMatrix(
        m_CurrentMatrixTransform->GetMatrix() );
      regRigid->GetTypedTransform()->SetOffset(
        m_CurrentMatrixTransform->GetOffset() );
      regRigid->SetInitialTransformParameters(
        regRigid->GetTypedTransform()->GetParameters() );
      regRigid->SetInitialTransformFixedParameters(
        regRigid->GetTypedTransform()->GetFixedParameters() );
      }

    regRigid->Update();

    m_RigidTransform = RigidTransformType::New();
    m_RigidTransform->SetFixedParameters(
      regRigid->GetTypedTransform()->GetFixedParameters() );
    // must call GetAffineTransform here because the typed transform
    // is a versor and has only 6 parameters (in this code the type
    // RigidTransform is a 12 parameter transform)
    m_RigidTransform->SetParametersByValue(
      regRigid->GetAffineTransform()->GetParameters() );
    m_CurrentMatrixTransform = regRigid->GetAffineTransform();

    m_FinalMetricValue = regRigid->GetFinalMetricValue();
    m_RigidMetricValue = m_FinalMetricValue;

    m_CompletedStage = RIGID_STAGE;
    m_CompletedResampling = false;
    }

  if( m_EnableAffineRegistration )
    {
    this->AffineRegND<ImageDimension>();
    }
}

template <class TSpatialObject, class TImage>
const TImage *
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::ResampleImage( InterpolationMethodEnumType interpolationMethod,
                 const ImageType * movingImage,
                 const MatrixTransformType * matrixTransform,
                 PixelType defaultPixelValue, double portion)
{
  typedef InterpolateImageFunction<TImage, double>
    InterpolatorType;
  typedef NearestNeighborInterpolateImageFunction<TImage, double>
    NearestNeighborInterpolatorType;
  typedef LinearInterpolateImageFunction<TImage, double>
    LinearInterpolatorType;
  typedef BSplineInterpolateImageFunction<TImage, double>
    BSplineInterpolatorType;
  typedef WindowedSincInterpolateImageFunction<TImage, 4,
    Function::HammingWindowFunction<4>, ConstantBoundaryCondition<TImage>,
    double>
    SincInterpolatorType;
  typedef ResampleImageFilter<TImage, TImage, double>
    ResampleImageFilterType;

  typename InterpolatorType::Pointer interpolator = nullptr;

  switch( interpolationMethod )
    {
    case OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION:
      interpolator = NearestNeighborInterpolatorType::New();
      break;
    case OptimizedRegistrationMethodType::LINEAR_INTERPOLATION:
      interpolator = LinearInterpolatorType::New();
      break;
    case OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION:
      interpolator = BSplineInterpolatorType::New();
      ( static_cast<BSplineInterpolatorType *>(
        interpolator.GetPointer() ) )->SetSplineOrder( 3 );
      break;
    case OptimizedRegistrationMethodType::SINC_INTERPOLATION:
      interpolator = SincInterpolatorType::New();
      break;
    default:
      std::cerr
        << "ERROR: Interpolation function not supported"
        << " in itk::SpatialObjectToImageRegistrationHelper::ResampleImage"
        << std::endl;
      interpolator = LinearInterpolatorType::New();
      break;
    }

  if( movingImage == NULL
      && matrixTransform == NULL
      && bsplineTransform == NULL
      && m_CompletedResampling )
    {
    return m_CurrentMovingImage;
    }

  bool doLoaded = false;
  bool doMatrix = false;

  switch( m_CompletedStage )
    {
    default:
    case PRE_STAGE:
      break;
    case LOAD_STAGE:
      doLoaded = true;
      break;
    case INIT_STAGE:
    case RIGID_STAGE:
    case AFFINE_STAGE:
      doMatrix = true;
      break;
    }

  bool resampled = false;
  bool passedImage = false;
  typename TImage::ConstPointer mImage = m_CurrentMovingImage;
  if( movingImage != NULL )
    {
    mImage = movingImage;

    //doLoaded = true;
    //doMatrix = true;
    }

  typename AffineTransformType::ConstPointer aTrans =
    m_CurrentMatrixTransform.GetPointer();
  if( matrixTransform != NULL )
    {
    doLoaded = false;
    doMatrix = false;

    if( matrixTransform != NULL )
      {
      aTrans = matrixTransform;
      doMatrix = true;
      }
    if( bsplineTransform != NULL )
      {
      bTrans = bsplineTransform;
      }
    }

  interpolator->SetInputImage( mImage );

  if( doMatrix && aTrans.IsNotNull() )
    {
    if( this->GetReportProgress() )
      {
      std::cout << "Resampling using matrix." << std::endl;
      }
    // Register using Matrix
    typename ResampleImageFilterType::Pointer resampler =
      ResampleImageFilterType::New();
    resampler->SetInput( mImage );
    resampler->SetInterpolator( interpolator.GetPointer() );
    resampler->SetReferenceImage( m_FixedImage );
    resampler->UseReferenceImageOn();
    typename MatrixTransformType::Pointer tmpTrans = MatrixTransformType::New();
    tmpTrans->SetIdentity();
    tmpTrans->SetFixedParameters( aTrans->GetFixedParameters() );
    if( portion != 1.0 )
      {
      typename MatrixTransformType::ParametersType aTransParams =
        aTrans->GetParameters();
      typename MatrixTransformType::ParametersType tmpParams =
        tmpTrans->GetParameters();
      std::cout << "Original params = " << aTransParams << std::endl;
      for( unsigned int p=0; p<tmpParams.size(); ++p )
        {
        tmpParams[p] = tmpParams[p] + portion * (aTransParams[p]-tmpParams[p]);
        }
      tmpTrans->SetParameters( tmpParams );
      std::cout << "portion params = " << tmpParams << std::endl;
      }
    else
      {
      tmpTrans->SetParameters( aTrans->GetParameters() );
      }
    resampler->SetTransform( tmpTrans );
    resampler->SetDefaultPixelValue( defaultPixelValue );
    resampler->Update();
    m_CurrentMovingImage = resampler->GetOutput();
    m_MatrixTransformResampledImage = m_CurrentMovingImage;

    resampled = true;
    mImage = resampler->GetOutput();
    interpolator->SetInputImage( mImage );
    }

  if( !resampled )
    {
    // Warning: No registrations computed
    if( this->GetReportProgress() )
      {
      std::cout << "Resampling using identity transform." << std::endl;
      }
    typename RigidTransformType::Pointer tmpTransform =
      RigidTransformType::New();
    tmpTransform->SetIdentity();

    interpolator->SetInputImage( mImage );
    typename ResampleImageFilterType::Pointer resampler =
      ResampleImageFilterType::New();
    resampler->SetInput( mImage );
    resampler->SetInterpolator( interpolator.GetPointer() );
    // We should not be casting away constness here, but
    // SetOutputParametersFromImage
    // Does not change the image.  This is needed to workaround fixes
    // to ITK
    //typename ImageType::Pointer tmp = const_cast<ImageType *>(
      //m_FixedImage.GetPointer() );
    //resampler->SetOutputParametersFromImage( tmp );
    resampler->SetReferenceImage( m_FixedImage );
    resampler->UseReferenceImageOn();
    resampler->SetTransform( tmpTransform );
    resampler->SetDefaultPixelValue( defaultPixelValue );
    resampler->Update();

    mImage = resampler->GetOutput();
    }
  else if( !passedImage )
    {
    m_CompletedResampling = true;
    }

  mImage->Register();
  return mImage.GetPointer();
}

template <class TSpatialObject, class TImage>
typename TGroupSpatialObjectType::ConstPointer
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::GetFinalMovingGroupSpatialObject( InterpolationMethodEnumType interpolationMethod,
  PixelType defaultPixelValue)
{
  return ResampleImage( interpolationMethod, nullptr, nullptr, nullptr,
    defaultPixelValue, 1.0 );
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::LoadBaselineImage( const std::string & filename )
{
  typedef ImageFileReader<TImage> ImageReaderType;

  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();

  imageReader->SetFileName( filename );

  imageReader->Update();

  SetBaselineImage( imageReader->GetOutput() );
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::LoadParameters( const std::string & itkNotUsed(filename) )
{
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SaveParameters( const std::string & itkNotUsed(filename) )
{
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::LoadTransform( const std::string & filename, bool invert )
{
  typedef TransformFileReader                    TransformReaderType;
  typedef TransformReaderType::TransformListType TransformListType;

  TransformReaderType::Pointer transformReader = TransformReaderType::New();
  transformReader->SetFileName( filename );

  transformReader->Update();

  const TransformListType *transforms = transformReader->GetTransformList();
  TransformListType::const_iterator transformIt = transforms->begin();
  while( transformIt != transforms->end() )
    {
    if( !strcmp( (*transformIt)->GetNameOfClass(), "AffineTransform") )
      {
      typename MatrixTransformType::Pointer affine_read =
        static_cast<MatrixTransformType *>( (*transformIt).GetPointer() );
      typename MatrixTransformType::ConstPointer affine =
        affine_read.GetPointer();
      SetLoadedMatrixTransform( *affine.GetPointer(), invert );
      }

    ++transformIt;
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SaveTransform( const std::string & filename )
{
  typedef TransformFileWriter TransformWriterType;

  TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetFileName( filename );

  if( m_CurrentMatrixTransform.IsNotNull() )
    {
    transformWriter->SetInput( m_CurrentMatrixTransform );
    transformWriter->Update();
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SaveDisplacementField( const std::string &filename )
{
  typedef itk::Vector<PixelType,TImage::ImageDimension>  VectorType;
  typedef itk::Image<VectorType,TImage::ImageDimension>  DisplacementFieldType;
  typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;

  typename TImage::RegionType fixedImageRegion =
    m_FixedImage->GetBufferedRegion();

  typename DisplacementFieldType::Pointer field =
    DisplacementFieldType::New();
  field->SetRegions( fixedImageRegion );
  field->SetOrigin( m_FixedImage->GetOrigin() );
  field->SetSpacing( m_FixedImage->GetSpacing() );
  field->SetDirection( m_FixedImage->GetDirection() );
  field->Allocate();

  typename DisplacementFieldType::IndexType index;

  VectorType dx;

  FieldIterator it( field, fixedImageRegion );
  it.GoToBegin();

  while( ! it.IsAtEnd() )
    {
    index = it.GetIndex();
    field->TransformIndexToPhysicalPoint( index, fixedPoint );
    movingPoint = fixedPoint;
    if(m_InitialTransform.IsNotNull())
      {
      movingPoint = m_InitialTransform->TransformPoint( movingPoint );
      }
    if(m_RigidTransform.IsNotNull())
      {
      movingPoint = m_RigidTransform->TransformPoint( movingPoint );
      }
    if(m_AffineTransform.IsNotNull())
      {
      movingPoint = m_AffineTransform->TransformPoint( movingPoint );
      }
    dx = movingPoint - fixedPoint;
    it.Set( dx );
    ++it;
    }

  typedef itk::ImageFileWriter< DisplacementFieldType >  FieldWriterType;
  typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

  fieldWriter->SetInput( field );
  fieldWriter->SetFileName( filename );

  fieldWriter->Update();
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetLoadedMatrixTransform( const MatrixTransformType & tfm, bool invert )
{
  m_LoadedMatrixTransform = MatrixTransformType::New();
  m_LoadedMatrixTransform->SetIdentity();
  m_LoadedMatrixTransform->SetFixedParameters( tfm.GetFixedParameters() );
  m_LoadedMatrixTransform->SetCenter( tfm.GetCenter() );
  m_LoadedMatrixTransform->SetMatrix( tfm.GetMatrix() );
  m_LoadedMatrixTransform->SetOffset( tfm.GetOffset() );
  if( invert )
    {
    std::cout << "GetInverseTransform" << std::endl;
    typename MatrixTransformType::Pointer invTfm = MatrixTransformType::New();
    m_LoadedMatrixTransform->GetInverse( invTfm );
    m_LoadedMatrixTransform = invTfm;
    }

  m_EnableLoadedRegistration = true;
  m_LoadedTransformResampledImage = 0;
  m_CurrentMovingImage = m_MovingImage;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetFixedImageRegionOfInterest( const PointType & point1,
                       const PointType & point2 )
{
  m_FixedImageRegionOfInterestPoint1 = point1;
  m_FixedImageRegionOfInterestPoint2 = point2;
  m_UseFixedImageRegionOfInterest = true;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetFixedImageRegionOfInterest( const std::vector<float> & points )
{
  if( points.size() != 2 * ImageDimension )
    {
    throw
      "Error: points to SetFixedImageRegionOfInterest is not twice image dimension";
    }
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_FixedImageRegionOfInterestPoint1[i] = points[i];
    m_FixedImageRegionOfInterestPoint2[i] = points[ImageDimension + i];
    }
  m_UseFixedImageRegionOfInterest = true;
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetFixedLandmarks( const std::vector<std::vector<float> > & fixedLandmarks )
{
  m_FixedLandmarks.clear();
  for( std::vector<std::vector<float> >::const_iterator i =
    fixedLandmarks.begin();
    i != fixedLandmarks.end(); ++i )
    {
    LandmarkPointType landmark;
    std::copy(i->begin(), i->end(), landmark.Begin() );
    m_FixedLandmarks.push_back(landmark);
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::SetMovingLandmarks( const std::vector<std::vector<float> > & movingLandmarks )
{
  m_MovingLandmarks.clear();
  for( std::vector<std::vector<float> >::const_iterator i =
    movingLandmarks.begin();
    i != movingLandmarks.end(); ++i )
    {
    LandmarkPointType landmark;
    std::copy(i->begin(), i->end(), landmark.Begin() );
    m_MovingLandmarks.push_back(landmark);
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::PrintSelfHelper( std::ostream & os, Indent indent,
                   const std::string & basename,
                   MetricMethodEnumType metric,
                   InterpolationMethodEnumType interpolation ) const
{
  switch( metric )
    {
    case OptimizedRegistrationMethodType::MATTES_MI_METRIC:
      os << indent << basename << " Metric Method = MATTES_MI_METRIC"
        << std::endl;
      break;
    case OptimizedRegistrationMethodType::NORMALIZED_CORRELATION_METRIC:
      os << indent << basename
        << " Metric Method = NORMALIZED_CORRELATION_METRIC" << std::endl;
      break;
    case OptimizedRegistrationMethodType::MEAN_SQUARED_ERROR_METRIC:
      os << indent << basename
        << " Metric Method = MEAN_SQUARED_ERROR_METRIC" << std::endl;
      break;
    default:
      os << indent << basename << " Metric Method = UNKNOWN" << std::endl;
      break;
    }
  os << indent << std::endl;

  switch( interpolation )
    {
    case OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = NEAREST_NEIGHBOR_INTERPOLATION"
        << std::endl;
      break;
    case OptimizedRegistrationMethodType::LINEAR_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = LINEAR_INTERPOLATION" << std::endl;
      break;
    case OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = BSPLINE_INTERPOLATION" << std::endl;
      break;
    case OptimizedRegistrationMethodType::SINC_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = SINC_INTERPOLATION" << std::endl;
      break;
    default:
      os << indent << basename
        << " Interpolation Method = UNKNOWN" << std::endl;
      break;
    }
}

template <class TSpatialObject, class TImage>
void
SpatialObjectToImageRegistrationHelper<TSpatialObject, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_FixedImage.IsNotNull() )
    {
    os << indent << "Fixed Image = " << m_FixedImage << std::endl;
    }
  if( m_MovingImage.IsNotNull() )
    {
    os << indent << "Moving Image = " << m_MovingImage << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Use region of interest = " << m_UseFixedImageRegionOfInterest
    << std::endl;
  os << indent << "Region of interest point1 = "
    << m_FixedImageRegionOfInterestPoint1 << std::endl;
  os << indent << "Region of interest point2 = "
    << m_FixedImageRegionOfInterestPoint2 << std::endl;
  os << indent << std::endl;
  os << indent << "Use Fixed Image Mask Object = "
    << m_UseFixedImageMaskObject << std::endl;
  os << indent << std::endl;
  if( m_FixedImageMaskObject.IsNotNull() )
    {
    os << indent << "Fixed Image Mask Object = " << m_FixedImageMaskObject
      << std::endl;
    }
  os << indent << "Use Moving Image Mask Object = "
    << m_UseMovingSpatialObjectMaskObject << std::endl;
  os << indent << std::endl;
  if( m_MovingSpatialObjectMaskObject.IsNotNull() )
    {
    os << indent << "Moving Image Mask Object = "
      << m_MovingSpatialObjectMaskObject << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Random Number Seed = " << m_RandomNumberSeed
    << std::endl;
  os << indent << std::endl;
  os << indent << "Enable Loaded Registration = "
    << m_EnableLoadedRegistration << std::endl;
  os << indent << "Enable Initial Registration = "
    << m_EnableInitialRegistration << std::endl;
  os << indent << "Enable Rigid Registration = "
    << m_EnableRigidRegistration << std::endl;
  os << indent << "Enable Affine Registration = "
    << m_EnableAffineRegistration << std::endl;
  os << indent << std::endl;
  os << indent << "Expected Offset (in Pixels) Magnitude = "
    << m_ExpectedOffsetMagnitude << std::endl;
  os << indent << "Expected Rotation Magnitude = "
    << m_ExpectedRotationMagnitude << std::endl;
  os << indent << "Expected Scale Magnitude = " << m_ExpectedScaleMagnitude
    << std::endl;
  os << indent << "Expected Skew Magnitude = " << m_ExpectedSkewMagnitude
    << std::endl;
  os << indent << "Completed Initialization = "
    << m_CompletedInitialization << std::endl;
  os << indent << "Completed Resampling = " << m_CompletedResampling
    << std::endl;
  os << indent << std::endl;
  os << indent << "Rigid Metric Value = " << m_RigidMetricValue
    << std::endl;
  os << indent << "Affine Metric Value = " << m_AffineMetricValue
    << std::endl;
  os << indent << "Final Metric Value = " << m_FinalMetricValue
    << std::endl;
  os << indent << std::endl;
  os << indent << "Report Progress = " << m_ReportProgress << std::endl;
  os << indent << std::endl;
  if( m_CurrentMovingSpatialObject.IsNotNull() )
    {
    os << indent << "Current Moving Image = " << m_CurrentMovingSpatialObject
      << std::endl;
    }
  else
    {
    os << indent << "Current Moving Image = NULL" << std::endl;
    }
  if( m_CurrentMatrixTransform.IsNotNull() )
    {
    os << indent << "Current Matrix Transform = "
      << m_CurrentMatrixTransform << std::endl;
    }
  else
    {
    os << indent << "Current Matrix Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  if( m_LoadedTransformResampledImage.IsNotNull() )
    {
    os << indent << "Loaded Transform Resampled Image = "
      << m_LoadedTransformResampledImage << std::endl;
    }
  else
    {
    os << indent << "Loaded Transform Resampled Image = NULL" << std::endl;
    }
  if( m_MatrixTransformResampledImage.IsNotNull() )
    {
    os << indent << "Matrix Transform Resampled Image = "
      << m_MatrixTransformResampledImage << std::endl;
    }
  else
    {
    os << indent << "Matrix Transform Resampled Image = NULL" << std::endl;
    }
  os << indent << std::endl;
  if( m_LoadedMatrixTransform.IsNotNull() )
    {
    os << indent << "Loaded Matrix Transform = " << m_LoadedMatrixTransform
      << std::endl;
    }
  else
    {
    os << indent << "Loaded Matrix Transform = NULL" << std::endl;
    }
  os << indent << std::endl;

  switch( m_InitialMethodEnum )
    {
    case INIT_WITH_NONE:
      os << indent << "Initial Registration Enum = INIT_WITH_NONE"
        << std::endl;
      break;
    case INIT_WITH_CURRENT_RESULTS:
      os << indent
        << "Initial Registration Enum = INIT_WITH_CURRENT_RESULTS"
        << std::endl;
      break;
    case INIT_WITH_IMAGE_CENTERS:
      os << indent
        << "Initial Registration Enum = INIT_WITH_IMAGE_CENTERS"
        << std::endl;
      break;
    case INIT_WITH_CENTERS_OF_MASS:
      os << indent
        << "Initial Registration Enum = INIT_WITH_CENTERS_OF_MASS"
        << std::endl;
      break;
    case INIT_WITH_SECOND_MOMENTS:
      os << indent
        << "Initial Registration Enum = INIT_WITH_SECOND_MOMENTS"
        << std::endl;
      break;
    default:
      os << indent << "Initial Registration Enum = UNKNOWN" << std::endl;
      break;
    }
  if( m_InitialTransform.IsNotNull() )
    {
    os << indent << "Initial Transform = " << m_InitialTransform
      << std::endl;
    }
  else
    {
    os << indent << "Initial Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Rigid Sampling Ratio = " << m_RigidSamplingRatio
    << std::endl;
  os << indent << "Rigid Target Error = " << m_RigidTargetError
    << std::endl;
  os << indent << "Rigid Max Iterations = " << m_RigidMaxIterations
    << std::endl;
  PrintSelfHelper( os, indent, "Rigid", m_RigidMetricMethodEnum,
                   m_RigidInterpolationMethodEnum );
  os << indent << std::endl;
  if( m_RigidTransform.IsNotNull() )
    {
    os << indent << "Rigid Transform = " << m_RigidTransform << std::endl;
    }
  else
    {
    os << indent << "Rigid Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Affine Sampling Ratio = " << m_AffineSamplingRatio
    << std::endl;
  os << indent << "Affine Target Error = " << m_AffineTargetError
    << std::endl;
  os << indent << "Affine Max Iterations = " << m_AffineMaxIterations
    << std::endl;
  PrintSelfHelper( os, indent, "Affine", m_AffineMetricMethodEnum,
    m_AffineInterpolationMethodEnum );
  os << indent << std::endl;
  if( m_AffineTransform.IsNotNull() )
    {
    os << indent << "Affine Transform = " << m_AffineTransform
      << std::endl;
    }
  else
    {
    os << indent << "Affine Transform = NULL" << std::endl;
    }
  os << indent << std::endl;

}

}

#endif
