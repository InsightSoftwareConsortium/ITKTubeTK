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

#ifndef __itktubeSpatialObjectToImageRegistrationHelper_txx
#define __itktubeSpatialObjectToImageRegistrationHelper_txx

#include "itktubeSpatialObjectToImageRegistrationHelper.h"

#include "itktubePointBasedSpatialObjectTransformFilter.h"

#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkVector.h"
#include "itkAffineTransform.h"

namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

  m_LoadedTransformResampledSpatialObject = NULL;
  m_MatrixTransformResampledSpatialObject = NULL;

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
    OptimizedRegistrationMethodType::IMAGE_INTENSITY_METRIC;
  m_RigidMetricValue = 0.0;

  // Affine
  m_AffineSamplingRatio = 0.1;
  m_AffineTargetError = 0.0001;
  m_AffineMaxIterations = 500;
  m_AffineTransform = NULL;
  m_AffineMetricMethodEnum = OptimizedRegistrationMethodType::IMAGE_INTENSITY_METRIC;
  m_AffineMetricValue = 0.0;

  m_Observer = nullptr;
}

template <unsigned int ObjectDimension, class TImage>
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::~SpatialObjectToImageRegistrationHelper()
{
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SetFixedImage( const TImage * fixedImage )
{
  if( this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;

    this->Modified();

    m_CompletedStage = PRE_STAGE;

    m_CompletedInitialization = false;
    m_CompletedResampling = false;
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SetMovingSpatialObject( const SpatialObjectType * movingSpatialObject )
{
  if( this->m_MovingSpatialObject.GetPointer() != movingSpatialObject )
    {
    this->m_MovingSpatialObject = movingSpatialObject;

    this->Modified();

    m_CompletedStage = PRE_STAGE;

    m_CompletedInitialization = false;
    m_CompletedResampling = false;
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SetFixedImageMaskObject( const ImageMaskObjectType * maskObject )
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SetMovingSpatialObjectMaskObject( const SpatialObjectMaskObjectType *
  maskObject )
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SetMetric( MetricMethodEnumType metric )
{
  this->SetRigidMetricMethodEnum( metric );
  this->SetAffineMetricMethodEnum( metric );
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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
    m_CurrentMovingSpatialObject = GetFinalMovingSpatialObject();
    }
  else
    {
    m_CurrentMovingSpatialObject = m_MovingSpatialObject;
    }

  m_LoadedTransformResampledSpatialObject = 0;
  m_MatrixTransformResampledSpatialObject = 0;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::AffineRegND( Image< double, 2 > * itkNotUsed( tmpImage ) )
{
  if( this->GetReportProgress() )
    {
    std::cout << "*** AFFINE REGISTRATION ***" << std::endl;
    }

  typename Affine2DRegistrationMethodType::Pointer regAff
    = Affine2DRegistrationMethodType::New();
  if( this->GetObserver() )
    {
    regAff->SetObserver( this->GetObserver() );
    }
  regAff->SetRandomNumberSeed( m_RandomNumberSeed );
  regAff->SetReportProgress( m_ReportProgress );
  regAff->SetMovingSpatialObject( m_CurrentMovingSpatialObject );
  regAff->SetFixedImage( m_FixedImage );
  regAff->SetSamplingRatio( m_AffineSamplingRatio );
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
  regAff->SetMetricMethodEnum( m_AffineMetricMethodEnum );

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
    regAff->GetTypedTransform()->SetMatrix(
      m_CurrentMatrixTransform->GetMatrix() );
    regAff->GetTypedTransform()->SetOffset(
      m_CurrentMatrixTransform->GetOffset() );
    regAff->SetInitialTransformFixedParameters(
      regAff->GetTypedTransform()->GetFixedParameters() );
    regAff->SetInitialTransformParameters(
      regAff->GetTypedTransform()->GetParameters() );
    }

  regAff->Update();

  m_AffineTransform = regAff->GetAffineTransform();
  m_CurrentMatrixTransform = m_AffineTransform;

  m_FinalMetricValue = regAff->GetFinalMetricValue();
  m_AffineMetricValue = m_FinalMetricValue;

  m_CompletedStage = AFFINE_STAGE;
  m_CompletedResampling = false;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::AffineRegND( Image< double, 3 > * itkNotUsed( tmpImage ) )
{
  if( this->GetReportProgress() )
    {
    std::cout << "*** AFFINE REGISTRATION ***" << std::endl;
    }

  typename Affine3DRegistrationMethodType::Pointer regAff =
    Affine3DRegistrationMethodType::New();
  if( this->GetObserver() )
    {
    regAff->SetObserver( this->GetObserver() );
    }
  regAff->SetRandomNumberSeed( m_RandomNumberSeed );
  regAff->SetReportProgress( m_ReportProgress );
  regAff->SetMovingSpatialObject( m_CurrentMovingSpatialObject );
  regAff->SetFixedImage( m_FixedImage );
  regAff->SetSamplingRatio( m_AffineSamplingRatio );
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
  regAff->SetMetricMethodEnum( m_AffineMetricMethodEnum );

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
    regAff->SetInitialTransformFixedParameters(
      regAff->GetTypedTransform()->GetFixedParameters() );
    regAff->SetInitialTransformParameters(
      regAff->GetTypedTransform()->GetParameters() );
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
template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::Update( void )
{
  if( !(this->m_CompletedInitialization) )
    {
    this->Initialize();
    }

  std::cout << "HERE" << std::endl;
  if( m_EnableLoadedRegistration
      && m_LoadedMatrixTransform.IsNotNull() )
    {
    if( m_LoadedTransformResampledSpatialObject.IsNotNull() )
      {
      m_CurrentMovingSpatialObject = m_LoadedTransformResampledSpatialObject;
      if( this->GetReportProgress() )
        {
        std::cout << "*** Using existing loaded transform ***" << std::endl;
        }
      }

    m_MatrixTransformResampledSpatialObject = 0;

    m_CompletedStage = LOAD_STAGE;
    m_CompletedResampling = true;

    m_CurrentMatrixTransform = 0;
    }

  std::cout << "HERE1" << std::endl;
  if( this->GetReportProgress() )
    {
    std::cout << "*** INITIAL REGISTRATION ***" << std::endl;
    }

  typename InitialRegistrationMethodType::Pointer regInit =
    InitialRegistrationMethodType::New();
  if( this->GetObserver() )
    {
    regInit->SetObserver( this->GetObserver() );
    }
  regInit->SetReportProgress( m_ReportProgress );
  regInit->SetMovingSpatialObject( m_CurrentMovingSpatialObject );
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
  std::cout << "HERE2" << std::endl;
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
      case INIT_WITH_LANDMARKS:
        regInit->SetUseLandmarks( true );
        regInit->SetFixedLandmarks( m_FixedLandmarks );
        regInit->SetMovingLandmarks( m_MovingLandmarks );
        break;
      case INIT_WITH_CURRENT_RESULTS:
      default:
        break;
      }
    std::cout << "RH: Update: InitialReg" << std::endl;
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

  std::cout << "RH: Update: Init done" << std::endl;
  m_CurrentMatrixTransform = m_InitialTransform;

  m_CompletedStage = INIT_STAGE;
  m_CompletedResampling = false;

  if( m_EnableRigidRegistration )
    {
    std::cout << "RH: Update: registration" << std::endl;
    if( this->GetReportProgress() )
      {
      std::cout << "*** RIGID REGISTRATION ***" << std::endl;
      }

    typename RigidRegistrationMethodType::Pointer regRigid =
     RigidRegistrationMethodType::New();
    if( this->GetObserver() )
      {
      regRigid->SetObserver( this->GetObserver() );
      }
    regRigid->SetRandomNumberSeed( m_RandomNumberSeed );
    if( !m_UseEvolutionaryOptimization )
      {
      regRigid->SetUseEvolutionaryOptimization( false );
      }
    regRigid->SetReportProgress( m_ReportProgress );
    regRigid->SetMovingSpatialObject( m_CurrentMovingSpatialObject );
    regRigid->SetFixedImage( m_FixedImage );
    regRigid->SetSamplingRatio( m_RigidSamplingRatio );
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
    regRigid->SetMetricMethodEnum( m_RigidMetricMethodEnum );
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

    std::cout << "RH: Update: set transform parameters scales" << std::endl;
    if( m_CurrentMatrixTransform.IsNotNull() )
      {
      regRigid->GetTypedTransform()->SetCenter(
        m_CurrentMatrixTransform->GetCenter() );
      regRigid->GetTypedTransform()->SetMatrix(
        m_CurrentMatrixTransform->GetMatrix() );
      regRigid->GetTypedTransform()->SetOffset(
        m_CurrentMatrixTransform->GetOffset() );
      regRigid->SetInitialTransformFixedParameters(
        regRigid->GetTypedTransform()->GetFixedParameters() );
      regRigid->SetInitialTransformParameters(
        regRigid->GetTypedTransform()->GetParameters() );
      }

    std::cout << "RH: Update: update" << std::endl;
    regRigid->Update();
    std::cout << "RH: Update: update done" << std::endl;

    m_RigidTransform = RigidTransformType::New();
    m_RigidTransform->SetFixedParameters(
      regRigid->GetTypedTransform()->GetFixedParameters() );
    // Rigid transform is stored as a matrix+offset (affine) transform 
    // because no rigid transform is templated over dimension.
    // See itktubeRigidSpatialObjectToImageRegistrationMethod class.
    m_RigidTransform->SetParametersByValue(
      regRigid->GetAffineTransform()->GetParameters() );
    m_CurrentMatrixTransform = regRigid->GetAffineTransform();

    m_FinalMetricValue = regRigid->GetFinalMetricValue();
    m_RigidMetricValue = m_FinalMetricValue;

    m_CompletedStage = RIGID_STAGE;
    m_CompletedResampling = false;
    }

  std::cout << "RH: Update: done" << std::endl;

  if( m_EnableAffineRegistration )
    {
    std::cout << "RH: Update: affine" << std::endl;
    this->AffineRegND<ImageDimension>();
    }
}

template <unsigned int ObjectDimension, class TImage>
const SpatialObject<ObjectDimension> *
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::ResampleSpatialObject( const SpatialObjectType * movingSpatialObject,
                 const MatrixTransformType * matrixTransform, double portion)
{
  typedef PointBasedSpatialObjectTransformFilter<MatrixTransformType,
    ObjectDimension> ResampleFilterType;

  if( movingSpatialObject == NULL
      && matrixTransform == NULL
      && portion == 1.0
      && m_CompletedResampling )
    {
    return m_CurrentMovingSpatialObject.GetPointer();
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

  typename SpatialObjectType::ConstPointer inputSO =
    m_CurrentMovingSpatialObject;
  if( movingSpatialObject != NULL )
    {
    inputSO = movingSpatialObject;

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
    }

  if( doMatrix && aTrans.IsNotNull() )
    {
    if( this->GetReportProgress() )
      {
      std::cout << "Resampling using matrix." << std::endl;
      }
    // Register using Matrix
    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput( inputSO );
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
    resampler->Update();
    m_MatrixTransformResampledSpatialObject = resampler->GetOutput();
    if( movingSpatialObject == NULL
        && matrixTransform == NULL
        && portion == 1.0 )
      {
      m_CurrentMovingSpatialObject = resampler->GetOutput();
      m_CompletedResampling = true;
      }
    }

  m_CurrentMovingSpatialObject->Register();
  return m_CurrentMovingSpatialObject.GetPointer();
}

template <unsigned int ObjectDimension, class TImage>
const SpatialObject<ObjectDimension> *
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::GetFinalMovingSpatialObject( void )
{
  return ResampleSpatialObject();
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::LoadParameters( const std::string & itkNotUsed(filename) )
{
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::SaveParameters( const std::string & itkNotUsed(filename) )
{
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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
  m_LoadedTransformResampledSpatialObject = 0;
  m_CurrentMovingSpatialObject = m_MovingSpatialObject;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::PrintSelfHelper( std::ostream & os, Indent indent,
                   const std::string & basename,
                   MetricMethodEnumType metric ) const
{
  switch( metric )
    {
    case OptimizedRegistrationMethodType::IMAGE_INTENSITY_METRIC:
      os << indent << basename << " Metric Method = IMAGE_INTENSITY_METRIC"
        << std::endl;
      break;
    default:
      os << indent << basename << " Metric Method = UNKNOWN" << std::endl;
      break;
    }
  os << indent << std::endl;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationHelper<ObjectDimension, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_FixedImage.IsNotNull() )
    {
    os << indent << "Fixed Image = " << m_FixedImage << std::endl;
    }
  if( m_MovingSpatialObject.IsNotNull() )
    {
    os << indent << "Moving SpatialObject = " << m_MovingSpatialObject << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Use Fixed Image Mask Object = "
    << m_UseFixedImageMaskObject << std::endl;
  os << indent << std::endl;
  if( m_FixedImageMaskObject.IsNotNull() )
    {
    os << indent << "Fixed Image Mask Object = " << m_FixedImageMaskObject
      << std::endl;
    }
  os << indent << "Use Moving SpatialObject Mask Object = "
    << m_UseMovingSpatialObjectMaskObject << std::endl;
  os << indent << std::endl;
  if( m_MovingSpatialObjectMaskObject.IsNotNull() )
    {
    os << indent << "Moving SpatialObject Mask Object = "
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
    os << indent << "Current Moving SpatialObject = " << m_CurrentMovingSpatialObject
      << std::endl;
    }
  else
    {
    os << indent << "Current Moving SpatialObject = NULL" << std::endl;
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
  if( m_LoadedTransformResampledSpatialObject.IsNotNull() )
    {
    os << indent << "Loaded Transform Resampled SpatialObject = "
      << m_LoadedTransformResampledSpatialObject << std::endl;
    }
  else
    {
    os << indent << "Loaded Transform Resampled SpatialObject = NULL" << std::endl;
    }
  if( m_MatrixTransformResampledSpatialObject.IsNotNull() )
    {
    os << indent << "Matrix Transform Resampled SpatialObject = "
      << m_MatrixTransformResampledSpatialObject << std::endl;
    }
  else
    {
    os << indent << "Matrix Transform Resampled SpatialObject = NULL" << std::endl;
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
  PrintSelfHelper( os, indent, "Rigid", m_RigidMetricMethodEnum );
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
  PrintSelfHelper( os, indent, "Affine", m_AffineMetricMethodEnum );
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

}; // tube

}; // itk

#endif
