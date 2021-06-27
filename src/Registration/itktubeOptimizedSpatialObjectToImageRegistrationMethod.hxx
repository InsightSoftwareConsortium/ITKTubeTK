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

#ifndef __itktubeOptimizedSpatialObjectToImageRegistrationMethod_txx
#define __itktubeOptimizedSpatialObjectToImageRegistrationMethod_txx

#include "itkOptimizedSpatialObjectToImageRegistrationMethod.h"

#include "itkPointBasedSpatialObjectToImageMetric.h"

#include "itkSpatialObjectToImageRegistrationMethod.h"

#include "itkRealTimeClock.h"

#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkFRPROptimizer.h"
#include "itkImageMaskSpatialObject.h"

#include "itkImage.h"
#include <itkConstantBoundaryCondition.h>

#include <sstream>

namespace itk
{

namespace tube
{

class SpatialObjectToImageRegistrationViewer
  : public Command
{
public:
  typedef SpatialObjectToImageRegistrationViewer Self;
  typedef Command                 Superclass;
  typedef SmartPointer<Self>      Pointer;

  itkTypeMacro( SpatialObjectToImageRegistrationViewer, Command );

  itkNewMacro( SpatialObjectToImageRegistrationViewer );

  typedef SingleValuedNonLinearOptimizer OptimizerType;

  itkSetMacro(DontShowParameters, bool);
  itkSetMacro(UpdateInterval, int);

  void Execute( Object * caller, const EventObject & event ) override
  {
    Execute( (const Object *)caller, event );
  }

  void Execute( const Object * object, const EventObject & event ) override
  {
    if( typeid( event ) != typeid( IterationEvent ) || object == NULL )
      {
      return;
      }

    const OptimizerType * opt = dynamic_cast<const OptimizerType *>(object);

    if( ++m_Iteration % m_UpdateInterval == 0 )
      {
      RealTimeClock::TimeStampType t = m_Clock->GetTimeInSeconds();
      if( !m_DontShowParameters )
        {
        std::cout << "   " << m_Iteration << " : "
                  << opt->GetCurrentPosition() << " = "
                  << opt->GetValue( opt->GetCurrentPosition() )
                  << "   (" << (t - m_LastTime) / m_UpdateInterval << "s)"
                  << std::endl;
        }
      else
        {
        std::cout << "   " << m_Iteration << " : "
                  << opt->GetValue( opt->GetCurrentPosition() )
                  << "   (" << (t - m_LastTime) / m_UpdateInterval << "s)"
                  << std::endl;
        }
      m_LastTime = t;
      }
  }

  void Update()
  {
    this->Execute( (const Object *)NULL, IterationEvent() );
  }

protected:

  RealTimeClock::Pointer       m_Clock;
  RealTimeClock::TimeStampType m_LastTime;

  int  m_Iteration;
  int  m_UpdateInterval;
  bool m_DontShowParameters;

  SpatialObjectToImageRegistrationViewer()
  {
    m_Clock = RealTimeClock::New();
    m_LastTime = m_Clock->GetTimeInSeconds();
    m_Iteration = 0;
    m_UpdateInterval = 1;
    m_DontShowParameters = false;
  };
  ~SpatialObjectToImageRegistrationViewer()
  {
  };

};

template <unsigned int ObjectDimension, class TImage>
OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::OptimizedSpatialObjectToImageRegistrationMethod( void )
{
  m_InitialTransformParameters = TransformParametersType(1);
  m_InitialTransformParameters.Fill( 0.0f ); \

  m_InitialTransformFixedParameters = TransformParametersType(1);
  m_InitialTransformFixedParameters.Fill( 0.0f ); \

  m_LastTransformParameters = TransformParametersType(1);
  m_LastTransformParameters.Fill( 0.0f );

  m_TransformParametersScales = TransformParametersScalesType(1);
  m_TransformParametersScales.Fill( 1.0f );

  // The following MaxIterations value is a good default for rigid
  //   registration.  Other derived registration methods should use
  //   their own default.
  m_MaxIterations = 100;

  m_UseEvolutionaryOptimization = true;

  m_SamplingRatio = 0.01;

  m_TargetError = 0.00001;

  m_RandomNumberSeed = 0;

  m_TransformMethodEnum = RIGID_TRANSFORM;

  m_MetricMethodEnum = POINTS_TO_IMAGE_METRIC;

  m_FinalMetricValue = 0;

}

template <unsigned int ObjectDimension, class TImage>
OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::~OptimizedSpatialObjectToImageRegistrationMethod( void )
{
}

template <unsigned int ObjectDimension, class TImage>
void
OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GenerateData( void )
{
  if( this->GetReportProgress() )
    {
    std::cout << "UPDATE START" << std::endl;
    }

  this->Initialize();

  this->GetTransform()->SetParametersByValue(
    this->GetInitialTransformParameters() );

  typename MetricType::Pointer metric;

  switch( this->GetMetricMethodEnum() )
    {
    default:
    case POINTS_TO_IMAGE_METRIC:
        {
        typedef PointBasedSpatialObjectToImageMetric<ObjectDimension, TImage>
          TypedMetricType;

        typename TypedMetricType::Pointer typedMetric = TypedMetricType::New();
        metric = typedMetric;
        }
      break;
    }
  if( m_RandomNumberSeed != 0 )
    {
    metric->ReinitializeSeed( m_RandomNumberSeed );
    }
  else
    {
    metric->ReinitializeSeed();
    }

  typename ImageType::ConstPointer fixedImage = this->GetFixedImage();
  typename ImageType::ConstPointer movingSpatialObject = this->GetMovingSpatialObject();

  metric->SetFixedImage( fixedImage );
  metric->SetMovingSpatialObject( movingSpatialObject );

  metric->SetSamplingRatio( m_SamplingRatio );

  metric->SetUseFixedImageRegionOfInterest( this->GetUseFixedImageRegionOfInterest() );
  if( this->GetUseFixedImageRegionOfInterest() )
    {
    metric->SetFixedImageRegionOfInterestPoint1( this->GetFixedImageRegionOfInterestPoint1() );
    metric->SetFixedImageRegionOfInterestPoint2( this->GetFixedImageRegionOfInterestPoint2() );
    }

  metric->SetUseFixedImageMaskObject( this->GetUseFixedImageMaskObject() );
  if( this->GetUseFixedImageMaskObject() )
    {
    metric->SetFixedImageMaskObject( this->GetFixedImageMaskObject() );
    }

  metric->SetUseMovingSpatialObjectMaskObject(
    this->GetUseMovingSpatialObjectMaskObject() );
  if( this->GetUseMovingSpatialObjectMaskObject() )
    {
    metric->SetMovingSpatialObjectMaskObject(
      this->GetMovingSpatialObjectMaskObject() );
    }

  try
    {
    this->Optimize(metric);
    }
  catch( itk::ExceptionObject& exception )
    {
    std::cerr << "Optimization threw an exception." << std::endl;
    std::cerr << exception << std::endl ;
    }

  if( this->GetReportProgress() )
    {
    std::cout << "UPDATE END" << std::endl;
    }
}

template <unsigned int ObjectDimension, class TImage>
void
OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::Optimize( MetricType * metric )
{
  typedef SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage> RegType;

  if( m_UseEvolutionaryOptimization )
    {
    if( this->GetReportProgress() )
      {
      std::cout << "EVOLUTIONARY START" << std::endl;
      }

    typedef OnePlusOneEvolutionaryOptimizer EvoOptimizerType;
    EvoOptimizerType::Pointer evoOpt = EvoOptimizerType::New();

    evoOpt->SetNormalVariateGenerator( Statistics::NormalVariateGenerator
                                       ::New() );
    evoOpt->SetEpsilon( this->GetTargetError() );
    evoOpt->Initialize( 0.1 );
    evoOpt->SetCatchGetValueException( true );
    evoOpt->SetMetricWorstPossibleValue( 100 );
    evoOpt->SetScales( this->GetTransformParametersScales() );
    evoOpt->SetMaximumIteration( this->GetMaxIterations() );

    if( this->GetObserver() )
      {
      evoOpt->AddObserver( IterationEvent(), this->GetObserver() );
      }

    if( this->GetReportProgress() )
      {
      typedef SpatialObjectToImageRegistrationViewer ViewerCommandType;
      typename ViewerCommandType::Pointer command = ViewerCommandType::New();
      if( this->GetTransform()->GetNumberOfParameters() > 16 )
        {
        command->SetDontShowParameters( true );
        }
      evoOpt->AddObserver( IterationEvent(), command );
      }

    typename RegType::Pointer reg = RegType::New();
    typename ImageType::ConstPointer fixedImage = this->GetFixedImage();
    typename SpatialObjectType::ConstPointer movingSpatialObject =
      this->GetMovingSpatialObject();
    reg->SetFixedImage( fixedImage );
    reg->SetMovingSpatialObject( movingSpatialObject );
    reg->SetFixedImageRegion( this->GetFixedImage()
                              ->GetLargestPossibleRegion() );
    reg->SetTransform( this->GetTransform() );
    reg->SetInitialTransformParameters(
      this->GetInitialTransformParameters() );
    reg->SetMetric( metric );
    reg->SetOptimizer( evoOpt );

    try
      {
      reg->Update();
      }
    catch( ... )
      {
      std::cout << "Exception caught in evolutionary registration."
                << std::endl;
      std::cout << "Continuing using best values..." << std::endl;
      }

    if( reg->GetLastTransformParameters().size() !=
      this->GetTransform()->GetNumberOfParameters() )
      {
      this->GetTransform()->SetParametersByValue(
        this->GetInitialTransformParameters() );
      }
    else
      {
      this->GetTransform()->SetParametersByValue(
        reg->GetLastTransformParameters() );
      }

    m_FinalMetricValue = reg->GetOptimizer()->GetValue(
      this->GetTransform()->GetParameters() );

    this->SetLastTransformParameters(
      this->GetTransform()->GetParameters() );

    if( this->GetReportProgress() )
      {
      std::cout << "EVOLUTIONARY END" << std::endl;
      }
    }
  else
    {
    this->GetTransform()->SetParametersByValue(
      this->GetInitialTransformParameters() );
    }

  if( this->GetReportProgress() )
    {
    std::cout << "GRADIENT START" << std::endl;
    }

  typedef FRPROptimizer GradOptimizerType;
  GradOptimizerType::Pointer gradOpt = GradOptimizerType::New();

  gradOpt->SetMaximize( false );
  gradOpt->SetCatchGetValueException( true );
  gradOpt->SetMetricWorstPossibleValue( 100 );
  gradOpt->SetStepLength( 0.25 );
  gradOpt->SetStepTolerance( this->GetTargetError() );
  gradOpt->SetMaximumIteration( this->GetMaxIterations() );
  gradOpt->SetMaximumLineIteration( 10 );
  gradOpt->SetScales( this->GetTransformParametersScales() );
  gradOpt->SetUseUnitLengthGradient(true);
  gradOpt->SetToFletchReeves();

  if( this->GetReportProgress() )
    {
    typedef SpatialObjectToImageRegistrationViewer ViewerCommandType;
    typename ViewerCommandType::Pointer command = ViewerCommandType::New();
    if( this->GetTransform()->GetNumberOfParameters() > 16 )
      {
      command->SetDontShowParameters( true );
      }
    gradOpt->AddObserver( IterationEvent(), command );
    }
  if( this->GetObserver() )
    {
    gradOpt->AddObserver( IterationEvent(), this->GetObserver() );
    }

  typename RegType::Pointer reg = RegType::New();
  typename ImageType::ConstPointer fixedImage = this->GetFixedImage();
  typename ImageType::ConstPointer movingSpatialObject =
    this->GetMovingSpatialObject();
  reg->SetFixedImage( fixedImage );
  reg->SetMovingSpatialObject( movingSpatialObject );
  reg->SetFixedImageRegion( this->GetFixedImage()
    ->GetLargestPossibleRegion() );
  reg->SetTransform( this->GetTransform() );
  reg->SetInitialTransformParameters(
    this->GetTransform()->GetParameters() );
  reg->SetMetric( metric );
  reg->SetOptimizer( gradOpt );

  bool failure = false;
  try
    {
    reg->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cout << "Exception caught during gradient registration."
              << excep << std::endl;
    std::cout << "Continuing using best values..." << std::endl;
    std::cout << "  Pos = " << reg->GetLastTransformParameters()
              << std::endl << std::endl;
    if( reg->GetLastTransformParameters().size()
        != reg->GetInitialTransformParameters().size() )
      {
      std::cout << "  Invalid position, using initial parameters." << std::endl;
      failure = true;
      }
    }
  catch( ... )
    {
    std::cout << "Exception caught in gradient registration."
              << std::endl;
    std::cout << "Continuing using best values..." << std::endl;
    std::cout << "  Pos = " << reg->GetLastTransformParameters()
              << std::endl << std::endl;
    if( reg->GetLastTransformParameters().size()
        != reg->GetInitialTransformParameters().size() )
      {
      std::cout << "  Invalid position, using initial parameters." << std::endl;
      failure = true;
      }
    }

  if( failure )
    {
    m_FinalMetricValue = reg->GetOptimizer()->GetValue(
        reg->GetInitialTransformParameters() );
    this->SetLastTransformParameters( reg->GetInitialTransformParameters() );
    this->GetTransform()->SetParametersByValue(
      this->GetInitialTransformParameters() );
    }
  else
    {
    m_FinalMetricValue = reg->GetOptimizer()->GetValue(
        reg->GetLastTransformParameters() );
    this->SetLastTransformParameters( reg->GetLastTransformParameters() );
    this->GetTransform()->SetParametersByValue(
      this->GetLastTransformParameters() );
    }

  if( this->GetReportProgress() )
    {
    std::cout << "GRADIENT END" << std::endl;
    }
}

template <unsigned int ObjectDimension, class TImage>
void
OptimizedSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Initial Transform Parameters = " <<
    m_InitialTransformParameters << std::endl;

  os << indent << "Initial Transform Fixed Parameters = " <<
    m_InitialTransformFixedParameters << std::endl;

  os << indent << "Last Transform Parameters = " <<
    m_LastTransformParameters << std::endl;

  os << indent << "Transform Parameter Scales = " <<
    m_TransformParametersScales << std::endl;

  os << indent << "Max Iterations = " << m_MaxIterations << std::endl;

  os << indent << "Use Evolutionary Optimization = " <<
    m_UseEvolutionaryOptimization << std::endl;

  os << indent << "Sampling Ratio = " << m_SamplingRatio << std::endl;

  os << indent << "Target Error = " << m_TargetError << std::endl;

  switch( m_MetricMethodEnum )
    {
    default:
    case POINTS_TO_IMAGE_METRIC:
      os << indent << "Metric method = Tube to image metric" << std::endl;
      break;
    }
}

}; // tube

}; // itk

#endif
