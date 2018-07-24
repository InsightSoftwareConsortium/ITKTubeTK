/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 (Tue, 10 Jul 2007) $
  Version:   $Revision: 0 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineImageToImageRegistrationMethod_txx
#define __itkBSplineImageToImageRegistrationMethod_txx

#include "itkBSplineImageToImageRegistrationMethod.h"

#include "itkBSplineTransformInitializer.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkLBFGSBOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"

#include "itkRealTimeClock.h"
#include "itkCommand.h"

namespace itk
{

class BSplineImageRegistrationViewer
  : public Command
{
public:
  typedef BSplineImageRegistrationViewer Self;
  typedef Command                        Superclass;
  typedef SmartPointer<Self>             Pointer;

  itkTypeMacro( BSplineImageRegistrationViewer, Command );

  itkNewMacro( BSplineImageRegistrationViewer );

  typedef SingleValuedNonLinearOptimizer OptimizerType;

  itkSetMacro(DontShowParameters, bool);
  itkSetMacro(UpdateInterval, int);

  void Execute( Object * caller, const EventObject & event )
    {
    Execute( (const Object *)caller, event );
    }

  void Execute( const Object * object, const EventObject & event )
    {
    if( typeid( event ) != typeid( IterationEvent ) || object == NULL )
      {
      return;
      }

    const OptimizerType * opt = dynamic_cast< const OptimizerType * >(
      object );

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

  BSplineImageRegistrationViewer()
    {
    m_Clock = RealTimeClock::New();
    m_LastTime = m_Clock->GetTimeInSeconds();
    m_Iteration = 0;
    m_UpdateInterval = 1;
    m_DontShowParameters = false;
    }

  ~BSplineImageRegistrationViewer()
    {
    }

};

template <class TImage>
BSplineImageToImageRegistrationMethod<TImage>
::BSplineImageToImageRegistrationMethod( void )
{
  m_NumberOfControlPoints = 10;
  m_NumberOfLevels = 3;
  m_ExpectedDeformationMagnitude = 10;
  m_GradientOptimizeOnly = false;
  this->SetTransformMethodEnum( Superclass::BSPLINE_TRANSFORM );

  // Override superclass defaults:
  this->SetMaxIterations( 40 );
  this->SetNumberOfSamples( 800000 );
  this->SetInterpolationMethodEnum( Superclass::BSPLINE_INTERPOLATION );
}

template <class TImage>
BSplineImageToImageRegistrationMethod<TImage>
::~BSplineImageToImageRegistrationMethod( void )
{
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::GenerateData( void )
{
  typename TransformType::Pointer tmpTrans = TransformType::New();
  this->SetTransform( tmpTrans );

  typename TransformType::MeshSizeType meshSize;
  meshSize.Fill( this->GetNumberOfControlPoints() - 3 );

  typedef itk::BSplineTransformInitializer< TransformType, ImageType > 
    TransformInitializerType;
  typename TransformInitializerType::Pointer initializer = 
    TransformInitializerType::New();
  initializer->SetTransform( tmpTrans );
  initializer->SetImage( this->GetFixedImage() );
  initializer->SetTransformDomainMeshSize( meshSize );
  initializer->InitializeTransform();
  tmpTrans->SetIdentity();

  const unsigned int numberOfParameters =
    this->GetTypedTransform()->GetNumberOfParameters();
  std::cout << "   numberOfParameters = "
    << numberOfParameters << std::endl;

  /* Remember the fixed parameters for this transform */
  this->SetInitialTransformFixedParameters(
    this->GetTypedTransform()->GetFixedParameters() );

  /* Make sure Initial Transform Parameters are valid */
  if( numberOfParameters !=
    this->GetInitialTransformParameters().GetSize() )
    {
    std::cout
      << "ERROR: numberOfParameters != InitialTransformParameters.size()"
      << std::endl << "   Using identity trasnform." << std::endl;
    std::cout << numberOfParameters << " != " 
      << this->GetInitialTransformParameters().GetSize() << std::endl;
    typename Superclass::TransformParametersType params(
      numberOfParameters );
    params.Fill( 0.0 );
    this->SetInitialTransformParameters( params );
    }

  /* Set scales = expected amount of movement of a control point */
  typename Superclass::TransformParametersType params(
    numberOfParameters );
  typename TransformType::SpacingType spacing = this->GetFixedImage()
    ->GetSpacing();
  for( unsigned int id = 0; id < ImageDimension; ++id )
    {
    params[id] = 1.0 / (m_ExpectedDeformationMagnitude * spacing[id]);
    }
  this->SetTransformParametersScales( params );

  this->Superclass::GenerateData();
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::Optimize( MetricType * metric, InterpolatorType * interpolator )
{
  if( this->GetGradientOptimizeOnly() )
    {
    this->GradientOptimize( metric, interpolator );
    }
  else
    {
    this->MultiResolutionOptimize( metric, interpolator );
    }
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::GradientOptimize( MetricType * metric,
                    InterpolatorType * interpolator )
{
  //if( this->GetReportProgress() )
    {
    std::cout << "BSpline GRADIENT START" << std::endl;
    }

  /* Setup FRPR - set params specific to this optimizer
  typedef FRPROptimizer GradOptimizerType;
  GradOptimizerType::Pointer gradOpt;
  gradOpt = GradOptimizerType::New();
  gradOpt->SetMaximize( false );
  gradOpt->SetCatchGetValueException( true );
  gradOpt->SetMetricWorstPossibleValue( 0 );
  gradOpt->SetStepLength( 0.25 );
  gradOpt->SetStepTolerance( this->GetTargetError() );
  gradOpt->SetMaximumIteration( this->GetMaxIterations() );
  gradOpt->SetMaximumLineIteration( 10 );
  gradOpt->SetScales( this->GetTransformParametersScales() );
  gradOpt->SetUseUnitLengthGradient(true);
  gradOpt->SetToFletchReeves();*/

  /* GradientDescent */
  typedef GradientDescentOptimizer         GradOptimizerType;
  GradOptimizerType::Pointer gradOpt;
  gradOpt = GradOptimizerType::New();
  gradOpt->SetLearningRate( 0.25 );
  gradOpt->SetMaximize( false );
  gradOpt->SetNumberOfIterations( this->GetMaxIterations() );

  /* LBFGSB */
  /*
  int numberOfParameters = this->GetTransform()->GetNumberOfParameters();
  typedef LBFGSBOptimizer                  GradOptimizerType;
  GradOptimizerType::Pointer gradOpt;
  gradOpt = GradOptimizerType::New();
  GradOptimizerType::Pointer tmpOpt = static_cast<GradOptimizerType *>(
    gradOpt.GetPointer() );

  GradOptimizerType::BoundSelectionType boundSelect( numberOfParameters );
  GradOptimizerType::BoundValueType upperBound( numberOfParameters );
  GradOptimizerType::BoundValueType lowerBound( numberOfParameters );
  boundSelect.Fill( 0 );
  upperBound.Fill( 0.0 );
  lowerBound.Fill( 0.0 );
  tmpOpt->SetBoundSelection( boundSelect );
  tmpOpt->SetUpperBound( upperBound );
  tmpOpt->SetLowerBound( lowerBound );
  tmpOpt->SetCostFunctionConvergenceFactor( 1000000 );
  tmpOpt->SetProjectedGradientTolerance( 1e-6 );
  tmpOpt->SetMaximumNumberOfIterations( 1 );//this->GetMaxIterations() );
  tmpOpt->SetMaximumNumberOfEvaluations( 1 );//this->GetMaxIterations() / 10 );
  tmpOpt->SetMaximumNumberOfCorrections( 1 );//this->GetMaxIterations() / 20 );
  */

  /**/
  /*  Create observers for reporting progress */
  /**/
  //if( this->GetReportProgress() )
    {
    typedef BSplineImageRegistrationViewer ViewerCommandType;
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

  /**/
  /* Setup a standard itk::ImageToImageRegistration method
   *  and plug-in optimizer, interpolator, etc. */
  /**/
  typedef ImageRegistrationMethod<TImage, TImage> ImageRegistrationType;
  typename ImageRegistrationType::Pointer reg =
    ImageRegistrationType::New();
  typename ImageType::ConstPointer fixedImage = this->GetFixedImage();
  typename ImageType::ConstPointer movingImage = this->GetMovingImage();
  reg->SetFixedImage( fixedImage );
  reg->SetMovingImage( movingImage );
  reg->SetFixedImageRegion( this->GetFixedImage()
                            ->GetLargestPossibleRegion() );
  reg->SetTransform( this->GetTransform() );
  reg->SetInitialTransformParameters(
    this->GetInitialTransformParameters() );
  reg->GetTransform()->SetParametersByValue(
    this->GetInitialTransformParameters() );
  reg->SetMetric( metric );
  reg->SetOptimizer( gradOpt );
  reg->SetInterpolator( interpolator );
  reg->SetDebug( true );

  //if( this->GetReportProgress() )
    {
    typename TransformType::OutputPointType p;
    for( unsigned int pdi = 0; pdi<ImageDimension; ++pdi )
      {
      p[pdi] =
        this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[pdi]
        / 2.0;
      }
    p = reg->GetTransform()->TransformPoint(p);
    std::cout << "Initial Point = " << p << std::endl;
    }

  /**/
  /*  Optimize! */
  /**/
  bool failure = false;
  try
    {
    std::cout << "  InitialParams = " << this->GetInitialTransformParameters()
      << std::endl;
    std::cout << "   reg->Update()" << std::endl;
    reg->Update();
    std::cout << "  FinalParams = " << reg->GetLastTransformParameters()
      << std::endl;
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cout
      << "Exception caught during gradient registration of BSpline."
      << excep << std::endl;
    std::cout << "Continuing using best values..." << std::endl;
    std::cout << "  Pos = " << reg->GetLastTransformParameters()
              << std::endl << std::endl;
    if( reg->GetLastTransformParameters().size()
        != reg->GetInitialTransformParameters().size() )
      {
      std::cout << "  Invalid position, using initial parameters."
                << std::endl;
      reg->GetTransform()->SetParametersByValue(
        reg->GetInitialTransformParameters() );
      failure = true;
      }
    }
  catch( ... )
    {
    std::cerr << "Error in gradient optimization of BSpline." << std::endl;
    std::cout << "Continuing using best values..." << std::endl;
    std::cout << "  Pos = " << reg->GetLastTransformParameters()
              << std::endl << std::endl;
    if( reg->GetLastTransformParameters().size()
        != reg->GetInitialTransformParameters().size() )
      {
      std::cout << "  Invalid position, using initial parameters."
        << std::endl;
      reg->GetTransform()->SetParametersByValue(
        reg->GetInitialTransformParameters() );
      failure = true;
      }
    }

  /**/
  /*  Record results */
  /**/
  if( failure )
    {
    this->SetFinalMetricValue( reg->GetOptimizer()
      ->GetValue( reg->GetInitialTransformParameters() ) );

    this->SetLastTransformParameters(
      reg->GetInitialTransformParameters() );
    this->GetTransform()->SetParametersByValue(
      this->GetInitialTransformParameters() );
    }
  else
    {
    this->SetFinalMetricValue( reg->GetOptimizer()
      ->GetValue( reg->GetLastTransformParameters() ) );

    this->SetLastTransformParameters( reg->GetLastTransformParameters() );
    this->GetTransform()->SetParametersByValue(
      this->GetLastTransformParameters() );
    }

  //if( this->GetReportProgress() )
    {
    typename TransformType::OutputPointType p;
    for( unsigned int pdi = 0; pdi<ImageDimension; ++pdi )
      {
      p[pdi] =
        this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[pdi]
        / 2.0;
      }
    p = reg->GetTransform()->TransformPoint(p);
    std::cout << "Resulting Point = " << p << std::endl;
    }

  if( this->GetReportProgress() )
    {
    std::cout << "BSpline GRADIENT END" << std::endl;
    }
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::MultiResolutionOptimize( MetricType * itkNotUsed(metric),
                           InterpolatorType * itkNotUsed(interpolator) )
{
  if( this->GetReportProgress() )
    {
    std::cout << "BSpline MULTIRESOLUTION START" << std::endl;
    }

  typedef RecursiveMultiResolutionPyramidImageFilter<ImageType,
                                                     ImageType>
  PyramidType;
  typename PyramidType::Pointer fixedPyramid = PyramidType::New();
  typename PyramidType::Pointer movingPyramid = PyramidType::New();

  /**/
  /* Determine the control points, samples, and scales to be used at
   * each level */
  /**/
  double       controlPointFactor = 2;
  unsigned int levelNumberOfControlPoints =
    this->GetNumberOfControlPoints();
  double levelScale = 1;
  unsigned int numberOfLevelsUsing m_NumberOfLevels;
  if( this->m_NumberOfLevels > 1 )
    {
    for( unsigned int level = 1; level < this->m_NumberOfLevels; level++ )
      {
      levelNumberOfControlPoints = (unsigned int)(
        levelNumberOfControlPoints / controlPointFactor );
      levelScale *= controlPointFactor;
      if( levelNumberOfControlPoints < 3 )  // splineOrder
        {
        levelNumberOfControlPoints = 3;
        numberOfLevelsUsing = level;
        break;
        }
      }
    }

  /**/
  /* Setup the multi-scale image pyramids */
  /**/
  fixedPyramid->SetNumberOfLevels( numberOfLevelsUsing );
  movingPyramid->SetNumberOfLevels( numberOfLevelsUsing );

  typename ImageType::SpacingType fixedSpacing =
    this->GetFixedImage()->GetSpacing();
  typename ImageType::SpacingType movingSpacing =
    this->GetFixedImage()->GetSpacing();

  typename PyramidType::ScheduleType fixedSchedule =
    fixedPyramid->GetSchedule();
  typename PyramidType::ScheduleType movingSchedule =
    movingPyramid->GetSchedule();

  /**/
  /*   First, determine the pyramid at level 0 */
  /**/
  unsigned int level = 0;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    fixedSchedule[0][i] = (unsigned int)(levelScale);
    movingSchedule[0][i] = (unsigned int)(levelScale);
    }

  /**/
  /*   Second, determine the pyramid at the remaining levels */
  /**/
  for( level = 1; level < numberOfLevelsUsing; level++ )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      fixedSchedule[level][i] = (int)(fixedSchedule[level - 1][i]
                                      / controlPointFactor);
      if( fixedSchedule[level][i] < 1 )
        {
        fixedSchedule[level][i] = 1;
        }
      movingSchedule[level][i] = (int)(movingSchedule[level - 1][i]
                                       / controlPointFactor);
      if( movingSchedule[level][i] < 1 )
        {
        movingSchedule[level][i] = 1;
        }
      }
    }

  /**/
  /*   Third, apply pyramid to fixed image */
  /**/
  fixedPyramid->SetSchedule( fixedSchedule );
  fixedPyramid->SetInput( this->GetFixedImage() );
  fixedPyramid->Update();

  /**/
  /*   Fourth, apply pyramid to moving image */
  /**/
  movingPyramid->SetSchedule( movingSchedule );
  movingPyramid->SetInput( this->GetMovingImage() );
  movingPyramid->Update();

  /**/
  /* Assign initial transform parameters at coarse level based on
   *   initial transform parameters - initial transform parameters are
   *   set in the GenerateData() method */
  /**/
  typename Superclass::TransformParametersType levelParameters;
  this->ResampleControlGrid( levelNumberOfControlPoints, levelParameters );
  /* Perform registration at each level */
  for( level = 0; level < numberOfLevelsUsing; level++ )
    {
    //if( this->GetReportProgress() )
      {
      std::cout << "MULTIRESOLUTION LEVEL = " << level << std::endl;
      std::cout << "   Number of control points = "
        << levelNumberOfControlPoints << std::endl;
      std::cout << "   Fixed image = "
        << fixedPyramid->GetOutput(level)
          ->GetLargestPossibleRegion().GetSize() << std::endl;
      std::cout << "   Moving image = "
        << movingPyramid->GetOutput(level)
          ->GetLargestPossibleRegion().GetSize() << std::endl;
      std::cout << "   Parameters.size = "
        << levelParameters.size()
        << std::endl;
      }

    /**/
    /* Get the fixed and moving images for this pyramid level */
    /**/
    typename ImageType::ConstPointer fixedImage =
      fixedPyramid->GetOutput(level);
    typename ImageType::ConstPointer movingImage =
      movingPyramid->GetOutput(level);

    /*
    typedef itk::ImageFileWriter< ImageType > FileWriterType;
    typename FileWriterType::Pointer writer = FileWriterType::New();
    std::stringstream ss;
    std::string name;

    if( this->GetReportProgress() )
      {
      writer->SetInput( fixedImage );
      ss << "level" << level << "Fixed.mha";
      ss >> name;
      writer->SetFileName( name );
      writer->Update();

      writer->SetInput( movingImage );
      ss.clear();
      ss << "level" << level << "Moving.mha";
      ss >> name;
      writer->SetFileName( name );
      writer->Update();
      }
    */

    double levelFactor = levelNumberOfControlPoints /
      (double)(this->GetNumberOfControlPoints() );

    double levelDeformationMagnitude =
      this->GetExpectedDeformationMagnitude(); // * levelFactor;
    // should be physical units, so no need to re-scale per level

    unsigned int levelNumberOfSamples = (unsigned int)(
      this->GetNumberOfSamples() / levelFactor);
    unsigned int fixedImageNumberOfSamples =
      fixedImage->GetLargestPossibleRegion().GetNumberOfPixels();
    if( levelNumberOfSamples > fixedImageNumberOfSamples )
      {
      levelNumberOfSamples = fixedImageNumberOfSamples;
      }

    //if( this->GetReportProgress() )
      {
      std::cout << "   Deformation magnitude = "
        << levelDeformationMagnitude << std::endl;
      std::cout << "   Number of samples = " << levelNumberOfSamples
        << std::endl;
      }

    /**/
    /* Create a BSpline registration module (an instance of this class!)
     * that
     * will perform gradient optimization (instead of pyramid
     * optimization). */
    /**/
    typedef BSplineImageToImageRegistrationMethod<ImageType>
      BSplineRegType;
    typename BSplineRegType::Pointer reg = BSplineRegType::New();
    reg->SetReportProgress( this->GetReportProgress() );
    reg->SetFixedImage( fixedImage );
    reg->SetMovingImage( movingImage );
    reg->SetNumberOfControlPoints( levelNumberOfControlPoints );
    reg->SetNumberOfSamples( levelNumberOfSamples );
    reg->SetExpectedDeformationMagnitude( levelDeformationMagnitude );
    reg->SetGradientOptimizeOnly( true );
    reg->SetTargetError( this->GetTargetError() );
    reg->SetSampleFromOverlap( this->GetSampleFromOverlap() );
    reg->SetFixedImageSamplesIntensityThreshold(
      this->GetFixedImageSamplesIntensityThreshold() );
    reg->SetUseFixedImageSamplesIntensityThreshold(
      this->GetUseFixedImageSamplesIntensityThreshold() );
    reg->SetMaxIterations( (unsigned int)(this->GetMaxIterations()
      * levelFactor) );
    reg->SetMetricMethodEnum( this->GetMetricMethodEnum() );
    reg->SetInterpolationMethodEnum( this->GetInterpolationMethodEnum() );
    std::cout << "pre levelParameters = " << levelParameters << std::endl;
    reg->SetInitialTransformParameters( levelParameters );
    // For the last two levels (the ones at the highest resolution, use
    //   user-specified values of MinimizeMemory, otherwise do not
    //   minimizeMemory so as to maximize speed.
    if( level >= numberOfLevelsUsing - 2 )
      {
      reg->SetMinimizeMemory( this->GetMinimizeMemory() );
      }
    else
      {
      reg->SetMinimizeMemory( false );
      }

    try
      {
      std::cout << "   reg->GetInitialTransformParameters() ="
        << reg->GetInitialTransformParameters() << std::endl;
      std::cout << "   reg->Update()" << std::endl;
      reg->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Exception caught during level registration."
                << excep << std::endl;
      std::cout << "Current Matrix Transform = " << std::endl;
      reg->GetTransform()->Print(std::cout, 2);
      }
    catch( ... )
      {
      std::cout << "Uncaught exception during helper class registration."
                << std::endl;
      }

    if( level < numberOfLevelsUsing - 1 )
      {
      std::cout << "post levelParameters = "
        << reg->GetLastTransformParameters() << std::endl;
      levelNumberOfControlPoints = (unsigned int)(
        levelNumberOfControlPoints * controlPointFactor);
      if( levelNumberOfControlPoints > this->GetNumberOfControlPoints() ||
        level == numberOfLevelsUsing - 2 )
        {
        levelNumberOfControlPoints = this->GetNumberOfControlPoints();
        }

      if( levelNumberOfControlPoints != reg->GetNumberOfControlPoints() )
        {
        //if( this->GetReportProgress() )
          {
          std::cout << "   Resampling grid..." << std::endl;
          }
        reg->ResampleControlGrid( levelNumberOfControlPoints,
          levelParameters );
        }
      else
        {
        levelParameters = reg->GetLastTransformParameters();
        }
      std::cout << "post resample levelParameters = " << levelParameters
        << std::endl;
      }
    else
      {
      /**/
      /*  Remember the results */
      /**/
      this->SetFinalMetricValue( reg->GetFinalMetricValue() );
      this->SetLastTransformParameters(
        reg->GetLastTransformParameters() );
      this->GetTransform()->SetParametersByValue(
        this->GetLastTransformParameters() );
      std::cout << "final levelParameters = "
        << this->GetLastTransformParameters() << std::endl;
      }

    if( this->GetReportProgress() )
      {
      std::cout << "   Level done." << std::endl;
      }
    }

  if( this->GetReportProgress() )
    {
    std::cout << "BSpline MULTIRESOLUTION END" << std::endl;
    }
}

template <class TImage>
typename BSplineImageToImageRegistrationMethod<TImage>::TransformType
* BSplineImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
const typename BSplineImageToImageRegistrationMethod<TImage>::TransformType
* BSplineImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
typename BSplineImageToImageRegistrationMethod<TImage>::BSplineTransformPointer
BSplineImageToImageRegistrationMethod<TImage>
::GetBSplineTransform( void ) const
{
  typename BSplineTransformType::Pointer trans =
    BSplineTransformType::New();

  trans->SetFixedParameters(
    this->GetTypedTransform()->GetFixedParameters() );
  trans->SetParametersByValue(
    this->GetTypedTransform()->GetParameters() );

  return trans;
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::ResampleControlGrid(int numberOfControlPoints,
  ParametersType & parameters )
{
  typename TransformType::Pointer tmpTrans = TransformType::New();

  typename TransformType::MeshSizeType meshSize;
  meshSize.Fill( numberOfControlPoints - 3 );

  typedef itk::BSplineTransformInitializer< TransformType, ImageType > 
    TransformInitializerType;
  typename TransformInitializerType::Pointer initializer = 
    TransformInitializerType::New();
  initializer->SetTransform( tmpTrans );
  initializer->SetImage( this->GetFixedImage() );
  initializer->SetTransformDomainMeshSize( meshSize );
  initializer->InitializeTransform();
  tmpTrans->SetIdentity();

  int numberOfParameters = tmpTrans->GetParameters().GetSize();
  std::cout << "Number of parameters = " << numberOfParameters << std::endl;

  parameters.SetSize( numberOfParameters );

  int parameterCounter = 0;
  //if( this->GetReportProgress() )
    {
    typename TransformType::OutputPointType p;
    for( unsigned int pdi = 0; pdi<ImageDimension; ++pdi )
      {
      p[pdi] =
        this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[pdi]
        / 2.0;
      }
    p = this->GetTransform()->TransformPoint(p);
    std::cout << "Pre upsample Point = " << p << std::endl;
    }

  typedef typename BSplineTransformType::ImageType
    ParametersImageType;
  typedef ResampleImageFilter<ParametersImageType, ParametersImageType>
    ResamplerType;
  typedef BSplineResampleImageFunction<ParametersImageType, double>
    ResamplerFunctionType;
  typedef IdentityTransform<double, ImageDimension>
    IdentityTransformType;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    typename ResamplerType::Pointer upsampler = ResamplerType::New();

    typename ResamplerFunctionType::Pointer resamplerFunction =
      ResamplerFunctionType::New();
    resamplerFunction->SetSplineOrder(3);

    typename IdentityTransformType::Pointer identity =
      IdentityTransformType::New();

    /*
    if( this->GetReportProgress() )
      {
      /typedef itk::ImageFileWriter<ParametersImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      std::stringstream ss;
      std::string name;
      ss << "inCoeImage" << k << ".mha";
      ss >> name;
      writer->SetInput(
        this->GetTypedTransform()->GetCoefficientImages()[k] );
      writer->SetFileName( name );
      try
        {
        writer->Update();
        }
      catch( ... )
        {
        std::cout << "Error writing coefficient image.  Ignoring."
          << std::endl;
        }
      }
    */

    upsampler->SetInput( this->GetTypedTransform()
      ->GetCoefficientImages()[k] );
    upsampler->SetInterpolator( resamplerFunction );
    upsampler->SetTransform( identity );
    upsampler->SetSize(
      tmpTrans->GetCoefficientImages()[k]->GetLargestPossibleRegion().GetSize() );
    upsampler->SetOutputSpacing(
      tmpTrans->GetCoefficientImages()[k]->GetSpacing() );
    upsampler->SetOutputOrigin(
      tmpTrans->GetCoefficientImages()[k]->GetOrigin() );
    upsampler->SetOutputDirection( 
      tmpTrans->GetCoefficientImages()[k]->GetDirection() );
    try
      {
      upsampler->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Exception in upsampler: " << excep << std::endl;
      }
    catch( ... )
      {
      std::cout << "Uncaught exception in upsampler" << std::endl;
      }

    typedef BSplineDecompositionImageFilter<ParametersImageType,
                                            ParametersImageType>
    DecompositionType;
    typename DecompositionType::Pointer decomposition =
      DecompositionType::New();

    decomposition->SetSplineOrder( 3 );
    decomposition->SetInput( upsampler->GetOutput() );
    try
      {
      decomposition->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cout << "Exception in decomposition: " << excep << std::endl;
      }
    catch( ... )
      {
      std::cout << "Uncaught exception in decomposition" << std::endl;
      }

    typename ParametersImageType::Pointer newCoefficients =
      decomposition->GetOutput();

    /*
    if( this->GetReportProgress() )
      {
      /typedef itk::ImageFileWriter<ParametersImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      std::stringstream ss;
      std::string name;
      ss << "outCoeImage" << k << ".mha";
      ss >> name;
      writer->SetInput( newCoefficients );
      writer->SetFileName( name );
      try
        {
        writer->Update();
        }
      catch( ... )
        {
        std::cout << "Error while writing coefficient image.  Ignoring."
          << std::endl;
        }
      }
    */

    std::cout << "Copying new coefficients for dim " << k << std::endl;
    // copy the coefficients into the parameter array
    typedef ImageRegionIterator<ParametersImageType> Iterator;
    Iterator it( newCoefficients,
                 newCoefficients->GetLargestPossibleRegion() );
    while( !it.IsAtEnd() )
      {
      parameters[parameterCounter++] = it.Get();
      ++it;
      }
    std::cout << "   Set parameters = " << parameterCounter << std::endl;
    }

  //if( this->GetReportProgress() )
    {
    typename TransformType::OutputPointType p;
    for( unsigned int pdi = 0; pdi<ImageDimension; ++pdi )
      {
      p[pdi] =
        this->GetMovingImage()->GetLargestPossibleRegion().GetSize()[pdi]
        / 2.0;
      }
    p = this->GetTransform()->TransformPoint(p);
    std::cout << "Post upsample Point = " << p << std::endl;
    }
}

template <class TImage>
void
BSplineImageToImageRegistrationMethod<TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif
