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

#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkFRPROptimizer.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itktubeSmoothingRecursiveGaussianImageFilter.h>
#include <itkTimeProbesCollectorBase.h>

#include "DeblendTomosynthesisSlicesUsingPriorCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

namespace itk
{

namespace tube
{

template< class TPixel, unsigned int VDimension >
class BlendCostFunction : public SingleValuedCostFunction
{
public:

  typedef BlendCostFunction                       Self;
  typedef SingleValuedCostFunction                Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkTypeMacro( BlendCostFunction, SingleValuedCostFunction );

  itkNewMacro( Self );

  typedef Superclass::MeasureType         MeasureType;
  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef itk::Image<TPixel, VDimension>  ImageType;

  unsigned int GetNumberOfParameters( void ) const
    {
    return 3;
    }

  void SetImageTop( typename ImageType::Pointer _top )
    {
    m_ImageTop = _top;
    }

  void SetImageMiddle( typename ImageType::Pointer _middle )
    {
    m_ImageMiddle = _middle;
    }

  void SetImageBottom( typename ImageType::Pointer _bottom )
    {
    m_ImageBottom = _bottom;
    }

  void SetImageMiddleTarget( typename ImageType::Pointer _targetMiddle )
    {
    m_ImageMiddleTarget = _targetMiddle;
    }

  void SetMetricMask( typename ImageType::Pointer _MetricMask )
    {
    m_MetricMask = _MetricMask;
    }

  void SetImageOutput( typename ImageType::Pointer _output )
    {
    m_ImageOutput = _output;
    }

  void SetScales( ParametersType & _scales )
    {
    m_Scales = _scales;
    }

  void Initialize( void )
    {
    m_CallsToGetValue = 0;
    }

  void GetDerivative( const ParametersType & params,
                      DerivativeType & deriv ) const
    {
    ParametersType tmpP = params;
    deriv = params;

    for( unsigned int i=0; i<this->GetNumberOfParameters(); i++ )
      {
      tmpP[i] = params[i] - 1.0 / m_Scales[i];
      double tf = this->GetValue( tmpP );
      tmpP[i] = params[i] + 1.0 / m_Scales[i];
      deriv[i] = this->GetValue( tmpP ) - tf;
      tmpP[i] = params[i];
      }
    }

  MeasureType GetValue( const ParametersType & params ) const
    {
    typedef itk::ImageRegionConstIterator< ImageType >
      ConstImageIteratorType;
    typedef itk::ImageRegionIterator< ImageType >
      ImageIteratorType;

    double result = 0;
    double sum255 = 0;
    double sumNot = 0;
    double sums255 = 0;
    double sumsNot = 0;
    unsigned int count255 = 0;
    unsigned int countNot = 0;

    ConstImageIteratorType iterBottomB( m_ImageBottom,
      m_ImageBottom->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMiddle( m_ImageMiddle,
      m_ImageMiddle->GetLargestPossibleRegion() );
    ConstImageIteratorType iterTopB( m_ImageTop,
      m_ImageTop->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMiddleTarget( m_ImageMiddleTarget,
      m_ImageMiddleTarget->GetLargestPossibleRegion() );
    ImageIteratorType iterOutput( m_ImageOutput,
      m_ImageOutput->GetLargestPossibleRegion() );

    typename ImageType::Pointer tmpMask = m_ImageOutput;
    if( m_MetricMask.IsNotNull() )
      {
      tmpMask = m_MetricMask;
      }
    ConstImageIteratorType iterMask( tmpMask,
      tmpMask->GetLargestPossibleRegion() );
    while( !iterMiddle.IsAtEnd() )
      {
      float tf = ( params[0] * iterBottomB.Get() +
        iterMiddle.Get() +
        params[1] * iterTopB.Get() )
        + params[2];

      iterOutput.Set( tf );

      if( m_MetricMask.IsNull() )
        {
        double diff = iterMiddleTarget.Get() - tf;
        result += diff * diff;
        }
      else if( iterMask.Get() != 0 )
        {
        double diff = iterMiddleTarget.Get() - tf;
        result += diff * diff;
        if( iterMask.Get() == 255 )
          {
          sum255 += tf;
          sums255 += tf * tf;
          ++count255;
          }
        else
          {
          sumNot += tf;
          sumsNot += tf * tf;
          ++countNot;
          }
        }

      ++iterBottomB;
      ++iterMiddle;
      ++iterTopB;
      ++iterMiddleTarget;
      ++iterMask;
      ++iterOutput;
      }

    if( count255 > 0 && countNot > 0 )
      {
      double mean255 = sum255/count255;
      double meanNot = sumNot/countNot;

      double stdDev255 = std::sqrt( sums255/count255 - mean255*mean255 );
      double stdDevNot = std::sqrt( sumsNot/countNot - meanNot*meanNot );

      result = - vnl_math_abs( mean255 - meanNot )
        / std::sqrt( stdDev255 * stdDevNot );
      }

    std::cout << ++m_CallsToGetValue << " : "
              << params[0] << ", "
              << params[1] << ", "
              << params[2] << ", ";
    std::cout << " : result = " << result << std::endl;

    return result;
    }

protected:

  BlendCostFunction( void )
    : m_Mode( 0 ), m_CallsToGetValue( 0 ) {}
  virtual ~BlendCostFunction( void ) {}

  void PrintSelf( std::ostream & os, Indent indent ) const
    {
    Superclass::PrintSelf( os, indent );
    }

private:

  BlendCostFunction( const Self & );
  void operator=( const Self & );

  unsigned int                        m_Mode;

  typename ImageType::Pointer         m_ImageTop;
  typename ImageType::Pointer         m_ImageMiddle;
  typename ImageType::Pointer         m_ImageBottom;
  typename ImageType::Pointer         m_ImageMiddleTarget;
  typename ImageType::Pointer         m_MetricMask;
  mutable typename ImageType::Pointer m_ImageOutput;

  ParametersType                      m_Scales;

  mutable unsigned int                m_CallsToGetValue;

}; // End class BlendCostFunction

template< class TPixel, unsigned int VDimension >
class BlendScaleCostFunction : public SingleValuedCostFunction
{
public:

  typedef BlendScaleCostFunction                  Self;
  typedef SingleValuedCostFunction                Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkTypeMacro( BlendScaleCostFunction, SingleValuedCostFunction );

  itkNewMacro( Self );

  typedef Superclass::MeasureType         MeasureType;
  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef itk::Image<TPixel, VDimension>  ImageType;

  typedef itk::tube::SmoothingRecursiveGaussianImageFilter< ImageType,
    ImageType >                           BlurFilterType;

  unsigned int GetNumberOfParameters( void ) const
    {
    return 4;
    }

  void SetImageTop( typename ImageType::Pointer _top )
    {
    m_ImageTop = _top;
    }

  void SetImageMiddle( typename ImageType::Pointer _middle )
    {
    m_ImageMiddle = _middle;
    }

  void SetImageBottom( typename ImageType::Pointer _bottom )
    {
    m_ImageBottom = _bottom;
    }

  void SetImageMiddleTarget( typename ImageType::Pointer _targetMiddle )
    {
    m_ImageMiddleTarget = _targetMiddle;
    }

  void SetMetricMask( typename ImageType::Pointer _MetricMask )
    {
    m_MetricMask = _MetricMask;
    }

  void SetImageOutput( typename ImageType::Pointer _output )
    {
    m_ImageOutput = _output;
    }

  void SetScales( ParametersType & _scales )
    {
    m_Scales = _scales;
    }

  void Initialize( void )
    {
    m_CallsToGetValue = 0;
    }

  void GetDerivative( const ParametersType & params,
                      DerivativeType & deriv ) const
    {
    ParametersType tmpP = params;
    deriv = params;

    for( unsigned int i=0; i<this->GetNumberOfParameters(); i++ )
      {
      tmpP[i] = params[i] - 1.0 / m_Scales[i];
      double tf = this->GetValue( tmpP );
      tmpP[i] = params[i] + 1.0 / m_Scales[i];
      deriv[i] = this->GetValue( tmpP ) - tf;
      tmpP[i] = params[i];
      }
    }

  MeasureType GetValue( const ParametersType & params ) const
    {
    typedef itk::ImageRegionConstIterator< ImageType >
      ConstImageIteratorType;
    typedef itk::ImageRegionIterator< ImageType >
      ImageIteratorType;

    typename BlurFilterType::Pointer filterBottom = BlurFilterType::New();
    filterBottom->SetInput( m_ImageBottom );
    typename BlurFilterType::Pointer filterTop = BlurFilterType::New();
    filterTop->SetInput( m_ImageTop );
    if( params[3] > 0.333 )
      {
      filterBottom->SetSigma( params[3] );
      filterTop->SetSigma( params[3] );
      }
    else
      {
      filterBottom->SetSigma( 0.333 );
      filterTop->SetSigma( 0.333 );
      }
    filterBottom->Update();
    typename ImageType::Pointer imageBottomB = filterBottom->GetOutput();
    filterTop->Update();
    typename ImageType::Pointer imageTopB = filterTop->GetOutput();

    double result = 0;
    double sum255 = 0;
    double sumNot = 0;
    double sums255 = 0;
    double sumsNot = 0;
    unsigned int count255 = 0;
    unsigned int countNot = 0;

    ConstImageIteratorType iterBottomB( imageBottomB,
      imageBottomB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMiddle( m_ImageMiddle,
      m_ImageMiddle->GetLargestPossibleRegion() );
    ConstImageIteratorType iterTopB( imageTopB,
      imageTopB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMiddleTarget( m_ImageMiddleTarget,
      m_ImageMiddleTarget->GetLargestPossibleRegion() );
    ImageIteratorType iterOutput( m_ImageOutput,
      m_ImageOutput->GetLargestPossibleRegion() );

    typename ImageType::Pointer tmpMask = m_ImageOutput;
    if( m_MetricMask.IsNotNull() )
      {
      tmpMask = m_MetricMask;
      }
    ConstImageIteratorType iterMask( tmpMask,
      tmpMask->GetLargestPossibleRegion() );
    while( !iterMiddle.IsAtEnd() )
      {
      float tf = ( params[0] * iterBottomB.Get() +
        iterMiddle.Get() +
        params[1] * iterTopB.Get() )
        + params[2];

      iterOutput.Set( tf );

      if( m_MetricMask.IsNull() )
        {
        double diff = iterMiddleTarget.Get() - tf;
        result += diff * diff;
        }
      else if( iterMask.Get() != 0 )
        {
        double diff = iterMiddleTarget.Get() - tf;
        result += diff * diff;
        if( iterMask.Get() == 255 )
          {
          sum255 += tf;
          sums255 += tf * tf;
          ++count255;
          }
        else
          {
          sumNot += tf;
          sumsNot += tf * tf;
          ++countNot;
          }
        }

      ++iterBottomB;
      ++iterMiddle;
      ++iterTopB;
      ++iterMiddleTarget;
      ++iterMask;
      ++iterOutput;
      }

    if( count255 > 0 && countNot > 0 )
      {
      double mean255 = sum255/count255;
      double meanNot = sumNot/countNot;

      double stdDev255 = std::sqrt( sums255/count255 - mean255*mean255 );
      double stdDevNot = std::sqrt( sumsNot/countNot - meanNot*meanNot );

      result = - vnl_math_abs( mean255 - meanNot )
        / std::sqrt( stdDev255 * stdDevNot );
      }

    std::cout << ++m_CallsToGetValue << " : "
              << params[0] << ", "
              << params[1] << ", "
              << params[2] << ", "
              << params[3];
    std::cout << " : result = " << result << std::endl;

    return result;
    }

protected:

  BlendScaleCostFunction( void )
    : m_Mode( 0 ), m_CallsToGetValue( 0 ) {}
  virtual ~BlendScaleCostFunction( void ) {}

  void PrintSelf( std::ostream & os, Indent indent ) const
    {
    Superclass::PrintSelf( os, indent );
    }

private:

  BlendScaleCostFunction( const Self & );
  void operator=( const Self & );

  unsigned int                        m_Mode;

  typename ImageType::Pointer         m_ImageTop;
  typename ImageType::Pointer         m_ImageMiddle;
  typename ImageType::Pointer         m_ImageBottom;
  typename ImageType::Pointer         m_ImageMiddleTarget;
  typename ImageType::Pointer         m_MetricMask;
  mutable typename ImageType::Pointer m_ImageOutput;

  ParametersType                      m_Scales;

  mutable unsigned int                m_CallsToGetValue;

}; // End class BlendScaleCostFunction

} // End namespace tube

} // End namespace itk

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter( "DeblendImages",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                PixelType;
  typedef itk::Image< PixelType, VDimension >  ImageType;

  /** Read input images */
  typename ImageType::Pointer imageBottom;
  typename ImageType::Pointer imageMiddle;
  typename ImageType::Pointer imageTop;
  typename ImageType::Pointer imageMiddleTarget;

  timeCollector.Start( "Read" );
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    typename ReaderType::Pointer readerBottom = ReaderType::New();
    typename ReaderType::Pointer readerMiddle = ReaderType::New();
    typename ReaderType::Pointer readerTop = ReaderType::New();
    typename ReaderType::Pointer readerMiddleTarget = ReaderType::New();

    //read input image
    readerBottom->SetFileName( inputBottom.c_str() );
    readerMiddle->SetFileName( inputMiddle.c_str() );
    readerTop->SetFileName( inputTop.c_str() );
    readerMiddleTarget->SetFileName( inputMiddleTarget.c_str() );

    try
      {
      readerBottom->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading bottom. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerMiddle->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading middle. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerTop->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading top. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerMiddleTarget->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading mask. Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    imageBottom = readerBottom->GetOutput();
    imageMiddle = readerMiddle->GetOutput();
    imageTop = readerTop->GetOutput();
    imageMiddleTarget = readerMiddleTarget->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop( "Read" );

  //
  // Generate output image
  //
  typename ImageType::Pointer imageOutput = ImageType::New();
  imageOutput->CopyInformation( imageMiddle );
  imageOutput->SetRegions( imageMiddle->GetLargestPossibleRegion() );
  imageOutput->Allocate();

  itk::Array<double> blendParams( 3 );
  blendParams[0] = alpha;
  blendParams[1] = gamma;
  blendParams[2] = offset;
  if( 1 )
    {
    typedef itk::tube::BlendCostFunction< PixelType, VDimension >
                                                  BlendCostFunctionType;
    typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
    typedef itk::FRPROptimizer                    OptimizerType;

    typename BlendCostFunctionType::Pointer costFunc =
      BlendCostFunctionType::New();
    costFunc->SetImageBottom( imageBottom );
    costFunc->SetImageMiddle( imageMiddle );
    costFunc->SetImageTop( imageTop );
    costFunc->SetImageMiddleTarget( imageMiddleTarget );

    if( metricMask.size() > 0 )
      {
      typedef itk::ImageFileReader< ImageType >   ReaderType;
      typename ReaderType::Pointer readerMask = ReaderType::New();
      readerMask->SetFileName( metricMask.c_str() );
      readerMask->Update();
      costFunc->SetMetricMask( readerMask->GetOutput() );
      }

    InitialOptimizerType::Pointer initOptimizer = InitialOptimizerType::New();
    itk::Statistics::NormalVariateGenerator::Pointer normGen =
      itk::Statistics::NormalVariateGenerator::New();
    if( seed != 0 )
      {
      normGen->Initialize( seed );
      }
    initOptimizer->SetNormalVariateGenerator( normGen );
    initOptimizer->Initialize( 0.1 );
    initOptimizer->SetMetricWorstPossibleValue( 100 );
    initOptimizer->SetMaximumIteration( iterations*0.5 );
    initOptimizer->SetMaximize( false );

    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetUseUnitLengthGradient( true );
    optimizer->SetMaximumIteration( iterations*0.4 );
    optimizer->SetMaximumLineIteration( iterations*0.2 );
    optimizer->SetStepLength( 0.1 );
    optimizer->SetStepTolerance( 0.001 );
    optimizer->SetValueTolerance( 1 );
    optimizer->SetMaximize( false );

    OptimizerType::ScalesType blendScales( 3 );
    blendScales[0] = 1.0 / 0.1;
    blendScales[1] = 1.0 / 0.1;
    blendScales[2] = 1.0 / 0.1;

    typename BlendCostFunctionType::ParametersType costFunctionScales( 3 );
    costFunctionScales[0] = blendScales[0];
    costFunctionScales[1] = blendScales[1];
    costFunctionScales[2] = blendScales[2];

    OptimizerType::ScalesType blendScales2( 3 );
    blendScales2[0] = blendScales[0] * blendScales[0];
    blendScales2[1] = blendScales[1] * blendScales[1];
    blendScales2[2] = blendScales[2] * blendScales[2];

    // OnePlusOne should be passed squared-scales
    initOptimizer->SetScales( blendScales2 );
    optimizer->SetScales( blendScales );
    costFunc->SetScales( costFunctionScales );

    initOptimizer->SetCostFunction( costFunc );
    optimizer->SetCostFunction( costFunc );

    costFunc->SetImageOutput( imageOutput );
    costFunc->Initialize();

    initOptimizer->SetInitialPosition( blendParams );
    if( iterations > 1 )
      {
      initOptimizer->StartOptimization();
      blendParams = initOptimizer->GetCurrentPosition();
      }
    double result = costFunc->GetValue( blendParams );
    std::cout << "Intermediate blendParams = " << blendParams
              << " Result = " << result << std::endl;

    optimizer->SetInitialPosition( blendParams );
    if( iterations > 1 )
      {
      optimizer->StartOptimization();
      blendParams = optimizer->GetCurrentPosition();
      }
    result = costFunc->GetValue( blendParams );
    std::cout << "Winning blendParams = " << blendParams
              << " Result = " << result << std::endl;
    }

  itk::Array<double> blendScaleParams( 4 );
  blendScaleParams[0] = blendParams[0];
  blendScaleParams[1] = blendParams[1];
  blendScaleParams[2] = blendParams[2];
  blendScaleParams[3] = sigma;
  if( 1 )
    {
    typedef itk::tube::BlendScaleCostFunction< PixelType, VDimension >
                                                  BlendScaleCostFunctionType;
    typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
    typedef itk::FRPROptimizer                    OptimizerType;

    typename BlendScaleCostFunctionType::Pointer costFunc =
      BlendScaleCostFunctionType::New();
    costFunc->SetImageBottom( imageBottom );
    costFunc->SetImageMiddle( imageMiddle );
    costFunc->SetImageTop( imageTop );
    costFunc->SetImageMiddleTarget( imageMiddleTarget );

    if( metricMask.size() > 0 )
      {
      typedef itk::ImageFileReader< ImageType >   ReaderType;
      typename ReaderType::Pointer readerMask = ReaderType::New();
      readerMask->SetFileName( metricMask.c_str() );
      readerMask->Update();
      costFunc->SetMetricMask( readerMask->GetOutput() );
      }

    InitialOptimizerType::Pointer initOptimizer = InitialOptimizerType::New();
    itk::Statistics::NormalVariateGenerator::Pointer normGen =
      itk::Statistics::NormalVariateGenerator::New();
    if( seed != 0 )
      {
      normGen->Initialize( seed );
      }
    initOptimizer->SetNormalVariateGenerator( normGen );
    initOptimizer->Initialize( 0.1 );
    initOptimizer->SetMetricWorstPossibleValue( 100 );
    initOptimizer->SetMaximumIteration( iterations*0.5 );
    initOptimizer->SetMaximize( false );

    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetUseUnitLengthGradient( true );
    optimizer->SetMaximumIteration( iterations*0.8 );
    optimizer->SetMaximumLineIteration( iterations*0.2 );
    optimizer->SetStepLength( 0.5 );
    optimizer->SetStepTolerance( 0.001 );
    optimizer->SetValueTolerance( 1 );
    optimizer->SetMaximize( false );

    OptimizerType::ScalesType blendScaleScales( 4 );
    blendScaleScales[0] = 1.0 / 0.01;
    blendScaleScales[1] = 1.0 / 0.01;
    blendScaleScales[2] = 1.0 / 0.01;
    blendScaleScales[3] = 1.0 / 0.5;

    typename BlendScaleCostFunctionType::ParametersType
      costFunctionScaleScales( 4 );
    costFunctionScaleScales[0] = blendScaleScales[0];
    costFunctionScaleScales[1] = blendScaleScales[1];
    costFunctionScaleScales[2] = blendScaleScales[2];
    costFunctionScaleScales[3] = blendScaleScales[3];

    OptimizerType::ScalesType blendScaleScales2( 4 );
    blendScaleScales2[0] = blendScaleScales[0] * blendScaleScales[0];
    blendScaleScales2[1] = blendScaleScales[1] * blendScaleScales[1];
    blendScaleScales2[2] = blendScaleScales[2] * blendScaleScales[2];
    blendScaleScales2[3] = blendScaleScales[3] * blendScaleScales[3];

    // OnePlusOne should be passed squared-scales
    initOptimizer->SetScales( blendScaleScales2 );
    optimizer->SetScales( blendScaleScales );
    costFunc->SetScales( costFunctionScaleScales );

    initOptimizer->SetCostFunction( costFunc );
    optimizer->SetCostFunction( costFunc );

    costFunc->SetImageOutput( imageOutput );
    costFunc->Initialize();

    //initOptimizer->SetInitialPosition( blendScaleParams );
    //initOptimizer->StartOptimization();

    //blendScaleParams = initOptimizer->GetCurrentPosition();
    double result = costFunc->GetValue( blendScaleParams );
    std::cout << "Intermediate blendScaleParams = " << blendScaleParams
              << " Result = " << result << std::endl;

    optimizer->SetInitialPosition( blendScaleParams );
    if( iterations > 1 )
      {
      optimizer->StartOptimization();
      blendScaleParams = optimizer->GetCurrentPosition();
      }
    result = costFunc->GetValue( blendScaleParams );
    std::cout << "Winning blendScaleParams = " << blendScaleParams
              << " Result = " << result << std::endl;
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputMiddle.c_str() );
  writer->SetInput( imageOutput );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume. Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  progressReporter.Report( 1.0 );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputMiddle, argc, argv );
}
