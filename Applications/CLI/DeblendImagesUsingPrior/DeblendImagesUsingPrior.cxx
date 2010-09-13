/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkOrientedImage.h"
#include "itkImageSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Application-specific includes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAmoebaOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNormalizeImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "DeblendImagesUsingPriorCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

namespace itk {

template< class pixelT, unsigned int dimensionT >
class BlendCostFunction
: public SingleValuedCostFunction
{
public:
  
  typedef BlendCostFunction                       Self;
  typedef SingleValuedCostFunction                Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkTypeMacro( BlendCostFunction, SingleValuedCostFunction );

  itkNewMacro( Self );

  typedef Superclass::MeasureType                 MeasureType;
  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef itk::OrientedImage<pixelT, dimensionT>  ImageType;

  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
                                                  BlurFilterType;
  
  unsigned int GetNumberOfParameters( void ) const
    { 
    return 6;
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
  
  void Initialize( void )
    {
    m_FilterTop = BlurFilterType::New();
    m_FilterTop->SetInput( m_ImageTop );

    m_FilterBottom = BlurFilterType::New();
    m_FilterBottom->SetInput( m_ImageBottom );

    m_CallsToGetValue = 0;
    }

  void GetDerivative( const ParametersType & params,
                      DerivativeType & deriv ) const
    {
    ParametersType unusedP = params;
    DerivativeType unusedD = deriv;
    return;
    }

  MeasureType GetValue( const ParametersType & params ) const
    {
    typedef itk::ImageRegionConstIterator< ImageType > 
      ConstImageIteratorType;
    typedef itk::ImageRegionIterator< ImageType > 
      ImageIteratorType;

    if( params[4] < 0.2 || params[5] < 0.2 )
      {
      return 100;
      }
    
    m_FilterTop->SetSigma( params[4] );
    m_FilterTop->Update();

    m_FilterBottom->SetSigma( params[5] );
    m_FilterBottom->Update();

    typename ImageType::Pointer imageTopB = m_FilterTop->GetOutput();
    typename ImageType::Pointer imageBottomB = m_FilterBottom->GetOutput();
  
    double result = 0;

    ConstImageIteratorType iterTopB( imageTopB,
      imageTopB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterBottomB( imageBottomB,
      imageBottomB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMiddle( m_ImageMiddle,
      m_ImageMiddle->GetLargestPossibleRegion() );
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
      float tf = ( ( params[0] * iterTopB.Get() +
        params[1] * iterBottomB.Get() +
        params[2] * iterMiddle.Get() ) /
        ( params[0] + params[1] + params[2] ) ) +
        params[3];

      iterOutput.Set( tf );

      if( iterMask.Get() )
        {
        double diff = iterMiddleTarget.Get() - tf;
        result += diff * diff;
        }

      ++iterMask;
      ++iterMiddle;
      ++iterMiddleTarget;
      ++iterTopB;
      ++iterBottomB;
      ++iterOutput;
      }

    std::cout << ++m_CallsToGetValue << ": "
              << params[0] << ", "
              << params[1] << ", "
              << params[2] << ", "
              << params[3] << ", "
              << params[4] << ", "
              << params[5] 
              << " : result =" << result << std::endl;

    return result;
    }

protected:
  
  BlendCostFunction() {};
  virtual ~BlendCostFunction() {};

  void PrintSelf( std::ostream & os, Indent indent ) const
    {
    Superclass::PrintSelf( os, indent );
    }

private:

  BlendCostFunction( const Self & );
  void operator=( const Self & );

  typename ImageType::Pointer         m_ImageTop;
  typename ImageType::Pointer         m_ImageMiddle;
  typename ImageType::Pointer         m_ImageBottom;
  typename ImageType::Pointer         m_ImageMiddleTarget;
  typename ImageType::Pointer         m_MetricMask;
  mutable typename ImageType::Pointer m_ImageOutput;

  typename BlurFilterType::Pointer    m_FilterTop;
  typename BlurFilterType::Pointer    m_FilterBottom;
  
  mutable unsigned int                m_CallsToGetValue;

};

}; //namespace itk


template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter( "DeblendImagesUsingPrior",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                         PixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >  ImageType;
  
  /** Read input images */
  typename ImageType::Pointer imageTop;
  typename ImageType::Pointer imageMiddle;
  typename ImageType::Pointer imageBottom;
  typename ImageType::Pointer imageMiddleTarget;

  timeCollector.Start("Read");
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;
  
    typename ReaderType::Pointer readerTop = ReaderType::New();
    typename ReaderType::Pointer readerMiddle = ReaderType::New();
    typename ReaderType::Pointer readerBottom = ReaderType::New();
    typename ReaderType::Pointer readerMiddleTarget = ReaderType::New();
  
    //read input image  
    readerTop->SetFileName( inputTop.c_str() );
    readerMiddle->SetFileName( inputMiddle.c_str() );
    readerBottom->SetFileName( inputBottom.c_str() );
    readerMiddleTarget->SetFileName( inputMiddleTarget.c_str() );
  
    try
      {
      readerTop->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading top. Exception caught: " 
                          + std::string(err.GetDescription()) );
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
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerBottom->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading bottom. Exception caught: " 
                          + std::string(err.GetDescription()) );
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
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
  
    imageTop = readerTop->GetOutput();
    imageMiddle = readerMiddle->GetOutput();
    imageBottom = readerBottom->GetOutput();
    imageMiddleTarget = readerMiddleTarget->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop("Read");

  typedef itk::BlendCostFunction< PixelType, dimensionT >   
                                                BlendCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
  typedef itk::AmoebaOptimizer                  OptimizerType;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

  itk::Array<double> params(6);
  params[0] = alpha;
  params[1] = beta;
  params[2] = gamma;
  params[3] = offset;
  params[4] = sigma;
  params[5] = sigma;

  typename BlendCostFunctionType::Pointer costFunc = 
    BlendCostFunctionType::New();
  costFunc->SetImageTop( imageTop );
  costFunc->SetImageMiddle( imageMiddle );
  costFunc->SetImageBottom( imageBottom );
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
  optimizer->SetMaximumNumberOfIterations( iterations*0.5 );
  optimizer->SetMaximize( false );

  OptimizerType::ScalesType scales( 6 );
  scales[0] = 1.0 / 0.5;
  scales[1] = 1.0 / 0.5;
  scales[2] = 1.0 / 0.5;
  scales[3] = 1.0 / 0.5;
  scales[4] = 1.0 / 2.0;
  scales[5] = 1.0 / 2.0;


  OptimizerType::ScalesType scales2( 6 );
  scales2[0] = scales[0] * scales[0];
  scales2[1] = scales[1] * scales[1];
  scales2[2] = scales[2] * scales[2];
  scales2[3] = scales[3] * scales[3];
  scales2[4] = scales[4] * scales[4];
  scales2[5] = scales[5] * scales[5];

  // OnePlusOne should be passed squared-scales
  initOptimizer->SetScales( scales2 );
  optimizer->SetScales( scales2 );

  initOptimizer->SetCostFunction( costFunc );

  optimizer->SetScales( scales );
  optimizer->SetCostFunction( costFunc );

  //
  // Generate output image
  //
  typename ImageType::Pointer imageOutput = ImageType::New();
  imageOutput->CopyInformation( imageMiddle );
  imageOutput->SetRegions( imageMiddle->GetLargestPossibleRegion() );
  imageOutput->Allocate();
  costFunc->SetImageOutput( imageOutput );
  costFunc->Initialize();

  initOptimizer->SetInitialPosition( params );
  initOptimizer->StartOptimization();

  params = initOptimizer->GetCurrentPosition();
  double result = costFunc->GetValue( params );

  std::cout << "Intermediate params = " << params 
            << " Result = " << result << std::endl;
  optimizer->SetInitialPosition( initOptimizer->GetCurrentPosition() );
  optimizer->StartOptimization();

  params = optimizer->GetCurrentPosition();
  result = costFunc->GetValue( params );
  std::cout << "Winning params = " << params 
            << " Result = " << result << std::endl;

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
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  
  progressReporter.Report( 1.0 );
  progressReporter.End( );

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char **argv )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputMiddle, argc, argv );
}
