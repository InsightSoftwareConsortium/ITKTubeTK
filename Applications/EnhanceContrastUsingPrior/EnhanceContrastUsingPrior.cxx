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

#include "itkImage.h"
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
#include "itkFRPROptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNormalizeImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "EnhanceContrastUsingPriorCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

namespace itk {

template< class pixelT, unsigned int dimensionT >
class ContrastCostFunction
: public SingleValuedCostFunction
{
public:

  typedef ContrastCostFunction                    Self;
  typedef SingleValuedCostFunction                Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkTypeMacro( ContrastCostFunction, SingleValuedCostFunction );

  itkNewMacro( Self );

  typedef Superclass::MeasureType         MeasureType;
  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef itk::Image<pixelT, dimensionT>  ImageType;

  typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType, ImageType >
                                                  BlurFilterType;

  unsigned int GetNumberOfParameters( void ) const
    {
    return 3;
    }

  void SetInputImage( typename ImageType::Pointer _inputImage )
    {
    m_InputImage = _inputImage;
    }

  void SetInputMask( typename ImageType::Pointer _maskImage )
    {
    m_InputMask = _maskImage;
    }

  void SetMaskObjectValue( int _maskObjectValue )
    {
    m_MaskObjectValue = _maskObjectValue;
    }

  void SetMaskBackgroundValue( int _maskBackgroundValue )
    {
    m_MaskBackgroundValue = _maskBackgroundValue;
    }

  void SetOutputImage( typename ImageType::Pointer _output )
    {
    m_OutputImage = _output;
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
      tmpP[i] = params[i] - 0.5 / m_Scales[i];
      double tf = this->GetValue( tmpP );
      tmpP[i] = params[i] + 0.5 / m_Scales[i];
      deriv[i] = this->GetValue( tmpP ) - tf;
      tmpP[i] = params[i];
      }
    }

  MeasureType GetValue( const ParametersType & params ) const
    {
    typename BlurFilterType::Pointer filterInputObj = BlurFilterType::New();
    filterInputObj->SetInput( m_InputImage );
    double sigmaObj = params[0];
    if( sigmaObj > 0.3 && sigmaObj < 100 )
      {
      filterInputObj->SetSigma( sigmaObj );
      }
    else
      {
      return 100;
      }
    filterInputObj->Update();
    typename ImageType::Pointer imgObj = filterInputObj->GetOutput();

    typename BlurFilterType::Pointer filterInputBkg =
      BlurFilterType::New();
    filterInputBkg->SetInput( m_InputImage );
    double sigmaBkg = params[1];
    if( sigmaBkg > sigmaObj && sigmaBkg < 100 )
      {
      filterInputBkg->SetSigma( sigmaBkg );
      }
    else
      {
      return 100;
      }
    filterInputBkg->Update();
    typename ImageType::Pointer imgBkg = filterInputBkg->GetOutput();

    double alpha = params[2];

    double meanObj = 0;
    double stdDevObj = 0;
    double countObj = 0;
    double meanBkg = 0;
    double stdDevBkg = 0;
    double countBkg = 0;

    double sumObj = 0;
    double sumsObj = 0;
    double sumBkg = 0;
    double sumsBkg = 0;

    typedef ImageRegionIterator< ImageType >       ImageIteratorType;
    typedef ImageRegionConstIterator< ImageType >  ConstImageIteratorType;

    ConstImageIteratorType iterObj( imgObj,
      imgObj->GetLargestPossibleRegion() );
    ConstImageIteratorType iterBkg( imgBkg,
      imgBkg->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMask( m_InputMask,
      m_InputMask->GetLargestPossibleRegion() );
    ImageIteratorType iterOut( m_OutputImage,
      m_OutputImage->GetLargestPossibleRegion() );

    double meanRawBkg = 0;
    double countRawBkg = 0;
    while( !iterBkg.IsAtEnd() )
      {
      meanRawBkg += iterBkg.Get();
      ++countRawBkg;
      ++iterBkg;
      }
    meanRawBkg /= countRawBkg;

    iterBkg.GoToBegin();
    while( !iterObj.IsAtEnd() )
      {
      double tf = iterObj.Get() * ( 1 + alpha * (iterBkg.Get()-meanRawBkg) );
      if( iterMask.Get() == m_MaskObjectValue )
        {
        sumObj += tf;
        sumsObj += tf * tf;
        ++countObj;
        }
      else if( iterMask.Get() == m_MaskBackgroundValue )
        {
        sumBkg += tf;
        sumsBkg += tf * tf;
        ++countBkg;
        }
      iterOut.Set( tf );
      ++iterObj;
      ++iterBkg;
      ++iterMask;
      ++iterOut;
      }

    if( countObj > 0 )
      {
      meanObj = sumObj/countObj;
      stdDevObj = vcl_sqrt( sumsObj/countObj - meanObj*meanObj );
      }
    if( countBkg > 0 )
      {
      meanBkg = sumBkg/countBkg;
      stdDevBkg = vcl_sqrt( sumsBkg/countBkg - stdDevBkg*stdDevBkg );
      }

    double dp = vnl_math_abs(meanObj - meanBkg) / (stdDevObj * stdDevBkg);

    std::cout << ++m_CallsToGetValue << " : "
              << params[0] << ", " << params[1] << ", "
              << params[2] << ": "
              << meanObj << " (" << stdDevObj << ") "
              << meanBkg << " (" << stdDevBkg << ") "
              << " : " << dp << std::endl;

    return dp;
    }

protected:

  ContrastCostFunction() : m_InputMean(0.0),
                           m_MaskObjectValue(0),
                           m_MaskBackgroundValue(0),
                           m_CallsToGetValue(0) {}
  virtual ~ContrastCostFunction() {}

  void PrintSelf( std::ostream & os, Indent indent ) const
    {
    Superclass::PrintSelf( os, indent );
    }

private:

  ContrastCostFunction( const Self & );
  void operator=( const Self & );

  typename ImageType::Pointer         m_InputImage;
  typename ImageType::Pointer         m_InputMask;
  mutable typename ImageType::Pointer m_OutputImage;

  double                              m_InputMean;

  unsigned int                        m_MaskObjectValue;
  unsigned int                        m_MaskBackgroundValue;

  ParametersType                      m_Scales;

  mutable unsigned int                m_CallsToGetValue;

};

} //namespace itk

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "ContrastImage", CLPProcessInformation );
  progressReporter.Start();

  typedef float                                 PixelType;
  typedef itk::Image< PixelType,  dimensionT >  ImageType;

  /** Read input images */
  typename ImageType::Pointer inputImage;
  typename ImageType::Pointer inputMask;

  timeCollector.Start("Read");
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;

    typename ReaderType::Pointer readerInputImage = ReaderType::New();
    typename ReaderType::Pointer readerInputMask = ReaderType::New();

    //read input image
    readerInputImage->SetFileName( inputImageName.c_str() );
    readerInputMask->SetFileName( inputMaskName.c_str() );

    try
      {
      readerInputImage->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input image. Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    try
      {
      readerInputMask->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading input mask. Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    inputImage = readerInputImage->GetOutput();
    inputMask = readerInputMask->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop("Read");

  //
  // Generate output image
  //
  typename ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->Allocate();
  progressReporter.Report( 0.2 );

  typedef itk::ContrastCostFunction< PixelType, dimensionT >
                                                ContrastCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
  typedef itk::FRPROptimizer                    OptimizerType;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

  ImageIteratorType iter( inputImage,
    inputImage->GetLargestPossibleRegion() );
  double inputMin = iter.Get();
  double inputMax = inputMin;
  while( !iter.IsAtEnd() )
    {
    double tf = iter.Get();
    if( tf < inputMin )
      {
      inputMin = tf;
      }
    else if( tf > inputMax )
      {
      inputMax = tf;
      }
    ++iter;
    }

  itk::Array<double> params(3);
  params[0] = objectScale;
  params[1] = backgroundScale;
  params[2] = 20*(inputMax - inputMin);

  typename ContrastCostFunctionType::Pointer costFunc =
    ContrastCostFunctionType::New();
  costFunc->SetInputImage( inputImage );
  costFunc->SetInputMask( inputMask );
  costFunc->SetOutputImage( outputImage );
  costFunc->SetMaskObjectValue( maskObjectValue );
  costFunc->SetMaskBackgroundValue( maskBackgroundValue );

  InitialOptimizerType::Pointer initOptimizer =
    InitialOptimizerType::New();
  itk::Statistics::NormalVariateGenerator::Pointer normGen =
    itk::Statistics::NormalVariateGenerator::New();
  if( seed > 0 )
    {
    normGen->Initialize( seed );
    }
  initOptimizer->SetNormalVariateGenerator( normGen );
  initOptimizer->Initialize( 1.0 );
  initOptimizer->SetMetricWorstPossibleValue( 101 );
  initOptimizer->SetMaximumIteration( iterations*0.5 );
  initOptimizer->SetMaximize( true );

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetUseUnitLengthGradient( true );
  optimizer->SetMaximumIteration( iterations*0.4 );
  optimizer->SetMaximumLineIteration( iterations*0.2 );
  optimizer->SetStepLength( 0.1 );
  optimizer->SetStepTolerance( 0.001 );
  optimizer->SetValueTolerance( 0.01 );
  optimizer->SetMaximize( true );

  OptimizerType::ScalesType scales( 3 );
  scales[0] = 1.0 / 0.1;
  scales[1] = 1.0 / 2.0;
  scales[2] = 1.0 / (params[2]/10);

  typename ContrastCostFunctionType::ParametersType costFunctionScales( 3 );
  costFunctionScales[0] = scales[0];
  costFunctionScales[1] = scales[1];
  costFunctionScales[2] = scales[2];

  OptimizerType::ScalesType scales2( 3 );
  scales2[0] = scales[0] * scales[0];
  scales2[1] = scales[1] * scales[1];
  scales2[2] = scales[2] * scales[2];

  // OnePlusOne should be passed squared-scales
  initOptimizer->SetScales( scales2 );
  optimizer->SetScales( scales );
  costFunc->SetScales( costFunctionScales );

  initOptimizer->SetCostFunction( costFunc );
  optimizer->SetCostFunction( costFunc );

  progressReporter.Report( 0.25 );

  costFunc->SetOutputImage( outputImage );
  costFunc->Initialize();

  initOptimizer->SetInitialPosition( params );

  progressReporter.Report( 0.3 );

  initOptimizer->StartOptimization();

  progressReporter.Report( 0.5 );

  params = initOptimizer->GetCurrentPosition();
  double result = costFunc->GetValue( params );
  std::cout << "Intermediate params = " << params
            << " Result = " << result << std::endl;

  optimizer->SetInitialPosition( params );
  optimizer->StartOptimization();

  progressReporter.Report( 0.8 );

  params = optimizer->GetCurrentPosition();
  result = costFunc->GetValue( params );
  std::cout << "Winning params = " << params
            << " Result = " << result << std::endl;

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputImageName.c_str() );
  writer->SetInput( outputImage );
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

  return tube::ParseArgsAndCallDoIt( inputImageName, argc, argv );
}
