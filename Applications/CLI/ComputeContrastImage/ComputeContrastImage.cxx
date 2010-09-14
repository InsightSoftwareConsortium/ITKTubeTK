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
#include "ComputeContrastImageCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

namespace itk {

template< class pixelT, unsigned int dimensionT >
class ContrastCostFunction
: public SingleValuedCostFunction
{
public:
  
  typedef ContrastCostFunction                       Self;
  typedef SingleValuedCostFunction                Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  itkTypeMacro( ContrastCostFunction, SingleValuedCostFunction );

  itkNewMacro( Self );

  typedef Superclass::MeasureType                 MeasureType;
  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef itk::OrientedImage<pixelT, dimensionT>  ImageType;

  typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType, ImageType >
                                                  BlurFilterType;

  unsigned int GetNumberOfParameters( void ) const
    { 
    return 2;
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

    typedef itk::ImageRegionIterator< ImageType > 
      ImageIteratorType;
    typedef itk::ImageRegionConstIterator< ImageType > 
      ConstImageIteratorType;

    m_MaskedImage = ImageType::New();
    m_MaskedImage->SetRegions( 
      m_InputImage->GetLargestPossibleRegion() );
    m_MaskedImage->CopyInformation( m_InputImage );
    m_MaskedImage->Allocate();

    double mean = 0;
    unsigned int count = 0;

    std::cout << "Computing background mean" << std::endl;

    ConstImageIteratorType iterInput( m_InputImage,
      m_InputImage->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMask( m_InputMask,
      m_InputMask->GetLargestPossibleRegion() );
    while( !iterInput.IsAtEnd() )
      {
      if( iterMask.Get() == m_MaskBackgroundValue )
        {
        mean += iterInput.Get();
        ++count;
        }
      ++iterInput;
      ++iterMask;
      }

    if( count > 0 )
      {
      mean /= count;
      }

    std::cout << "  Background mean = " << mean << std::endl;

    ImageIteratorType iterMaskedImage( m_MaskedImage,
      m_MaskedImage->GetLargestPossibleRegion() );
    iterMask.GoToBegin();
    while( !iterMask.IsAtEnd() )
      {
      if( iterMask.Get() != m_MaskBackgroundValue )
        {
        iterMaskedImage.Set( mean );
        }
      ++iterMaskedImage;
      ++iterMask;
      }
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
    typename BlurFilterType::Pointer filterInput = BlurFilterType::New();
    filterInput->SetInput( m_InputImage );
    if( params[0] > 0.333 )
      {
      if( params[0] < 10 )
        {
        filterInput->SetSigma( params[0] );
        }
      else
        {
        return 0;
        //filterInput->SetSigma( 60 );
        }
      }
    else
      {
      return 0;
      //filterInput->SetSigma( 0.333 );
      }
    filterInput->Update();

    typename BlurFilterType::Pointer filterMaskedImage =
      BlurFilterType::New();
    filterMaskedImage->SetInput( m_MaskedImage );
    if( params[1] > 0.333 )
      {
      if( params[1] < 60 )
        {
        filterMaskedImage->SetSigma( params[1] );
        }
      else
        {
        return 0;
        //filterMaskedImage->SetSigma( 60 );
        }
      }
    else
      {
      return 0;
      //filterMaskedImage->SetSigma( 0.333 );
      }
    filterMaskedImage->Update();

    typename ImageType::Pointer maskedInputB = 
      filterMaskedImage->GetOutput();
    typename ImageType::Pointer inputB = 
      filterInput->GetOutput();

    typedef itk::ImageRegionConstIterator< ImageType > 
      ConstImageIteratorType;
    typedef itk::ImageRegionIterator< ImageType > 
      ImageIteratorType;

    double result = 0;
    double sumObj = 0;
    double sumBkg = 0;
    double sumsObj = 0;
    double sumsBkg = 0;
    unsigned int countObj = 0;
    unsigned int countBkg = 0;

    ConstImageIteratorType iterInputB( inputB,
      inputB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMaskedImageB( maskedInputB,
      maskedInputB->GetLargestPossibleRegion() );
    ConstImageIteratorType iterMask( m_InputMask,
      m_InputMask->GetLargestPossibleRegion() );
    ImageIteratorType iterOutput( m_OutputImage,
      m_OutputImage->GetLargestPossibleRegion() );
    while( !iterInputB.IsAtEnd() )
      {
      double tf = iterInputB.Get() - iterMaskedImageB.Get();
      iterOutput.Set( tf );
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

      ++iterInputB;
      ++iterMaskedImageB;
      ++iterMask;
      ++iterOutput;
      }

    if( countObj > 0 && countBkg > 0 )
      {
      double meanObj = sumObj/countObj;
      double meanBkg = sumBkg/countBkg;
  
      double stdDevObj = vcl_sqrt( sumsObj/countObj - meanObj*meanObj );
      double stdDevBkg = vcl_sqrt( sumsBkg/countBkg - meanBkg*meanBkg );
  
      double denom = stdDevObj * stdDevBkg;
      if( denom > 0 )
        {
        result = - vnl_math_abs(meanObj - meanBkg) 
          / vcl_sqrt( denom );
        }
      else
        {
        std::cout << "Warning: 0 std dev." << std::endl;
        }
      }

    std::cout << ++m_CallsToGetValue << " : "
              << params[0] << ", "
              << params[1] << ", ";
    std::cout << " : result = " << result << std::endl;

    return result;
    }

protected:
  
  ContrastCostFunction() {};
  virtual ~ContrastCostFunction() {};

  void PrintSelf( std::ostream & os, Indent indent ) const
    {
    Superclass::PrintSelf( os, indent );
    }

private:

  ContrastCostFunction( const Self & );
  void operator=( const Self & );

  typename ImageType::Pointer         m_InputImage;
  typename ImageType::Pointer         m_InputMask;
  typename ImageType::Pointer         m_MaskedImage;
  mutable typename ImageType::Pointer m_OutputImage;

  unsigned int                        m_MaskObjectValue;
  unsigned int                        m_MaskBackgroundValue;

  ParametersType                      m_Scales;

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
  tube::CLIProgressReporter    progressReporter( "ComputeContrastImage",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                         PixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >  ImageType;
  
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

  itk::Array<double> params(2);
  params[0] = objectScale;
  params[1] = backgroundScale;

  typedef itk::ContrastCostFunction< PixelType, dimensionT >   
                                                ContrastCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
  typedef itk::FRPROptimizer                    OptimizerType;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;
 
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
 
  OptimizerType::ScalesType scales( 2 );
  scales[0] = 1.0 / 0.5;
  scales[1] = 1.0 / 1.0;
 
  OptimizerType::ScalesType scales2( 2 );
  scales2[0] = scales[0] * scales[0];
  scales2[1] = scales[1] * scales[1];
 
  // OnePlusOne should be passed squared-scales
  initOptimizer->SetScales( scales2 );
  optimizer->SetScales( scales );
  costFunc->SetScales( scales );
 
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
