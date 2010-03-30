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
#include "itkImageRegionIterator.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
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

    typedef BlendCostFunction           Self;
    typedef SingleValuedCostFunction    Superclass;
    typedef SmartPointer< Self >        Pointer;
    typedef SmartPointer< const Self >  ConstPointer;

    itkTypeMacro( BlendCostFunction, SingleValuedCostFunction );

    itkNewMacro( Self );

    typedef Superclass::MeasureType     MeasureType;
    typedef Superclass::ParametersType  ParametersType;
    typedef Superclass::DerivativeType  DerivativeType;

    typedef itk::OrientedImage<pixelT, dimensionT>  ImageType;

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
    void SetMaskMiddle( typename ImageType::Pointer _maskMiddle )
      {
      m_MaskMiddle = _maskMiddle;
      }
    void SetImageOutput( typename ImageType::Pointer _output )
      {
      m_ImageOutput = _output;
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
      static unsigned int calls = 0;

      typedef itk::ImageRegionIterator< ImageType >
        ImageIteratorType;
      typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
        FilterType;

      if( params[3] < 0.2 )
        {
        return 10000;
        }

      typename FilterType::Pointer filterTop = FilterType::New();
      filterTop->SetInput( m_ImageTop );
      filterTop->SetSigma( params[3] );
      filterTop->Update();
      typename ImageType::Pointer imageTopB = filterTop->GetOutput();

      typename FilterType::Pointer filterBottom = FilterType::New();
      filterBottom->SetInput( m_ImageBottom );
      filterBottom->SetSigma( params[3] );
      filterBottom->Update();
      typename ImageType::Pointer imageBottomB = filterBottom->GetOutput();

      ImageIteratorType iterTopB( imageTopB,
                               imageTopB->GetLargestPossibleRegion() );
      ImageIteratorType iterBottomB( imageBottomB,
                               imageBottomB->GetLargestPossibleRegion() );
      ImageIteratorType iterMiddle( m_ImageMiddle,
                               m_ImageMiddle->GetLargestPossibleRegion() );
      ImageIteratorType iterOutput( m_ImageOutput,
                               m_ImageOutput->GetLargestPossibleRegion() );
      while( !iterMiddle.IsAtEnd() )
        {
        float tf = params[0] * iterTopB.Get() +
                   params[1] * iterBottomB.Get() +
                   params[2];

        iterOutput.Set( iterMiddle.Get() - tf );

        ++iterMiddle;
        ++iterTopB;
        ++iterBottomB;
        ++iterOutput;
        }

      typedef itk::IdentityTransform< double, dimensionT > 
        TransformType;
      typename TransformType::Pointer transform = TransformType::New();

      typedef itk::LinearInterpolateImageFunction< ImageType, double > 
        InterpolatorType;
      typename InterpolatorType::Pointer interpolator = 
        InterpolatorType::New();
      interpolator->SetInputImage( m_ImageOutput );

      typename ImageType::SizeType size;
      size = m_ImageMiddle->GetLargestPossibleRegion().GetSize();

      typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType >
                                                             MetricType;
      typename MetricType::Pointer metric = MetricType::New();

      double samplingRate = 0.5;
      metric->SetFixedImage( m_MaskMiddle );
      metric->SetMovingImage( m_ImageOutput );
      metric->SetFixedImageRegion( m_MaskMiddle->
                                     GetLargestPossibleRegion() );
      metric->SetTransform( transform );
      metric->SetInterpolator( interpolator );
      metric->SetNumberOfSpatialSamples( size[0]*size[1]*samplingRate );
      metric->Initialize();
      metric->MultiThreadingInitialize();
    
      double result = metric->GetValue( transform->GetParameters() ) ;

      std::cout << ++calls << ": "
                << params[0] << ", "
                << params[1] << ", "
                << params[2] << ", "
                << params[3] 
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
    typename ImageType::Pointer         m_MaskMiddle;
    mutable typename ImageType::Pointer m_ImageOutput;

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
  typename ImageType::Pointer maskMiddle;

  timeCollector.Start("Read");
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;
  
    typename ReaderType::Pointer readerTop = ReaderType::New();
    typename ReaderType::Pointer readerMiddle = ReaderType::New();
    typename ReaderType::Pointer readerBottom = ReaderType::New();
    typename ReaderType::Pointer readerMask = ReaderType::New();
  
    //read input image  
    readerTop->SetFileName( inputTop.c_str() );
    readerMiddle->SetFileName( inputMiddle.c_str() );
    readerBottom->SetFileName( inputBottom.c_str() );
    readerMask->SetFileName( inputMask.c_str() );
  
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
      readerMask->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading mask. Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
  
    typedef itk::NormalizeImageFilter< ImageType, ImageType > NormFilterType;
    typename NormFilterType::Pointer normTop = NormFilterType::New();
    normTop->SetInput( readerTop->GetOutput() );
    normTop->Update();
    imageTop = normTop->GetOutput();

    typename NormFilterType::Pointer normMiddle = NormFilterType::New();
    normMiddle->SetInput( readerMiddle->GetOutput() );
    normMiddle->Update();
    imageMiddle = normMiddle->GetOutput();

    typename NormFilterType::Pointer normBottom = NormFilterType::New();
    normBottom->SetInput( readerBottom->GetOutput() );
    normBottom->Update();
    imageBottom = normBottom->GetOutput();

    typename NormFilterType::Pointer normMask = NormFilterType::New();
    normMask->SetInput( readerMask->GetOutput() );
    normMask->Update();
    maskMiddle = normMask->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop("Read");

  typedef itk::BlendCostFunction< PixelType, dimensionT >   
                                                BlendCostFunctionType;
  typedef itk::AmoebaOptimizer                    OptimizerType;;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

   itk::Array<double> params(4);
  params[0] = alpha;
  params[1] = beta;
  params[2] = gamma;
  params[3] = sigma;

  typename BlendCostFunctionType::Pointer costFunc = 
    BlendCostFunctionType::New();
  costFunc->SetImageTop( imageTop );
  costFunc->SetImageMiddle( imageMiddle );
  costFunc->SetImageBottom( imageBottom );
  costFunc->SetMaskMiddle( maskMiddle );

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfIterations( iterations );
  optimizer->SetMaximize( false );

  OptimizerType::ScalesType scales( 4 );
  scales[0] = 1.0 / 0.05;
  scales[1] = 1.0 / 0.05;
  scales[2] = 1.0 / 1;
  scales[3] = 1.0 / 0.2;
  optimizer->SetScales( scales );

  optimizer->SetCostFunction( costFunc );
  optimizer->SetInitialPosition( params );

  //
  // Generate output image
  //
  typename ImageType::Pointer imageOutput = ImageType::New();
  imageOutput->CopyInformation( imageMiddle );
  imageOutput->SetRegions( imageMiddle->GetLargestPossibleRegion() );
  imageOutput->Allocate();
  costFunc->SetImageOutput( imageOutput );

  optimizer->StartOptimization();

  params = optimizer->GetCurrentPosition();
  double result = costFunc->GetValue( params );
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
