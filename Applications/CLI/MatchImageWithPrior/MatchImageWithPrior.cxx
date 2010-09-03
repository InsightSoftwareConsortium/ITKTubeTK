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
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Application-specific includes
#include "tubeCompareImageWithPrior.h"
#include "itkOptBrent1D.h"
#include "itkSplineApproximation1D.h"
#include "itkSplineND.h"

#include <map>

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

template< class pixelT, unsigned int dimensionT >
class MyMIWPFunc:
  public itk::UserFunc< vnl_vector<int>, double > 
  {
  public:

    typedef tube::CompareImageWithPrior< pixelT, dimensionT > ImageEvalType;  
    MyMIWPFunc( ImageEvalType & eval )
      {
      cEval = eval;
      cGof = 0;
      };

    const double & value( const vnl_vector<int> & x )
      {
      cEval.SetErode( x[0] );
      cEval.SetDilate( x[1] );
      cEval.SetGaussianBlur( x[2] );
      cEval.Update();
      cGof = cEval.GetGoodnessOfFit();
      return cGof;
      };

  private:

    ImageEvalType cEval;
    double cGof;

  };

// Must include CLP before including tubeCLIHleperFunctions
#include "MatchImageWithPriorCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with Slicer GUI
  tube::CLIProgressReporter    progressReporter( "MatchImageWithPrior",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                         PixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >  ImageType;
  
  /** Read input images */
  typename ImageType::Pointer curVolume;
  typename ImageType::Pointer curMask;

  timeCollector.Start("Read");
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;
  
    typename ReaderType::Pointer readerVolume = ReaderType::New();
    typename ReaderType::Pointer readerMask = ReaderType::New();
  
    //read input image  
    readerVolume->SetFileName( inputVolume.c_str() );
    readerMask->SetFileName( inputMask.c_str() );
  
    try
      {
      readerVolume->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume. Exception caught: " 
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
  
    curVolume = readerVolume->GetOutput();
    curMask = readerMask->GetOutput();
    }
  progressReporter.Report( 0.1 );
  timeCollector.Stop("Read");


  typename ImageType::SizeType inputSize = curVolume
                                           ->GetLargestPossibleRegion()
                                           .GetSize();
  typename ImageType::Pointer orgMask = curMask;

  typename ImageType::Pointer metricMaskImage = NULL;
  if( metricMask.size() != 0 )
    {
    typedef itk::ImageFileReader< ImageType >   ReaderType;
    typename ReaderType::Pointer readerMetricMask = ReaderType::New();
    readerMetricMask->SetFileName( metricMask.c_str() );
    readerMetricMask->Update();
    metricMaskImage = readerMetricMask->GetOutput();
    }
  progressReporter.Report( 0.2 );

  if( foreground != 1 || background != 0 )
    {
    timeCollector.Start("Fg/Bg");

    typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > 
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( curMask );
    if( foreground != 1 )
      {
      filter->SetLowerThreshold( foreground );
      filter->SetUpperThreshold( foreground );
      filter->SetInsideValue( 1 );
      filter->SetOutsideValue( 0 );
      }
    else
      {
      filter->SetLowerThreshold( background );
      filter->SetUpperThreshold( background );
      filter->SetInsideValue( 0 );
      filter->SetOutsideValue( 1 );
      }
    filter->Update();
    curMask = filter->GetOutput();

    timeCollector.Stop("Fg/Bg");
    }
  progressReporter.Report( 0.3 );
  
  int erodeBest = erode;
  int dilateBest = dilate;
  float gaussianBlurBest = gaussianBlur;

  typedef tube::CompareImageWithPrior< pixelT, dimensionT > ImageEvalType;  
  ImageEvalType eval;
  eval.SetVolumeImage( curVolume );
  eval.SetMaskImage( curMask );
  eval.SetOriginalMaskImage( orgMask );
  eval.SetForeground( foreground );
  eval.SetBackground( background );
  eval.SetBoundarySize( outputBoundary );
  eval.SetTimeCollector( &timeCollector );
  eval.SetNormalize( true );
  eval.SetProgressReporter( &progressReporter, 0.3, 0.1 );

  if( metricMaskImage.IsNotNull() )
    {
    eval.SetMetricMask( metricMaskImage );
    }

  typedef typename ImageEvalType::RegistrationMethodType::TransformType  
    TransformType;
  typename TransformType::Pointer regTfm;
  if( loadTransform.size() == 0 )
    {
    eval.SetUseRegistrationTransform( false );
    }
  else
    {
    itk::TransformFileReader::Pointer treader = 
      itk::TransformFileReader::New();
    treader->SetFileName(loadTransform);
    treader->Update();  
    typename TransformType::Pointer transform;
    transform = static_cast< TransformType * >(
      treader->GetTransformList()->front().GetPointer() );

    eval.SetRegistrationTransform( transform );
    eval.SetUseRegistrationTransform( true );
    }

  if( disableRegistrationOptimization )
    {
    eval.SetUseRegistrationOptimization( false );
    if( !eval.GetUseRegistrationTransform() )
      {
      eval.SetUseRegistration( false );
      }
    }

  eval.SetErode( erode );
  eval.SetDilate( dilate );
  eval.SetGaussianBlur( gaussianBlur );

  eval.SetSamplingRate( samplingRate );

  eval.Update();

  if( loadTransform.size() == 0 )
    {
    regTfm = eval.GetRegistrationTransform();
    }

  if( saveTransform.size()>0 )
    {
    itk::TransformFileWriter::Pointer twriter = 
      itk::TransformFileWriter::New();
    twriter->SetInput(regTfm);
    twriter->SetFileName(saveTransform);
    twriter->Update();  
    }

  double gof = eval.GetGoodnessOfFit();

  double gofBest = gof;
  if( !disableParameterOptimization )
    {
    double dilateBase = dilate/2;
    double erodeBase = erode/2;
    eval.SetDilate( dilateBase );
    eval.SetErode( erodeBase );
    eval.SetGaussianBlur( 0 );
    eval.SetUseRegistrationTransform( true );
    eval.SetRegistrationTransform( regTfm );
    eval.SetUseRegistrationOptimization( false );
    eval.SetOriginalMaskImage( orgMask );
    eval.SetBoundarySize( outputBoundary );
    eval.SetNormalize( false );
    eval.Update();

    curVolume = eval.GetVolumeImage();
    curMask = eval.GetMaskImage();

    eval.SetVolumeImage( curVolume );
    eval.SetMaskImage( curMask );

    eval.SetUseRegistration( false );
    eval.SetProgressReporter( &progressReporter, 0.4, 0.5 );

    MyMIWPFunc< pixelT, dimensionT > * myFunc = new 
      MyMIWPFunc< pixelT, dimensionT >( eval );
    itk::SplineApproximation1D * spline1D = new 
      itk::SplineApproximation1D();
    itk::OptBrent1D * opt = new itk::OptBrent1D( );
    itk::SplineND spline( 3, myFunc, spline1D, opt );

    vnl_vector< int > xMin(3);
    xMin.fill( 1 );
    vnl_vector< int > xMax(3);
    xMax.fill( 12 );
    spline.xMin( xMin );
    spline.xMax( xMax );

    vnl_vector< double > x(3);
    x[0] = erodeBase;
    x[1] = dilateBase;
    x[2] = gaussianBlur;

    spline.extreme( x, &gofBest );

    erodeBest = x[0] + erodeBase;
    dilateBest = x[1] + dilateBase;
    gaussianBlurBest = x[2];
    }

  curVolume = eval.GetVolumeImage();
  curMask = eval.GetMaskImage();

  std::cout << "Erode = " << erodeBest << std::endl;
  std::cout << "Dilate = " << dilateBest << std::endl;
  std::cout << "Blur = " << gaussianBlurBest << std::endl;

  progressReporter.Report( 0.9 );

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writerVolume = ImageWriterType::New();
  typename ImageWriterType::Pointer writerMask = ImageWriterType::New();

  writerVolume->SetFileName( outputVolume.c_str() );
  writerVolume->SetInput( curVolume );
  writerVolume->SetUseCompression( true );

  writerMask->SetFileName( outputMask.c_str() );
  writerMask->SetInput( curMask );
  writerMask->SetUseCompression( true );

  try
    {
    writerVolume->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume. Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  
  try
    {
    writerMask->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing mask. Exception caught: " 
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

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
