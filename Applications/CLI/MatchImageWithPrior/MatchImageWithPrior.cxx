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
#include "tubeCompareImageWithPrior.h"

#include <map>

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

class triple 
{
public:
  int erode;
  int dilate;
  float gaussianBlur;
  bool operator<(const triple &other) const 
    {
    return (this->erode < other.erode 
            && this->dilate < other.dilate 
            && this->gaussianBlur < other.gaussianBlur);
    }
  bool operator==(const triple &other) const 
    {
    return (this->erode == other.erode 
            && this->dilate == other.dilate 
            && this->gaussianBlur == other.gaussianBlur);
    }
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
  
  if( true )
    {
    typedef std::map< triple, float >     TrialListType;

    TrialListType trials;
    triple        curTrial;

    int erodeStep = 1;
    int dilateStep = 1;
    float gaussianBlurStep = 1;
    int erodeBest = erode;
    int dilateBest = dilate;
    float gaussianBlurBest = gaussianBlur;

    typedef tube::CompareImageWithPrior< pixelT, dimensionT > ImageEvalType;  
    ImageEvalType eval;
    eval.SetVolumeImage( curVolume );
    eval.SetMaskImage( curMask );
    eval.SetOriginalMaskImage( orgMask );
    eval.SetForeground( foreground );
    eval.SetBoundarySize( outputBoundary );
    eval.SetTimeCollector( &timeCollector );
    eval.SetProgressReporter( &progressReporter, 0.3, 0.1 );
    eval.SetUseRegistration( true );
    eval.SetUseRegistrationTransform( false );
    eval.SetErode( erode );
    eval.SetDilate( dilate );
    eval.SetGaussianBlur( gaussianBlur );
    eval.Update();

    typedef typename ImageEvalType::RegistrationMethodType::TransformType  
      TransformType;
    typename TransformType::Pointer regTfm;
    regTfm = eval.GetRegistrationTransform();

    double gof = eval.GetGoodnessOfFit();

    curTrial.erode = erode;
    curTrial.dilate = dilate;
    curTrial.gaussianBlur = gaussianBlur;
    trials[curTrial] = gof;

    double gofBest = gof;
    bool erodeFlip = false;
    bool dilateFlip = false;
    bool gaussianBlurFlip = false;
    bool done = false;
    std::cout << "params = " << erode << ", " 
                             << dilate << ", " 
                             << gaussianBlur
              << "  val = " << gof << std::endl;
    while( !done )
      {
      if( erode+erodeStep < 0 )
        {
        erodeStep = 0;
        erode = 0;
        erodeFlip = true;
        }
      if( dilate+dilateStep < 0 )
        {
        dilateStep = 0;
        dilate = 0;
        dilateFlip = true;
        }
      if( gaussianBlur+gaussianBlurStep < 0 )
        {
        done = true;
        continue;
        }

      curTrial.erode = erode+erodeStep;
      curTrial.dilate = dilate+dilateStep;
      curTrial.gaussianBlur = gaussianBlur+gaussianBlurStep;
      typename TrialListType::iterator iter = trials.begin();
      while( iter != trials.end() )
        {
        if( iter->first == curTrial )
          {
          break;
          }
        ++iter;
        }
      if( iter == trials.end() )
        {
        eval.SetVolumeImage( curVolume );
        eval.SetMaskImage( curMask );
        eval.SetOriginalMaskImage( orgMask );
        eval.SetErode( erode+erodeStep );
        eval.SetDilate( dilate+dilateStep );
        eval.SetBoundarySize( outputBoundary );
        eval.SetGaussianBlur( gaussianBlur+gaussianBlurStep );
        eval.SetUseRegistration( true );
        eval.SetUseRegistrationTransform( true );
        eval.SetRegistrationTransform( regTfm );
        eval.Update();
        gof = eval.GetGoodnessOfFit();

        trials[curTrial] = gof;
        }
      else
        {
        std::cout << "*** Repeated point ***" << std::endl;
        gof = iter->second;
        }

      std::cout << "params = " << erode+erodeStep << ", " 
                               << dilate+dilateStep << ", " 
                               << gaussianBlur+gaussianBlurStep
                << "  val = " << gof << std::endl;

      if( gof > gofBest )
        {
        erode += erodeStep;
        dilate += dilateStep;
        gaussianBlur += gaussianBlurStep;

        gofBest = gof;
        erodeBest = erode;
        dilateBest = dilate;
        gaussianBlurBest = gaussianBlur;

        erodeFlip = false;
        if( erodeStep == 0 )
          {
          erodeStep = 1;
          }
        dilateFlip = false;
        if( dilateStep == 0 )
          {
          dilateStep = 1;
          }
        gaussianBlurFlip = false;
        }
      else
        {
        if( !erodeFlip )
          {
          erodeFlip = true;
          erodeStep *= -1;
          }
        else 
          {
          erodeStep = 0;
          if( !dilateFlip )
            {
            dilateFlip = true;
            dilateStep *= -1;
            }
          else
            {
            dilateStep = 0;
            if( !gaussianBlurFlip )
              {
              gaussianBlurFlip = true;
              gaussianBlurStep *= -1;
              }
            else
              {
              done = true;
              }
            }
          }
        }
      }

    eval.SetVolumeImage( curVolume );
    eval.SetMaskImage( curMask );
    eval.SetOriginalMaskImage( orgMask );
    eval.SetErode( erodeBest );
    eval.SetDilate( dilateBest );
    eval.SetGaussianBlur( gaussianBlurBest );
    eval.SetBoundarySize( outputBoundary );
    eval.Update();
    eval.SetProgressReporter( &progressReporter, 0.4, 0.5 );
    gof = eval.GetGoodnessOfFit();

    curVolume = eval.GetVolumeImage();
    curMask = eval.GetMaskImage();

    std::cout << "Best Params = " << erodeBest << ", " 
                             << dilateBest << ", " 
                             << gaussianBlurBest
              << "   Best Value = " << gof << std::endl;
    }
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
