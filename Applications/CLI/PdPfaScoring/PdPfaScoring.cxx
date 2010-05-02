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

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "tubePdPfaScorer.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "PdPfaScoringCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "SampleCLIApplication",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef unsigned char                                    PixelType;
  typedef itk::OrientedImage< PixelType, dimensionT >      ImageType;
  typedef itk::ImageFileReader< ImageType >                ReaderType;
  
  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  typename ReaderType::Pointer maskReader = ReaderType::New();
  maskReader->SetFileName( inputPrior.c_str() );
  typename ReaderType::Pointer modifiedMaskReader = ReaderType::New();
  modifiedMaskReader->SetFileName( modifiedPrior.c_str() );

  // read the input volume while handling exceptions
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading input volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  // read the input prior while handling exceptions
  try
    {
    maskReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading prior: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  // read the modified prior while handling exceptions
  try
    {
    modifiedMaskReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading modfied prior: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop("Load data");

  double progress = 0.1;
  progressReporter.Report( progress );

  typename ImageType::Pointer input = reader->GetOutput();
  typename ImageType::Pointer prior = maskReader->GetOutput();
  typename ImageType::Pointer modPrior = modifiedMaskReader->GetOutput();

  timeCollector.Start( "Compute Change Statistics" );
  int totalNumChanges;
  int totalNumChangesFound;
  int totalNumFalsePositives;
  tube::PdPfaScorer< PixelType, dimensionT > scorer;
  scorer.ComputeChangeStatistics( input, modPrior, prior, 
                                  static_cast<int>( featureSize ),
                                  totalNumChanges, totalNumChangesFound,
                                  totalNumFalsePositives );

  timeCollector.Stop( "Compute Change Statistics" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End( );
  
  timeCollector.Report();

  // Print at the very end
  std::cout << std::endl;
  std::cout << "***RESULTS***" << std::endl;
  std::cout << "Total Changes: " << totalNumChanges 
            << std::endl;
  std::cout << "Total Changes Found: " << totalNumChangesFound 
            << std::endl;
  std::cout << "Total False Positives: " << totalNumFalsePositives 
            << std::endl;


  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
