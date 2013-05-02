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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// Must declare a forward declaration of DoIt before tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char **argv );

// Must include CLP file before tubeCLIHelperFunctions
#include "EnhanceTubesUsingDiffusionCLP.h"

#include "tubeCLIHelperFunctions.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTubeEnhancingDiffusion2DImageFilter.h"

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char **argv )
{
  PARSE_ARGS;

  typedef pixelT                                           PixelType;
  typedef itk::Image< PixelType,  dimensionT  >            ImageType;
  typedef itk::ImageFileReader< ImageType >                ReaderType;
  typedef itk::TubeEnhancingDiffusion2DImageFilter< PixelType,
                                                    dimensionT  >
                                                           FilterType;

  typename ReaderType::Pointer reader = ReaderType::New();

  //read input image
  reader->SetFileName( inputVolume.c_str() );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  typename FilterType::Pointer filter = NULL;
  filter = FilterType::New();

  // set input image
  filter->SetInput( reader->GetOutput() );
  filter->SetDefaultPars();

  // set paramters
  filter->SetTimeStep( timeStep );
  filter->SetIterations( numIterations );
  filter->SetRecalculateTubeness( recalculateTubeness );

  filter->SetBeta( beta );
  filter->SetGamma( gamma );

  filter->SetEpsilon( epsilon );
  filter->SetOmega( omega );
  filter->SetSensitivity( sensitivity );

  // Compute scales and then set them
  std::vector< float > scales( numSigmaSteps );
  double deltaSigma = maxSigma - minSigma;
  for( int i = 0; i < numSigmaSteps; i++ )
    {
    scales[i] = minSigma + i * ( deltaSigma / numSigmaSteps );
    }
  filter->SetScales( scales );

  // compute vesselness image
  filter->Update();

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput ( filter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
