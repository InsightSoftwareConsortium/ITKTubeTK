/*=========================================================================

Program:   Maverick
Module:    $RCSfile: config.h,v $

Copyright (c) Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. No part of this software may be reproduced, distributed,
or modified, in any form or by any means, without permission in writing from
Kitware Inc.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF,
EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN
"AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

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

#include "itkVesselEnhancingDiffusion2DImageFilter.h"

#include "HessianVesselness2DCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef float                                                   PixelType;
  typedef itk::Image< PixelType,  2 >                             ImageType;
  typedef itk::ImageFileReader< ImageType >                       ReaderType;
  typedef itk::VesselEnhancingDiffusion2DImageFilter< PixelType > FilterType;
  
  ReaderType::Pointer reader = ReaderType::New();

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
  
  FilterType::Pointer filter = NULL;
  filter = FilterType::New();
  
  // set input image
  filter->SetInput( reader->GetOutput() );
  filter->SetDefaultPars();

  // set paramters
  filter->SetTimeStep( timeStep );
  filter->SetIterations( numIterations );
  filter->SetRecalculateVesselness( recalculateVesselness );
  
  filter->SetBeta( beta );
  filter->SetGamma( gamma );
 
  filter->SetEpsilon( epsilon );
  filter->SetOmega( omega );
  filter->SetSensitivity( sensitivity );
 
  // Compute scales and then set them
  std::vector< float > scales( numSigmaSteps );
  double deltaSigma = maxSigma - minSigma;
  for( int i = 0; i <= numSigmaSteps; i++ )
    {
    scales[i] = minSigma + i * ( deltaSigma / numSigmaSteps );
    }
  filter->SetScales( scales );

  // compute vesselness image
  filter->Update();
  
  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();

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

