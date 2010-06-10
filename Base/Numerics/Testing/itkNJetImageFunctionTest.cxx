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

#include <itkImage.h>
#include "itkFilterWatcher.h"
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkNJetImageFunction.h>

int itkNJetImageFunctionTest(int argc, char* argv [] ) 
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " function inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
    }
  
  // Define the dimension of the images
  const unsigned int Dimension = 2;

  // Define the pixel type
  typedef float PixelType;
  
  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
 
  // Declare the type for the Filter
  typedef itk::NJetImageFunction< ImageType > FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during read:\n"  << e;
    return EXIT_FAILURE;
    }

  ImageType::Pointer inputImage = reader->GetOutput();

  typedef itk::NJetImageFunction< ImageType >   FunctionType;
  FunctionType::Pointer func = FunctionType::New();
  func->SetInputImage( inputImage );

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->Allocate();

  bool error = false;

  func->ComputeStatistics();
  func->UseProjection( false );
  if( func->GetMin() != 100 )
    {
    error = true;
    std::cerr << "Min = " << func->GetMin() << " != 100" << std::endl;
    }
  if( func->GetMax() != 1000 )
    {
    error = true;
    std::cerr << "Max = " << func->GetMax() << " != 1000" << std::endl;
    }

  FunctionType::VectorType v1, v2;
  v1[0] = 1;
  v2[1] = 1;

  int function = atoi( argv[1] );

  itk::ImageRegionIteratorWithIndex< ImageType > outIter( outputImage,
    outputImage->GetLargestPossibleRegion() );
  ImageType::PointType pnt;
  outIter.GoToBegin();
  while( !outIter.IsAtEnd() )
    {
    inputImage->TransformIndexToPhysicalPoint( outIter.GetIndex(), pnt);
    switch( function )
      {
      case 0:
        {
        outIter.Set( func->Evaluate( pnt, 2.0 ) );
        break;
        }
      case 1:
        {
        outIter.Set( func->Evaluate( pnt, v1, 2.0 ) );
        break;
        }
      case 2:
        {
        outIter.Set( func->Evaluate( pnt, v1, v2, 2.0 ) );
        break;
        }
      case 3:
        {
        outIter.Set( func->EvaluateAtIndex( outIter.GetIndex(), 2.0 ) );
        break;
        }
      case 4:
        {
        outIter.Set( func->EvaluateAtIndex( outIter.GetIndex(), v1, 2.0 ) );
        break;
        }
      case 5:
        {
        outIter.Set( func->EvaluateAtIndex( outIter.GetIndex(),
            v1, v2, 2.0 ) );
        break;
        }
      case 6:
        {
        outIter.Set( func->Derivative( pnt, 2.0 )[0] );
        break;
        }
      case 7:
        {
        outIter.Set( func->Derivative( pnt, v1, 2.0 )[0] );
        break;
        }
      case 8:
        {
        outIter.Set( func->Derivative( pnt, v1, v2, 2.0 )[0] );
        break;
        }
      case 9:
        {
        outIter.Set( func->DerivativeAtIndex( outIter.GetIndex(), 2.0 )[0] );
        break;
        }
      case 10:
        {
        outIter.Set( func->DerivativeAtIndex( outIter.GetIndex(), v1, 2.0 )[0] );
        break;
        }
      case 11:
        {
        outIter.Set( func->DerivativeAtIndex( outIter.GetIndex(), v1, v2, 2.0 )[0] );
        break;
        }
      case 12:
        {
        double val;
        func->ValueAndDerivative(pnt, val, 2.0 );
        outIter.Set( val );
        break;
        }
      case 13:
        {
        double val;
        func->ValueAndDerivative(pnt, val, v1, 2.0 );
        outIter.Set( val );
        break;
        }
      case 14:
        {
        double val;
        func->ValueAndDerivative(pnt, val, v1, v2, 2.0 );
        outIter.Set( val );
        break;
        }
      case 15:
        {
        double val;
        func->ValueAndDerivativeAtIndex( outIter.GetIndex(), val, 2.0 );
        outIter.Set( val );
        break;
        }
      case 16:
        {
        double val;
        func->ValueAndDerivativeAtIndex( outIter.GetIndex(), val, v1, 2.0 );
        outIter.Set( val );
        break;
        }
      case 17:
        {
        double val;
        func->ValueAndDerivativeAtIndex( outIter.GetIndex(), val, v1, v2, 2.0 );
        outIter.Set( val );
        break;
        }
      }
    ++outIter;
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( outputImage );
  
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;

}

