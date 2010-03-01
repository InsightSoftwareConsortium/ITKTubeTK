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

#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkShiftScaleImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "SimulateTomosynthesisCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef float                                                   PixelType;
  typedef itk::Image< PixelType,  dimensionT >                    ImageType;
  typedef itk::ImageFileReader< ImageType >                       ReaderType;
  
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

  typename ImageType::Pointer curImage = reader->GetOutput();

  if( foreground != 1 || background != 0 )
    {
    typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( curImage );
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
    curImage = filter->GetOutput();
    }
  
  typedef itk::BinaryBallStructuringElement<PixelType, dimensionT >  BallType;
  BallType ball;
  ball.SetRadius( 1 );
  ball.CreateStructuringElement();
  if( erodeDilate && erodeRadius > 0 )
    {
    typedef itk::BinaryErodeImageFilter
                 <ImageType, ImageType, BallType>       ErodeFilterType;

    for(int r=0; r<erodeRadius; r++)
      {
      typename ErodeFilterType::Pointer filter = ErodeFilterType::New();
      filter->SetBackgroundValue( 0 );
      filter->SetErodeValue( 1 );
      filter->SetKernel( ball );
      filter->SetInput( curImage );
      filter->Update();
      curImage = filter->GetOutput();
      }
    }
  if( dilateRadius > 0 )
    {
    typedef itk::BinaryDilateImageFilter
                 <ImageType, ImageType, BallType>       DilateFilterType;
    for(int r=0; r<dilateRadius; r++)
      {
      typename DilateFilterType::Pointer filter = DilateFilterType::New();
      filter->SetKernel( ball );
      filter->SetDilateValue( 1 );
      filter->SetInput( curImage );
      filter->Update();
      curImage = filter->GetOutput();
      }
    }
  if( !erodeDilate && erodeRadius > 0 )
    {
    typedef itk::BinaryErodeImageFilter
                 <ImageType, ImageType, BallType>       ErodeFilterType;

    for(int r=0; r<erodeRadius; r++)
      {
      typename ErodeFilterType::Pointer filter = ErodeFilterType::New();
      filter->SetBackgroundValue( 0 );
      filter->SetErodeValue( 1 );
      filter->SetKernel( ball );
      filter->SetInput( curImage );
      filter->Update();
      curImage = filter->GetOutput();
      }
    }
  
  if( gaussianBlurStdDev > 0 )
    {
    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    for(unsigned int i=0; i<dimensionT; i++)
      {
      filter = FilterType::New();
      filter->SetInput( curImage );
      // filter->SetNormalizeAcrossScale( true );
      filter->SetSigma( gaussianBlurStdDev );

      filter->SetOrder( 
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );

      filter->Update();
      curImage = filter->GetOutput();
      }
    }

  if( outOffset != 0 || outScale != 1 )
    {
    typedef itk::ShiftScaleImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter = FilterType::New();
    filter->SetInput( curImage );
    filter->SetShift( outOffset );
    filter->SetScale( outScale );

    filter->Update();
    curImage = filter->GetOutput();
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput ( curImage );

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

int main( int argc, char **argv )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}

