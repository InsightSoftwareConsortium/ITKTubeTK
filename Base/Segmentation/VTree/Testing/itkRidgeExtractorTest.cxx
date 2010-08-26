/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 ( the "License" );
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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "../itkRidgeExtractor.h"

int itkRidgeExtractorTest( int argc, char * argv[] )
  {
  if( argc != 3 )
    {
    std::cout 
      << "itkRidgeExtractorTest <inputImage> <outputImage>" 
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  ImageType::SizeType size = im->GetLargestPossibleRegion().GetSize();

  typedef itk::RidgeExtractor<ImageType> RidgeOpType;
  RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  ridgeOp->SetInputImage( im );
  ridgeOp->SetStepX( 0.75 );

  double dataMin = ridgeOp->GetDataMin();
  std::cout << "Data min = " << dataMin << std::endl;
  double dataMax = ridgeOp->GetDataMax();
  std::cout << "Data max = " << dataMax << std::endl;
  double threshT = ridgeOp->GetThreshT();
  std::cout << "Theta threshold = " << threshT << std::endl;
  double threshX = ridgeOp->GetThreshX();
  std::cout << "X threshold = " << threshX << std::endl;
  double threshRidgeness = ridgeOp->GetThreshRidgeness();
  std::cout << "Ridgeness threshold = " << threshRidgeness << std::endl;
  double threshRidgenessStart = ridgeOp->GetThreshRidgenessStart();
  std::cout << "RidgenessStart threshold = " << threshRidgenessStart 
    << std::endl;
  double threshCurvature = ridgeOp->GetThreshCurvature();
  std::cout << "Curvature threshold = " << threshCurvature << std::endl;
  double threshCurvatureStart = ridgeOp->GetThreshCurvatureStart();
  std::cout << "CurvatureStart threshold = " << threshCurvatureStart 
    << std::endl;
  double threshRoundness = ridgeOp->GetThreshRoundness();
  std::cout << "Roundness threshold = " << threshRoundness << std::endl;
  double threshRoundnessStart = ridgeOp->GetThreshRoundnessStart();
  std::cout << "RoundnessStart threshold = " << threshRoundnessStart 
    << std::endl;
  itk::Index<3> extractBoundMin = ridgeOp->GetExtractBoundMin();
  std::cout << "Extract bound min = " << extractBoundMin << std::endl;
  itk::Index<3> extractBoundMax = ridgeOp->GetExtractBoundMax();
  std::cout << "Extract bound max = " << extractBoundMax << std::endl;
  
  ridgeOp->SetScale( 2.0 );
  ridgeOp->SetExtent( 3.0 );
  ridgeOp->SetDynamicScale( true );

  double recoveryMax = ridgeOp->GetRecoveryMax();
  std::cout << "Recovery max = " << recoveryMax << std::endl;

  ImageType::Pointer imOut = ImageType::New();
  imOut->SetRegions( im->GetLargestPossibleRegion() );
  imOut->CopyInformation( im );
  imOut->Allocate();
  imOut->FillBuffer( 0 );

  double roundness;
  double curvature;
  unsigned int skip = 0;
  itk::ImageRegionIteratorWithIndex<ImageType> itOut( imOut,
    imOut->GetLargestPossibleRegion() );
  ImageType::IndexType indx;
  itk::ContinuousIndex<double, 3> contIndx;
  ImageType::PointType pnt;
  itOut.GoToBegin();
  ImageType::IndexType startIndx = itOut.GetIndex();
  std::cout << "Start..." << std::endl;
  bool firstSearch = false;
  while( !itOut.IsAtEnd() )
    {
    contIndx = itOut.GetIndex();
    switch( ( itOut.GetIndex()[2] - startIndx[2] ) % 5 )
      {
      default:
      case 0:
        std::cout << "Intensity: " << contIndx << std::endl;
        itOut.Set( ridgeOp->Intensity( itOut.GetIndex() ) );
        firstSearch = true;
        break;
      case 1:
        std::cout << "Ridgeness: " << contIndx << std::endl;
        itOut.Set( ridgeOp->Ridgeness( contIndx, roundness, curvature ) );
        firstSearch = true;
        break;
      case 2:
        std::cout << "Roundness: " << contIndx << std::endl;
        ridgeOp->Ridgeness( contIndx, roundness, curvature );
        itOut.Set( roundness );
        firstSearch = true;
        break;
      case 3:
        std::cout << "Curvature: " << contIndx << std::endl;
        ridgeOp->Ridgeness( contIndx, roundness, curvature );
        itOut.Set( curvature );
        firstSearch = true;
        break;
      case 4:
        if( firstSearch )
          {
          skip = (int)(contIndx[1]-startIndx[1])%4 * 3;
          firstSearch = false;
          }
        if( ++skip > 12 )
          {
          skip = 0;
          std::cout << "Local ridge: " << contIndx << std::endl;
          contIndx[1] = size[1]/2;
          if( ridgeOp->LocalRidge( contIndx ) )
            {
            std::cout << "   leads to: " << contIndx << std::endl;
            for( unsigned int i=0; i<3; i++ )
              {
              indx[i] = contIndx[i];
              }
            if( vnl_math_abs( indx[2] - itOut.GetIndex()[2] ) < 4 )
              {
              indx[2] = itOut.GetIndex()[2];
              int v = imOut->GetPixel( indx );
              imOut->SetPixel( indx, v+1 );
              }
            }
          }
        break;
      }
    ++itOut;
    }
  std::cout << "...end" << std::endl;

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imWriter = ImageWriterType::New();
  imWriter->SetFileName( argv[2] );
  imWriter->SetInput( imOut );
  imWriter->Update();

  return EXIT_SUCCESS;
  }

