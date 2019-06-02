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

#include "itktubeRidgeExtractor.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

int itktubeRidgeExtractorTest( int argc, char * argv[] )
{
  if( argc != 3 )
    {
    std::cout
      << "itktubeRidgeExtractorTest <inputImage> <outputImage>"
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

  typedef itk::tube::RidgeExtractor<ImageType> RidgeOpType;
  RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  ridgeOp->SetInputImage( im );
  ridgeOp->SetStepX( 0.75 );

  double dataMin = ridgeOp->GetDataMin();
  std::cout << "Data min = " << dataMin << std::endl;
  double dataMax = ridgeOp->GetDataMax();
  std::cout << "Data max = " << dataMax << std::endl;
  double maxT = ridgeOp->GetMaxTangentChange();
  std::cout << "Tangent change threshold = " << maxT << std::endl;
  double maxX = ridgeOp->GetMaxXChange();
  std::cout << "X threshold = " << maxX << std::endl;
  double minRidgeness = ridgeOp->GetMinRidgeness();
  std::cout << "Ridgeness threshold = " << minRidgeness << std::endl;
  double minRidgenessStart = ridgeOp->GetMinRidgenessStart();
  std::cout << "RidgenessStart threshold = " << minRidgenessStart
    << std::endl;
  double minCurvature = ridgeOp->GetMinCurvature();
  std::cout << "Curvature threshold = " << minCurvature << std::endl;
  double minCurvatureStart = ridgeOp->GetMinCurvatureStart();
  std::cout << "CurvatureStart threshold = " << minCurvatureStart
    << std::endl;
  double minRoundness = ridgeOp->GetMinRoundness();
  std::cout << "Roundness threshold = " << minRoundness << std::endl;
  double minRoundnessStart = ridgeOp->GetMinRoundnessStart();
  std::cout << "RoundnessStart threshold = " << minRoundnessStart
    << std::endl;
  itk::Index<3> extractBoundMin = ridgeOp->GetExtractBoundMinInIndexSpace();
  std::cout << "Extract bound min = " << extractBoundMin << std::endl;
  itk::Index<3> extractBoundMax = ridgeOp->GetExtractBoundMaxInIndexSpace();
  std::cout << "Extract bound max = " << extractBoundMax << std::endl;

  ridgeOp->SetScale( 2.0 );
  ridgeOp->SetDynamicScale( false );

  double recoveryMax = ridgeOp->GetMaxRecoveryAttempts();
  std::cout << "Recovery max = " << recoveryMax << std::endl;

  ImageType::Pointer imOut = ImageType::New();
  imOut->SetRegions( im->GetLargestPossibleRegion() );
  imOut->CopyInformation( im );
  imOut->Allocate();
  imOut->FillBuffer( 0 );

  double intensity;
  double roundness;
  double curvature;
  double levelness;
  unsigned int skip = 0;
  itk::ImageRegionIteratorWithIndex<ImageType> itOut( imOut,
    imOut->GetLargestPossibleRegion() );
  ImageType::IndexType indx;
  itk::ContinuousIndex<double, 3> contIndx;
  itOut.GoToBegin();
  ImageType::IndexType startIndx = itOut.GetIndex();
  std::cout << "Start..." << std::endl;
  while( !itOut.IsAtEnd() )
    {
    contIndx = itOut.GetIndex();
    int taskHash = ( int )( ( contIndx[2] - startIndx[2] ) / 2 );
    if( taskHash > 5 )
      {
      if( taskHash > 6 )
        {
        break;
        }
      taskHash = 5;
      }
    switch( taskHash )
      {
      default:
      case 0:
        {
        itOut.Set( ridgeOp->IntensityInIndexSpace( itOut.GetIndex() ) );
        break;
        }
      case 1:
        {
        itOut.Set( ridgeOp->Ridgeness( contIndx, intensity, roundness,
          curvature, levelness ) );
        break;
        }
      case 2:
        {
        ridgeOp->Ridgeness( contIndx, intensity, roundness,
          curvature, levelness );
        itOut.Set( roundness );
        break;
        }
      case 3:
        {
        ridgeOp->Ridgeness( contIndx, intensity, roundness,
          curvature, levelness );
        itOut.Set( curvature );
        break;
        }
      case 4:
        {
        ridgeOp->Ridgeness( contIndx,  intensity,roundness,
          curvature, levelness );
        itOut.Set( levelness );
        break;
        }
      case 5:
        {
        if( ++skip > 37 )
          {
          skip = 0;
          std::cout << "Local ridge: " << contIndx << std::endl;
          contIndx[1] = size[1]/2;
          if( ridgeOp->LocalRidge( contIndx ) == RidgeOpType::SUCCESS )
            {
            std::cout << "   leads to: " << contIndx << std::endl;
            for( unsigned int i=0; i<3; i++ )
              {
              indx[i] = contIndx[i];
              }
            imOut->SetPixel( indx, imOut->GetPixel( indx ) + 1 );
            std::cout << "      ridgeness = " <<
              ridgeOp->GetCurrentRidgeness() << std::endl;
            std::cout << "      roundness = " <<
              ridgeOp->GetCurrentRoundness() << std::endl;
            std::cout << "      curvature = " <<
              ridgeOp->GetCurrentCurvature() << std::endl;
            std::cout << "      levelness = " <<
              ridgeOp->GetCurrentLevelness() << std::endl;
            }
          }
        break;
        }
      }
    ++itOut;
    }
  std::cout << "...end" << std::endl;

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imWriter = ImageWriterType::New();
  imWriter->SetFileName( argv[2] );
  imWriter->SetInput( imOut );
  imWriter->SetUseCompression( true );
  imWriter->Update();

  return EXIT_SUCCESS;
}
