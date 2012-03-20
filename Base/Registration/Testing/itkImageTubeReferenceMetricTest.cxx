/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2012/03/19
  Version:   $Revision: 1.0 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageToTubeRigidMetric.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTubeSpatialObjectPoint.h"

/**
 *  This test is a base to generate image and spatial objects for the
 *  registration/metric testing process.
 */

int itkImageTubeReferenceMetricTest(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Output_Image " << "Output_Tube "
              << "Output_TubeAsImage "
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<double, 3>                             Image3DType;
  typedef itk::ImageRegionIteratorWithIndex< Image3DType >  Image3DIteratorType;
  typedef itk::TubeSpatialObject<3>                         TubeType;
  typedef itk::TubeSpatialObjectPoint<3>                    TubePointType;
  typedef itk::GroupSpatialObject<3>                        TubeNetType;

  typedef itk::SpatialObjectToImageFilter< TubeNetType, Image3DType >
    SpatialObjectToImageFilterType;
  typedef itk::ImageFileWriter<Image3DType>                 ImageWriterType;
  typedef itk::SpatialObjectWriter<3>                       TubeWriterType;

  Image3DType::SizeType imageSize;
  imageSize[0] = 128;
  imageSize[1] = 128;
  imageSize[2] = 128;

  // Create fixed image
  std::cout << "Create fixed image..." << std::endl;
  Image3DType::Pointer fixedImage = Image3DType::New();
  fixedImage->SetRegions(imageSize);
  fixedImage->Allocate();
  fixedImage->FillBuffer(0);
  fixedImage->Update();

  std::cout << "Start Filling Images (Square Tube)..." << std::endl;
  Image3DIteratorType fixedIt( fixedImage, fixedImage->GetBufferedRegion() );
  int pixelIndex = 1;
  for ( fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt, ++pixelIndex )
    {
    Image3DType::IndexType index = fixedIt.GetIndex();
    if ((index[0]>=55)&&(index[0]<=65)&&(index[1]>=55)&&(index[1]<=65))
      {
        fixedIt.Set(255 - 20 * (pixelIndex % 5)); // Center brighter
      }
    }

  // guassian blur the images to increase the likelihood of vessel
  // spatial object overlapping.
  std::cout << "Apply gaussian blur..." << std::endl;
  typedef itk::RecursiveGaussianImageFilter<Image3DType, Image3DType>
    GaussianBlurFilterType;
  GaussianBlurFilterType::Pointer blurFilters[3];
  for (int i = 0; i < 3; i++)
    {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma(3.0);
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection(i);
    }
  blurFilters[0]->SetInput(fixedImage);
  blurFilters[1]->SetInput(blurFilters[0]->GetOutput());
  blurFilters[2]->SetInput(blurFilters[1]->GetOutput());
  try
    {
    blurFilters[0]->Update();
    blurFilters[1]->Update();
    blurFilters[2]->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // write image
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( argv[1] );
  std::cout << "Write imageFile: " << argv[1] << std::endl;
  imageWriter->SetInput( blurFilters[2]->GetOutput() );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Create tube
  std::cout << "Create moving image..." << std::endl;

  TubeType::Pointer tube = TubeType::New();
  TubePointType point;
  point.SetRadius( 2.0 );

  for (int i = -750; i < 750; ++i)
    {
    point.SetPosition( 15, 15, i / 100.);
    tube->GetPoints().push_back(point);
    }

  TubeNetType::Pointer group = TubeNetType::New();
  group->AddSpatialObject( tube );

  std::cout << "Write tubeFile: " << argv[2] << std::endl;
  TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
  tubeWriter->SetFileName( argv[2] );
  tubeWriter->SetInput( group );
  try
    {
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Transform the tube into an image
  SpatialObjectToImageFilterType::Pointer imageFilter =
  SpatialObjectToImageFilterType::New();
  imageFilter->SetInput( group );
  imageFilter->SetSize( imageSize );

  Image3DType::PointType origin;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  imageFilter->SetOrigin( origin );

  // write image
  ImageWriterType::Pointer imageTubeWriter = ImageWriterType::New();
  imageTubeWriter->SetFileName( argv[3] );
  std::cout << "Write tubeAsImageFile: " << argv[3] << std::endl;
  imageTubeWriter->SetInput( imageFilter->GetOutput() );
  try
    {
    imageTubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
