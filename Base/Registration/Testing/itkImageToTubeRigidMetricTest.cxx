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

#include "itkImageFileReader.h"


/**
 *  This test exercised the metric evaluation methods in the
 *  itkImageToTubeRigidMetric class. Two 3D binary images are
 *  created for testing purposes -- one of a rectangle (tube) and another of the
 *  same rectangle translated in both x, y, z (first then rotate).
 *
 */

int itkImageToTubeRigidMetricTest(int, char* [] )
{
  typedef itk::Image<double, 3>                             Image3DType;
  typedef itk::ImageRegionIteratorWithIndex< Image3DType >  Image3DIteratorType;
  typedef itk::TubeSpatialObject<3>                         TubeType;
  typedef itk::TubeSpatialObjectPoint<3>                    TubePointType;
  typedef itk::GroupSpatialObject<3>                        TubeNetType;

  typedef itk::SpatialObjectToImageFilter< TubeNetType, Image3DType >
    SpatialObjectToImageFilterType;
  typedef itk::ImageFileWriter<Image3DType>                 Image3DWriterType;
  typedef itk::SpatialObjectWriter<3>                       TubeWriterType;
  typedef itk::SpatialObjectReader<3>                       TubeNetReaderType;

  typedef itk::ImageToTubeRigidMetric<Image3DType, TubeNetType>   MetricType;
  typedef itk::Array<double>                                      ParametersType;
  typedef MetricType::InterpolatorType                            InterpolatorType;
  typedef MetricType::TransformType                               TransformType;

  // Create fixed image
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

  std::cout << "Start Filling Images..." << std::endl;
  Image3DIteratorType fixedIt( fixedImage, fixedImage->GetBufferedRegion() );
  for ( fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt )
    {
    Image3DType::IndexType index = fixedIt.GetIndex();
    if ((index[0]>=48)&&(index[0]<=80)&&(index[1]>=48)&&(index[1]<=80))
      {
      fixedIt.Set(255);
      }
    }

  // guassian blur the images to increase the likelihood of vessel
  // spatial object overlapping.
  // ----> TODO with a proper model
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

  std::cout << "End." << std::endl;

  // Create moving image
  std::cout << "Create moving image..." << std::endl;
  TubeType::Pointer movingImage = TubeType::New();

  //Create some tube points
  TubePointType point1_1;
  point1_1.SetPosition(10.5, 10.0, 15.0);
  point1_1.SetRadius(4.0);

  TubePointType point1_2;
  point1_2.SetPosition(50, 30.5, 30.0);
  point1_2.SetRadius(2.0);

  TubePointType point1_3;
  point1_3.SetPosition(65.0, 45.5, 30.0);
  point1_3.SetRadius(3.0);

  TubePointType point1_4;
  point1_4.SetPosition(100.0, 115.0, 60.0);
  point1_4.SetRadius(2.0);

  TubePointType point1_5;
  point1_5.SetPosition(170.0, 60.0, 200.0);
  point1_5.SetRadius(1.55);

  TubePointType point1_1_1;
  point1_1_1.SetPosition(100.0, 120.0, 60.0);
  point1_1_1.SetRadius(2.0);

  TubePointType point1_1_2;
  point1_1_2.SetPosition(65.0, 130.0, 60.0);
  point1_1_2.SetRadius(1.0);

  TubePointType point1_1_3;
  point1_1_3.SetPosition(35.0, 130.0, 90.0);
  point1_1_3.SetRadius(3.5);

  //Create "root" tube
  movingImage->GetPoints().push_back(point1_1);
  movingImage->GetPoints().push_back(point1_2);
  movingImage->GetPoints().push_back(point1_3);
  movingImage->GetPoints().push_back(point1_4);
  movingImage->GetPoints().push_back(point1_5);
  movingImage->SetRoot(true);

  std::cout << "End." << std::endl;

  TubeNetType::Pointer group = TubeNetType::New();
  group->AddSpatialObject(movingImage);

  //----Image
  /*SpatialObjectToImageFilterType::Pointer imageFilter =
    SpatialObjectToImageFilterType::New();
  imageFilter->SetInput( group );
  imageFilter->Update();*/

  std::cout << "End." << std::endl;

  std::cout << "Outputing both image ... "  << std::endl;
  Image3DWriterType::Pointer fixedImagetWriter = Image3DWriterType::New();
  fixedImagetWriter->SetFileName("fixedObject.mha");
  fixedImagetWriter->SetInput(fixedImage);
  fixedImagetWriter->Update();

  /*Image3DWriterType::Pointer movingImagetWriter = Image3DWriterType::New();
  movingImagetWriter->SetFileName("spatialObject.mha");
  movingImagetWriter->SetInput(imageFilter->GetOutput());
  movingImagetWriter->Update();
  TubeWriterType::Pointer movingImagetWriter = TubeWriterType::New();
  movingImagetWriter->SetFileName("spatialObject.mha");
  movingImagetWriter->SetInput(tube1);
  movingImagetWriter->Update();*/

  std::cout << "end." << std::endl;

  // From File
  TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName("DataTest/Branch.tre");
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // Initialize the metric
  MetricType::Pointer metric = MetricType::New();
  metric->SetExtent( 3 );
  metric->SetSampling( 1 );
  metric->SetVerbose( true );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  TransformType::Pointer transform = TransformType::New();
  TransformType::ParametersType parameters = transform->GetParameters();

  metric->SetFixedImage( fixedImage );
 // metric->SetFixedImage( imageFilter->GetOutput() );
  metric->SetMovingSpatialObject ( group );
  //metric->SetMovingSpatialObject ( vesselReader->GetGroup() );
  metric->SetInterpolator( interpolator );
  metric->SetTransform( transform );
  try
    {
    metric->Initialize();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while initializing metric." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  MetricType::MeasureType value = metric->GetValue( parameters );
  std::cout << "Metric Measure: " << value << std::endl;

  return EXIT_SUCCESS;
}
