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
#include "itkEuler3DTransform.h"
#include "itkImageToTubeRigidMetric.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkTubeSpatialObjectPoint.h"

/**
 *  This test exercised the metric evaluation methods in the
 *  itkImageToTubeRigidMetric class. Two 3D binary images are
 *  created for testing purposes -- one of a rectangle (tube) and another of the
 *  same rectangle translated in both x, y, z (first then rotate).
 */

int itkImageToTubeRigidMetricTest(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image " << "Input_Vessel "
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<double, 3>                             Image3DType;
  typedef itk::ImageRegionIteratorWithIndex< Image3DType >  Image3DIteratorType;
  typedef itk::TubeSpatialObject<3>                         TubeType;
  typedef itk::TubeSpatialObjectPoint<3>                    TubePointType;
  typedef itk::GroupSpatialObject<3>                        TubeNetType;

  typedef itk::ImageFileReader<Image3DType>                 ImageReaderType;
  typedef itk::ImageFileWriter<Image3DType>                 Image3DWriterType;
  typedef itk::SpatialObjectWriter<3>                       TubeWriterType;
  typedef itk::SpatialObjectReader<3>                       TubeNetReaderType;

  typedef itk::ImageToTubeRigidMetric<Image3DType, TubeNetType>   MetricType;
  typedef itk::Array<double>                                      ParametersType;
  typedef MetricType::InterpolatorType                            InterpolatorType;
  typedef MetricType::TransformType                               TransformType;

  // read image (fixedImage)
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName("TubeImageTransformed.mha"); //argv[1]);
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // read tube (spatialObject)
  TubeNetReaderType::Pointer tubeReader = TubeNetReaderType::New();
  tubeReader->SetFileName("TubeOutM.tre"); //argv[2]);
  try
    {
    tubeReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Initialize the metric
  MetricType::Pointer metric = MetricType::New();
  metric->SetExtent( 3 );
  metric->SetSampling( 20 );
  metric->SetVerbose( true );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  TransformType::Pointer transform = TransformType::New();
  TransformType::ParametersType parameters = transform->GetParameters();

  metric->SetFixedImage( imageReader->GetOutput() );
  metric->SetMovingSpatialObject ( tubeReader->GetGroup() );
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

  metric->GetValue( parameters );

  return EXIT_SUCCESS;
}
