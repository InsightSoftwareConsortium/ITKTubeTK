/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageToTubeRigidRegistration.h"
#include "itkSpatialObjectReader.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransform.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTubeToTubeTransformFilter.h"
#include "itkEuler3DTransform.h"
#include "itkMath.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

int itkImageToTubeRigidRegistrationTest(int argc, char* argv [] )
{

  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Image " << "Input_Vessel " << "Output_Image " 
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::GroupSpatialObject<3>                      TubeNetType;
  typedef itk::SpatialObjectReader<3>                     TubeNetReaderType;
  typedef itk::Image<double, 3>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                 ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>                 ImageWriterType;
  typedef itk::ImageToTubeRigidRegistration<ImageType, TubeNetType>  
                                                          RegistrationFilterType;
  typedef itk::Euler3DTransform<double>                   TransformType; 
  typedef itk::TubeToTubeTransformFilter<TransformType,3> TubeTransformFilterType;

  // read image
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(argv[1]);
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // guassian blur the original input image to increase the likelihood of vessel
  // spatial object overlapping with the vessel image at their initial alignment.
  // this enlarges the convergence zone.
  typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType>
                                                           GaussianBlurFilterType;
  GaussianBlurFilterType::Pointer blurFilters[3];
  for (int i = 0; i < 3; i++)
    {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma(3.0);
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection(i);
    }
  blurFilters[0]->SetInput(imageReader->GetOutput());
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

  // read vessel  
  TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName(argv[2]);
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // register the vessel and the image
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize( 137593424 );
  
  double parameterScales[6] = {30.0, 30.0, 30.0, 1.0, 1.0, 1.0};
  double initialPose[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (argc > 9)
    {
    for (int i = 0; i < 6; i++)
      {
      initialPose[i] = atof(argv[4+i]);
      }
    }
  
  RegistrationFilterType::Pointer  registrationFilter = RegistrationFilterType::New();
  
  registrationFilter->SetFixedImage( blurFilters[2]->GetOutput() );
  registrationFilter->SetMovingSpatialObject( vesselReader->GetGroup() );
  registrationFilter->SetNumberOfIteration( 1000 );
  registrationFilter->SetLearningRate(0.1);
  registrationFilter->SetInitialPosition(initialPose);
  registrationFilter->SetParametersScale(parameterScales);
  registrationFilter->SetVerbose(0);
  try
    {
    registrationFilter->Initialize();
    registrationFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // validate the registration result
  TransformType* outputTransform = 
    dynamic_cast<TransformType *>(registrationFilter->GetTransform());
  
  TransformType::Pointer invTransform = TransformType::New();
  outputTransform->GetInverse(invTransform);
  
  std::cout << "Registration result: ";
  for (int i = 0; i < 6; i++)
    {
    std::cout << outputTransform->GetParameters().GetElement(i) << " ";
    }
  std::cout << std::endl;

  // create transform filter
  TubeTransformFilterType::Pointer transformFilter = TubeTransformFilterType::New();
  transformFilter->SetInput(vesselReader->GetGroup());
  transformFilter->SetScale(1.0);
  transformFilter->SetTransform(outputTransform);
  
  try
    {
    transformFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  
  typedef itk::SpatialObjectToImageFilter<TubeNetType, ImageType> 
                                              SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer vesselToImageFilter = 
    SpatialObjectToImageFilterType::New();

  std::cout << "Converting transformed vessel model into a binary image ... ";
  vesselToImageFilter->SetInput(transformFilter->GetOutput());
  vesselToImageFilter->SetSize(
    imageReader->GetOutput()->GetLargestPossibleRegion().GetSize());
  vesselToImageFilter->SetOrigin(imageReader->GetOutput()->GetOrigin());
  vesselToImageFilter->SetSpacing(imageReader->GetOutput()->GetSpacing());
  vesselToImageFilter->SetInsideValue(1.0);
  vesselToImageFilter->SetOutsideValue(0.0);
  vesselToImageFilter->Update();
  std::cout << "end." << std::endl;
  
  std::cout << "Outputing result image ... ";
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName(argv[3]);
  imageWriter->SetInput(vesselToImageFilter->GetOutput());
  imageWriter->Update();
  std::cout << "end." << std::endl;
  
  return EXIT_SUCCESS;

}

