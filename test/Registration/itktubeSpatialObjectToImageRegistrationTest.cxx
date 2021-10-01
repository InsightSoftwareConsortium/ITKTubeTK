/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "itktubeSpatialObjectToImageRegistrationHelper.h"

#include "itktubePointBasedSpatialObjectTransformFilter.h"

#include "itktubeSubSampleSpatialObjectFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkSpatialObjectWriter.h>
#include <itkComposeScaleSkewVersor3DTransform.h>

int itktubeSpatialObjectToImageRegistrationTest( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image "
              << " Input_Group"
              << " Output_Tubes"
              << " Output_Image"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];
  const char * inputGroup = argv[2];
  const char * outputTubes = argv[3];
  const char * outputImage = argv[4];

  enum { ObjectDimension = 3 };
  typedef double            FloatType;

  typedef itk::TubeSpatialObject< ObjectDimension >            TubeType;
  typedef itk::GroupSpatialObject< ObjectDimension >           GroupType;
  typedef itk::SpatialObjectReader< ObjectDimension >          SOReaderType;
  typedef itk::Image< FloatType, ObjectDimension >             ImageType;
  typedef itk::ImageFileReader< ImageType >              ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >              ImageWriterType;
  typedef itk::tube::SpatialObjectToImageRegistrationHelper< 3, ImageType >
                                                         RegistrationHelperType;
  typedef RegistrationHelperType::MatrixTransformType    TransformType;
  typedef itk::tube::PointBasedSpatialObjectTransformFilter< TransformType, ObjectDimension >
                                                         TubeTransformFilterType;

  std::cout << "start" << std::endl;
  // read image
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( inputImage );

  // Gaussian blur the original input image to increase the likelihood of vessel
  // spatial object overlapping with the vessel image at their initial alignment.
  // this enlarges the convergence zone.
  typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType>
                                                           GaussianBlurFilterType;
  GaussianBlurFilterType::Pointer blurFilters[ObjectDimension];
  for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
    {
    blurFilters[ii] = GaussianBlurFilterType::New();
    blurFilters[ii]->SetSigma( 2.0 );
    blurFilters[ii]->SetZeroOrder();
    blurFilters[ii]->SetDirection( ii );
    }
  blurFilters[0]->SetInput( imageReader->GetOutput() );
  blurFilters[1]->SetInput( blurFilters[0]->GetOutput() );
  blurFilters[2]->SetInput( blurFilters[1]->GetOutput() );
  try
    {
    blurFilters[2]->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // read group
  SOReaderType::Pointer groupReader = SOReaderType::New();
  groupReader->SetFileName( inputGroup );
  try
    {
    groupReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "subsample" << std::endl;
  // subsample points in vessel
  typedef itk::tube::SubSampleSpatialObjectFilter< ObjectDimension > SubSampleFilterType;
  SubSampleFilterType::Pointer subSampleFilter =
    SubSampleFilterType::New();
  subSampleFilter->SetInput( groupReader->GetGroup() );
  subSampleFilter->SetSampling( 20 );
  try
    {
    subSampleFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }


  RegistrationHelperType::Pointer registrationHelper =
    RegistrationHelperType::New();

  registrationHelper->SetFixedImage( blurFilters[2]->GetOutput() );
  registrationHelper->SetMovingSpatialObject( subSampleFilter->GetOutput() );
  registrationHelper->SetRegistration(RegistrationHelperType::RegistrationMethodEnumType::RIGID);

  try
    {
    registrationHelper->Initialize();
    registrationHelper->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // validate the registration result
  TransformType::ConstPointer outputTransform =
    registrationHelper->GetCurrentMatrixTransform();
  TransformType::ParametersType lastParameters =
    outputTransform->GetParameters();

  TransformType::Pointer inverseTransform = TransformType::New();
  outputTransform->GetInverse( inverseTransform );

  itk::Matrix< double, ObjectDimension, ObjectDimension > rotationMatrix;
  itk::Vector< double, ObjectDimension > translation;
  rotationMatrix = outputTransform->GetMatrix();
  translation = outputTransform->GetTranslation();

  std::cout << "Rotation matrix: " << std::endl;
  std::string indent( "    " );
  std::cout << indent << rotationMatrix( 0,0 ) << " "
    << rotationMatrix( 0,1 ) << " " << rotationMatrix( 0,2 ) << std::endl;
  std::cout << indent << rotationMatrix( 1,0 ) << " "
    << rotationMatrix( 1,1 ) << " " << rotationMatrix( 1,2 ) << std::endl;
  std::cout << indent << rotationMatrix( 2,0 ) << " "
    << rotationMatrix( 2,1 ) << " " << rotationMatrix( 2,2 ) << std::endl;

  std::cout << "Translations: " << std::endl;
  std::cout << indent << translation[0] << " "
    << translation[1] << " " << translation[2] << std::endl;

  double knownResult[] = { 0.00697,
    -0.0049,
    0.0052,
    -1.8605,
    -0.2859,
    -0.9287 };
  int result = EXIT_SUCCESS;
  std::cout << "Parameters: " << std::endl;
  for( unsigned int ii = 0; ii < 6; ++ii )
    {
    std::cout << "Obtained: " << lastParameters[ii] << std::endl;
    std::cout << "Known:    " << knownResult[ii] << std::endl;
    if( std::abs( lastParameters[ii] - knownResult[ii] ) > 0.01 )
      {
      std::cerr << "Registration did not converge to correct parameter!" << std::endl;
      result = EXIT_FAILURE;
      }
    }

  // Transform the input tubes.
  TubeTransformFilterType::Pointer transformFilter = TubeTransformFilterType::New();
  transformFilter->SetInput( groupReader->GetGroup() );
  transformFilter->SetTransform( outputTransform.GetPointer() );
  std::cout << "Outputting transformed tubes...";
  try
    {
    transformFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    result = EXIT_FAILURE;
    }
  std::cout << " done.\n";

  // Write the transformed tube to file.
  typedef itk::SpatialObjectWriter< ObjectDimension > SOWriterType;
  SOWriterType::Pointer tubesWriter = SOWriterType::New();
  tubesWriter->SetInput( transformFilter->GetOutput() );
  tubesWriter->SetFileName( outputTubes );
  try
    {
    tubesWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject: " << err << std::endl;
    result = EXIT_FAILURE;
    }

  typedef itk::SpatialObjectToImageFilter< GroupType, ImageType>
                                              SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer vesselToImageFilter =
    SpatialObjectToImageFilterType::New();

  ImageType::Pointer img = imageReader->GetOutput();

  const double decimationFactor = 1.0;
  ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
  size[0] = size[0] / decimationFactor;
  size[1] = size[1] / decimationFactor;
  size[2] = size[2] / decimationFactor;

  ImageType::SpacingType spacing = img->GetSpacing();
  spacing[0] = spacing[0] * decimationFactor;
  spacing[1] = spacing[1] * decimationFactor;
  spacing[2] = spacing[2] * decimationFactor;

  std::cout << "Converting transformed vessel model into a binary image ... ";
  vesselToImageFilter->SetInput( dynamic_cast< const GroupType *>(transformFilter->GetOutput()) );
  vesselToImageFilter->SetSize( size );
  vesselToImageFilter->SetSpacing( spacing );
  vesselToImageFilter->SetOrigin( img->GetOrigin() );
  vesselToImageFilter->SetInsideValue( 1.0 );
  vesselToImageFilter->SetOutsideValue( 0.0 );
  vesselToImageFilter->Update();
  std::cout << "done." << std::endl;

  std::cout << "Outputting tube sampled as an image ... ";
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( outputImage );
  imageWriter->SetInput( vesselToImageFilter->GetOutput() );
  imageWriter->Update();
  std::cout << "done." << std::endl;

  return result;
}
