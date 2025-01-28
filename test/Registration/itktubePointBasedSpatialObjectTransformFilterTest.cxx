/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#include "itktubePointBasedSpatialObjectTransformFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkEuler3DTransform.h>

#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkSpatialObjectWriter.h>

int
itktubePointBasedSpatialObjectTransformFilterTest(int argc, char * argv[])
{

  if (argc < 12)
  {
    std::cerr << "Missing Parameters: " << argv[0] << " Input_Vessel " << "Output_Vessel "
              << "Example_Image " << "Output_Image "
              << "R1 R2 R3 T1 T2 T3 Write_Image_Flag" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::GroupSpatialObject<3>    GroupType;
  typedef itk::SpatialObjectReader<3>   SOReaderType;
  typedef itk::SpatialObjectWriter<3>   SOWriterType;
  typedef itk::Euler3DTransform<double> TransformType;

  typedef itk::tube::PointBasedSpatialObjectTransformFilter<TransformType, 3> TransformFilterType;

  // read in vessel
  SOReaderType::Pointer reader = SOReaderType::New();
  reader->SetFileName(argv[1]);

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // generate transform
  TransformType::Pointer transform = TransformType::New();
  itk::Vector<double, 3> rotation;
  rotation[0] = std::atof(argv[5]); //-0.5/itk::Math::one_over_pi;
  rotation[1] = std::atof(argv[6]); //-0.5/itk::Math::one_over_pi;
  rotation[2] = std::atof(argv[7]); // 0.0

  std::cout << "input rotation: " << rotation[0] << " " << rotation[1] << " " << rotation[2] << std::endl;

  itk::Vector<double, 3> translation;
  translation[0] = std::atof(argv[8]);
  translation[1] = std::atof(argv[9]);
  translation[2] = std::atof(argv[10]);

  std::cout << "input translation: " << translation[0] << " " << translation[1] << " " << translation[2] << std::endl;

  double ca = cos(rotation[0]);
  double sa = sin(rotation[0]);
  double cb = cos(rotation[1]);
  double sb = sin(rotation[1]);
  double cg = cos(rotation[2]);
  double sg = sin(rotation[2]);

  itk::Matrix<double, 3, 3> rotationMatrix;
  rotationMatrix[0][0] = ca * cb;
  rotationMatrix[0][1] = ca * sb * sg - sa * cg;
  rotationMatrix[0][2] = ca * sb * cg + sa * sg;
  rotationMatrix[1][0] = sa * cb;
  rotationMatrix[1][1] = sa * sb * sg + ca * cg;
  rotationMatrix[1][2] = sa * sb * cg - ca * sg;
  rotationMatrix[2][0] = -sb;
  rotationMatrix[2][1] = cb * sg;
  rotationMatrix[2][2] = cb * cg;

  transform->SetMatrix(rotationMatrix);

  std::cout << rotationMatrix(0, 0) << " " << rotationMatrix(0, 1) << " " << rotationMatrix(0, 2) << std::endl;
  std::cout << rotationMatrix(1, 0) << " " << rotationMatrix(1, 1) << " " << rotationMatrix(1, 2) << std::endl;
  std::cout << rotationMatrix(2, 0) << " " << rotationMatrix(2, 1) << " " << rotationMatrix(2, 2) << std::endl;

  // transform->SetRotation( rotation[0], rotation[1], rotation[2] );

  transform->Translate(translation);

  std::cout << translation[0] << " " << translation[1] << " " << translation[2] << std::endl;

  // create transform filter
  TransformFilterType::Pointer transformFilter = TransformFilterType::New();
  transformFilter->SetInput(reader->GetGroup());
  transformFilter->SetTransform(transform);

  try
  {
    transformFilter->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // write vessel
  SOWriterType::Pointer writer = SOWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(transformFilter->GetOutput());

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  if (std::atoi(argv[11]))
  {
    // Write vessel as an image
    typedef itk::Image<double, 3>           ImageType;
    typedef itk::ImageFileReader<ImageType> ImageReaderType;
    typedef itk::ImageFileWriter<ImageType> ImageWriterType;

    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader->SetFileName(argv[3]);
    imageReader->Update();

    typedef itk::SpatialObjectToImageFilter<GroupType, ImageType> SpatialObjectToImageFilterType;
    SpatialObjectToImageFilterType::Pointer vesselToImageFilter = SpatialObjectToImageFilterType::New();

    vesselToImageFilter->SetInput(dynamic_cast<GroupType *>(transformFilter->GetOutput()));
    vesselToImageFilter->SetSize(imageReader->GetOutput()->GetLargestPossibleRegion().GetSize());
    vesselToImageFilter->SetOrigin(imageReader->GetOutput()->GetOrigin());
    vesselToImageFilter->SetSpacing(imageReader->GetOutput()->GetSpacing());
    vesselToImageFilter->SetInsideValue(1.0);
    vesselToImageFilter->SetOutsideValue(0.0);
    vesselToImageFilter->Update();

    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetFileName(argv[4]);
    imageWriter->SetUseCompression(true);
    imageWriter->SetInput(vesselToImageFilter->GetOutput());
    imageWriter->Update();
  }

  return EXIT_SUCCESS;
}
