/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "itktubeBasisFeatureVectorGenerator.h"
#include "itktubeRidgeFFTFeatureVectorGenerator.h"

int
itktubeRidgeBasisFeatureVectorGeneratorTest(int argc, char * argv[])
{
  if (argc != 7)
  {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage maskImage objId bkgId outputLDA0Image outputLDA1Image" << std::endl;
    return EXIT_FAILURE;
  }

  // Define the dimension of the images
  enum
  {
    Dimension = 2
  };

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension> ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typedef itk::Image<unsigned char, Dimension> LabelMapType;
  typedef itk::ImageFileReader<LabelMapType>   LabelMapReaderType;


  // Declare the type for the Filter
  typedef itk::tube::RidgeFFTFeatureVectorGenerator<ImageType>            FilterType;
  typedef itk::tube::BasisFeatureVectorGenerator<ImageType, LabelMapType> BasisFilterType;

  // Create the reader
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during input read:" << std::endl << e;
    return EXIT_FAILURE;
  }
  ImageType::Pointer inputImage = reader->GetOutput();

  // Create the mask reader
  LabelMapReaderType::Pointer mReader = LabelMapReaderType::New();
  mReader->SetFileName(argv[2]);
  try
  {
    mReader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during input mask read:" << std::endl << e;
    return EXIT_FAILURE;
  }
  LabelMapType::Pointer maskImage = mReader->GetOutput();

  FilterType::RidgeScalesType scales(3);
  scales[0] = 0.4;
  scales[1] = 0.8;
  scales[2] = 1.6;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(inputImage);
  filter->SetScales(scales);
  filter->SetUpdateWhitenStatisticsOnUpdate(true);
  filter->Update();
  std::cout << filter << std::endl;

  BasisFilterType::Pointer basisFilter = BasisFilterType::New();
  basisFilter->SetInputFeatureVectorGenerator(filter.GetPointer());
  basisFilter->SetInput(inputImage);
  basisFilter->SetLabelMap(maskImage);
  int objId = atoi(argv[3]);
  int bkgId = atoi(argv[4]);
  basisFilter->SetObjectId(objId);
  basisFilter->AddObjectId(bkgId);
  basisFilter->AddObjectId(0);
  std::cout << "Start" << std::endl;
  std::cout << basisFilter << std::endl;
  basisFilter->SetNumberOfLDABasisToUseAsFeatures(2);
  basisFilter->SetNumberOfPCABasisToUseAsFeatures(2);
  basisFilter->SetUpdateWhitenStatisticsOnUpdate(true);
  basisFilter->Update();
  std::cout << "Stop" << std::endl;
  std::cout << basisFilter << std::endl;

  basisFilter->SetLabelMap(NULL);

  WriterType::Pointer featureImage0Writer = WriterType::New();
  featureImage0Writer->SetFileName(argv[5]);
  featureImage0Writer->SetUseCompression(true);
  featureImage0Writer->SetInput(basisFilter->GetFeatureImage(0));
  try
  {
    featureImage0Writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
  }

  WriterType::Pointer featureImage1Writer = WriterType::New();
  featureImage1Writer->SetFileName(argv[6]);
  featureImage1Writer->SetUseCompression(true);
  featureImage1Writer->SetInput(basisFilter->GetFeatureImage(1));
  try
  {
    featureImage1Writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during write:" << std::endl << e;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
