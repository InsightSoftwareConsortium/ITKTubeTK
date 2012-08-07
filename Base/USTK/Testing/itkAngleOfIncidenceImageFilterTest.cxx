/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSheetnessMeasureImageFilter.h"
#include "itkAngleOfIncidenceImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkThresholdImageFilter.h"

int itkAngleOfIncidenceImageFilterTest(int argc ,char* argv [] )
{
  // Argument parsing.
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " Input_Ultrasound_Image Output_Sheetness_Image Output_Angle_Of_Incidence_Image" << std::endl;
    return EXIT_FAILURE;
    }

  // Types.
  static const unsigned int Dimension = 3;

  typedef unsigned short                              UltrasoundPixelType;
  typedef itk::Image< UltrasoundPixelType, Dimension > UltrasoundImageType;

  typedef float
    AngleOfIncidencePixelType;
  typedef itk::Image< AngleOfIncidencePixelType, Dimension >
    AngleOfIncidenceImageType;

  // Reader.
  typedef itk::ImageFileReader< UltrasoundImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  // Declare the type for the Hessian filter
  typedef itk::HessianRecursiveGaussianImageFilter<
                                            ImageType >  HessianFilterType;

  typedef HessianFilterType::OutputImageType HessianImageType;

  // Declare the type for the sheetness measure filter
  typedef itk::SheetnessMeasureImageFilter< float >  AngleOfIncidenceImageFilter;

  typedef AngleOfIncidenceImageFilter::OutputImageType SheetnessImageType;


  // Create a Hessian Filter
  HessianFilterType::Pointer filterHessian = HessianFilterType::New();

  // Create a sheetness Filter
  AngleOfIncidenceImageFilter::Pointer filterSheetness = AngleOfIncidenceImageFilter::New();


  // Connect the input images
  filterHessian->SetInput( inputImage );
  filterSheetness->SetInput( filterHessian->GetOutput() );

  // Select the value of Sigma
  filterHessian->SetSigma( 0.5 );


  // Execute the filter
  std::cout << "Generate sheetness measure" << std::endl;
  filterSheetness->Update();

  //Write out the sheetness image
  typedef AngleOfIncidenceImageFilter::OutputImageType SheetnessImageType;

  typedef itk::ImageFileWriter<SheetnessImageType>     SheetnessImageWriterType;
  SheetnessImageWriterType::Pointer writer= SheetnessImageWriterType::New();
  std::cout<< "Writing out sheetness measure image" << std::endl;
  writer->SetFileName(argv[2]);
  writer->SetInput(filterSheetness->GetOutput());
  writer->Update();

  //Generate a binary image by threshlding the sheetness measure
  typedef itk::ThresholdImageFilter< UltrasoundImageType
 >  ThresholdFilterType;

  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput( filterSheetness->GetOutput();
  thresholdFilter->SetOutsideValue( 0 );
  thresholdFilter->ThresholdBelow ( 180 );
  thresholdFilter->Update();

  //Compute the angle of incidence measure
  typedef itk::AngleOfIncidenceImageFilter
        < UltrasoundImageType, AngleOfIncidenceImageType >   AngleOfIncidenceImageFilterType;

  // Create a sheetness Filter
  AngleOfIncidenceImageFilter::Pointer filterAngleOfIncidence = AngleOfIncidenceImageFilter::New();
  filterAngleOfIncidence->SetInput( thresholdFilter->GetOutput() );

  //Write out the Angle of Incidence image
  typedef AngleOfIncidencesFilterType::OutputImageType AngleOfIncidencesImageType;
  typedef itk::ImageFileWriter<AngleOfIncidencesImageType>     AngleOfIncidencesImageWriterType;
  AngleOfIncidencesImageWriterType::Pointer writer= AngleOfIncidencesImageWriterType::New();
  std::cout<< "Writing out angle of incidence measure image" << std::endl;
  writer->SetFileName(argv[3]);
  writer->SetInput(filterAngleOfIncidence->GetOutput());
  writer->Update();


  return EXIT_SUCCESS;

}
