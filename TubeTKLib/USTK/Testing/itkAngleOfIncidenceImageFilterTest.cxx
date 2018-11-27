/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
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

#include "itkAngleOfIncidenceImageFilter.h"
#include "itktubeSheetnessMeasureImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkThresholdImageFilter.h>

int itkAngleOfIncidenceImageFilterTest( int argc, char * argv[] )
{
  // Argument parsing.
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " Input_Ultrasound_Image Output_Sheetness_Image Output_Angle_Of_Incidence_Image UltrasoundOriginX UltrasoundOrignY UltrasoundOriginZ [Sheetness_Threshold]" << std::endl;
    return EXIT_FAILURE;
    }

  // Types.
  enum { Dimension = 3 };

  typedef unsigned short                              UltrasoundPixelType;
  typedef itk::Image< UltrasoundPixelType, Dimension > UltrasoundImageType;

  typedef float                                              AngleOfIncidencePixelType;
  typedef itk::Image< AngleOfIncidencePixelType, Dimension > AngleOfIncidenceImageType;

  // Reader,
  typedef itk::ImageFileReader< UltrasoundImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  // Declare the type for the Hessian filter
  typedef itk::HessianRecursiveGaussianImageFilter< UltrasoundImageType >  HessianFilterType;

  // Declare the type for the sheetness measure filter
  typedef itk::tube::SheetnessMeasureImageFilter< float >  SheetnessMeasureImageFilterType;


  // Create a Hessian Filter
  HessianFilterType::Pointer filterHessian = HessianFilterType::New();

  // Create a sheetness Filter
  SheetnessMeasureImageFilterType::Pointer filterSheetness = SheetnessMeasureImageFilterType::New();


  // Connect the input images
  filterHessian->SetInput( reader->GetOutput() );
  filterSheetness->SetInput( filterHessian->GetOutput() );

  // Select the value of Sigma
  filterHessian->SetSigma( 0.5 );

  // Execute the filter
  std::cout << "Generate sheetness measure" << std::endl;
  filterSheetness->Update();

  //Write out the sheetness image
  typedef SheetnessMeasureImageFilterType::OutputImageType SheetnessImageType;

  typedef itk::ImageFileWriter<SheetnessImageType>     SheetnessImageWriterType;
  SheetnessImageWriterType::Pointer writer= SheetnessImageWriterType::New();
  std::cout<< "Writing out sheetness measure image" << std::endl;
  writer->SetFileName( argv[2] );
  writer->SetInput( filterSheetness->GetOutput() );
  writer->Update();

  //Generate a binary image by thresholding the sheetness measure
  typedef itk::ThresholdImageFilter< SheetnessImageType >  ThresholdFilterType;

  ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput( filterSheetness->GetOutput() );
  thresholdFilter->SetOutsideValue( 0 );

  double sheetnessThresholdValue = 0.1;
  if( argc > 7 )
    {
    sheetnessThresholdValue = std::atof( argv[7] );
    }
  thresholdFilter->ThresholdBelow ( sheetnessThresholdValue );
  thresholdFilter->Update();

  //Compute the angle of incidence measure
  typedef itk::AngleOfIncidenceImageFilter
        < SheetnessImageType, AngleOfIncidenceImageType >   AngleOfIncidenceImageFilterType;

  // Create a sheetness Filter
  AngleOfIncidenceImageFilterType::Pointer filterAngleOfIncidence = AngleOfIncidenceImageFilterType::New();

  //Read in the ultrasound probe origin ( X,Y,Z )
  double UltrasoundProbeOriginX = std::atof( argv[4] );
  double UltrasoundProbeOriginY = std::atof( argv[5] );
  double UltrasoundProbeOriginZ = std::atof( argv[6] );

  itk::Vector< double, 3 > UltrasoundProbeOriginVector;

  UltrasoundProbeOriginVector[0] = UltrasoundProbeOriginX;
  UltrasoundProbeOriginVector[1] = UltrasoundProbeOriginY;
  UltrasoundProbeOriginVector[2] = UltrasoundProbeOriginZ;

  filterAngleOfIncidence->SetUltrasoundProbeOrigin( UltrasoundProbeOriginVector );
  filterAngleOfIncidence->SetInput( thresholdFilter->GetOutput() );

  //Write out the Angle of Incidence image
  typedef AngleOfIncidenceImageFilterType::OutputImageType AngleOfIncidencesImageType;
  typedef itk::ImageFileWriter<AngleOfIncidencesImageType>     AngleOfIncidencesImageWriterType;
  AngleOfIncidencesImageWriterType::Pointer angleOfIncidenceWriter= AngleOfIncidencesImageWriterType::New();
  std::cout<< "Writing out angle of incidence measure image" << std::endl;
  angleOfIncidenceWriter->SetFileName( argv[3] );
  angleOfIncidenceWriter->SetInput( filterAngleOfIncidence->GetOutput() );
  angleOfIncidenceWriter->Update();


  return EXIT_SUCCESS;
}
