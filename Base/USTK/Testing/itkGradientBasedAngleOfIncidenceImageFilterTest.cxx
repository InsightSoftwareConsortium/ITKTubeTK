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

#include "itkGradientBasedAngleOfIncidenceImageFilter.h"

#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int itkGradientBasedAngleOfIncidenceImageFilterTest( int argc, char * argv[] )
{
  // Argument parsing.
  if( argc < 7 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: "
              << argv[0]
              << " inputImage"
              << " outputImage"
              << " outputImageWithGradientRecursiveGaussian"
              << " probeType( CURVILINEAR|LINEAR )"
              << " ( originX|beamDirectionX )"
              << " ( originY|beamDirectionY )"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];
  const char * outputImage = argv[2];
  const char * outputImageWithGradientRecursiveGaussian = argv[3];
  const char * probeType = argv[4];
  const char * originX = argv[5];
  const char * originY = argv[6];
  const char * beamDirectionX = argv[5];
  const char * beamDirectionY = argv[6];

  // Types
  enum { Dimension = 2 };

  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // Reader
  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );

  // Calculate the angle of incidence
  typedef itk::GradientBasedAngleOfIncidenceImageFilter< ImageType, ImageType >
    AngleOfIncidenceFilterType;
  AngleOfIncidenceFilterType::Pointer angleOfIncidenceFilter =
    AngleOfIncidenceFilterType::New();
  angleOfIncidenceFilter->SetInput( reader->GetOutput() );
  AngleOfIncidenceFilterType::OriginType probeOrigin;
  std::istringstream istrm;
  istrm.str( originX );
  istrm >> probeOrigin[0];
  istrm.clear();
  istrm.str( originY );
  istrm >> probeOrigin[1];
  istrm.clear();
  AngleOfIncidenceFilterType::BeamDirectionType beamDirection;
  istrm.str( beamDirectionX );
  istrm >> beamDirection[0];
  istrm.clear();
  istrm.str( beamDirectionY );
  istrm >> beamDirection[1];
  const std::string probeTypeStr( probeType );
  if( probeTypeStr.compare( "CURVILINEAR" ) == 0 )
    {
    angleOfIncidenceFilter->SetUltrasoundProbeType( AngleOfIncidenceFilterType::CURVILINEAR );
    angleOfIncidenceFilter->SetUltrasoundProbeOrigin( probeOrigin );
    }
  else if( probeTypeStr.compare( "LINEAR" ) == 0 )
    {
    angleOfIncidenceFilter->SetUltrasoundProbeType( AngleOfIncidenceFilterType::LINEAR );
    angleOfIncidenceFilter->SetUltrasoundProbeBeamDirection( beamDirection );
    }
  else
    {
    std::cerr << "Unknown probe type" << std::endl;
    return EXIT_FAILURE;
    }

  // Writer
  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputImage );
  writer->SetInput( angleOfIncidenceFilter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  // Use a different gradient filter.
  typedef itk::GradientRecursiveGaussianImageFilter<
    AngleOfIncidenceFilterType::OperatorImageType,
    AngleOfIncidenceFilterType::GradientOutputImageType
      >
      GradientRecursiveGaussianFilterType;
  GradientRecursiveGaussianFilterType::Pointer gradientRecursiveGaussianFilter =
    GradientRecursiveGaussianFilterType::New();
  gradientRecursiveGaussianFilter->SetSigma( 0.5 );
  angleOfIncidenceFilter->SetGradientFilter( gradientRecursiveGaussianFilter );
  angleOfIncidenceFilter->SetGradientMagnitudeTolerance( 0.5e-3 );

  writer->SetFileName( outputImageWithGradientRecursiveGaussian );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
