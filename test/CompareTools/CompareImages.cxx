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
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"

#include "itktubeDifferenceImageFilter.h"

#include "metaCommand.h"
#include "metaUtils.h"

#include <iostream>
#include <fstream>
#include <sstream>

#define ITK_TEST_DIMENSION_MAX 6

int RegressionTestImage( const char *, const char *, int, bool, double, int, 
  int, double );

int main( int argc, char **argv )
{
  if( argc < 5 )
    {
    std::cerr << "Usage:" << std::endl;
    std::cerr << "-t testImage" << std::endl;
    std::cerr << "-b baselineImage1[,baselineImage2[,baselineImage3[, ...]"
      << std::endl;
    std::cerr << "    If more than one baselineImage, test will pass if any"
      << std::endl;
    std::cerr << "      of them match the testImage." << std::endl;
    std::cerr << "-i toleranceIntensity" << std::endl;
    std::cerr << "-r toleranceRadius" << std::endl;
    std::cerr << "-n toleranceNumberOfPixels" << std::endl;
    std::cerr << "-c toleranceCoordinates" << std::endl;
    return -1;
    }
  int bestBaselineStatus = 2001;

  // Process some command-line arguments intended for BatchMake
  MetaCommand command;

  // Option for setting the tolerable difference in intensity values
  // between the two images.
  command.SetOption( "toleranceIntensity", "i", false, 
    "Acceptable differences in pixels intensity" );
  command.AddOptionField( "toleranceIntensity", "value", MetaCommand::FLOAT,
    true );

  // Option for setting the radius of the neighborhood around a pixel
  // to search for similar intensity values.
  command.SetOption( "toleranceRadius", "r", false, 
    "Neighbor pixels to look for similar values" );
  command.AddOptionField( "toleranceRadius", "value", MetaCommand::INT, true );

  // Option for setting the number of pixel that can be tolerated to
  // have different intensities.
  command.SetOption( "toleranceNumberOfPixels", "n", false, 
    "Number of Pixels that are acceptable to have intensity differences" );
  command.AddOptionField( "toleranceNumberOfPixels", "value", MetaCommand::INT,
    true );

  // Option for setting the number of pixel that can be tolerated to
  // have different intensities.
  command.SetOption( "toleranceCoordinates", "c", false,
    "Acceptable differences in image origin/spacing" );
  command.AddOptionField( "toleranceCoordinates", "value", MetaCommand::FLOAT,
    true );

  // Option for setting the filename of the test image.
  command.SetOption( "testImage", "t", true, 
    "Filename of the image to be tested against the baseline images" );
  command.AddOptionField( "testImage", "filename", MetaCommand::STRING, true );

  // Option for setting the filename of a single baseline image.
  command.SetOption( "baselineImages", "b", true, 
    "Baseline images filename" );
  command.AddOptionField( "baselineImages", "filenames", MetaCommand::STRING,
    true );


  command.Parse( argc, argv );


  double toleranceIntensity = 0.0;
  unsigned int  toleranceRadius = 0;
  unsigned long toleranceNumberOfPixels = 0;
  double toleranceCoordinates = 1e-6;
  std::string testImageFilename;

  // If a value of intensity tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceIntensity" ) )
    {
    toleranceIntensity =
      command.GetValueAsFloat( "toleranceIntensity", "value" );
    }

  // If a value of neighborhood radius tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceRadius" ) )
    {
    toleranceRadius =
      command.GetValueAsInt( "toleranceRadius", "value" );
    }

  // If a value of number of pixels tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceNumberOfPixels" ) )
    {
    toleranceNumberOfPixels =
      command.GetValueAsInt( "toleranceNumberOfPixels", "value" );
    }

  // If a value of coordinates tolerance was given in the command line
  if( command.GetOptionWasSet( "toleranceCoordinates" ) )
    {
    toleranceCoordinates =
      command.GetValueAsFloat( "toleranceCoordinates", "value" );
    }

  // Get the filename of the image to be tested
  if( command.GetOptionWasSet( "testImage" ) )
    {
    testImageFilename =
      command.GetValueAsString( "testImage", "filename" );
    }

  std::vector< std::string > baselineImageFilenames;
  baselineImageFilenames.clear();

  // Get the filename of the base line image
  std::string filenames = command.GetValueAsString( "baselineImages", "filenames" );
  MET_StringToVector( filenames, baselineImageFilenames );

  std::string bestBaselineFilename;

  try
    {
    typedef std::vector< std::string >::const_iterator  NameIteratorType;
    NameIteratorType baselineImageItr = baselineImageFilenames.begin();
    while( baselineImageItr != baselineImageFilenames.end() )
      {
      const int currentStatus =
        RegressionTestImage( 
            testImageFilename.c_str(), baselineImageItr->c_str(), 
            0, false, toleranceIntensity, toleranceRadius, 
            toleranceNumberOfPixels, toleranceCoordinates );
      if( currentStatus < bestBaselineStatus )
        {
        bestBaselineStatus = currentStatus;
        bestBaselineFilename = *baselineImageItr;
        }
      if( bestBaselineStatus == 0 )
        {
        break;
        }
      ++baselineImageItr;
      }
    // generate images of our closest match
    if( bestBaselineStatus == 0 )
      {
      RegressionTestImage( 
        testImageFilename.c_str(), 
        bestBaselineFilename.c_str(), 1, false, 
        toleranceIntensity, toleranceRadius, toleranceNumberOfPixels,
        toleranceCoordinates );
      }
    else
      {
      RegressionTestImage( 
        testImageFilename.c_str(), 
        bestBaselineFilename.c_str(), 1, true, 
        toleranceIntensity, toleranceRadius, toleranceNumberOfPixels,
        toleranceCoordinates );
      }
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "ITK test driver caught an ITK exception:\n";
    std::cerr << e << "\n";
    bestBaselineStatus = -1;
    }
  catch( const std::exception& e )
    {
    std::cerr << "ITK test driver caught an exception:\n";
    std::cerr << e.what() << "\n";
    bestBaselineStatus = -1;
    }
  catch( ... )
    {
    std::cerr << "ITK test driver caught an unknown exception!!!\n";
    bestBaselineStatus = -1;
    }
  std::cout << bestBaselineStatus << std::endl;
  return bestBaselineStatus;
}

// Regression Testing Code
int RegressionTestImage ( const char *testImageFilename,
  const char *baselineImageFilename, int reportErrors,
  bool createDifferenceImage, double intensityTolerance, 
  int radiusTolerance, int numberOfPixelsTolerance,
  double coordinatesTolerance )
{
  // Use the factory mechanism to read the test and baseline files and
  //   convert them to double
  typedef itk::Image<double, ITK_TEST_DIMENSION_MAX>         ImageType;
  typedef itk::Image<unsigned char, ITK_TEST_DIMENSION_MAX>  OutputType;
  typedef itk::Image<unsigned char, 2>                       DiffOutputType;
  typedef itk::ImageFileReader<ImageType>                    ReaderType;

  // Read the baseline file
  ReaderType::Pointer baselineReader = ReaderType::New();
    baselineReader->SetFileName( baselineImageFilename );
  try
    {
    baselineReader->UpdateLargestPossibleRegion();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Exception detected while reading " << baselineImageFilename
      << " : "  << e;
    return 1000;
    }

  // Read the file generated by the test
  ReaderType::Pointer testReader = ReaderType::New();
    testReader->SetFileName( testImageFilename );
  try
    {
    testReader->UpdateLargestPossibleRegion();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Exception detected while reading " << testImageFilename
      << " : "  << e << std::endl;
    return 1000;
    }

  // The sizes of the baseline and test image must match
  ImageType::SizeType baselineSize;
    baselineSize = baselineReader->GetOutput()->
      GetLargestPossibleRegion().GetSize();
  ImageType::SizeType testSize;
    testSize = testReader->GetOutput()->GetLargestPossibleRegion().GetSize();

  if ( baselineSize != testSize )
    {
    std::cerr << "The size of the Baseline image and Test image do not match!"
      << std::endl;
    std::cerr << "Baseline image: " << baselineImageFilename
              << " has size " << baselineSize << std::endl;
    std::cerr << "Test image:     " << testImageFilename
              << " has size " << testSize << std::endl;
    return 1;
    }

  // Now compare the two images
  typedef itk::tube::DifferenceImageFilter<ImageType, ImageType> DiffType;
  DiffType::Pointer diff = DiffType::New();
    diff->SetValidInput( baselineReader->GetOutput() );
    diff->SetTestInput( testReader->GetOutput() );

    diff->SetDifferenceThreshold( intensityTolerance );
    diff->SetToleranceRadius( radiusTolerance );
    diff->SetCoordinateTolerance( coordinatesTolerance );

    diff->UpdateLargestPossibleRegion();

  bool differenceFailed = false;

  double averageIntensityDifference = diff->GetTotalDifference();

  unsigned long numberOfPixelsWithDifferences =
    diff->GetNumberOfPixelsWithDifferences();

  if( averageIntensityDifference > 0.0 )
    {
    if( static_cast<int>( numberOfPixelsWithDifferences ) >
      numberOfPixelsTolerance )
      {
      differenceFailed = true;
      }
    else
      {
      differenceFailed = false;
      }
    }
  else
    {
    differenceFailed = false;
    }

  if ( reportErrors )
    {
    typedef itk::RescaleIntensityImageFilter<ImageType, OutputType>
      RescaleType;
    typedef itk::ExtractImageFilter<OutputType, DiffOutputType>
      ExtractType;
    typedef itk::ImageFileWriter<DiffOutputType>
      WriterType;
    typedef itk::ImageRegion<ITK_TEST_DIMENSION_MAX>
      RegionType;

    OutputType::IndexType index; index.Fill( 0 );
    OutputType::SizeType size; size.Fill( 0 );

    RescaleType::Pointer rescale = RescaleType::New();

    rescale->SetOutputMinimum(
      itk::NumericTraits<unsigned char>::NonpositiveMin() );
    rescale->SetOutputMaximum( itk::NumericTraits<unsigned char>::max() );
    rescale->SetInput( diff->GetOutput() );
    rescale->UpdateLargestPossibleRegion();

    RegionType region;
    region.SetIndex( index );

    size = rescale->GetOutput()->GetLargestPossibleRegion().GetSize();
    for ( unsigned int i = 2; i < ITK_TEST_DIMENSION_MAX; i++ )
      {
      size[i] = 0;
      }
    region.SetSize( size );

    ExtractType::Pointer extract = ExtractType::New();
    extract->SetInput( rescale->GetOutput() );
    extract->SetExtractionRegion( region );
    extract->SetDirectionCollapseToIdentity();

    WriterType::Pointer writer = WriterType::New();
      writer->SetInput( extract->GetOutput() );
    if( createDifferenceImage )
      {
      // if there are discrepencies, create an diff image
      std::cout
        << "<DartMeasurement name=\"ImageError\" type=\"numeric/double\">";
      std::cout << averageIntensityDifference;
      std::cout <<  "</DartMeasurement>" << std::endl;

      std::cout
        << "<DartMeasurement name=\"NumberOfPixelsError\" type=\"numeric/int\">";
      std::cout << numberOfPixelsWithDifferences;
      std::cout <<  "</DartMeasurement>" << std::endl;

      std::ostringstream diffName;
      diffName << testImageFilename << ".diff.png";
      try
        {
        rescale->SetInput( diff->GetOutput() );
        rescale->Update();
        }
      catch( const std::exception& e )
        {
        std::cerr << "Error during rescale of " << diffName.str() << std::endl;
        std::cerr << e.what() << "\n";
        }
      catch ( ... )
        {
        std::cerr << "Error during rescale of " << diffName.str() << std::endl;
        }
      writer->SetFileName( diffName.str().c_str() );
      try
        {
        writer->Update();
        }
      catch( const std::exception& e )
        {
        std::cerr << "Error during write of " << diffName.str() << std::endl;
        std::cerr << e.what() << "\n";
        }
      catch ( ... )
        {
        std::cerr << "Error during write of " << diffName.str() << std::endl;
        }

      std::cout << "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
      std::cout << diffName.str();
      std::cout << "</DartMeasurementFile>" << std::endl;
      }
    std::ostringstream baseName;
    baseName << testImageFilename << ".base.png";
    try
      {
      rescale->SetInput( baselineReader->GetOutput() );
      rescale->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during rescale of " << baseName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during rescale of " << baseName.str() << std::endl;
      }
    try
      {
      writer->SetFileName( baseName.str().c_str() );
      writer->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during write of " << baseName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during write of " << baseName.str() << std::endl;
      }

    std::cout
      << "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
    std::cout << baseName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;

    std::ostringstream testName;
    testName << testImageFilename << ".test.png";
    try
      {
      rescale->SetInput( testReader->GetOutput() );
      rescale->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during rescale of " << testName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during rescale of " << testName.str() << std::endl;
      }
    try
      {
      writer->SetFileName( testName.str().c_str() );
      writer->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during write of " << testName.str() << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch ( ... )
      {
      std::cerr << "Error during write of " << testName.str() << std::endl;
      }

    std::cout << "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
    std::cout << testName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;


    }
  return differenceFailed;
}
