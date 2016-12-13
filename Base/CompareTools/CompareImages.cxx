/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: CompareImages.cxx, v $
  Language:  C++
  Date:      $Date: 2008-11-09 18:18:52 $
  Version:   $Revision: 1.12 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itktubeDifferenceImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>

#include "CompareImagesCLP.h"

// Description:
// Get the ComponentType and dimension of the image
void GetImageInformation( std::string fileName,
                          itk::ImageIOBase::IOComponentType &componentType,
                          unsigned int & dimension )
{
  // Find out the component type of the image in file
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( fileName.c_str(),
                                        itk::ImageIOFactory::ReadMode );
  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to
  // read the meta data from the image file.
  imageIO->SetFileName( fileName.c_str() );
  imageIO->ReadImageInformation();

  componentType = imageIO->GetComponentType();
  dimension = imageIO->GetNumberOfDimensions();
}

template< unsigned int VDimension >
int RegressionTestImage( const char *, const char *, bool, bool, double,
  int, int, double, bool, const char * );

int main( int argc, char * argv[] )
{
  int bestBaselineStatus = 2001;
  PARSE_ARGS;

  std::string bestBaselineFilename;
  itk::ImageIOBase::IOComponentType cType = itk::ImageIOBase::UCHAR;
  unsigned int dims = 0;

  GetImageInformation( testImageFilename, cType, dims );

  try
    {
    typedef std::vector< std::string >::const_iterator  nameIterator;
    nameIterator baselineImageItr = baselineImageFilenames.begin();
    while( baselineImageItr != baselineImageFilenames.end() )
      {
      int currentStatus = 0;
      if( dims == 2 )
        {
        currentStatus = RegressionTestImage<2>( 
          testImageFilename.c_str(), baselineImageItr->c_str(),
          false, false, toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          false, outputImageFilename.c_str() );
        }
      else
        {
        currentStatus = RegressionTestImage<3>( 
          testImageFilename.c_str(), baselineImageItr->c_str(),
          false, false, toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          false, outputImageFilename.c_str() );
        }
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

    bool writeOutputImage = false;
    // generate images of our closest match
    if( bestBaselineStatus == 0 )
      {
      if( dims == 2 )
        {
        RegressionTestImage<2>( 
          testImageFilename.c_str(),
          bestBaselineFilename.c_str(), true, false,
          toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          writeOutputImage, outputImageFilename.c_str() );
        }
      else
        {
        RegressionTestImage<3>( 
          testImageFilename.c_str(),
          bestBaselineFilename.c_str(), true, false,
          toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          writeOutputImage, outputImageFilename.c_str() );
        }
      }
    else
      {
      if( dims == 2 )
        {
        RegressionTestImage<2>( 
          testImageFilename.c_str(),
          bestBaselineFilename.c_str(), true, true,
          toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          writeOutputImage, outputImageFilename.c_str() );
        }
      else
        {
        RegressionTestImage<3>( 
          testImageFilename.c_str(),
          bestBaselineFilename.c_str(), true, true,
          toleranceIntensity, toleranceRadius,
          toleranceNumberOfPixels, toleranceCoordinates,
          writeOutputImage, outputImageFilename.c_str() );
        }
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
template< unsigned int VDimension >
int
RegressionTestImage( const char *testImageFilename,
  const char *baselineImageFilename,
  bool reportErrors,
  bool createDifferenceImage,
  double intensityTolerance,
  int radiusTolerance,
  int numberOfPixelsTolerance,
  double coordinateTolerance,
  bool writeOutputImage,
  const char *outputImageFilename )
{
  // Use the factory mechanism to read the test and baseline files and
  //  convert them to double
  typedef itk::Image< double, VDimension >         ImageType;
  typedef itk::Image< unsigned char, VDimension >  OutputType;
  typedef itk::Image< unsigned char, 2 >           DiffOutputType;
  typedef itk::ImageFileReader< ImageType >        ReaderType;

  // Read the baseline file
  typename ReaderType::Pointer baselineReader = ReaderType::New();
  baselineReader->SetFileName( baselineImageFilename );
  try
    {
    baselineReader->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception detected while reading "
      << baselineImageFilename << " : "  << e;
    return 1000;
    }

  // Read the file generated by the test
  typename ReaderType::Pointer testReader = ReaderType::New();
  testReader->SetFileName( testImageFilename );
  try
    {
    testReader->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception detected while reading " << testImageFilename
      << " : "  << e << std::endl;
    return 1000;
    }

  // The sizes of the baseline and test image must match
  typename ImageType::SizeType baselineSize;
  baselineSize = baselineReader->GetOutput()->GetLargestPossibleRegion().
    GetSize();
  typename ImageType::SizeType testSize;
  testSize = testReader->GetOutput()->GetLargestPossibleRegion().GetSize();

  if( baselineSize != testSize )
    {
    std::cerr <<
      "The size of the Baseline image and Test image do not match!"
      << std::endl;
    std::cerr << "Baseline image: " << baselineImageFilename
      << " has size " << baselineSize << std::endl;
    std::cerr << "Test image:     " << testImageFilename
      << " has size " << testSize << std::endl;
    return 1;
    }

  // Now compare the two images
  typedef itk::tube::DifferenceImageFilter<ImageType, ImageType> DiffType;
  typename DiffType::Pointer diff = DiffType::New();
  diff->SetValidInput( baselineReader->GetOutput() );
  diff->SetTestInput( testReader->GetOutput() );

  diff->SetDifferenceThreshold( intensityTolerance );
  diff->SetToleranceRadius( radiusTolerance );
  diff->SetCoordinateTolerance( coordinateTolerance );

  diff->UpdateLargestPossibleRegion();

  bool differenceFailed = false;
  double averageIntensityDifference = diff->GetMeanDifference();
  unsigned long numberOfPixelsWithDifferences = diff->
    GetNumberOfPixelsWithDifferences();


  if( static_cast<int>( numberOfPixelsWithDifferences ) >
    numberOfPixelsTolerance )
    {
    differenceFailed = true;
    }

  if( writeOutputImage )
    {
    typedef itk::ImageFileWriter< ImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputImageFilename );
    writer->SetInput( diff->GetOutput() );
    writer->Update();
    }

  if( reportErrors )
    {
    typedef itk::RescaleIntensityImageFilter< ImageType, OutputType >
      RescaleType;
    typedef itk::ExtractImageFilter< OutputType, DiffOutputType >
      ExtractType;
    typedef itk::ImageFileWriter< DiffOutputType >
      WriterType;
    typedef itk::ImageRegion< VDimension >
      RegionType;

    typename OutputType::IndexType index;
    index.Fill( 0 );

    typename OutputType::SizeType size;
    size.Fill( 0 );

    typename RescaleType::Pointer rescale = RescaleType::New();

    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 255 );
    rescale->SetInput( diff->GetOutput() );
    rescale->Update();

    RegionType region;
    region.SetIndex( index );

    size = rescale->GetOutput()->GetLargestPossibleRegion().GetSize();
    for( unsigned int i = 2; i < VDimension; i++ )
      {
      size[i] = 0;
      }
    region.SetSize( size );

    typename ExtractType::Pointer extract = ExtractType::New();
    extract->SetDirectionCollapseToIdentity();

    extract->SetInput( rescale->GetOutput() );
    extract->SetExtractionRegion( region );

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( extract->GetOutput() );
    if( createDifferenceImage )
      {
      // if there are discrepencies, create an diff image
      std::cout <<
        "<DartMeasurement name=\"ImageError\" type=\"numeric/double\">";
      std::cout << averageIntensityDifference;
      std::cout <<  "</DartMeasurement>" << std::endl;

      std::cout <<
        "<DartMeasurement name=\"NumberOfPixelsError\" type=\"numeric/int\">";
      std::cout << numberOfPixelsWithDifferences;
      std::cout <<  "</DartMeasurement>" << std::endl;

      std::ostringstream diffName;
      diffName << testImageFilename << ".diff.png";
      try
        {
        rescale->SetOutputMinimum( 0 );
        rescale->SetOutputMaximum( 255 );
        rescale->SetInput( diff->GetOutput() );
        rescale->Update();
        }
      catch( const std::exception& e )
        {
        std::cerr << "Error during rescale of " << diffName.str()
          << std::endl;
        std::cerr << e.what() << "\n";
        }
      catch( ... )
        {
        std::cerr << "Error during rescale of " << diffName.str()
          << std::endl;
        }
      writer->SetFileName( diffName.str().c_str() );
      try
        {
        writer->Update();
        }
      catch( const std::exception& e )
        {
        std::cerr << "Error during write of " << diffName.str()
          << std::endl;
        std::cerr << e.what() << "\n";
        }
      catch( ... )
        {
        std::cerr << "Error during write of " << diffName.str()
          << std::endl;
        }

      std::cout <<
        "<DartMeasurementFile name=\"DifferenceImage\" type=\"image/png\">";
      std::cout << diffName.str();
      std::cout << "</DartMeasurementFile>" << std::endl;
      }
    std::ostringstream baseName;
    baseName << testImageFilename << ".base.png";
    try
      {
      rescale->SetOutputMinimum( 0 );
      rescale->SetOutputMaximum( 255 );
      rescale->SetInput( baselineReader->GetOutput() );
      rescale->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during rescale of " << baseName.str()
        << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch( ... )
      {
      std::cerr << "Error during rescale of " << baseName.str()
        << std::endl;
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
    catch( ... )
      {
      std::cerr << "Error during write of " << baseName.str() << std::endl;
      }

    std::cout <<
      "<DartMeasurementFile name=\"BaselineImage\" type=\"image/png\">";
    std::cout << baseName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;

    std::ostringstream testName;
    testName << testImageFilename << ".test.png";
    try
      {
      rescale->SetOutputMinimum( 0 );
      rescale->SetOutputMaximum( 255 );
      rescale->SetInput( testReader->GetOutput() );
      rescale->Update();
      }
    catch( const std::exception& e )
      {
      std::cerr << "Error during rescale of " << testName.str()
        << std::endl;
      std::cerr << e.what() << "\n";
      }
    catch( ... )
      {
      std::cerr << "Error during rescale of " << testName.str()
        << std::endl;
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
    catch( ... )
      {
      std::cerr << "Error during write of " << testName.str() << std::endl;
      }

    std::cout <<
      "<DartMeasurementFile name=\"TestImage\" type=\"image/png\">";
    std::cout << testName.str();
    std::cout << "</DartMeasurementFile>" << std::endl;


    }
  return differenceFailed;
}
