/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkTimeProbesCollectorBase.h>

#include <itkExtractImageFilter.h>

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <metaUtils.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "Convert4DImageTo3DImagesCLP.h"

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "Convert4DImageTo3DImages",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                    InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;

  typedef TPixel                                    OutputPixelType;
  typedef itk::Image< OutputPixelType, 3 >          OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType  >  WriterType;

  typedef itk::ExtractImageFilter< InputImageType, OutputImageType  >
                                                    FilterType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );

  double progress = 0.1;
  progressReporter.Report( progress );

  typename InputImageType::Pointer inIm = reader->GetOutput();

  typename InputImageType::RegionType inImRegion =
    inIm->GetLargestPossibleRegion();
  typename InputImageType::SizeType inImSize = inImRegion.GetSize();
  typename InputImageType::SizeType inImToOutImSize = inImSize;
  inImToOutImSize[3] = 0;
  typename InputImageType::IndexType inImToOutImIndex = inImRegion.GetIndex();
  std::cout << "Input size = " << inImSize << std::endl;
  std::cout << "Input To Output size = " << inImToOutImSize << std::endl;

  typename InputImageType::RegionType inImToOutImRegion =
    inIm->GetLargestPossibleRegion();
  inImToOutImRegion.SetSize( inImToOutImSize );

  progress = 0.2;
  progressReporter.Report( progress );
  for( unsigned int d4 = 0; d4 < inImSize[3]; ++d4 )
    {
    timeCollector.Start( "Save data" );

    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( inIm );
    inImToOutImRegion.SetIndex( inImToOutImIndex );
    ++inImToOutImIndex[3];
    filter->SetExtractionRegion( inImToOutImRegion );
    filter->SetDirectionCollapseToSubmatrix();
    filter->Update();
    typename OutputImageType::Pointer outIm = filter->GetOutput();

    progress = 0.8 * ( d4 / (double)inImSize[3] ) + 0.2;
    progressReporter.Report( progress );

    std::string outputImageFileName = inputImageFileName;
    std::ostringstream zeroPadNum;
    zeroPadNum << std::setw(3) << std::setfill( '0' ) << d4;
    std::string newSuffix;
    newSuffix = zeroPadNum.str() + ".mha";
    MET_SetFileSuffix( outputImageFileName, newSuffix );

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputImageFileName );
    writer->SetInput( filter->GetOutput() );
    writer->SetUseCompression( true );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
        + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Save data" );
    }

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
  }

namespace tube
{

// Get the component type and dimension of the image.
void GetImageInformation( const std::string & fileName,
                          itk::ImageIOBase::IOComponentEnum & componentType,
                          unsigned int & dimension )
{
  typedef itk::ImageIOBase     ImageIOType;
  typedef itk::ImageIOFactory  ImageIOFactoryType;

  ImageIOType::Pointer imageIO =
    ImageIOFactoryType::CreateImageIO( fileName.c_str(),
      itk::IOFileModeEnum::ReadMode );

  if( imageIO )
    {
    // Read the metadata from the image file.
    imageIO->SetFileName( fileName.c_str() );
    imageIO->ReadImageInformation();

    componentType = imageIO->GetComponentType();
    dimension = imageIO->GetNumberOfDimensions();
    }
  else
    {
    tubeErrorMacro( << "No ImageIO was found." );
    }
}

}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef itk::ImageIOBase              ImageIOType;
  typedef ImageIOType::IOComponentEnum  IOComponentType;

  IOComponentType componentType = IOComponentType::UNKNOWNCOMPONENTTYPE;
  unsigned int dimension = 0;
  try
    {
    tube::GetImageInformation( inputImageFileName, componentType, dimension );

    if( dimension == 4 )
      {
      switch( componentType )
        {
        case IOComponentType::UCHAR:
          return DoIt< unsigned char, 4 >( argc, argv );
        case IOComponentType::CHAR:
          return DoIt< char, 4 >( argc, argv );
        case IOComponentType::USHORT:
          return DoIt< unsigned short, 4 >( argc, argv );
        case IOComponentType::SHORT:
          return DoIt< short, 4 >( argc, argv );
        case IOComponentType::FLOAT:
          return DoIt< float, 4 >( argc, argv );
        case IOComponentType::DOUBLE:
          return DoIt< double, 4 >( argc, argv );
        case IOComponentType::INT:
          return DoIt< int, 4 >( argc, argv );
        case IOComponentType::UINT:
          return DoIt< unsigned int, 4 >( argc, argv );
        case IOComponentType::UNKNOWNCOMPONENTTYPE:
        default:
          tubeErrorMacro( << "Unknown component type." );
          return EXIT_FAILURE;
        }
      }
    tubeErrorMacro( << "Dimension size of " << dimension << " not supported" );
    return EXIT_FAILURE;
    }
  catch( itk::ExceptionObject & ex )
    {
    tubeErrorMacro( << "ITK exception caught. " << ex );
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    tubeErrorMacro( << "Exception caught." );
    return EXIT_FAILURE;
    }

  return EXIT_FAILURE;
}
