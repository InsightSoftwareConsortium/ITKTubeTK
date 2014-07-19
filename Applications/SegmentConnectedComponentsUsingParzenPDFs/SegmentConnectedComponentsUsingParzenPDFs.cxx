/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "itktubePDFSegmenter.h"
#include "tubeMacro.h"

#include "SegmentConnectedComponentsUsingParzenPDFsCLP.h"

// Get the component type and dimension of the image.
void GetImageInformation( const std::string & fileName,
                          itk::ImageIOBase::IOComponentType & componentType,
                          unsigned int & dimension )
{
  typedef itk::ImageIOBase     ImageIOType;
  typedef itk::ImageIOFactory  ImageIOFactoryType;

  ImageIOType::Pointer imageIO =
    ImageIOFactoryType::CreateImageIO( fileName.c_str(),
                                       ImageIOFactoryType::ReadMode );

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

template< class TInputImage, class TLabelMap >
bool
CheckImageAttributes( const TInputImage * input,
                     const TLabelMap * mask )
{
  assert( input );
  assert( mask );
  const typename TInputImage::PointType inputOrigin = input->GetOrigin();
  const typename TLabelMap::PointType maskOrigin = mask->GetOrigin();
  const typename TInputImage::RegionType inputRegion =
    input->GetLargestPossibleRegion();
  const typename TLabelMap::RegionType maskRegion =
    mask->GetLargestPossibleRegion();
  return (
    inputOrigin.GetVnlVector().is_equal( maskOrigin.GetVnlVector(), 0.01 )
    && inputRegion.GetIndex() == maskRegion.GetIndex()
    && inputRegion.GetSize() == maskRegion.GetSize() );
}

template< class T, unsigned int N, unsigned int dimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;

  typedef T                                        InputPixelType;
  typedef itk::Image< InputPixelType, dimension >  InputImageType;
  typedef itk::Image< float, dimension >           ProbImageType;
  typedef itk::Image< unsigned short, dimension >  LabelMapType;

  typedef itk::Image< float, N >                   PDFImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileReader< LabelMapType >     LabelMapReaderType;
  typedef itk::ImageFileWriter< LabelMapType >     LabelMapWriterType;
  typedef itk::ImageFileWriter< ProbImageType >    ProbImageWriterType;
  typedef itk::ImageFileWriter< PDFImageType >     PDFImageWriterType;
  typedef itk::ImageFileReader< PDFImageType >     PDFImageReaderType;

  typedef itk::tube::PDFSegmenter< InputImageType, N, LabelMapType >
    PDFSegmenterType;
  typename PDFSegmenterType::Pointer pdfSegmenter = PDFSegmenterType::New();

  timeCollector.Start( "LoadData" );

  typename LabelMapReaderType::Pointer  inLabelMapReader =
    LabelMapReaderType::New();
  inLabelMapReader->SetFileName( labelmap.c_str() );
  inLabelMapReader->Update();
  pdfSegmenter->SetLabelMap( inLabelMapReader->GetOutput() );

  typename ImageReaderType::Pointer reader;
  for( unsigned int i = 0; i < N; i++ )
    {
    reader = ImageReaderType::New();
    if( i == 0 )
      {
      reader->SetFileName( inputVolume1.c_str() );
      }
    else if( i == 1 )
      {
      reader->SetFileName( inputVolume2.c_str() );
      }
    else if( i == 2 )
      {
      reader->SetFileName( inputVolume3.c_str() );
      }
    else if( i == 3 )
      {
      reader->SetFileName( inputVolume4.c_str() );
      }
    else
      {
      std::cout << "ERROR: current command line xml file limits"
                << " this filter to 4 input images" << std::endl;
      return 1;
      }
    reader->Update();
    if( !CheckImageAttributes( reader->GetOutput(),
        inLabelMapReader->GetOutput() ) )
      {
      std::cout << "Image attributes of inputVolume" << i+1
        << " and label map do not match.  Please check size, spacing, "
        << "origin." << std::endl;
      return EXIT_FAILURE;
      }
    pdfSegmenter->SetInput( i, reader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  for( unsigned int o = 0; o < objectId.size(); o++ )
    {
    pdfSegmenter->AddObjectId( objectId[o] );
    }
  pdfSegmenter->SetVoidId( voidId );
  pdfSegmenter->SetErodeRadius( erodeRadius );
  pdfSegmenter->SetHoleFillIterations( holeFillIterations );
  pdfSegmenter->SetDilateFirst( dilateFirst );
  if( objectPDFWeight.size() == objectId.size() )
    {
    for( unsigned int i = 0; i < objectPDFWeight.size(); i++ )
      {
      pdfSegmenter->SetObjectPDFWeight( i, objectPDFWeight[i] );
      }
    }
  pdfSegmenter->SetProbabilityImageSmoothingStandardDeviation(
    probSmoothingStdDev );
  pdfSegmenter->SetDraft( draft );
  pdfSegmenter->SetReclassifyNotObjectLabels( reclassifyNotObjectLabels );
  pdfSegmenter->SetReclassifyObjectLabels( reclassifyObjectLabels );
  if( forceClassification )
    {
    pdfSegmenter->SetForceClassification( true );
    }

  if( loadClassPDFBase.size() > 0 )
    {
    unsigned int numClasses = pdfSegmenter->GetNumberOfClasses();
    std::cout << "loading classes" << std::endl;
    for( unsigned int i = 0; i < numClasses; i++ )
      {
      std::string fname = loadClassPDFBase;
      char c[80];
      std::sprintf(c, ".c%u.mha", i );
      fname += std::string( c );
      typename PDFImageReaderType::Pointer pdfImageReader =
        PDFImageReaderType::New();
      pdfImageReader->SetFileName( fname.c_str() );
      pdfImageReader->Update();
      typename PDFImageType::Pointer img = pdfImageReader->GetOutput();
      if( i == 0 )
        {
        std::vector< double > tmpOrigin;
        tmpOrigin.resize( dimension );
        std::vector< double > tmpSpacing;
        tmpSpacing.resize( dimension );
        typename PDFImageType::PointType origin = img->GetOrigin();
        typename PDFImageType::SpacingType spacing = img->GetSpacing();
        for( unsigned int d=0; d<N; ++d )
          {
          tmpOrigin[d] = origin[d];
          tmpSpacing[d] = spacing[d];
          }
        pdfSegmenter->SetBinMin( tmpOrigin );
        pdfSegmenter->SetBinSize( tmpSpacing );
        }
      pdfSegmenter->SetClassPDFImage( i, img );
      }
    pdfSegmenter->ClassifyImages();
    }
  else
    {
    pdfSegmenter->Update();
    pdfSegmenter->ClassifyImages();
    }

  timeCollector.Start( "Save" );

  if( saveClassProbabilityVolumeBase.size() > 0 )
    {
    unsigned int numClasses = pdfSegmenter->GetNumberOfClasses();
    for( unsigned int i = 0; i < numClasses; i++ )
      {
      std::string fname = saveClassProbabilityVolumeBase;
      char c[80];
      std::sprintf(c, ".c%u.mha", i );
      fname += std::string( c );
      typename ProbImageWriterType::Pointer probImageWriter =
        ProbImageWriterType::New();
      probImageWriter->SetFileName( fname.c_str() );
      probImageWriter->SetInput( pdfSegmenter->
        GetClassProbabilityForInput( i ) );
      probImageWriter->Update();
      }
    }

  typename LabelMapWriterType::Pointer writer = LabelMapWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( pdfSegmenter->GetLabelMap() );
  writer->Update();

  if( saveClassPDFBase.size() > 0 )
    {
    unsigned int numClasses = pdfSegmenter->GetNumberOfClasses();
    for( unsigned int i = 0; i < numClasses; i++ )
      {
      std::string fname = saveClassPDFBase;
      char c[80];
      std::sprintf(c, ".c%u.mha", i );
      fname += std::string( c );
      typename PDFImageWriterType::Pointer pdfImageWriter =
        PDFImageWriterType::New();
      pdfImageWriter->SetFileName( fname.c_str() );
      pdfImageWriter->SetInput( pdfSegmenter->GetClassPDFImage( i ) );
      pdfImageWriter->Update();
      }
    }

  timeCollector.Stop( "Save" );

  timeCollector.Report();

  return 0;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  if( objectId.size() < 2)
    {
    std::cout << "Please specify a foreground and a background object."
      << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    unsigned int imageDimension = 3;
    itk::ImageIOBase::IOComponentType imageType;

    GetImageInformation( inputVolume1, imageType, imageDimension );

    int N = 1;
    if( !inputVolume2.empty() )
      {
      ++N;
      if( !inputVolume3.empty() )
        {
        ++N;
        if( !inputVolume4.empty() )
          {
          ++N;
          }
        }
      }

    if( imageDimension == 2 )
      {
      switch( imageType )
        {
        case itk::ImageIOBase::UCHAR:
          if( N == 1 )
            {
            return DoIt<unsigned char, 1, 2>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<unsigned char, 2, 2>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<unsigned char, 3, 2>( argc, argv );
            }
          else
            {
            return DoIt<unsigned char, 4, 2>( argc, argv );
            }
          break;
        case itk::ImageIOBase::USHORT:
          if( N == 1 )
            {
            return DoIt<unsigned short, 1, 2>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<unsigned short, 2, 2>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<unsigned short, 3, 2>( argc, argv );
            }
          else
            {
            return DoIt<unsigned short, 4, 2>( argc, argv );
            }
          break;
        case itk::ImageIOBase::CHAR:
        case itk::ImageIOBase::SHORT:
          if( N == 1 )
            {
            return DoIt<short, 1, 2>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<short, 2, 2>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<short, 3, 2>( argc, argv );
            }
          else
            {
            return DoIt<short, 4, 2>( argc, argv );
            }
          break;
        case itk::ImageIOBase::UINT:
        case itk::ImageIOBase::INT:
        case itk::ImageIOBase::ULONG:
        case itk::ImageIOBase::LONG:
        case itk::ImageIOBase::FLOAT:
        case itk::ImageIOBase::DOUBLE:
          if( N == 1 )
            {
            return DoIt<float, 1, 2>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<float, 2, 2>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<float, 3, 2>( argc, argv );
            }
          else
            {
            return DoIt<float, 4, 2>( argc, argv );
            }
          break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown image type" << std::endl;
          break;
        }
      }
    else if( imageDimension == 3 )
      {
      switch( imageType )
        {
        case itk::ImageIOBase::UCHAR:
          if( N == 1 )
            {
            return DoIt<unsigned char, 1, 3>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<unsigned char, 2, 3>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<unsigned char, 3, 3>( argc, argv );
            }
          else
            {
            return DoIt<unsigned char, 4, 3>( argc, argv );
            }
          break;
        case itk::ImageIOBase::USHORT:
          if( N == 1 )
            {
            return DoIt<unsigned short, 1, 3>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<unsigned short, 2, 3>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<unsigned short, 3, 3>( argc, argv );
            }
          else
            {
            return DoIt<unsigned short, 4, 3>( argc, argv );
            }
          break;
        case itk::ImageIOBase::CHAR:
        case itk::ImageIOBase::SHORT:
          if( N == 1 )
            {
            return DoIt<short, 1, 3>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<short, 2, 3>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<short, 3, 3>( argc, argv );
            }
          else
            {
            return DoIt<short, 4, 3>( argc, argv );
            }
          break;
        case itk::ImageIOBase::UINT:
        case itk::ImageIOBase::INT:
        case itk::ImageIOBase::ULONG:
        case itk::ImageIOBase::LONG:
        case itk::ImageIOBase::FLOAT:
        case itk::ImageIOBase::DOUBLE:
          if( N == 1 )
            {
            return DoIt<float, 1, 3>( argc, argv );
            }
          else if( N == 2 )
            {
            return DoIt<float, 2, 3>( argc, argv );
            }
          else if( N == 3 )
            {
            return DoIt<float, 3, 3>( argc, argv );
            }
          else
            {
            return DoIt<float, 4, 3>( argc, argv );
            }
          break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown image type" << std::endl;
          return EXIT_FAILURE;
          break;
        }
      }
    else
      {
      std::cout << "Only 2 and 3 dimensional images supported."
        << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
