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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkPDFSegmenter.h"

#include "PDFSegmenterCLP.h"

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbesCollectorBase.h"

// Description:
// Get the PixelType and ComponentType from fileName
void GetImageType (std::string fileName,
  itk::ImageIOBase::IOPixelType &pixelType,
  itk::ImageIOBase::IOComponentType &componentType)
{
  typedef itk::Image<short, 3> ImageType;
  itk::ImageFileReader<ImageType>::Pointer imageReader =
        itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(fileName.c_str());
  imageReader->UpdateOutputInformation();

  pixelType = imageReader->GetImageIO()->GetPixelType();
  componentType = imageReader->GetImageIO()->GetComponentType();
}

// Description:
// Get the PixelTypes and ComponentTypes from fileNames
void GetImageTypes (std::vector<std::string> fileNames,
  std::vector<itk::ImageIOBase::IOPixelType> &pixelTypes,
  std::vector<itk::ImageIOBase::IOComponentType> &componentTypes)
{
  pixelTypes.clear();
  componentTypes.clear();

  // For each file, find the pixel and component type
  for (std::vector<std::string>::size_type i = 0; i < fileNames.size(); i++)
    {
    itk::ImageIOBase::IOPixelType pixelType;
    itk::ImageIOBase::IOComponentType componentType;
    GetImageType (fileNames[i],
                  pixelType,
                  componentType);
    pixelTypes.push_back(pixelType);
    componentTypes.push_back(componentType);
    }
}

template <class T, unsigned int N>
int DoIt( int argc, char *argv[] )
{
  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;

  typedef T                                        InputPixelType;
  typedef itk::OrientedImage< InputPixelType, 3 >  InputImageType;
  typedef itk::OrientedImage< unsigned short, 3 >  MaskImageType;
  typedef itk::OrientedImage< float, 3 >           ProbImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< MaskImageType >    MaskWriterType;
  typedef itk::ImageFileWriter< ProbImageType >    ProbImageWriterType;

  typedef itk::PDFSegmenter< InputImageType, N, MaskImageType >
    PDFSegmenterType;
  typename PDFSegmenterType::Pointer pdfSegmenter = PDFSegmenterType::New();

  timeCollector.Start("LoadData");

  typename ImageReaderType::Pointer reader;
  int j = N;
  if( useTexture )
    {
    --j;
    }
  for(int i=0; i<j; i++)
    {
    reader = ImageReaderType::New();
    if(i == 0)
      {
      reader->SetFileName( inputVolume1.c_str() );
      reader->Update();
      pdfSegmenter->SetInputVolume1( reader->GetOutput() );
      }
    else if(i == 1)
      {
      reader->SetFileName( inputVolume2.c_str() );
      reader->Update();
      pdfSegmenter->SetInputVolume2( reader->GetOutput() );
      }
    else if(i == 2)
      {
      reader->SetFileName( inputVolume3.c_str() );
      reader->Update();
      pdfSegmenter->SetInputVolume3( reader->GetOutput() );
      }
    else
      {
      std::cout << "ERROR: current command line xml file limits"
                << " this filter to 3 input images" << std::endl;
      return 1;
      }
    }

  MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
  inMaskReader->SetFileName( labelmap.c_str() );
  inMaskReader->Update();
  pdfSegmenter->SetLabelmap( inMaskReader->GetOutput() );

  timeCollector.Stop("LoadData");

  pdfSegmenter->SetObjectId( objectId[0] );
  if( objectId.size() > 1 )
    {
    for( unsigned int o=1; o<objectId.size(); o++ )
      {
      pdfSegmenter->AddObjectId( objectId[o] );
      }
    }
  pdfSegmenter->SetVoidId( voidId );
  pdfSegmenter->SetUseTexture( useTexture );
  pdfSegmenter->SetErodeRadius( erodeRadius );
  pdfSegmenter->SetHoleFillIterations( holeFillIterations );
  pdfSegmenter->SetFprWeight( fprWeight );
  pdfSegmenter->SetProbabilitySmoothingStandardDeviation(
    probSmoothingStdDev );
  pdfSegmenter->SetDraft( draft );
  pdfSegmenter->SetReclassifyNotObjectMask( reclassifyNotObjectMask );
  pdfSegmenter->SetReclassifyObjectMask( reclassifyObjectMask );

  pdfSegmenter->Update();

  timeCollector.Start("Save");

  if( pdfSegmenter->GetProbabilityImage(0) != NULL
    && probabilityVolume0.size() > 2 )
    {
    ProbImageWriterType::Pointer probImageWriter =
      ProbImageWriterType::New();
    probImageWriter->SetFileName( probabilityVolume0.c_str() );
    probImageWriter->SetInput( *(pdfSegmenter->GetProbabilityImage(0)) );
    probImageWriter->Update();
    }
  if( objectId.size() > 1 && pdfSegmenter->GetProbabilityImage(1) != NULL
    && probabilityVolume1.size() > 2 )
    {
    ProbImageWriterType::Pointer probImageWriter =
      ProbImageWriterType::New();
    probImageWriter->SetFileName( probabilityVolume1.c_str() );
    probImageWriter->SetInput( *(pdfSegmenter->GetProbabilityImage(1)) );
    probImageWriter->Update();
    }
  if( objectId.size() > 2 && pdfSegmenter->GetProbabilityImage(2) != NULL
    && probabilityVolume2.size() > 2 )
    {
    ProbImageWriterType::Pointer probImageWriter =
      ProbImageWriterType::New();
    probImageWriter->SetFileName( probabilityVolume2.c_str() );
    probImageWriter->SetInput( *(pdfSegmenter->GetProbabilityImage(2)) );
    probImageWriter->Update();
    }

  MaskWriterType::Pointer writer = MaskWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( pdfSegmenter->GetLabelmap() );
  writer->Update();

  timeCollector.Stop("Save");

  timeCollector.Report();

  return 0;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    GetImageType (inputVolume1, pixelType, componentType);

    int N = 1;
    if(inputVolume2.length() > 1)
      {
      ++N;
      if(inputVolume3.length() > 1)
        {
        ++N;
        }
      }
    if(useTexture)
      {
      ++N;
      }

    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        if(N == 1)
          {
          return DoIt<unsigned char, 1>( argc, argv );
          }
        else if(N == 2)
          {
          return DoIt<unsigned char, 2>( argc, argv );
          }
        else if(N == 3)
          {
          return DoIt<unsigned char, 3>( argc, argv );
          }
        else
          {
          return DoIt<unsigned char, 4>( argc, argv );
          }
        break;
      case itk::ImageIOBase::USHORT:
        if(N == 1)
          {
          return DoIt<unsigned short, 1>( argc, argv );
          }
        else if(N == 2)
          {
          return DoIt<unsigned short, 2>( argc, argv );
          }
        else if(N == 3)
          {
          return DoIt<unsigned short, 3>( argc, argv );
          }
        else
          {
          return DoIt<unsigned short, 4>( argc, argv );
          }
        break;
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::SHORT:
        if(N == 1)
          {
          return DoIt<short, 1>( argc, argv );
          }
        else if(N == 2)
          {
          return DoIt<short, 2>( argc, argv );
          }
        else if(N == 3)
          {
          return DoIt<short, 3>( argc, argv );
          }
        else
          {
          return DoIt<short, 4>( argc, argv );
          }
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::INT:
      case itk::ImageIOBase::ULONG:
      case itk::ImageIOBase::LONG:
      case itk::ImageIOBase::FLOAT:
      case itk::ImageIOBase::DOUBLE:
        if(N == 1)
          {
          return DoIt<float, 1>( argc, argv );
          }
        else if(N == 2)
          {
          return DoIt<float, 2>( argc, argv );
          }
        else if(N == 3)
          {
          return DoIt<float, 3>( argc, argv );
          }
        else
          {
          return DoIt<float, 4>( argc, argv );
          }
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
