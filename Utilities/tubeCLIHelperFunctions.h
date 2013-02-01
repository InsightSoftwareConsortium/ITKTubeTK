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

#ifndef __tubeCLIHelperFunctions_h
#define __tubeCLIHelperFunctions_h

#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"

namespace tube {

// Description:
// Get the ComponentType and dimension of the image
void GetImageInformation( const std::string & fileName,
                          itk::ImageIOBase::IOComponentType &componentType,
                          unsigned int & dimension )
  {
  // Find out the component type of the image in file
  typedef itk::ImageIOBase::IOComponentType  PixelType;

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

#ifndef PARSE_ARGS_FLOAT_ONLY

int ParseArgsAndCallDoIt( std::string inputImage,
                          int argc,
                          char **argv )
  {
  itk::ImageIOBase::IOComponentType componentType;
  unsigned int dimension;

  try
    {
    GetImageInformation(inputImage, componentType, dimension );
    if( dimension == 2 )
      {
      switch( componentType )
        {
        case itk::ImageIOBase::UCHAR:
          return DoIt<unsigned char, 2>( argc, argv );
        case itk::ImageIOBase::CHAR:
          return DoIt<char, 2>( argc, argv );
        case itk::ImageIOBase::USHORT:
          return DoIt<unsigned short, 2>( argc, argv );
        case itk::ImageIOBase::SHORT:
          return DoIt<short, 2>( argc, argv );
        case itk::ImageIOBase::UINT:
          return DoIt<unsigned int, 2>( argc, argv );
        case itk::ImageIOBase::INT:
          return DoIt<int, 2>( argc, argv );
        case itk::ImageIOBase::ULONG:
          return DoIt<unsigned long, 2>( argc, argv );
        case itk::ImageIOBase::LONG:
          return DoIt<long, 2>( argc, argv );
        case itk::ImageIOBase::FLOAT:
          return DoIt<float, 2>( argc, argv );
        case itk::ImageIOBase::DOUBLE:
          return DoIt<double, 2>( argc, argv );
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown component type" << std::endl;
          return EXIT_FAILURE;
        }
      }
    else if( dimension == 3 )
      {
      switch( componentType )
        {
        case itk::ImageIOBase::UCHAR:
          return DoIt<unsigned char, 3>( argc, argv );
        case itk::ImageIOBase::CHAR:
          return DoIt<char, 3>( argc, argv );
        case itk::ImageIOBase::USHORT:
          return DoIt<unsigned short, 3>( argc, argv );
        case itk::ImageIOBase::SHORT:
          return DoIt<short, 3>( argc, argv );
        case itk::ImageIOBase::UINT:
          return DoIt<unsigned int, 3>( argc, argv );
        case itk::ImageIOBase::INT:
          return DoIt<int, 3>( argc, argv );
        case itk::ImageIOBase::ULONG:
          return DoIt<unsigned long, 3>( argc, argv );
        case itk::ImageIOBase::LONG:
          return DoIt<long, 3>( argc, argv );
        case itk::ImageIOBase::FLOAT:
          return DoIt<float, 3>( argc, argv );
        case itk::ImageIOBase::DOUBLE:
          return DoIt<double, 3>( argc, argv );
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
          std::cout << "unknown component type" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << argv[0] << ": ITK exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }

#else // PARSE_ARGS_FLOAT_ONLY

int ParseArgsAndCallDoIt( std::string inputImage,
                          int argc,
                          char **argv )
  {
  itk::ImageIOBase::IOComponentType componentType;
  unsigned int dimension;

  try
    {
    GetImageInformation(inputImage, componentType, dimension );
    if( dimension == 2 )
      {
      return DoIt<float, 2>( argc, argv );
      }
    else if( dimension == 3 )
      {
      return DoIt<float, 3>( argc, argv );
      }
    }
  catch( itk::ExceptionObject &excep )
    {
    std::cerr << argv[0] << ": itk exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }

#endif

}; // namespace tube

#endif

