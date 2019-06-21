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

#ifndef __tubeCLIHelperFunctions_h
#define __tubeCLIHelperFunctions_h

#include "tubeMacro.h"

#include <itkImage.h>
#include <itkImageIOFactory.h>

namespace tube
{

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

int ParseArgsAndCallDoIt( const std::string & inputImage, int argc,
                          char * argv[] )
{
  typedef itk::ImageIOBase              ImageIOType;
  typedef ImageIOType::IOComponentType  IOComponentType;

  IOComponentType componentType = ImageIOType::UNKNOWNCOMPONENTTYPE;
  unsigned int dimension = 0;
  try
    {
    GetImageInformation( inputImage, componentType, dimension );

#ifndef PARSE_ARGS_FLOAT_ONLY

#ifndef PARSE_ARGS_3D_ONLY
    if( dimension == 2 )
      {
      switch( componentType )
        {
        case ImageIOType::UCHAR:
          return DoIt< unsigned char, 2 >( argc, argv );
        case ImageIOType::USHORT:
          return DoIt< unsigned short, 2 >( argc, argv );
        case ImageIOType::SHORT:
          return DoIt< short, 2 >( argc, argv );
#ifndef PARSE_ARGS_INT_ONLY
        case ImageIOType::FLOAT:
          return DoIt< float, 2 >( argc, argv );
#endif
        case ImageIOType::INT:
          return DoIt< int, 2 >( argc, argv );
        case ImageIOType::UNKNOWNCOMPONENTTYPE:
        default:
          tubeErrorMacro( << "Unknown component type." );
          return EXIT_FAILURE;
        }
      }
#endif

    if( dimension == 3 )
      {
      switch( componentType )
        {
        case ImageIOType::UCHAR:
          return DoIt < unsigned char, 3 >( argc, argv );
        case ImageIOType::USHORT:
          return DoIt < unsigned short, 3 >( argc, argv );
        case ImageIOType::SHORT:
          return DoIt< short, 3 >( argc, argv );
#ifndef PARSE_ARGS_INT_ONLY
        case ImageIOType::FLOAT:
          return DoIt < float, 3 >( argc, argv );
#endif
        case ImageIOType::INT:
          return DoIt< int, 3 >( argc, argv );
        case ImageIOType::UNKNOWNCOMPONENTTYPE:
        default:
          tubeErrorMacro( << "Unknown component type." );
          return EXIT_FAILURE;
        }
      }

#else
#ifndef PARSE_ARGS_3D_ONLY
    if( dimension == 2 )
      {
      return DoIt< float, 2 >( argc, argv );
      }
#endif
    if( dimension == 3 )
      {
      return DoIt< float, 3 >( argc, argv );
      }

#endif // End !defined( PARSE_ARGS_FLOAT_ONLY )

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

  return EXIT_SUCCESS;
}

} // End namespace tube

#endif // End !defined( __tubeCLIHelperFunctions_h )
