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

#include "tubeCLIFilterWatcher.h"

#include <itkImageFileReader.h>
#include <itkRecursiveGaussianImageFilter.h>

void testCallbackFunction( void * data )
{
  if( data )
    {
    std::cout << "testCallbackFunction with pointer: " << data << std::endl;
    }
  else
    {
    std::cout << "testCallbackFunction with pointer: NULL" << std::endl;
    }
}

int tubeCLIFilterWatcherTest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "inputImage" << std::endl;
    return EXIT_FAILURE;
    }


  typedef itk::Image< float, 2 >              ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
  FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );

  ModuleProcessInformation CLPProcessInfo;
  CLPProcessInfo.Initialize();
  CLPProcessInfo.SetProgressCallback( testCallbackFunction, NULL );
  tube::CLIFilterWatcher watcher( filter, "RecursiveGaussianFilter",
                                  &CLPProcessInfo );

  filter->Update();

  return EXIT_SUCCESS;
}
