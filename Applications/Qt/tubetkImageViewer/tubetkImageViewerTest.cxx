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
#include <qapplication.h>
#include "QtGlSliceView.h"
#include "tubetkImageViewer.h"
#include <qfiledialog.h>
#include <qslider.h>

#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char* argv[] )
{
  QApplication myApp( argc, argv );

  tubetkImageViewer myGUI( 0, 0, TRUE );

  typedef float                             PixelType;
  typedef itk::Image<PixelType, 3>          ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

  QString caption = "Open";
  QString directory = ".";
  QString filter = "Images (*.*)";

  QString filename = QFileDialog::getOpenFileName( 0, caption, directory, filter );

  if (filename.isNull())
    return 1;

  reader->SetFileName( filename.toLatin1() );

  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Done!" << std::endl;
  myGUI.SetInputImage( reader->GetOutput() );

  try
    {
    myGUI.show();
    myApp.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  return 0;

}
