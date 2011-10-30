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
#include "tubeImageViewer.h"
#include <qfiledialog.h>
#include <qslider.h>

#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageRegionIteratorWithIndex.h"

int Usage( int argc, char *argv[] )
{
  std::cout << argv[0] << " [ <inputImage> [ <inputOverlay> ] ]"
    << std::endl;
  return EXIT_FAILURE;
}

int main( int argc, char* argv[] )
{
  QApplication myApp( argc, argv );

  tubeImageViewer myGUI( 0, 0, TRUE );

  typedef float                             ImagePixelType;
  typedef itk::Image< ImagePixelType, 3 >   ImageType;

  typedef unsigned char                     OverlayPixelType;
  typedef itk::Image< OverlayPixelType, 3 > OverlayType;

  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileReader< OverlayType > OverlayReaderType;


  QString caption = "Open";
  QString directory = ".";
  QString filter = "Images (*.*)";

  QString imageFilename;
  QString overlayFilename;
  if( argc == 1 )
    {
    imageFilename = QFileDialog::getOpenFileName( 0, caption, directory,
      filter );
    }
  else
    {
    if( argc == 2 )
      {
      imageFilename = argv[1];
      }
    else
      {
      if( argc == 3 )
        {
        imageFilename = argv[1];
        overlayFilename = argv[2];
        }
      else
        {
        return Usage( argc, argv );
        }
      }
    }

  if( imageFilename.isNull() )
    {
    std::cout << "Filename not specified." << std::endl;
    return EXIT_FAILURE;
    }

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( imageFilename.toLatin1() );

  try
    {
    std::cout << "Loading image..." << std::endl;
    imageReader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file imageReader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  myGUI.SetInputImage( imageReader->GetOutput() );

  if( !overlayFilename.isNull() )
    {
    OverlayReaderType::Pointer overlayReader = OverlayReaderType::New();
    overlayReader->SetFileName( overlayFilename.toLatin1() );

    try
      {
      std::cout << "Loading overlay..." << std::endl;
      overlayReader->Update();
      }
    catch (itk::ExceptionObject & e)
      {
      std::cerr << "Exception in file overlayReader " << std::endl;
      std::cerr << e << std::endl;
      return EXIT_FAILURE;
      }
    myGUI.SetInputOverlay( overlayReader->GetOutput() );
    }

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
