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
#include "ImageEditorConfigure.h"
#include "ImageEditorCLP.h"

//Qt includes
#include <QApplication>
#include <QDebug>
#include <QFileDialog>

//QtImageViewer includes
#include "QtGlSliceView.h"
#include "QtSlicer.h"

//QtImageEditor includes
#include "QtImageEditor.h"


using namespace tube;

int execImageEditor(int argc, char* argv[])
{
  QApplication myApp( argc, argv );

  QtImageEditor qtSlicerWindow(0,0);

  qtSlicerWindow.setWindowTitle("ImageEditor");
  qtSlicerWindow.loadImage();
  qtSlicerWindow.show();

  int execReturn;
  try
    {
    execReturn = myApp.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  return execReturn;
}

int parseAndExecImageEditor(int argc, char* argv[])
{
  PARSE_ARGS;
  QApplication myApp( argc, argv );

  QtImageEditor qtSlicerWindow(0);

  qtSlicerWindow.setWindowTitle("ImageEditor");

  qtSlicerWindow.loadImage(QString::fromStdString(inputImage));
  if(!overlayImage.empty())
    {
    qtSlicerWindow.loadOverlay(QString::fromStdString(overlayImage));
    }

  qtSlicerWindow.OpenGlWindow->setOrientation(orientation);
  if(sliceOffset != -1)
    {
    qtSlicerWindow.OpenGlWindow->setSliceNum(sliceOffset);
    }
  if(maxIntensityArg.isSet())
    {
    qtSlicerWindow.OpenGlWindow->setMaxIntensity(maxIntensity);
    }
  if(minIntensityArg.isSet())
    {
    qtSlicerWindow.OpenGlWindow->setMinIntensity(minIntensity);
    }
  qtSlicerWindow.OpenGlWindow->setZoom(zoom);
  qtSlicerWindow.OpenGlWindow->transpose(transpose);
  qtSlicerWindow.OpenGlWindow->flipZ(zFlipped);
  qtSlicerWindow.OpenGlWindow->flipY(yFlipped);
  qtSlicerWindow.OpenGlWindow->flipX(xFlipped);
  qtSlicerWindow.OpenGlWindow->setOverlayOpacity(overlayOpacity);
  qtSlicerWindow.OpenGlWindow->setViewCrosshairs(crosshairs);
  qtSlicerWindow.OpenGlWindow->setViewDetails(details);
  qtSlicerWindow.OpenGlWindow->setViewValuePhysicalUnits(physicalUnits);
  qtSlicerWindow.OpenGlWindow->setViewValue(value);
  qtSlicerWindow.OpenGlWindow->setViewAxisLabel(axisLabel);
  qtSlicerWindow.OpenGlWindow->setViewClickedPoints(clickedPoints);
  qtSlicerWindow.OpenGlWindow->setImageMode(imageMode.c_str());
  qtSlicerWindow.OpenGlWindow->setIWModeMax(iwModeMax.c_str());
  qtSlicerWindow.OpenGlWindow->setIWModeMin(iwModeMin.c_str());
  if(sigmaArg.isSet())
    {
    qtSlicerWindow.setDisplaySigma(QString::number(sigma, 'f', 2));
    qtSlicerWindow.applyFilter();
    }
  qtSlicerWindow.OpenGlWindow->update();

  qtSlicerWindow.show();
  int execReturn;
  try
    {
    execReturn = myApp.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  return execReturn;
}

int main( int argc, char* argv[] )
{
#if !defined(BUILD_SHARED_LIBS)
  Q_INIT_RESOURCE(qtImageViewerResources);
#endif
  int res;
  if(argc == 1)
    {
    res = execImageEditor(argc, argv);
    }
  else
    {
    res = parseAndExecImageEditor(argc, argv);
    }
  return res;
}
