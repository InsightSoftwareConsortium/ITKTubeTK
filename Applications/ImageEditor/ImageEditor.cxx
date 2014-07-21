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

// Qt includes
#include <QApplication>
#include <QDebug>
#include <QFileDialog>

// QtImageViewer includes
#include "QtGlSliceView.h"
#include "QtImageViewer.h"

// QtImageEditor includes
#include "QtImageEditor.h"

int execImageEditor(int argc, char* argv[])
{
  QApplication myApp( argc, argv );

  tube::QtImageEditor imageEditor(0,0);
  //QtImageViewer imageEditor(0,0);
  imageEditor.setWindowTitle("ImageEditor");
  imageEditor.loadInputImage();
  imageEditor.show();
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

  tube::QtImageEditor imageEditor(0);
  //QtImageViewer imageEditor(0);

  imageEditor.setWindowTitle("ImageEditor");

  imageEditor.loadInputImage(QString::fromStdString(inputImage));
  if(!overlayImage.empty())
    {
    imageEditor.loadOverlayImage(QString::fromStdString(overlayImage));
    }

  imageEditor.sliceView()->setOrientation(orientation);
  if(sliceOffset != -1)
    {
    imageEditor.sliceView()->setSliceNum(sliceOffset);
    }
  if(maxIntensityArg.isSet())
    {
    imageEditor.sliceView()->setIWMax(maxIntensity);
    }
  if(minIntensityArg.isSet())
    {
    imageEditor.sliceView()->setIWMin(minIntensity);
    }
  imageEditor.sliceView()->setZoom(zoom);
  imageEditor.sliceView()->transpose(transpose);
  imageEditor.sliceView()->flipZ(zFlipped);
  imageEditor.sliceView()->flipY(yFlipped);
  imageEditor.sliceView()->flipX(xFlipped);
  imageEditor.sliceView()->setOverlayOpacity(overlayOpacity);
  imageEditor.sliceView()->setViewCrosshairs(crosshairs);
  imageEditor.sliceView()->setDisplayState(details);
  imageEditor.sliceView()->setViewValuePhysicalUnits(physicalUnits);
  imageEditor.sliceView()->setViewValue(value);
  imageEditor.sliceView()->setViewAxisLabel(axisLabel);
  imageEditor.sliceView()->setViewClickedPoints(clickedPoints);
  imageEditor.sliceView()->setImageMode(imageMode.c_str());
  imageEditor.sliceView()->setIWModeMax(iwModeMax.c_str());
  imageEditor.sliceView()->setIWModeMin(iwModeMin.c_str());
  if(sigmaArg.isSet())
    {
    imageEditor.setDisplaySigma(sigma);
    imageEditor.applyFilter();
    }
  imageEditor.show();
  imageEditor.sliceView()->update();
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
