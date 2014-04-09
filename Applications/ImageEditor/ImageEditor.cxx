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
//Qt includes
#include <QApplication>
#include <QDebug>
#include <QFileDialog>
#include <QPlastiqueStyle>

//QtImageViewer includes
#include "QtGlSliceView.h"
#include "QtSlicer.h"

//QtImageEditor includes
#include "QtImageEditor.h"

// itk includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMetaImageIOFactory.h"

#include "ImageEditorCLP.h"

using namespace tube;

int main( int argc, char* argv[] )
{
  PARSE_ARGS;

  QApplication myApp( argc, argv );

  QtImageEditor qtSlicerWindow(0,0);

  qtSlicerWindow.setWindowTitle("Insight Qt Slicer" );
  myApp.setStyle(new QPlastiqueStyle );
  QPalette p( QColor( 239, 239, 239 ) );
  myApp.setPalette( p );
  qtSlicerWindow.loadImage(inputImage);
  qtSlicerWindow.loadOverlay(overlayImage);

  qtSlicerWindow.OpenGlWindow->setOrientation(orientation);
  qtSlicerWindow.OpenGlWindow->setSliceNum(sliceOffset);
  qtSlicerWindow.OpenGlWindow->setMaxIntensity(maxIntensity);
  qtSlicerWindow.OpenGlWindow->setMinIntensity(minIntensity);
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
  double sig = sigma;
  qtSlicerWindow.setDisplaySigma(QString::number(sigma, 'f', 2));
  if(sig != 0) qtSlicerWindow.applyFilter();

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
