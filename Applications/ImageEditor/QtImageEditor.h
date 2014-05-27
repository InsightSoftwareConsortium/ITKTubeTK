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

#ifndef __QtImageEditor_h
#define __QtImageEditor_h

// Qt includes
#include <QDialog>
#include <QLineEdit>

// ITK includes
#include <itkImage.h>

// QtImageViewer includes
#include "ui_QtSlicerGUI.h"

// TubeTK includes
class QtOverlayControlsWidget;

namespace tube
{

class QtImageEditor
  : public QDialog, public Ui::GuiDialogBase
{
  Q_OBJECT
public:

  QtImageEditor( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
  ~QtImageEditor();

  typedef double                      PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;


public slots:
  void hideHelp();
  void showHelp(bool checked);
  void setMaximumSlice();
  void setDisplaySigma(QString value);
  void setInputImage(ImageType *newImData);
  void setDisplaySliceNumber(int number);
  bool loadImage(QString filePathToLoad = QString());
  void loadOverlay(QString overlayImagePath = QString());
  void applyFFT();
  void applyInverseFFT();
  void applyFilter();
  void displayFFT();
  void blurFilter();

private:
  class Internals;
  Internals* m_Internals;
};

} // End namespace tube

#endif
