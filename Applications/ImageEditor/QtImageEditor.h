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

//Qt includes
#include <QDialog>
#include <QDialogButtonBox>
#include <QLineEdit>
#include <QString>

//itk includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"

//QtImageViewer includes
#include "ui_QtSlicerGUI.h"
#include "ui_QtSlicerHelpGUI.h"

//TubeTK includes
#include "QtOverlayControlsWidget.h"
#include "ui_QtOverlayControlsWidgetGUI.h"

namespace tube
{

class QtImageEditor : public QDialog, public Ui::GuiDialogBase
{
  Q_OBJECT
public:

  QtImageEditor( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
  ~QtImageEditor();

  typedef itk::Image< double, 3 > ImageType;

  typedef double InputPixelType;
  typedef double OutputPixelType;

  typedef itk::Image< InputPixelType, 3 >  InputImageType;
  typedef itk::Image< OutputPixelType, 3 > OutputImageType;
  typedef itk::RecursiveGaussianImageFilter<
               InputImageType, OutputImageType > FilterType;

  typedef itk::ImageFileReader<ImageType>   ReaderType;

public:

public slots:
  //void toggleTextEdit(int viewDetail);
  void setMaximumSlice();
  void setDisplaySigma(QString value);
  void setInputImage(ImageType *newImData);
  void setDisplaySliceNumber(int number);
  int loadImage(std::string path);
  int loadImage();
  int loadOverlay(std::string path);
  void applyFilter();

private:
  QLineEdit               *m_SigmaLineEdit;
  ImageType               *m_ImageData;
  QtOverlayControlsWidget *m_OverlayWidget;
};

} // End namespace tube

#endif
