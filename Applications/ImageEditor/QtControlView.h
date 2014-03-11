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

#ifndef __QtControlView_h
#define __QtControlView_h

//Qt includes
#include <QDialog>

//itk includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"

//QtImageViewer includes
#include "ui_QtSlicerGUI.h"
#include "ui_QtSlicerHelpGUI.h"

namespace tube
{

class QtControlView : public QDialog, public Ui::GuiDialogBase
{
  Q_OBJECT
public:

  QtControlView( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
  ~QtControlView();

  typedef itk::Image< double, 3 > ImageType;

  typedef double InputPixelType;
  typedef double OutputPixelType;

  typedef itk::Image< InputPixelType, 3 >  InputImageType;
  typedef itk::Image< OutputPixelType, 3 > OutputImageType;
  typedef itk::RecursiveGaussianImageFilter<
               InputImageType, OutputImageType > FilterType;

public:

public slots:
  void setDisplaySigma(QString value);
  void setDisplayPosition(int x, int y , int z, double value);
  void setInputImage(ImageType *newImData);
  void setDisplaySliceNumber(int number);
  void setDisplayIMin(int value);
  void setDisplayIMax(int value);

private slots:
  void setFilter();

private:
  void setTab();

  QTabWidget   *m_TabWidget;
  QWidget      *m_FilterControlWidget;
  QGridLayout  *m_FilterGridLayout;
  QLineEdit    *m_SigmaLineEdit;
  QPushButton  *m_ApplyButton;
  ImageType    *m_ImageData;
};

} // End namespace tube

#endif
