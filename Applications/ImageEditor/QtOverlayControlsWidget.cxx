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

#include "QtOverlayControlsWidget.h"

//Qt includes
#include <QDebug>
#include <QDir>
#include <QFileDialog>

namespace tube{


QtOverlayControlsWidget::QtOverlayControlsWidget(QWidget* parent)
  : QWidget(parent)
{
  this->m_UI = new Ui::Overlay;
  m_UI->setupUi(this);
}

QtOverlayControlsWidget::~QtOverlayControlsWidget()
{
}

void QtOverlayControlsWidget::onOpacityChanged(int newOpacity)
{
  double opacity = (static_cast<double>(newOpacity)/100.);
  emit opacityChanged(opacity);
}


double QtOverlayControlsWidget::getOpacity() const
{
  return (this->m_UI->OverlayOpacity->value() / 100 );
}


void QtOverlayControlsWidget::setOpacity(double value)
{
  int valueSlider = value*100;
  this->m_UI->OverlayOpacity->setValue(valueSlider);
  this->m_UI->OverlayOpacityLineEdit->setText(QString::number(value, 'f', 3));
}


void QtOverlayControlsWidget::setSliceView(QtGlSliceView* sliceView)
{
  this->m_SliceView = sliceView;
  QObject::connect(this->m_UI->LoadOverlayButton, SIGNAL(clicked()),
                   this, SLOT(loadOverlay()));
  QObject::connect(this->m_UI->OverlayOpacity, SIGNAL(sliderMoved(int)), this,
                   SLOT(onOpacityChanged(int)));
  QObject::connect(this, SIGNAL(opacityChanged(double)), m_SliceView,
                   SLOT(setOverlayOpacity(double)));
  QObject::connect(m_SliceView, SIGNAL(overlayOpacityChanged(double)), this,
                   SLOT(setOpacity(double)));
  QObject::connect(m_SliceView, SIGNAL(validOverlayDataChanged(bool)),
                   this->m_UI->OverlayCheckBox, SLOT(setChecked(bool)));
  QObject::connect(this->m_UI->OverlayCheckBox, SIGNAL(toggled(bool)), this,
                   SLOT(setOverlayVisibility(bool)));
}


int QtOverlayControlsWidget::loadOverlay()
{
  OverlayReaderType::Pointer overlayReader = OverlayReaderType::New();
  QString pathOverlay = QFileDialog::getOpenFileName(
        0,"", QDir::currentPath());

  if(pathOverlay.isEmpty())
    {
    return 0;
    }
  overlayReader->SetFileName( pathOverlay.toLatin1().data() );

  qDebug() << "loading image " << pathOverlay << " ... ";
  try
    {
    overlayReader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Done!" << std::endl;
  setInputOverlay( overlayReader->GetOutput() );

  show();
}

int QtOverlayControlsWidget::loadOverlay(std::string path)
{
  OverlayReaderType::Pointer overlayReader = OverlayReaderType::New();
  QString pathOverlay = QString::fromStdString(path);

  if(pathOverlay.isNull())
    {
    pathOverlay = QFileDialog::getOpenFileName(
            0,"", QDir::currentPath());
    return 0;
    }

  if(pathOverlay.isEmpty())
    {
    return 0;
    }
  overlayReader->SetFileName( pathOverlay.toLatin1().data() );

  qDebug() << "loading image " << pathOverlay << " ... ";
  try
    {
    overlayReader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Done!" << std::endl;
  setInputOverlay( overlayReader->GetOutput() );

  show();
}

void QtOverlayControlsWidget::setInputOverlay(OverlayType* overlayImage)
{
  this->m_SliceView->setInputOverlay(overlayImage);
  this->m_UI->OverlayOpacity->setMaximum(this->m_SliceView->overlayOpacity()*100);
  this->setOpacity(this->m_SliceView->overlayOpacity());

  this->m_SliceView->setFocus();
  this->m_SliceView->update();
  this->m_UI->OverlayOpacityLineEdit->setText(QString::number
                              (this->m_SliceView->overlayOpacity(),'f',3));
}

void QtOverlayControlsWidget::setOverlayVisibility(bool show)
{
  if(show)
    {
    if(!(this->m_SliceView->validOverlayData()))
      {
      loadOverlay();
      }
    else
      {
      this->m_SliceView->setValidOverlayData(show);
      }
    }
  else
    {
    this->m_SliceView->setValidOverlayData(show);
    }
  this->m_SliceView->update();
}

}
