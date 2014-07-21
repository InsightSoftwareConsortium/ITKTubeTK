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

// Qt includes
#include <QDebug>
#include <QDir>
#include <QFileDialog>

// ImageEditor includes
#include "QtOverlayControlsWidget.h"
#include "ui_QtOverlayControlsWidget.h"

namespace tube
{


class QtOverlayControlsWidgetPrivate: public Ui_QtOverlayControlsWidget
{
  Q_DECLARE_PUBLIC(QtOverlayControlsWidget);
public:
  typedef Ui_QtOverlayControlsWidget Superclass;
  QtOverlayControlsWidgetPrivate(QtOverlayControlsWidget& obj);

  virtual void setupUi(QWidget* widgetToSetup);

  QtGlSliceView  *m_SliceView;

protected:
  QtOverlayControlsWidget* const q_ptr;
};


QtOverlayControlsWidgetPrivate::QtOverlayControlsWidgetPrivate(QtOverlayControlsWidget& obj)
  : q_ptr(&obj)
  , m_SliceView(0)
{
}

void QtOverlayControlsWidgetPrivate::setupUi(QWidget* widgetToSetup)
{
  Q_Q(QtOverlayControlsWidget);
  this->Superclass::setupUi(widgetToSetup);

  QObject::connect(this->LoadOverlayButton, SIGNAL(clicked()),
                   q, SLOT(loadOverlay()));
  QObject::connect(this->OverlayOpacitySlider, SIGNAL(sliderMoved(int)),
                   q, SLOT(onOpacityChanged(int)));

}

QtOverlayControlsWidget::QtOverlayControlsWidget(QWidget* parent)
  : QWidget(parent)
  , d_ptr(new QtOverlayControlsWidgetPrivate(*this))
{
  Q_D(QtOverlayControlsWidget);
  d->setupUi(this);
}


QtOverlayControlsWidget::~QtOverlayControlsWidget()
{
}


void QtOverlayControlsWidget::onOpacityChanged(int newOpacity)
{
  Q_D(QtOverlayControlsWidget);
  double opacity = (static_cast<double>(newOpacity)/100.);
  if (d->m_SliceView)
    {
    d->m_SliceView->setOverlayOpacity(opacity);
    }
  emit opacityChanged(opacity);
}


double QtOverlayControlsWidget::opacity() const
{
  Q_D(const QtOverlayControlsWidget);
  return (d->OverlayOpacitySpinBox->value() / 100 );
}


void QtOverlayControlsWidget::setOpacity(double value)
{
  Q_D(QtOverlayControlsWidget);
  int valueSlider = value*100;
  d->OverlayOpacitySlider->setValue(valueSlider);
  d->OverlayOpacitySpinBox->setValue(value);
}


void QtOverlayControlsWidget::setSliceView(QtGlSliceView* sliceView)
{
  Q_D(QtOverlayControlsWidget);
  d->m_SliceView = sliceView;
  QObject::connect(d->m_SliceView, SIGNAL(overlayOpacityChanged(double)),
                   this, SLOT(setOpacity(double)));
}


bool QtOverlayControlsWidget::loadOverlay(QString pathOverlay)
{
  Q_D(QtOverlayControlsWidget);
  OverlayReaderType::Pointer overlayReader = OverlayReaderType::New();

  if(pathOverlay.isEmpty())
    {
    pathOverlay = QFileDialog::getOpenFileName(
      0,"", QDir::currentPath());
    }
  if(pathOverlay.isEmpty())
    {
    return false;
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
  this->setInputOverlay( overlayReader->GetOutput() );
  return true;
}

void QtOverlayControlsWidget::setInputOverlay(OverlayType* overlayImage)
{
  Q_D(QtOverlayControlsWidget);
  d->m_SliceView->setInputOverlay(overlayImage);

  d->OverlayOpacitySpinBox->setEnabled(d->m_SliceView->validOverlayData());
  d->OverlayOpacitySlider->setEnabled(d->m_SliceView->validOverlayData());

  this->setOpacity(d->m_SliceView->overlayOpacity());
}

}
