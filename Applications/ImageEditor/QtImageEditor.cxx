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
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QTabWidget>

// QtImageViewer includes
#include "QtGlSliceView.h"

// ImageEditor includes
#include "QtImageEditor.h"
#include "QtOverlayControlsWidget.h"

// TubeTK includes
#include "itktubeGaussianDerivativeImageSource.h"
#include "itktubeFFTGaussianDerivativeIFFTFilter.h"

// ITK includes
#include <itkComplexToImaginaryImageFilter.h>
#include <itkComplexToModulusImageFilter.h>
#include <itkComplexToRealImageFilter.h>
#include <itkFFTShiftImageFilter.h>
#include <itkForwardFFTImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkInverseFFTImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkTimeProbe.h>

namespace tube
{

/**
 * Helper function for forceUpdate(). Not self-sufficient!
 */
void invalidateLayout(QLayout *layout) {
    // Recompute the given layout and all its child layouts.
    for (int i = 0; i < layout->count(); i++) {
        QLayoutItem *item = layout->itemAt(i);
        if (item->layout()) {
            invalidateLayout(item->layout());
        } else {
            item->invalidate();
        }
    }
    layout->invalidate();
    layout->update();
    layout->activate();
}

void forceUpdate(QWidget *widget) {
    // Update all child widgets.
    for (int i = 0; i < widget->children().size(); i++) {
        QObject *child = widget->children()[i];
        if (child->isWidgetType()) {
            forceUpdate((QWidget *)child);
        }
    }

    // Invalidate the layout of the widget.
    if (widget->layout()) {
        invalidateLayout(widget->layout());
    }
}


class QtImageEditorPrivate
{
public:
  typedef QtImageEditor::ImageType                     ImageType;
  typedef itk::Image<unsigned char, 3>                 UnsignedCharImageType;

  typedef itk::ImageFileReader<ImageType>              ReaderType;
  typedef itk::ImageFileWriter<UnsignedCharImageType>  WriterType;

  typedef itk::ForwardFFTImageFilter<ImageType>        FFTType;
  typedef FFTType::OutputImageType                     ComplexImageType;
  typedef itk::InverseFFTImageFilter<ComplexImageType, ImageType>
                                                       InverseFFTType;
  typedef itk::tube::GaussianDerivativeImageSource<ImageType>
                                                       GaussianDerivativeSourceType;
  typedef itk::tube::FFTGaussianDerivativeIFFTFilter<ImageType,ImageType>
                                                       FFTGaussianDerivativeIFFTType;
  typedef itk::FFTShiftImageFilter<ImageType, ImageType>
                                                       FFTShiftFilterType;
  typedef itk::MultiplyImageFilter<ComplexImageType, ImageType, ComplexImageType>
                                                       MultiplyFilterType;

  QtImageEditorPrivate(QtImageEditor& obj);
  void setupUi(QDialog* widgetToSetup);

protected:
  QtImageEditor              *q_ptr;

  ImageType                  *m_CurrentImageData;
  ImageType                  *m_ImageData;
  QtOverlayControlsWidget    *m_OverlayWidget;
  FFTType::Pointer            m_FFTFilter;
  QDoubleSpinBox             *m_SigmaSpinBox;
  MultiplyFilterType::Pointer m_MultiplyFilter;
  FFTShiftFilterType::Pointer m_FFTShiftFilter;
  InverseFFTType::Pointer     m_InverseFFTFilter;
  QTabWidget                 *m_TabWidget;
  QSpinBox                   *m_XOrderSpinBox;
  QSpinBox                   *m_YOrderSpinBox;
  QSpinBox                   *m_ZOrderSpinBox;

  GaussianDerivativeSourceType::Pointer createGaussianDerivative(
    GaussianDerivativeSourceType::OrdersType order);/*,
    ImageType::Pointer image);*/

  void setupFFTPipeline(ImageType::Pointer image);
private:
  Q_DECLARE_PUBLIC(QtImageEditor);
};

QtImageEditorPrivate::QtImageEditorPrivate(QtImageEditor& obj)
  : q_ptr(&obj)
  , m_CurrentImageData(0)
  , m_ImageData(0)
  , m_OverlayWidget(0)
  , m_SigmaSpinBox(0)
  , m_TabWidget(0)
  , m_XOrderSpinBox(0)
  , m_YOrderSpinBox(0)
  , m_ZOrderSpinBox(0)
{
  this->m_ImageData = 0;
  this->m_FFTFilter = FFTType::New();
  this->m_FFTShiftFilter = FFTShiftFilterType::New();
  this->m_MultiplyFilter = MultiplyFilterType::New();
  this->m_InverseFFTFilter = InverseFFTType::New();
}

void QtImageEditorPrivate::setupUi(QDialog* widgetToSetup)
{
  Q_UNUSED(widgetToSetup);
  Q_Q(QtImageEditor);

  QGridLayout* gridLayout = qobject_cast<QGridLayout*>(q->layout());
  QWidget* imageControls = gridLayout->itemAtPosition(1,0)->widget();

  this->m_TabWidget = new QTabWidget(q);
  this->m_TabWidget->insertTab(0, imageControls, "Controls");
  this->m_TabWidget->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);

  QWidget *filterControlWidget = new QWidget(this->m_TabWidget);
  this->m_TabWidget->insertTab(1, filterControlWidget, "FFT Filter");

  /// \todo create ui file
  QGridLayout *filterGridLayout = new QGridLayout(filterControlWidget);
  filterGridLayout->setContentsMargins(0, 0, 0, 0);

  QLabel *orderLabel = new QLabel(filterControlWidget);
  orderLabel->setText("Order:");
  filterGridLayout->addWidget(orderLabel, 0, 0);

  this->m_XOrderSpinBox = new QSpinBox(filterControlWidget);
  filterGridLayout->addWidget(this->m_XOrderSpinBox, 0, 1);

  this->m_YOrderSpinBox = new QSpinBox(filterControlWidget);
  filterGridLayout->addWidget(this->m_YOrderSpinBox, 0, 2);

  this->m_ZOrderSpinBox = new QSpinBox(filterControlWidget);
  filterGridLayout->addWidget(this->m_ZOrderSpinBox, 0, 3);

  QLabel *sigmaLabel = new QLabel(filterControlWidget);
  sigmaLabel->setText("Sigma:");
  filterGridLayout->addWidget(sigmaLabel, 2, 0);

  this->m_SigmaSpinBox = new QDoubleSpinBox();
  this->m_SigmaSpinBox->setMinimum(1);
  filterGridLayout->addWidget(this->m_SigmaSpinBox, 2, 1, 1, 3);

  QPushButton *fftButton = new QPushButton(filterControlWidget);
  fftButton->setText("FFT Gaussian Derivative");
  filterGridLayout->addWidget(fftButton, 0, 4);

  QPushButton *gaussianButton = new QPushButton(filterControlWidget);
  gaussianButton->setText("Gaussian Derivative in Frequency Domain");
  filterGridLayout->addWidget(gaussianButton, 2, 4);

  this->m_OverlayWidget = new QtOverlayControlsWidget(this->m_TabWidget);
  this->m_OverlayWidget->setSliceView(q->sliceView());
  this->m_TabWidget->insertTab(2, this->m_OverlayWidget, "Overlay");

  QDialogButtonBox* buttons = q->findChild<QDialogButtonBox*>();
  QPushButton *loadImageButton = new QPushButton(buttons);
  loadImageButton->setText("Load");
  buttons->addButton(loadImageButton, QDialogButtonBox::ActionRole);

  gridLayout->addWidget(this->m_TabWidget, 1, 0, 1, 2);

  QObject::connect(fftButton, SIGNAL(clicked()),
                   q, SLOT(useNewFilter()));
  QObject::connect(gaussianButton, SIGNAL(clicked()),
                   q, SLOT(applyFilter()));
  QObject::connect(loadImageButton, SIGNAL(clicked()),
                   q, SLOT(loadInputImage()));
}

void QtImageEditorPrivate::setupFFTPipeline(ImageType::Pointer image)
{
  this->m_FFTFilter->SetInput( image );
  this->m_MultiplyFilter->SetInput1( this->m_FFTFilter->GetOutput() );
  this->m_MultiplyFilter->SetInput2( this->m_FFTShiftFilter->GetOutput() );
  this->m_InverseFFTFilter->SetInput( this->m_MultiplyFilter->GetOutput() );
}


QtImageEditorPrivate::GaussianDerivativeSourceType::Pointer
QtImageEditorPrivate::createGaussianDerivative(
  GaussianDerivativeSourceType::OrdersType order)
{
  GaussianDerivativeSourceType::Pointer gaussianDerivativeSource =
    GaussianDerivativeSourceType::New();

  ComplexImageType::ConstPointer transformedInput =
    this->m_FFTFilter->GetOutput();
  const ComplexImageType::RegionType inputRegion(
    transformedInput->GetLargestPossibleRegion() );
  const ComplexImageType::SizeType inputSize =
    inputRegion.GetSize();
  const ComplexImageType::SpacingType inputSpacing =
    transformedInput->GetSpacing();
  const ComplexImageType::PointType inputOrigin =
    transformedInput->GetOrigin();
  const ComplexImageType::DirectionType inputDirection =
    transformedInput->GetDirection();

  gaussianDerivativeSource->SetSize( inputSize );
  gaussianDerivativeSource->SetSpacing( inputSpacing );
  gaussianDerivativeSource->SetOrigin( inputOrigin );
  gaussianDerivativeSource->SetDirection( inputDirection );

  GaussianDerivativeSourceType::SigmasType sigmas;
  GaussianDerivativeSourceType::PointType mean;
  const double sigma = this->m_SigmaSpinBox->value();

  for( unsigned int ii = 0; ii < 3; ++ii )
    {
    const double halfLength = inputSize[ii]  / 2.0;
    sigmas[ii] = sigma;
    mean[ii] = inputOrigin[ii] + halfLength;
    }
  mean = inputDirection * mean;
  gaussianDerivativeSource->SetSigmas( sigmas );
  gaussianDerivativeSource->SetMean( mean );
  gaussianDerivativeSource->SetOrders( order );

  gaussianDerivativeSource->Update();

  return gaussianDerivativeSource;
}


QtImageEditor::QtImageEditor(QWidget* _parent, Qt::WindowFlags fl )
  : Superclass( _parent, fl )
  , d_ptr(new QtImageEditorPrivate(*this))
{
  Q_D(QtImageEditor);
  this->setAttribute(Qt::WA_WState_ExplicitShowHide, true);
  this->setAttribute(Qt::WA_WState_Hidden, false);
  //this->setAttribute(Qt::WA_DontShowOnScreen, true);
  d->setupUi(this);
  this->layout()->invalidate();
  this->layout()->update();
  this->layout()->activate();
  //qDebug() << "1" << sizeHint() << this->geometry().size() << this->layout()->minimumSize() << this->minimumSizeHint();
  forceUpdate(this);
  //qDebug() << "2" << sizeHint() << this->geometry().size() << d->m_TabWidget->sizeHint();
  this->resize(this->sizeHint());
  //d->m_TabWidget->layout()->invalidate();
  //d->m_TabWidget->layout()->activate();
  //qDebug() << "3" << sizeHint() << this->geometry().size() << d->m_TabWidget->sizeHint();
  this->resize(this->sizeHint());
  //this->adjustSize();
  //qDebug() << "4" << sizeHint() << this->geometry().size();
  this->setAttribute(Qt::WA_WState_ExplicitShowHide, true);
  this->setAttribute(Qt::WA_WState_Hidden, true);
//  this->setAttribute(Qt::WA_DontShowOnScreen, false);
}


QtImageEditor::~QtImageEditor()
{
}


void QtImageEditor::setInputImage(ImageType* newImageData)
{
  Q_D(QtImageEditor);
  this->Superclass::setInputImage(newImageData);

  d->m_ImageData = newImageData;
  d->m_CurrentImageData = newImageData;
  d->setupFFTPipeline(d->m_ImageData);
}


void QtImageEditor::setOverlayImage(OverlayImageType* newOverlayImage)
{
  Q_D(QtImageEditor);
  this->Superclass::setOverlayImage(newOverlayImage);
  d->m_OverlayWidget->setInputOverlay(newOverlayImage);
}


void QtImageEditor::setDisplaySigma(double sigma)
{
  Q_D(QtImageEditor);
  d->m_SigmaSpinBox->setValue(sigma);
}


void QtImageEditor::applyFFT()
{
  Q_D(QtImageEditor);
  if (d->m_ImageData == NULL)
    {
    qDebug() << "No image to transform";
    return;
    }

  qDebug() << "Start FFT";

  TimeProbe clockFFT;
  clockFFT.Start();
  d->m_FFTFilter->Update();
  clockFFT.Stop();

  qDebug() << "FFT total time:" << clockFFT.GetTotal();
}


void QtImageEditor::applyFilter()
{
  Q_D(QtImageEditor);
  TimeProbe clockMultiply;
  this->applyFFT();
  clockMultiply.Start();
  QtImageEditorPrivate::GaussianDerivativeSourceType::OrdersType order;
  order[0] = d->m_XOrderSpinBox->value();
  order[1] = d->m_YOrderSpinBox->value();
  order[2] = d->m_ZOrderSpinBox->value();
  QtImageEditorPrivate::GaussianDerivativeSourceType::Pointer gaussianFilter =
    d->createGaussianDerivative( order );
  d->m_FFTShiftFilter->SetInput(gaussianFilter->GetOutput());
  d->m_FFTShiftFilter->Update();
  d->m_MultiplyFilter->Update();
  clockMultiply.Stop();

  qDebug() << "Multiply total time:" << clockMultiply.GetTotal();
  this->applyInverseFFT();
}


void QtImageEditor::applyInverseFFT()
{
  Q_D(QtImageEditor);
  if (d->m_ImageData == NULL)
    {
    return;
    }

  qDebug()<<"Start IFFT";

  TimeProbe clockIFFT;
  clockIFFT.Start();
  d->m_InverseFFTFilter->Update();
  clockIFFT.Stop();

  qDebug() << "IFFT total time:" << clockIFFT.GetTotal();

  this->Superclass::setInputImage(d->m_InverseFFTFilter->GetOutput());
}


void QtImageEditor::useNewFilter()
{
  Q_D(QtImageEditor);
  QtImageEditorPrivate::FFTGaussianDerivativeIFFTType::Pointer FFTgaussianDerivativeIFFT =
    QtImageEditorPrivate::FFTGaussianDerivativeIFFTType::New();
  FFTgaussianDerivativeIFFT->SetInput(d->m_ImageData);
  QtImageEditorPrivate::GaussianDerivativeSourceType::OrdersType order;
  order[0] = d->m_XOrderSpinBox->value();
  order[1] = d->m_YOrderSpinBox->value();
  order[2] = d->m_ZOrderSpinBox->value();
  FFTgaussianDerivativeIFFT->SetOrders(order);
  QtImageEditorPrivate::GaussianDerivativeSourceType::SigmasType sigmas;
  sigmas.Fill(d->m_SigmaSpinBox->value());
  FFTgaussianDerivativeIFFT->SetSigmas(sigmas);
  FFTgaussianDerivativeIFFT->Update();
  this->Superclass::setInputImage(FFTgaussianDerivativeIFFT->GetOutput());
}

void QtImageEditor::setControlsVisible(bool controlsVisible)
{
  Q_D(QtImageEditor);
  d->m_TabWidget->setVisible(controlsVisible);
  this->Superclass::setControlsVisible(controlsVisible);
}

} // End namespace tube
