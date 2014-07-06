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
#include <QFileDialog>

// QtImageViewer includes
#include "QtGlSliceView.h"
#include "ui_QtSlicerHelpGUI.h"

// ImageEditor includes
#include "QtImageEditor.h"
#include "QtOverlayControlsWidget.h"
#include "ui_QtOverlayControlsWidgetGUI.h"

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

class QtImageEditor::Internals
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
protected:
  ImageType                  *m_CurrentImageData;
  ImageType                  *m_ImageData;
  QtOverlayControlsWidget    *m_OverlayWidget;
  FFTType::Pointer            m_FFTFilter;
  QLineEdit                  *m_SigmaLineEdit;
  MultiplyFilterType::Pointer m_MultiplyFilter;
  FFTShiftFilterType::Pointer m_FFTShiftFilter;
  InverseFFTType::Pointer     m_InverseFFTFilter;
  QDialog                    *m_HelpDialog;
  QLineEdit                  *m_Order_x;
  QLineEdit                  *m_Order_y;
  QLineEdit                  *m_Order_z;


  GaussianDerivativeSourceType::Pointer createGaussianDerivative(
    GaussianDerivativeSourceType::OrdersType order);/*,
    ImageType::Pointer image);*/

  void setupFFTPipeline(ImageType::Pointer image);

  friend class QtImageEditor;
};


void QtImageEditor::Internals::setupFFTPipeline(ImageType::Pointer image)
{
  this->m_FFTFilter->SetInput( image );
  this->m_MultiplyFilter->SetInput1( this->m_FFTFilter->GetOutput() );
  this->m_MultiplyFilter->SetInput2( this->m_FFTShiftFilter->GetOutput() );

  this->m_InverseFFTFilter->SetInput( this->m_MultiplyFilter->GetOutput() );
}


QtImageEditor::Internals::GaussianDerivativeSourceType::Pointer
QtImageEditor::Internals::createGaussianDerivative(
  GaussianDerivativeSourceType::OrdersType order)/*,
  ImageType::Pointer image)*/
{
  GaussianDerivativeSourceType::Pointer gaussianDerivativeSource =
  GaussianDerivativeSourceType::New();

  ComplexImageType::ConstPointer transformedInput
    = this->m_FFTFilter->GetOutput();
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

  GaussianDerivativeSourceType::SigmasType sigma;
  GaussianDerivativeSourceType::PointType mean;
  const double sigmaLineEdit = this->m_SigmaLineEdit->text().toDouble();

  for( unsigned int ii = 0; ii < 3; ++ii )
    {
    const double halfLength = inputSize[ii]  / 2.0;
    sigma[ii] = sigmaLineEdit;
    mean[ii] = inputOrigin[ii] + halfLength;
    }
  mean = inputDirection * mean;
  gaussianDerivativeSource->SetSigmas( sigma );
  gaussianDerivativeSource->SetMean( mean );
  gaussianDerivativeSource->SetOrders(order);

  gaussianDerivativeSource->Update();

  return gaussianDerivativeSource;
}


QtImageEditor::QtImageEditor(QWidget* _parent, Qt::WindowFlags fl ) :
  QDialog( _parent, fl )
{
  this->m_Internals = new Internals;
  this->m_Internals->m_ImageData = 0;
  this->m_Internals->m_FFTFilter = Internals::FFTType::New();
  this->m_Internals->m_FFTShiftFilter = Internals::FFTShiftFilterType::New();
  this->m_Internals->m_MultiplyFilter = Internals::MultiplyFilterType::New();
  this->m_Internals->m_InverseFFTFilter = Internals::InverseFFTType::New();

  this->setupUi(this);
  this->Controls->setSliceView(this->OpenGlWindow);

  QTabWidget *tabWidget = new QTabWidget(this);
  tabWidget->setMaximumHeight(300);

  tabWidget->insertTab(0, this->Controls, "Controls");

  QWidget *filterControlWidget = new QWidget(tabWidget);
  tabWidget->insertTab(1, filterControlWidget, "FFT Filter");

  QGridLayout *filterGridLayout = new QGridLayout(filterControlWidget);
  this->m_Internals->m_SigmaLineEdit = new QLineEdit();
  this->m_Internals->m_SigmaLineEdit->setMaximumWidth(80);
  this->m_Internals->m_SigmaLineEdit->setText("1");
  filterGridLayout->addWidget(this->m_Internals->m_SigmaLineEdit, 2, 1);

  QLabel *orderLabel = new QLabel(filterControlWidget);
  orderLabel->setText("Order Vector:");
  filterGridLayout->addWidget(orderLabel, 0, 0);

  this->m_Internals->m_Order_x = new QLineEdit();
  filterGridLayout->addWidget(this->m_Internals->m_Order_x, 0, 1);

  this->m_Internals->m_Order_y = new QLineEdit();
  filterGridLayout->addWidget(this->m_Internals->m_Order_y, 0, 2);

  this->m_Internals->m_Order_z = new QLineEdit();
  filterGridLayout->addWidget(this->m_Internals->m_Order_z, 0, 3);

  QLabel *sigmaLabel = new QLabel(filterControlWidget);
  sigmaLabel->setText("Sigma:");
  filterGridLayout->addWidget(sigmaLabel, 2, 0);

  QPushButton *fftButton = new QPushButton();
  fftButton->setText("FFT");
  filterGridLayout->addWidget(fftButton, 0, 4);

  QPushButton *inverseFFTButton = new QPushButton();
  inverseFFTButton->setText("Inverse FFT");
  filterGridLayout->addWidget(inverseFFTButton, 1, 4);

  QPushButton *gaussianButton = new QPushButton();
  gaussianButton->setText("FFT-BLUR-INVERSE FFT");
  filterGridLayout->addWidget(gaussianButton, 2, 4);

  this->m_Internals->m_OverlayWidget = new QtOverlayControlsWidget(tabWidget);
  this->m_Internals->m_OverlayWidget->setSliceView(this->OpenGlWindow);
  tabWidget->insertTab(2, this->m_Internals->m_OverlayWidget, "Overlay");

  QDialogButtonBox *buttons = new QDialogButtonBox(Qt::Horizontal);
  buttons->addButton(this->ButtonOk, QDialogButtonBox::AcceptRole);
  QPushButton *loadImageButton = new QPushButton();
  loadImageButton->setText("Load");

  buttons->addButton(loadImageButton, QDialogButtonBox::ActionRole);
  buttons->addButton(this->ButtonHelp, QDialogButtonBox::HelpRole);


  this->gridLayout->addWidget(buttons, 2, 0,1,2);
  this->gridLayout->addWidget(tabWidget, 1, 0,1,2);

  QObject::connect(ButtonOk, SIGNAL(clicked()), this, SLOT(accept()));
  QObject::connect(ButtonHelp, SIGNAL(toggled(bool)), this, SLOT(showHelp(bool)));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), OpenGlWindow,
                   SLOT(changeSlice(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), SliceNumSlider,
                   SLOT(setValue(int)));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(fftButton, SIGNAL(clicked()), this, SLOT(useNewFilter()));
  QObject::connect(inverseFFTButton, SIGNAL(clicked()), this, SLOT(applyInverseFFT()));
  QObject::connect(gaussianButton, SIGNAL(clicked()), this, SLOT(applyFilter()));
  QObject::connect(this->m_Internals->m_SigmaLineEdit, SIGNAL(textChanged(QString)), this,
                   SLOT(setDisplaySigma(QString)));
  QObject::connect(OpenGlWindow, SIGNAL(orientationChanged(int)), this,
                   SLOT(setMaximumSlice()));
//  QObject::connect(OpenGlWindow, SIGNAL(viewDetailsChanged(int)), this,
//                   SLOT(toggleTextEdit(int)));
  QObject::connect(loadImageButton, SIGNAL(clicked()), this, SLOT(loadImage()));
}

QtImageEditor::~QtImageEditor()
{
}


void QtImageEditor::showHelp(bool checked)
{
  if(!checked && this->m_Internals->m_HelpDialog != 0)
    {
    this->m_Internals->m_HelpDialog->reject();
    }
  else
    {
    this->OpenGlWindow->showHelp();
    this->m_Internals->m_HelpDialog = this->OpenGlWindow->helpWindow();
    if(this->m_Internals->m_HelpDialog != 0)
      {
      QObject::connect(this->m_Internals->m_HelpDialog, SIGNAL(rejected()),
                       this, SLOT(hideHelp()),Qt::UniqueConnection);
      }
    }
}


void QtImageEditor::hideHelp()
{
  this->ButtonHelp->setChecked(false);
}


void QtImageEditor::setInputImage(ImageType* newImageData)
{
  if (this->m_Internals->m_ImageData == 0)
    {
    this->m_Internals->m_ImageData = newImageData;
    }
  this->m_Internals->m_CurrentImageData = newImageData;
  this->OpenGlWindow->setInputImage(newImageData);
  this->setMaximumSlice();
  this->OpenGlWindow->changeSlice(((this->OpenGlWindow->maxSliceNum() -1)/2));
  this->setDisplaySliceNumber(static_cast<int>
                              (this->OpenGlWindow->sliceNum()));
  this->Controls->setInputImage();
  this->OpenGlWindow->update();

  this->m_Internals->setupFFTPipeline(this->m_Internals->m_ImageData);
}


void QtImageEditor::setDisplaySliceNumber(int number)
{
  QString tempchar = QString::number(number);
  this->SliceValue->setText(tempchar);
}


bool QtImageEditor::loadImage(QString filePathToLoad)
{
  Internals::ReaderType::Pointer reader = Internals::ReaderType::New();
  if( filePathToLoad.isEmpty() )
    {
    filePathToLoad = QFileDialog::getOpenFileName(
        0,"", QDir::currentPath());
    }

  if(filePathToLoad.isEmpty())
    {
    return false;
    }
  reader->SetFileName( filePathToLoad.toLatin1().data() );
  QFileInfo filePath(filePathToLoad);
  setWindowTitle(filePath.fileName());

  std::cout << "Loading image " << filePathToLoad.toStdString() << "... ";
  try
    {
    reader->Update();
    }
  catch (ExceptionObject & e)
    {
    std::cerr << "exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "done!" << std::endl;
  this->setInputImage( reader->GetOutput() );

  return true;
}


void QtImageEditor::loadOverlay(QString overlayImagePath)
{
  this->m_Internals->m_OverlayWidget->loadOverlay(overlayImagePath);
}


void QtImageEditor::setDisplaySigma(QString value)
{
  this->m_Internals->m_SigmaLineEdit->setText(value);
}


void QtImageEditor::applyFFT()
{
  if(this->m_Internals->m_ImageData == NULL)
    {
    qDebug() << "No image to transform";
    return;
    }

  qDebug() << "Start FFT";

  TimeProbe clockFFT;
  clockFFT.Start();
  this->m_Internals->m_FFTFilter->Update();
  clockFFT.Stop();

  qDebug() << "FFT total time:" << clockFFT.GetTotal();
}


void QtImageEditor::applyFilter()
{
  TimeProbe clockMultiply;
  applyFFT();
  clockMultiply.Start();
  Internals::GaussianDerivativeSourceType::OrdersType order;
  order[0] = this->m_Internals->m_Order_x->text().toInt();
  order[1] = this->m_Internals->m_Order_y->text().toInt();
  order[2] = this->m_Internals->m_Order_z->text().toInt();
  Internals::GaussianDerivativeSourceType::Pointer gaussianFilter =
  this->m_Internals->createGaussianDerivative( order );
  this->m_Internals->m_FFTShiftFilter->SetInput(gaussianFilter->GetOutput());
  this->m_Internals->m_FFTShiftFilter->Update();
  this->m_Internals->m_MultiplyFilter->Update();
  clockMultiply.Stop();

  qDebug() << "Multiply total time:" << clockMultiply.GetTotal();
  applyInverseFFT();
}


void QtImageEditor::applyInverseFFT()
{
  if(this->m_Internals->m_ImageData == NULL)
    {
    return;
    }

  qDebug()<<"Start IFFT";

  TimeProbe clockIFFT;
  clockIFFT.Start();
  this->m_Internals->m_InverseFFTFilter->Update();
  clockIFFT.Stop();

  qDebug() << "IFFT total time:" << clockIFFT.GetTotal();

  this->setInputImage(
    this->m_Internals->m_InverseFFTFilter->GetOutput());
  this->OpenGlWindow->update();
}


void QtImageEditor::useNewFilter()
{
  Internals::FFTGaussianDerivativeIFFTType::Pointer FFTgaussianDerivativeIFFT =
  Internals::FFTGaussianDerivativeIFFTType::New();
  FFTgaussianDerivativeIFFT->SetInput(this->m_Internals->m_ImageData);
  Internals::GaussianDerivativeSourceType::OrdersType order;
  Internals::GaussianDerivativeSourceType::SigmasType sigma;
  order[0] = this->m_Internals->m_Order_x->text().toInt();
  order[1] = this->m_Internals->m_Order_y->text().toInt();
  order[2] = this->m_Internals->m_Order_z->text().toInt();
  FFTgaussianDerivativeIFFT->SetOrders(order);
  sigma.Fill(this->m_Internals->m_SigmaLineEdit->text().toDouble());
  FFTgaussianDerivativeIFFT->SetSigmas(sigma);
  FFTgaussianDerivativeIFFT->Update();
  this->setInputImage(FFTgaussianDerivativeIFFT->GetOutput());
  this->OpenGlWindow->update();
}

void QtImageEditor::setMaximumSlice()
{
  this->SliceNumSlider->setMaximum(static_cast<int>
                                   (this->OpenGlWindow->maxSliceNum() -1));
  this->SliceNumSlider->setValue(static_cast<int>
                                 (this->SliceValue->text().toInt()));
}

} // End namespace tube
