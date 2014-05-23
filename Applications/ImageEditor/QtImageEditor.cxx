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
#include <QDebug>
#include <QDialogButtonBox>
#include <QFileDialog>

//QtImageViewer includes
#include "QtGlSliceView.h"
#include "ui_QtSlicerHelpGUI.h"

//ImageEditor includes
#include "QtImageEditor.h"
#include "QtOverlayControlsWidget.h"
#include "ui_QtOverlayControlsWidgetGUI.h"


namespace tube {


QtImageEditor::QtImageEditor(QWidget* parent, Qt::WindowFlags fl ) :
  QDialog( parent, fl )
{
  this->m_ImageData = 0;
  this->m_FFTFilter = FFTType::New();
  this->m_ComplexeImage = 0;
  this->setupUi(this);
  this->Controls->setSliceView(this->OpenGlWindow);

  QTabWidget *tabWidget = new QTabWidget(this);
  tabWidget->setMaximumHeight(300);

  tabWidget->insertTab(0, this->Controls, "Controls");

  QWidget *filterControlWidget = new QWidget(tabWidget);
  tabWidget->insertTab(1, filterControlWidget, "FFT Filter");

  QGridLayout *filterGridLayout = new QGridLayout(filterControlWidget);
  this->m_SigmaLineEdit = new QLineEdit();
  this->m_SigmaLineEdit->setMaximumWidth(80);
  this->m_SigmaLineEdit->setText("1");
  filterGridLayout->addWidget(this->m_SigmaLineEdit, 2, 1);

  QLabel *sigmaLabel = new QLabel(filterControlWidget);
  sigmaLabel->setText("Sigma:");
  filterGridLayout->addWidget(sigmaLabel, 2, 0);

  QPushButton *fftButton = new QPushButton();
  fftButton->setText("FFT");
  filterGridLayout->addWidget(fftButton, 0, 2);

  QPushButton *inverseFFTButton = new QPushButton();
  inverseFFTButton->setText("Inverse FFT");
  filterGridLayout->addWidget(inverseFFTButton, 1, 2);

  QPushButton *gaussianButton = new QPushButton();
  gaussianButton->setText("FFT-BLUR-INVERSE FFT");
  filterGridLayout->addWidget(gaussianButton, 2, 2);

  this->m_OverlayWidget = new QtOverlayControlsWidget(tabWidget);
  this->m_OverlayWidget->setSliceView(this->OpenGlWindow);
  tabWidget->insertTab(2, this->m_OverlayWidget, "Overlay");

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
  QObject::connect(fftButton, SIGNAL(clicked()), this, SLOT(applyFFT()));
  QObject::connect(inverseFFTButton, SIGNAL(clicked()), this, SLOT(applyInverseFFT()));
  QObject::connect(gaussianButton, SIGNAL(clicked()), this, SLOT(applyFilter()));
  QObject::connect(this->m_SigmaLineEdit, SIGNAL(textChanged(QString)), this,
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
  if(!checked && this->m_HelpDialog != 0)
    {
    this->m_HelpDialog->reject();
    }
  else
    {
    this->OpenGlWindow->showHelp();
    this->m_HelpDialog = this->OpenGlWindow->helpWindow();
    if(this->m_HelpDialog != 0)
      {
      QObject::connect(m_HelpDialog, SIGNAL(rejected()), this,
                       SLOT(hideHelp()),Qt::UniqueConnection);
      }
    }
}


void QtImageEditor::hideHelp()
{
  this->ButtonHelp->setChecked(false);
}


void QtImageEditor::setInputImage(ImageType * newImData)
{
  this->m_ImageData = newImData;
  this->OpenGlWindow->setInputImage(newImData);
  setMaximumSlice();
  this->OpenGlWindow->changeSlice(((this->OpenGlWindow->maxSliceNum() -1)/2));
  this->setDisplaySliceNumber(static_cast<int>
                                (this->OpenGlWindow->sliceNum()));
  this->Controls->setInputImage();
  this->OpenGlWindow->update();
}


void QtImageEditor::setDisplaySliceNumber(int number)
{
  QString tempchar = QString::number(number);
  this->SliceValue->setText(tempchar);
}


bool QtImageEditor::loadImage(QString filePathToLoad)
{
  ReaderType::Pointer reader = ReaderType::New();
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

  qDebug() << "loading image " << filePathToLoad << " ... ";
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Done!" << std::endl;
  setInputImage( reader->GetOutput() );

  show();
  return true;
}


void QtImageEditor::loadOverlay(QString overlayImagePath)
{
  this->m_OverlayWidget->loadOverlay(overlayImagePath);
}


void QtImageEditor::setDisplaySigma(QString value)
{
  this->m_SigmaLineEdit->setText(value);
}


void QtImageEditor::blurFilter()
{
  this->m_FilterX = FilterType::New();
  this->m_FilterY = FilterType::New();
  this->m_FilterZ = FilterType::New();

  this->m_FilterX->SetDirection( 0 );   // 0 --> X direction
  this->m_FilterY->SetDirection( 1 );   // 1 --> Y direction
  this->m_FilterZ->SetDirection( 2 );   // 2 --> Z direction

  this->m_FilterX->SetOrder( FilterType::ZeroOrder );
  this->m_FilterY->SetOrder( FilterType::ZeroOrder );
  this->m_FilterZ->SetOrder( FilterType::ZeroOrder );

  this->m_FilterX->SetNormalizeAcrossScale( true );
  this->m_FilterY->SetNormalizeAcrossScale( true );
  this->m_FilterZ->SetNormalizeAcrossScale( true );

  this->m_FilterX->SetInput( this->m_FFTFilter->GetOutput() );
  this->m_FilterY->SetInput( this->m_FilterX->GetOutput() );
  this->m_FilterY->SetInput( this->m_FilterY->GetOutput() );

  const double sigma = m_SigmaLineEdit->text().toDouble();

  this->m_FilterX->SetSigma( sigma );
  this->m_FilterY->SetSigma( sigma );
  this->m_FilterZ->SetSigma( sigma );

  this->m_FilterX->Update();
  this->m_FilterY->Update();
  this->m_FilterZ->Update();
}


void QtImageEditor::applyFFT()
{
  if(this->m_ImageData == NULL)
    {
    return;
    }
  qDebug()<<"Start FFT";

  //Compute the FFT
  this->m_FFTFilter->SetInput(this->m_ImageData);
  this->m_FFTFilter->Update();
}


void QtImageEditor::applyInverseFFT()
{
  if(this->m_ImageData == NULL)
    {
    return;
    }
  qDebug()<<"Start IFFT";
  this->m_InverseFFTFilter = InverseFFTType::New();
  this->m_InverseFFTFilter->SetInput(this->m_FFTFilter->GetOutput());
  this->m_InverseFFTFilter->Update();
  setInputImage(this->m_InverseFFTFilter->GetOutput());
}


void QtImageEditor::applyFilter()
{
  itk::TimeProbe clockFFT;
  itk::TimeProbe clockGaussian;
  itk::TimeProbe clockIFFT;

  clockFFT.Start();
  applyFFT();
  clockFFT.Stop();

  qDebug()<<"Output FFT mean"<<clockFFT.GetMean();
  qDebug()<<"Output FFT total"<<clockFFT.GetTotal();
  // A Gaussian is used here to create a low-pass filter.
  clockGaussian.Start();
  qDebug()<<"Start gaussian filter";
  GaussianDerivativeSourceType::VectorType order;
  order[0] = 0;
  order[1] = 1;
  order[2] = 1;

  GaussianDerivativeSourceType::Pointer gaussian =
    this->createGaussianDerivative(order);
  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
  fftShiftFilter->SetInput( gaussian->GetOutput() );

  clockGaussian.Stop();
  qDebug()<<"Output Shift Filter mean"<<clockGaussian.GetMean();
  qDebug()<<"Output Shift Filter meantime"<<QString::fromStdString(
  clockGaussian.GetUnit());
  qDebug()<<"Output Shift Filter"<<clockGaussian.GetTotal();
  clockGaussian.Start();

  MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
  multiplyFilter->SetInput1( this->m_FFTFilter->GetOutput() );
  multiplyFilter->SetInput2( fftShiftFilter->GetOutput() );

  clockGaussian.Stop();
  qDebug()<<"Output Multiply Filter"<<clockGaussian.GetMean();
  qDebug()<<"Output Multiply Filter"<<clockGaussian.GetTotal();

  clockIFFT.Start();
  qDebug()<<"Start IFFT";
  this->m_InverseFFTFilter = InverseFFTType::New();
  this->m_InverseFFTFilter->SetInput(multiplyFilter->GetOutput());
  this->m_InverseFFTFilter->Update();

  clockIFFT.Stop();
  qDebug()<<"Mean Time"<<clockIFFT.GetMean();
  qDebug()<<"Total Time"<<clockIFFT.GetTotal();
  qDebug()<<"Total"<<clockIFFT.GetTotal() + clockFFT.GetTotal() +
  clockGaussian.GetTotal();

  setInputImage(m_InverseFFTFilter->GetOutput());
  this->OpenGlWindow->update();
  show();
}

QtImageEditor::GaussianDerivativeSourceType::Pointer
QtImageEditor::createGaussianDerivative(
    GaussianDerivativeSourceType::VectorType order)
{
  GaussianDerivativeSourceType::Pointer gaussianDerivativeSource =
  GaussianDerivativeSourceType::New();
  gaussianDerivativeSource->SetNormalized( true );
  this->m_FFTFilter->Update();

  ComplexImageType::ConstPointer transformedInput
    = this->m_FFTFilter->GetOutput();
  const ComplexImageType::RegionType inputRegion(
    transformedInput->GetLargestPossibleRegion() );
  const ComplexImageType::SizeType inputSize
    = inputRegion.GetSize();
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

  GaussianDerivativeSourceType::ArrayType sigma;
  GaussianDerivativeSourceType::PointType mean;
  const double sigmaLineEdit = m_SigmaLineEdit->text().toDouble();

  for( unsigned int ii = 0; ii < 3; ++ii )
    {
    const double halfLength = inputSize[ii]  / 2.0;
    sigma[ii] = sigmaLineEdit;
    mean[ii] = inputOrigin[ii] + halfLength;
    }
  mean = inputDirection * mean;
  gaussianDerivativeSource->SetSigma( sigma );
  gaussianDerivativeSource->SetMean( mean );
  gaussianDerivativeSource->SetOrdersVector(order);
  gaussianDerivativeSource->Update();
  return gaussianDerivativeSource;
}


void QtImageEditor::displayFFT()
{
  //Extract the real part
  RealFilterType::Pointer realFilter = RealFilterType::New();
  realFilter->SetInput(this->m_FFTFilter->GetOutput());
  realFilter->Update();

  RescaleFilterType::Pointer realRescaleFilter = RescaleFilterType::New();
  realRescaleFilter->SetInput(realFilter->GetOutput());
  realRescaleFilter->SetOutputMinimum(0);
  realRescaleFilter->SetOutputMaximum(255);
  realRescaleFilter->Update();

  //Extract the imaginary part
  ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
  imaginaryFilter->SetInput(this->m_FFTFilter->GetOutput());
  imaginaryFilter->Update();

  RescaleFilterType::Pointer imaginaryRescaleFilter = RescaleFilterType::New();
  imaginaryRescaleFilter->SetInput(imaginaryFilter->GetOutput());
  imaginaryRescaleFilter->SetOutputMinimum(0);
  imaginaryRescaleFilter->SetOutputMaximum(255);
  imaginaryRescaleFilter->Update();

  // Compute the magnitude
  ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
  modulusFilter->SetInput(this->m_FFTFilter->GetOutput());
  modulusFilter->Update();

  RescaleFilterType::Pointer magnitudeRescaleFilter = RescaleFilterType::New();
  magnitudeRescaleFilter->SetInput(modulusFilter->GetOutput());
  magnitudeRescaleFilter->SetOutputMinimum(0);
  magnitudeRescaleFilter->SetOutputMaximum(255);
  magnitudeRescaleFilter->Update();

  // Write the images
  WriterType::Pointer realWriter = WriterType::New();
  realWriter->SetFileName("real.png");
  realWriter->SetInput(realRescaleFilter->GetOutput());
  realWriter->Update();

  WriterType::Pointer imaginaryWriter = WriterType::New();
  imaginaryWriter->SetFileName("imaginary.png");
  imaginaryWriter->SetInput(imaginaryRescaleFilter->GetOutput());
  imaginaryWriter->Update();

  WriterType::Pointer magnitudeWriter = WriterType::New();
  magnitudeWriter->SetFileName("magnitude.png");
  magnitudeWriter->SetInput(magnitudeRescaleFilter->GetOutput());
  magnitudeWriter->Update();
}


void QtImageEditor::setMaximumSlice()
{
  this->SliceNumSlider->setMaximum(static_cast<int>
                                   (this->OpenGlWindow->maxSliceNum() -1));
  this->SliceNumSlider->setValue(static_cast<int>
                                 (this->SliceValue->text().toInt()));
}

}
