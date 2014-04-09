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
#include <QFileDialog>

//QtImageViewer includes
#include "QtGlSliceView.h"

//ImageEditor includes
#include "QtImageEditor.h"
#include "QtOverlayControlsWidget.h"

namespace tube{

QtImageEditor::QtImageEditor(QWidget* parent, Qt::WindowFlags fl ) :
  QDialog( parent, fl )
{
  this->m_ImageData = 0;
  this->setupUi(this);
  this->Controls->setSliceView(this->OpenGlWindow);

  QTabWidget *tabWidget = new QTabWidget(this);
  tabWidget->setMaximumHeight(300);

  tabWidget->insertTab(0, this->Controls, "Controls");

  QWidget *filterControlWidget = new QWidget(tabWidget);
  tabWidget->insertTab(1, filterControlWidget, "Filter");

  QGridLayout *filterGridLayout = new QGridLayout(filterControlWidget);
  this->m_SigmaLineEdit = new QLineEdit();
  this->m_SigmaLineEdit->setMaximumWidth(80);
  this->m_SigmaLineEdit->setText("0.2");
  filterGridLayout->addWidget(this->m_SigmaLineEdit, 0, 0);

  QPushButton *applyButton = new QPushButton();
  applyButton->setText("Apply");
  filterGridLayout->addWidget(applyButton, 0, 1);

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
  QObject::connect(ButtonHelp, SIGNAL(clicked()), OpenGlWindow, SLOT(showHelp()));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), OpenGlWindow,
                   SLOT(changeSlice(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), SliceNumSlider,
                   SLOT(setValue(int)));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(applyButton, SIGNAL(clicked()), this, SLOT(applyFilter()));
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

//void QtImageEditor::toggleTextEdit(int viewDetail)
//{
//  bool hidden;
//  viewDetail = VD_TEXTBOX ? hidden = false : hidden  = true;
//  this->Details->setHidden(hidden);
//}


void QtImageEditor::setInputImage(ImageType * newImData)
{
  this->m_ImageData = newImData;
  this->OpenGlWindow->setInputImage(newImData);
//  this->SliceNumSlider->setMaximum( static_cast<int>
//                                    (this->OpenGlWindow->maxSliceNum() -1));
  setMaximumSlice();
  this->OpenGlWindow->changeSlice(((this->OpenGlWindow->maxSliceNum() -1)/2));
  //this->SliceNumSlider->setValue(static_cast<int>
//                                   (this->OpenGlWindow->sliceNum()));
  this->setDisplaySliceNumber(static_cast<int>
                                (this->OpenGlWindow->sliceNum()));
  this->Controls->setInputImage();
  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}


void QtImageEditor::setDisplaySliceNumber(int number)
{
  QString tempchar = QString::number(number);
  this->SliceValue->setText(tempchar);
}


int QtImageEditor::loadImage(std::string path)
{
  ReaderType::Pointer reader = ReaderType::New();
  QString filePathToLoad = QString::fromStdString(path);
  if( filePathToLoad.isNull() )
    {
    filePathToLoad = QFileDialog::getOpenFileName(
        0,"", QDir::currentPath());
    }

  if(filePathToLoad.isEmpty())
    {
    return 0;
    }
  reader->SetFileName( filePathToLoad.toLatin1().data() );

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
}

int QtImageEditor::loadImage()
{
  ReaderType::Pointer reader = ReaderType::New();
  QString filePathToLoad = QFileDialog::getOpenFileName(
        0,"", QDir::currentPath());

  if(filePathToLoad.isEmpty())
    {
    return 0;
    }
  reader->SetFileName( filePathToLoad.toLatin1().data() );

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
}


int QtImageEditor::loadOverlay(std::string path)
{
  this->m_OverlayWidget->loadOverlay(path);
}


void QtImageEditor::setDisplaySigma(QString value)
{
  this->m_SigmaLineEdit->setText(value);
}


void QtImageEditor::applyFilter()
{
  if(this->m_ImageData == NULL)
    {
    return;
    }
  FilterType::Pointer filterX;
  FilterType::Pointer filterY;

  filterX = FilterType::New();
  filterY = FilterType::New();

  filterX->SetDirection( 0 );   // 0 --> X direction
  filterY->SetDirection( 1 );   // 1 --> Y direction

  filterX->SetOrder( FilterType::ZeroOrder );
  filterY->SetOrder( FilterType::ZeroOrder );

  filterX->SetNormalizeAcrossScale( true );
  filterY->SetNormalizeAcrossScale( true );

  filterX->SetInput( this->m_ImageData );
  filterY->SetInput( filterX->GetOutput() );

  const double sigma = m_SigmaLineEdit->text().toDouble();

  filterX->SetSigma( sigma );
  filterY->SetSigma( sigma );

  filterX->Update();
  filterY->Update();

  this->OpenGlWindow->update();
  setInputImage(filterY->GetOutput());
  show();
}


void QtImageEditor::setMaximumSlice()
{
  this->SliceNumSlider->setMaximum(static_cast<int>
                                   (this->OpenGlWindow->maxSliceNum() -1));
  this->SliceNumSlider->setValue(static_cast<int>
                                 (this->SliceValue->text().toInt()));
}

}
