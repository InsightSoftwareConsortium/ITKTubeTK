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

//#ifndef __QtControlView_cxx
//#define __QtControlView_cxx

//Qt includes
#include <QDebug>
#include <QLineEdit>
#include <QSlider>

//QtImageViewer includes
#include "QtGlSliceView.h"

//ImageEditor includes
#include "QtControlView.h"

namespace tube{

QtControlView::QtControlView(QWidget* parent, Qt::WindowFlags fl ) :
  QDialog( parent, fl )
{
  this->setupUi(this);
  setTab();
  QObject::connect(ButtonOk, SIGNAL(clicked()), this, SLOT(accept()));
  QObject::connect(ButtonHelp, SIGNAL(clicked()), OpenGlWindow, SLOT(showHelp()));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), OpenGlWindow,
                   SLOT(changeSlice(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), SliceNumSlider,
                   SLOT(setValue(int)));
  QObject::connect(OpenGlWindow, SIGNAL(positionChanged(int,int,int,double)), this,
                   SLOT(setDisplayPosition(int,int,int,double)));
  QObject::connect(IntensityMax, SIGNAL(sliderMoved(int)), OpenGlWindow,
                   SLOT(setMaxIntensity(int)));
  QObject::connect(IntensityMin, SIGNAL(sliderMoved(int)), OpenGlWindow,
                   SLOT(setMinIntensity(int)));
  QObject::connect(ZoomIn, SIGNAL(clicked()), OpenGlWindow, SLOT(zoomIn()));
  QObject::connect(ZoomOut, SIGNAL(clicked()), OpenGlWindow, SLOT(zoomOut()));
  QObject::connect(SliceNumSlider, SIGNAL(sliderMoved(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(OpenGlWindow, SIGNAL(sliceNumChanged(int)), this,
                   SLOT(setDisplaySliceNumber(int)));
  QObject::connect(OpenGlWindow, SIGNAL(maxIntensityChanged(int)), this,
                   SLOT(setDisplayIMax(int)));
  QObject::connect(OpenGlWindow, SIGNAL(minIntensityChanged(int)), this,
                   SLOT(setDisplayIMin(int)));
  QObject::connect(OpenGlWindow, SIGNAL(maxIntensityChanged(int)), IntensityMax,
                   SLOT(setValue(int)));
  QObject::connect(OpenGlWindow, SIGNAL(minIntensityChanged(int)), IntensityMin,
                   SLOT(setValue(int)));
  QObject::connect(OpenGlWindow, SIGNAL(updateDetails(QString)), Details,
                   SLOT(setText(QString)));
  QObject::connect(m_ApplyButton, SIGNAL(clicked()), this, SLOT(setFilter()));
  QObject::connect(m_SigmaLineEdit, SIGNAL(textChanged(QString)), this,
                   SLOT(setDisplaySigma(QString)));
  QObject::connect(OpenGlWindow, SIGNAL(orientationChanged(int)), this,
                   SLOT(setMaximumSlice()));

  this->m_ImageData = 0;

}

QtControlView::~QtControlView()
{
}


void QtControlView::setMaximumSlice()
{
  this->SliceNumSlider->setMaximum(static_cast<int>
                                   (this->OpenGlWindow->maxSliceNum() -1));
  this->SliceNumSlider->setValue(static_cast<int>
                                 (this->SliceValue->text().toInt()));
}


void QtControlView::setTab()
{
  this->m_TabWidget = new QTabWidget(this);
  this->m_TabWidget->setMaximumHeight(300);

  this->Controls->setParent(this->m_TabWidget);
  this->m_TabWidget->insertTab(0, this->Controls, "Control Widgets");

  this->m_FilterControlWidget = new QWidget(this->m_TabWidget);
  this->m_TabWidget->insertTab(1, this->m_FilterControlWidget, "Filter");

  this->m_FilterGridLayout = new QGridLayout(this->m_FilterControlWidget);
  this->m_SigmaLineEdit = new QLineEdit();
  this->m_SigmaLineEdit->setMaximumWidth(80);
  this->m_SigmaLineEdit->setText("0.2");
  this->m_FilterGridLayout->addWidget(this->m_SigmaLineEdit, 0, 0);

  this->m_ApplyButton = new QPushButton();
  this->m_ApplyButton->setText("Apply");
  this->m_FilterGridLayout->addWidget(this->m_ApplyButton, 0, 1);

  this->m_Buttons = new QWidget(this);
  this->m_ButtonsLayout = new QHBoxLayout(this->m_Buttons);
  this->m_ButtonsLayout->insertSpacing(0,250);
  this->m_ButtonsLayout->addWidget(this->ButtonHelp);
  this->m_ButtonsLayout->addWidget(this->ButtonOk);

  this->gridLayout->addWidget(this->m_Buttons, 2, 0,1,2);

  this->gridLayout->addWidget(m_TabWidget, 1, 0,1,2);
}

void QtControlView::setInputImage(ImageType * newImData)
{
  this->m_ImageData = newImData;
  this->OpenGlWindow->setInputImage(newImData);
  this->SliceNumSlider->setMaximum( static_cast<int>
                                    (this->OpenGlWindow->maxSliceNum() -1));
  this->OpenGlWindow->changeSlice(((this->OpenGlWindow->maxSliceNum() -1)/2));
  this->SliceNumSlider->setValue(static_cast<int>
                                   ((this->OpenGlWindow->maxSliceNum() -1)/2));
  this->setDisplaySliceNumber(static_cast<int>
                                ((this->OpenGlWindow->maxSliceNum() -1)/2));
  this->IntensityMin->setMinimum( static_cast<int>
                                  ( this->OpenGlWindow->minIntensity() ));
  this->IntensityMin->setMaximum( static_cast<int>
                                  ( this->OpenGlWindow->maxIntensity() ));
  this->IntensityMin->setValue( static_cast<int>
                                ( this->OpenGlWindow->minIntensity() ));
  this->IntensityMax->setMinimum( static_cast<int>
                                  ( this->OpenGlWindow->minIntensity() ));
  this->IntensityMax->setMaximum( static_cast<int>
                                  ( this->OpenGlWindow->maxIntensity() ));
  this->IntensityMax->setValue( static_cast<int>
                                ( this->OpenGlWindow->maxIntensity() ));

  char* tempchar = new char[20];
  sprintf(tempchar,"%.0f",this->OpenGlWindow->minIntensity());
  this->IntensityMinDisplay->setText(tempchar);
  sprintf(tempchar,"%.0f",this->OpenGlWindow->maxIntensity());
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;

  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}

void QtControlView::setDisplayPosition(int x,int y ,int z,double value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",x);
  this->PositionX->setText(tr(tempchar));
  sprintf(tempchar,"%d",y);
  this->PositionY->setText(tr(tempchar));
  sprintf(tempchar,"%d",z);
  this->PositionZ->setText(tr(tempchar));
  sprintf(tempchar,"%3.1f",value);
  this->PixelValue->setText(tr(tempchar));
  this->OpenGlWindow->update();
  delete tempchar;
}

void QtControlView::setDisplaySliceNumber(int number)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",number);
  number = this->SliceNumSlider->value();
  this->SliceValue->setText(tempchar);
  delete tempchar;
}

void QtControlView::setDisplayIMin(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMinDisplay->setText(tempchar);
  delete tempchar;
}

void QtControlView::setDisplayIMax(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;
}

void QtControlView::setDisplaySigma(QString value)
{
  this->m_SigmaLineEdit->setText(value);
}

void QtControlView::setFilter()
{
  if(this->m_ImageData == NULL)
    {
    return;
    }
  FilterType::Pointer FilterX;
  FilterType::Pointer FilterY;

  FilterX = FilterType::New();
  FilterY = FilterType::New();

  FilterX->SetDirection( 0 );   // 0 --> X direction
  FilterY->SetDirection( 1 );   // 1 --> Y direction

  FilterX->SetOrder( FilterType::ZeroOrder );
  FilterY->SetOrder( FilterType::ZeroOrder );

  FilterX->SetNormalizeAcrossScale( true );
  FilterY->SetNormalizeAcrossScale( true );

  FilterX->SetInput( this->m_ImageData );
  FilterY->SetInput( FilterX->GetOutput() );

  const double sigma = m_SigmaLineEdit->text().toDouble();

  FilterX->SetSigma( sigma );
  FilterY->SetSigma( sigma );

  FilterX->Update();
  FilterY->Update();

  this->OpenGlWindow->update();
  setInputImage(FilterY->GetOutput());
  show();
}

}
