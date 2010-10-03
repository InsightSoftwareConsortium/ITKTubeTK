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
#include "QtSlicer.h"
#include <iostream>
#include "QtGlSliceView.h"
#include <qlineedit.h>
#include <qslider.h>


/**
 *
 */
QtSlicer::QtSlicer( QWidget* parent,  const char* name, bool modal, Qt::WFlags fl )
:QDialog(parent)
{
    setupUi(this);
    this->IntensityMinLabel->setMinimumSize( this->IntensityMaxLabel->sizeHint() );
    this->SliceNum->setMaximumWidth(40);
}

/**  
 *  Destroys the object and frees any allocated resources
 */
QtSlicer::~QtSlicer()
{
}

void QtSlicer::DisplayPosition(int x,int y ,int z,float value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",x);
  PositionX->setText(QString(tempchar));
  sprintf(tempchar,"%d",y);
  PositionY->setText(QString(tempchar));
  sprintf(tempchar,"%d",z);
  PositionZ->setText(QString(tempchar));
  sprintf(tempchar,"%3.1f",value);
  PixelValue->setText(QString(tempchar));
  delete tempchar;
}

void QtSlicer::Help()
{
  Ui::HelpWindow * helpWindow = new Ui::HelpWindow();

  // FIXME: Used in Qt3, not needed in Qt4:
  // helpWindow->show();
}

void QtSlicer::SetInputImage(ImageType * newImData)
{
  this->OpenGlWindow->SetInputImage(newImData);
  this->OpenGlWindow->SetFlipY(true);
  
  int maxSlice = newImData->GetLargestPossibleRegion().GetSize()[2]-1;
  this->Slice->setMaximum(maxSlice);
  this->SliceNum->setMaximum(maxSlice);
  
  this->Slice->setValue(0);

  this->IntensityMin->setMinimum( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMin->setMaximum( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMinDisplay->setMinimum( this->IntensityMin->minimum() );
  this->IntensityMinDisplay->setMaximum( this->IntensityMin->maximum() );
  this->IntensityMin->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  
  this->IntensityMax->setMinimum( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMax->setMaximum( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMaxDisplay->setMinimum( this->IntensityMax->minimum() );
  this->IntensityMaxDisplay->setMaximum( this->IntensityMax->maximum() );
  this->IntensityMax->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  
  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}
