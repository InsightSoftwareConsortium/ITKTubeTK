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
  
  int maxSlice = newImData->GetLargestPossibleRegion().GetSize()[2]-1;
  this->Slice->setMaximum(maxSlice);
  this->SliceNum->setMaximum(maxSlice);
  
  this->Slice->setValue(newImData->GetLargestPossibleRegion().GetSize()[2]/2);

  this->IntensityMin->setMinimum( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMin->setMaximum( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMin->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMax->setMinimum( static_cast<int>( this->OpenGlWindow->GetIntensityMin() ));
  this->IntensityMax->setMaximum( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  this->IntensityMax->setValue( static_cast<int>( this->OpenGlWindow->GetIntensityMax() ));
  
  char* tempchar = new char[20];
  sprintf(tempchar,"%.0f",this->OpenGlWindow->GetIntensityMin());
  this->IntensityMinDisplay->setText(tempchar);
  sprintf(tempchar,"%.0f",this->OpenGlWindow->GetIntensityMax());
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;

  this->OpenGlWindow->show();
  this->OpenGlWindow->update();
}

void QtSlicer::DisplayIMin(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMinDisplay->setText(tempchar);
  delete tempchar;
}

void QtSlicer::DisplayIMax(int value)
{
  char* tempchar = new char[20];
  sprintf(tempchar,"%d",value);
  this->IntensityMaxDisplay->setText(tempchar);
  delete tempchar;
}
