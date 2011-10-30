/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#include "tubeImageViewer.h"
#include <iostream>
#include "QtGlSliceView.h"
#include <qlineedit.h>
#include <qslider.h>


/**
 *
 */
tubeImageViewer::tubeImageViewer(
  QWidget* _parent,  const char* itkNotUsed( _name ),
  bool itkNotUsed( _modal ), Qt::WFlags itkNotUsed( _fl ) )
:QDialog( _parent )
{
  setupUi( this );

  QObject::connect( this->Window2D, SIGNAL( PixelValueChanged( double ) ),
                    this, SLOT( this->DisplayClick( double ) ) );
  QObject::connect( this->HelpButton, SIGNAL( clicked() ),
                    this, SLOT( this->DisplayHelp() ) );
}

/**
 *  Destroys the object and frees any allocated resources
 */
tubeImageViewer::~tubeImageViewer()
{
}

void tubeImageViewer::DisplayClick( double v )
{
  std::cout << "Here" << std::endl;

  std::cout << "  x = " << this->Window2D->GetClickSelectX() << std::endl;
}

void tubeImageViewer::DisplayHelp()
{
  Ui::HelpWindow * helpWindow = new Ui::HelpWindow();
}

void tubeImageViewer::SetInputImage( ImageType * newImData )
{
  this->Window2D->SetInputImage( newImData );
  this->Window2D->SetFlipY( true );

  int maxSlice = newImData->GetLargestPossibleRegion().GetSize()[2]-1;
  this->Slice->setMaximum( maxSlice );
  this->SliceNum->setMaximum( maxSlice );

  this->Slice->setValue( 0 );

  this->IntensityMin->setMinimum( static_cast<int>(
    this->Window2D->GetIntensityMin() ) );
  this->IntensityMin->setMaximum( static_cast<int>(
    this->Window2D->GetIntensityMax() ) );
  this->IntensityMinDisplay->setMinimum( this->IntensityMin->minimum() );
  this->IntensityMinDisplay->setMaximum( this->IntensityMin->maximum() );
  this->IntensityMin->setValue( static_cast<int>(
    this->Window2D->GetIntensityMin() ) );

  this->IntensityMax->setMinimum( static_cast<int>(
    this->Window2D->GetIntensityMin() ) );
  this->IntensityMax->setMaximum( static_cast<int>(
    this->Window2D->GetIntensityMax() ) );
  this->IntensityMaxDisplay->setMinimum( this->IntensityMax->minimum() );
  this->IntensityMaxDisplay->setMaximum( this->IntensityMax->maximum() );
  this->IntensityMax->setValue( static_cast<int>(
    this->Window2D->GetIntensityMax() ) );

  this->Window2D->show();
  this->Window2D->update();
}

void tubeImageViewer::SetInputOverlay( OverlayType * newOverlayData )
{
  this->Window2D->SetInputOverlay( newOverlayData );

  this->Window2D->show();
  this->Window2D->update();
}
