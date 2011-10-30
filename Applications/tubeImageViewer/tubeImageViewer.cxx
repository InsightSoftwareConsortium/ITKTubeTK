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
tubeImageViewer::
tubeImageViewer( QWidget* _parent )
:QWidget( _parent )
{
  setupUi( this );
}

/**
 *  Destroys the object and frees any allocated resources
 */
tubeImageViewer::
~tubeImageViewer()
{
}

void
tubeImageViewer::
on_Window2D_PixelValueChanged( double v )
{
  std::cout << "Pixel" << std::endl;
}

void
tubeImageViewer::
on_HelpButton_clicked()
{
  std::cout << "Help" << std::endl;
  //Ui::HelpWindow * helpWindow = new Ui::HelpWindow();
}

void
tubeImageViewer::
on_QuitButton_clicked()
{
  this->close();
}

void
tubeImageViewer::
SetInputImage( ImageType * newImData )
{
  m_ImageData = newImData;

  this->Window2D->SetInputImage( newImData );
  this->Window2D->SetFlipY( true );

  ImageType::SizeType dimSize;
  dimSize = newImData->GetLargestPossibleRegion().GetSize();

  ImageType::SpacingType spacing;
  spacing = newImData->GetSpacing();

  ImageType::IndexType index0;
  index0 = newImData->GetLargestPossibleRegion().GetIndex();

  ImageType::SpacingType origin;
  for( unsigned int i=0; i<3; i++ )
    {
    origin[i] = newImData->GetOrigin()[i];
    }

  int maxSlice = dimSize[2]-1;
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

  QString infoText;
  QString tmpText;
  m_ImageInfoText.clear();
  infoText = "NDims = 3";
  m_ImageInfoText.push_back( QListWidgetItem( infoText ) );

  infoText = "DimSize = ";
  tmpText.setNum( dimSize[0] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( dimSize[1] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( dimSize[2] );
  infoText.push_back( tmpText );
  m_ImageInfoText.push_back( QListWidgetItem( infoText ) );

  infoText = "Index = ";
  tmpText.setNum( index0[0] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( index0[1] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( index0[2] );
  infoText.push_back( tmpText );
  m_ImageInfoText.push_back( QListWidgetItem( infoText ) );

  infoText = "Origin = ";
  tmpText.setNum( origin[0] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( origin[1] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( origin[2] );
  infoText.push_back( tmpText );
  m_ImageInfoText.push_back( QListWidgetItem( infoText ) );

  infoText = "Spacing = ";
  tmpText.setNum( spacing[0] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( spacing[1] );
  infoText.push_back( tmpText );
  infoText.push_back( ", " );
  tmpText.setNum( spacing[2] );
  infoText.push_back( tmpText );
  m_ImageInfoText.push_back( QListWidgetItem( infoText ) );

  for( unsigned int i=0; i< m_ImageInfoText.size(); i++ )
    {
    this->ImageInfoTextList->addItem( & m_ImageInfoText[i] );
    }

  this->Window2D->show();
  this->Window2D->update();
}

void
tubeImageViewer::
SetInputOverlay( OverlayType * newOverlayData )
{
  m_OverlayData = newOverlayData;

  this->Window2D->SetInputOverlay( newOverlayData );

  this->Window2D->show();
  this->Window2D->update();
}
