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
#ifndef tubeImageViewer_h
#define tubeImageViewer_h

#include <QObject>

#include "ui_tubeImageViewerGUI.h"
#include "ui_tubeImageViewerHelpGUI.h"

#include "itkImage.h"

class tubeImageViewer : public QDialog, Ui::Gui
{
public:

  typedef itk::Image< float, 3 >         ImageType;
  typedef itk::Image< unsigned char, 3 > OverlayType;

  tubeImageViewer( QWidget* parent = 0, const char* name = 0,
    bool modal = FALSE, Qt::WFlags fl = 0 );

  ~tubeImageViewer();

public slots:
  void DisplayClick( double v );
  void DisplayHelp( void );

public:

  void SetInputImage(ImageType * newImData);
  void SetInputOverlay(OverlayType * newOverlayData);
};

#endif
