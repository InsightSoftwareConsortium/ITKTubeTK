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

#ifndef __QtOverlayControlsWidget_h
#define __QtOverlayControlsWidget_h

//Qt includes
#include <QWidget>

//itk includes
#include "itkImage.h"
#include "itkImageFileReader.h"

//ImageViewer includes
#include "QtGlSliceView.h"

//ImageEditor includes
#include "ui_QtOverlayControlsWidgetGUI.h"


namespace tube
{

class QtOverlayControlsWidget : public QWidget
{
  Q_OBJECT

  Q_PROPERTY(double opacity READ getOpacity WRITE setOpacity NOTIFY opacityChanged);
public:
  QtOverlayControlsWidget(QWidget* parent);
  ~QtOverlayControlsWidget();
  typedef itk::Image< double, 3 >             ImageType;

  typedef unsigned char                       OverlayPixelType;
  typedef itk::Image<OverlayPixelType,3>      OverlayType;
  typedef itk::ImageFileReader<OverlayType>   OverlayReaderType;

  void setSliceView(QtGlSliceView *sliceView);
  double getOpacity() const;

public slots:
  /// Set the new opacity changed in the SliceView to the Slider
  void setOpacity(double value);
  bool loadOverlay(QString pathOverlay = QString());
  void setInputOverlay(OverlayType* overlayImage);
  /// Load Overlay if the check box is checked by a user, Disable overlay if
  /// the check box is unchecked by the user.
  void setOverlayVisibility(bool show);

protected slots:
  void onOpacityChanged(int newOpacity);

signals:
  void opacityChanged(double getOpacity);

private:
  QtGlSliceView  *m_SliceView;
  Ui::Overlay    *m_UI;

  Q_DISABLE_COPY(QtOverlayControlsWidget)
};
}

#endif
