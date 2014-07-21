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

// ImageViewer includes
#include "QtGlSliceView.h"


namespace tube
{
// ImageEditor includes
class QtOverlayControlsWidgetPrivate;

class QtOverlayControlsWidget: public QWidget
{
  Q_OBJECT

  Q_PROPERTY(double opacity READ opacity WRITE setOpacity NOTIFY opacityChanged);
public:
  QtOverlayControlsWidget(QWidget* parent);
  virtual ~QtOverlayControlsWidget();
  typedef itk::Image< double, 3 >             ImageType;

  typedef unsigned char                       OverlayPixelType;
  typedef itk::Image<OverlayPixelType,3>      OverlayType;
  typedef itk::ImageFileReader<OverlayType>   OverlayReaderType;

  void setSliceView(QtGlSliceView *sliceView);
  /// Return the opacity property value.
  /// \sa opacity
  double opacity() const;

public slots:
  /// Set the new opacity changed in the SliceView to the Slider
  void setOpacity(double value);

  bool loadOverlay(QString pathOverlay = QString());
  void setInputOverlay(OverlayType* overlayImage);

protected slots:
  void onOpacityChanged(int newOpacity);

signals:
  void opacityChanged(double opacity);

protected:
  QScopedPointer<tube::QtOverlayControlsWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(QtOverlayControlsWidget);

  Q_DISABLE_COPY(QtOverlayControlsWidget)
};

}

#endif
