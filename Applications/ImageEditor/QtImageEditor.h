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

#ifndef __QtImageEditor_h
#define __QtImageEditor_h

// QtImageViewer includes
#include "QtImageViewer.h"

namespace tube
{
// TubeTK includes
class QtImageEditorPrivate;

class QtImageEditor: public QtImageViewer
{
  Q_OBJECT
public:
  typedef QtImageViewer Superclass;

  QtImageEditor( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
  virtual ~QtImageEditor();

public:
  virtual void setInputImage(ImageType *newImData);
  virtual void setOverlayImage(OverlayImageType* newImData);

public slots:
  void setDisplaySigma(double value);
  void applyFFT();
  void applyInverseFFT();
  void applyFilter();
  void displayFFT();
  void blurFilter();
  void useNewFilter();

protected:
  QScopedPointer<QtImageEditorPrivate> d_ptr;
  virtual void setControlsVisible(bool visible);

private:
  Q_DECLARE_PRIVATE(QtImageEditor);
  Q_DISABLE_COPY(QtImageEditor);
};

} // End namespace tube

#endif
