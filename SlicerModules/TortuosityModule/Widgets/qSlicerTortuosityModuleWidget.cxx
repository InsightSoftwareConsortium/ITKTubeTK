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

// Qt includes
#include "qSlicerTortuosityModuleWidget.h"
#include "ui_qSlicerTortuosityModule.h"

// MRML includes


//------------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Tortuosity
class qSlicerTortuosityModuleWidgetPrivate :
 public Ui_qSlicerTortuosityModule
{
  Q_DECLARE_PUBLIC(qSlicerTortuosityModuleWidget);

protected:
  qSlicerTortuosityModuleWidget* const q_ptr;

public:
  qSlicerTortuosityModuleWidgetPrivate(
    qSlicerTortuosityModuleWidget& object);
  void init();
};

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidgetPrivate
::qSlicerTortuosityModuleWidgetPrivate(
  qSlicerTortuosityModuleWidget& object)
  : q_ptr(&object)
{
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidgetPrivate::init()
{
  Q_Q(qSlicerTortuosityModuleWidget);

  this->setupUi(q);
}

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidget::
qSlicerTortuosityModuleWidget(QWidget* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerTortuosityModuleWidgetPrivate(*this))
{
}

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidget::~qSlicerTortuosityModuleWidget()
{}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::setup()
{
  Q_D(qSlicerTortuosityModuleWidget);
  d->init();
}
