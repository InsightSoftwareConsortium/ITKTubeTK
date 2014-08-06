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

// Module includes
#include "qSlicerTortuosityModuleWidget.h"
#include "ui_qSlicerTortuosityModule.h"

// Qt includes
#include <QDebug>
#include <QFileDialog>

// MRML includes
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkSlicerTortuosityLogic.h"

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
  vtkSlicerTortuosityLogic* logic() const;

  vtkMRMLSpatialObjectsNode* currentSpatialObject;
};

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidgetPrivate
::qSlicerTortuosityModuleWidgetPrivate(
  qSlicerTortuosityModuleWidget& object)
  : q_ptr(&object)
{
  this->currentSpatialObject = 0;
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidgetPrivate::init()
{
  Q_Q(qSlicerTortuosityModuleWidget);

  this->setupUi(q);

  QObject::connect(
    this->CurrentSpatialObjectComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
    q, SLOT(setCurrentSpatialObjectsNode(vtkMRMLNode*)));

  QObject::connect(
    this->RunPushButton, SIGNAL(toggled(bool)),
    q, SLOT(runSelectedMetrics(bool)));

  QObject::connect(
    this->SaveCSVPushButton, SIGNAL(toggled(bool)),
    q, SLOT(saveCurrentSpatialObjectAsCSV(bool)));
}

//-----------------------------------------------------------------------------
vtkSlicerTortuosityLogic* qSlicerTortuosityModuleWidgetPrivate::logic() const
{
  Q_Q(const qSlicerTortuosityModuleWidget);
  return vtkSlicerTortuosityLogic::SafeDownCast(q->logic());
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
{
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::setup()
{
  Q_D(qSlicerTortuosityModuleWidget);
  d->init();
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget
::setCurrentSpatialObjectsNode(vtkMRMLNode* node)
{
  this->setCurrentSpatialObjectsNode(
    vtkMRMLSpatialObjectsNode::SafeDownCast(node));
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget
::setCurrentSpatialObjectsNode(vtkMRMLSpatialObjectsNode* node)
{
  Q_D(qSlicerTortuosityModuleWidget);
  if (d->currentSpatialObject == node)
    {
    return;
    }

  d->currentSpatialObject = node;

  d->RunPushButton->setEnabled(node != 0);
  d->SaveCSVPushButton->setEnabled(node != 0);
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::runMetrics(int flag)
{
  Q_D(qSlicerTortuosityModuleWidget);
  d->RunPushButton->setEnabled(false);
  if (!d->logic()->RunMetrics(d->currentSpatialObject, flag))
    {
    qCritical("Error while running metrics !");
    }
  d->RunPushButton->setChecked(false);
  d->RunPushButton->setEnabled(true);
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::runSelectedMetrics(bool run)
{
  Q_D(qSlicerTortuosityModuleWidget);
  if (!run)
    {
    return;
    }

  int flag = d->DistanceMetricCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::DistanceMetric : 0;

  flag |= d->InflectionCountMetricCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::InflectionCountMetric : 0;
  
  flag |= d->InflectionPointsCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::InflectionPoints : 0;

  flag |= d->SumOfAnglesMetricCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::SumOfAnglesMetric : 0;

  this->runMetrics(flag);
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::saveCurrentSpatialObjectAsCSV(bool save)
{
  Q_D(qSlicerTortuosityModuleWidget);

  if (!d->currentSpatialObject || !save)
    {
    return;
    }

  d->SaveCSVPushButton->setEnabled(false);

  QString filename =
    QFileDialog::getSaveFileName(
      this, "Save tortuosity as csv...", "", "Comma Separated Value (*.csv)");
  int flag = vtkSlicerTortuosityLogic::DistanceMetric
    | vtkSlicerTortuosityLogic::InflectionCountMetric
    | vtkSlicerTortuosityLogic::SumOfAnglesMetric;
  if (!d->logic()->SaveAsCSV(d->currentSpatialObject, filename.toLatin1(), flag))
    {
    QString msg = "Failed to write CSV at %1. Make sure you have run the "
      "tortuosity metrics at least once.";
    qCritical(msg.arg(filename).toLatin1());
    }

  d->SaveCSVPushButton->setChecked(false);
  d->SaveCSVPushButton->setEnabled(true);
}
