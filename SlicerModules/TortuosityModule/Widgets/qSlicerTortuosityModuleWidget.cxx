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

// TubeTK includes
#include "tubeTubeMath.h"

// MRML includes
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkSlicerTortuosityLogic.h"

//------------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_Tortuosity
class qSlicerTortuosityModuleWidgetPrivate :
 public Ui_qSlicerTortuosityModule
{
  Q_DECLARE_PUBLIC( qSlicerTortuosityModuleWidget );

public:
  qSlicerTortuosityModuleWidgetPrivate(
    qSlicerTortuosityModuleWidget& object );
  void init();
  vtkSlicerTortuosityLogic* logic() const;
  int getFlagFromCheckBoxes();

  vtkMRMLSpatialObjectsNode* currentSpatialObject;

protected:
  qSlicerTortuosityModuleWidget* const q_ptr;
};

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidgetPrivate
::qSlicerTortuosityModuleWidgetPrivate(
  qSlicerTortuosityModuleWidget& object )
  : q_ptr( &object )
{
  this->currentSpatialObject = 0;
}

//------------------------------------------------------------------------------
int qSlicerTortuosityModuleWidgetPrivate::getFlagFromCheckBoxes()
{
  int flag = this->BasicMetricsCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::BasicMetricsGroup : 0;

  flag |= this->OldMetricsCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::OldMetricsGroup : 0;

  flag |= this->CurvatureMetricsCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::CurvatureMetricsGroup : 0;

  flag |= this->HistogramMetricsCheckBox->isChecked() ?
    vtkSlicerTortuosityLogic::HistogramMetricsGroup : 0;

  return flag;
}


//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidgetPrivate::init()
{
  Q_Q(qSlicerTortuosityModuleWidget);

  this->setupUi( q );

  QObject::connect(
    this->CurrentSpatialObjectComboBox, SIGNAL( currentNodeChanged( vtkMRMLNode* ) ),
    q, SLOT( setCurrentSpatialObjectsNode( vtkMRMLNode* ) ) );

  QObject::connect(
    this->RunPushButton, SIGNAL( toggled( bool ) ),
    q, SLOT( runSelectedMetrics( bool ) ) );

  QObject::connect(
    this->SaveCSVPushButton, SIGNAL( toggled( bool ) ),
    q, SLOT( saveCurrentSpatialObjectAsCSV( bool ) ) );

  QObject::connect(
    this->LoadColorsFromCSVPushButton, SIGNAL( clicked() ),
    q, SLOT( loadColorsFromCSV() ) );

  QObject::connect(
    this->SmoothingMethodComboBox, SIGNAL( currentIndexChanged( int ) ),
    q, SLOT( smoothingMethodChanged( int ) ) );

  this->SmoothingMethodComboBox->addItem( "Average on Index",
                                         tube::SMOOTH_TUBE_USING_INDEX_AVERAGE );
  this->SmoothingMethodComboBox->addItem( "Gaussian on Index",
                                         tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN );
  this->SmoothingMethodComboBox->addItem( "Gaussian on Distance",
                                         tube::SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN );
  this->SmoothingMethodComboBox->setCurrentIndex(
    this->SmoothingMethodComboBox->findData( tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN ) );


  this->BasicMetricsLabel->setToolTip( "AverageRadiusMetric, "
                                         "ChordLengthMetric, "
                                         "PathLengthMetric" );
  this->OldMetricsLabel->setToolTip( "DistanceMetric, "
                                       "InflectionCountMetric, "
                                       "InflectionPointsMetric, "
                                       "SumOfAnglesMetric,"
                                       "SumOfTorsionMetric");
  this->CurvatureMetricsLabel->setToolTip( "CurvatureScalarMetric, "
                                             "InflectionCount1Metric, "
                                             "InflectionCount2Metric, "
                                             "Percentile95Metric, "
                                             "Tau4Metric, "
                                             "TotalCurvatureMetric, "
                                             "TotalSquared, "
                                             "CurvatureMetric" );
  this->HistogramMetricsLabel->setToolTip( "CurvatureHistogramMetrics" );
}

//-----------------------------------------------------------------------------
vtkSlicerTortuosityLogic* qSlicerTortuosityModuleWidgetPrivate::logic() const
{
  Q_Q( const qSlicerTortuosityModuleWidget );
  return vtkSlicerTortuosityLogic::SafeDownCast( q->logic() );
}

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidget::
qSlicerTortuosityModuleWidget( QWidget* _parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerTortuosityModuleWidgetPrivate( *this ) )
{
}

//------------------------------------------------------------------------------
qSlicerTortuosityModuleWidget::~qSlicerTortuosityModuleWidget()
{
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::setup()
{
  Q_D( qSlicerTortuosityModuleWidget );
  d->init();
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget
::setCurrentSpatialObjectsNode( vtkMRMLNode* node )
{
  this->setCurrentSpatialObjectsNode(
    vtkMRMLSpatialObjectsNode::SafeDownCast( node ) );
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget
::setCurrentSpatialObjectsNode( vtkMRMLSpatialObjectsNode* node )
{
  Q_D( qSlicerTortuosityModuleWidget );
  if ( d->currentSpatialObject == node )
    {
    return;
    }

  d->currentSpatialObject = node;

  d->RunPushButton->setEnabled( node != 0 );
  d->SaveCSVPushButton->setEnabled( node != 0 );
  d->LoadColorsFromCSVPushButton->setEnabled( node != 0 );
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::runMetrics( int flag )
{
  Q_D( qSlicerTortuosityModuleWidget );

  // Smoothing Scale
  double smoothingScale = d->SmoothingScaleSliderWidget->value();

  // Smoothing Method
  tube::SmoothTubeFunctionEnum smoothMethod =
    static_cast<tube::SmoothTubeFunctionEnum>( d->SmoothingMethodComboBox->itemData(
    d->SmoothingMethodComboBox->currentIndex() ).toInt() );

  //Subsampling Factor
  int subsampling = d->SubsamplingSliderWidget->value();

  d->RunPushButton->setEnabled( false );
  if ( !d->logic()->RunMetrics( d->currentSpatialObject, flag,
                              smoothMethod, smoothingScale, subsampling ) )
    {
    qCritical( "Error while running metrics !" );
    }
  d->RunPushButton->setChecked( false );
  d->RunPushButton->setEnabled( true );
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::loadColorsFromCSV()
{
  Q_D( qSlicerTortuosityModuleWidget );

  if ( !d->currentSpatialObject )
    {
    return;
    }

  d->LoadColorsFromCSVPushButton->setEnabled( false );

  QString filename = QFileDialog::getOpenFileName(
    this, "Load colors as csv...", "", "Comma Separated Value (*.csv)" );
  if ( !d->logic()->LoadColorsFromCSV( d->currentSpatialObject, filename.toLatin1() ) )
    {
    qCritical( "Error while loading colors !" );
    }

  d->LoadColorsFromCSVPushButton->setEnabled( true );
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::runSelectedMetrics( bool run )
{
  Q_D( qSlicerTortuosityModuleWidget );
  if ( !run )
    {
    return;
    }
  this->runMetrics( d->getFlagFromCheckBoxes() );
}

//------------------------------------------------------------------------------
void qSlicerTortuosityModuleWidget::saveCurrentSpatialObjectAsCSV( bool save )
{
  Q_D( qSlicerTortuosityModuleWidget );

  if ( !d->currentSpatialObject || !save )
    {
    return;
    }

  d->SaveCSVPushButton->setEnabled( false );

  QString filename =
    QFileDialog::getSaveFileName(
      this, "Save tortuosity as csv...", "", "Comma Separated Value (*.csv)" );

  int flag = d->getFlagFromCheckBoxes();

  if ( !d->logic()->SaveAsCSV( d->currentSpatialObject, filename.toLatin1(), flag ) )
    {
    qCritical( "Failed to write CSV at '%s'. Make sure you have run the "
              "tortuosity metrics at least once. ", qPrintable( filename ) );
    }

  d->SaveCSVPushButton->setChecked( false );
  d->SaveCSVPushButton->setEnabled( true );
}

void qSlicerTortuosityModuleWidget::smoothingMethodChanged( int index )
{
  Q_D( qSlicerTortuosityModuleWidget );

  switch ( index )
    {
    case tube::SMOOTH_TUBE_USING_INDEX_AVERAGE:
      d->SmoothingScaleSliderWidget->setDecimals( 0 );
      d->SmoothingScaleSliderWidget->setSingleStep( 1 );
      d->SmoothingScaleSliderWidget->setMinimum( 0 );
      d->SmoothingScaleSliderWidget->setMaximum( 100 );
      d->SmoothingScaleSliderWidget->setValue( 3 );
      d->SmoothingScaleSliderWidget->setToolTip( "Half the window size" );
      d->SmoothingScaleLabel->setToolTip( "Half the window size" );
      break;
    case tube::SMOOTH_TUBE_USING_INDEX_GAUSSIAN:
      d->SmoothingScaleSliderWidget->setDecimals( 2 );
      d->SmoothingScaleSliderWidget->setSingleStep( 0.1 );
      d->SmoothingScaleSliderWidget->setMinimum( 0 );
      d->SmoothingScaleSliderWidget->setMaximum( 50 );
      d->SmoothingScaleSliderWidget->setValue( 5 );
      d->SmoothingScaleSliderWidget->setToolTip( "Standard deviation" );
      d->SmoothingScaleLabel->setToolTip( "Standard deviation" );
      break;
    case tube::SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN:
      d->SmoothingScaleSliderWidget->setDecimals( 2 );
      d->SmoothingScaleSliderWidget->setSingleStep( 0.01 );
      d->SmoothingScaleSliderWidget->setMinimum( 0 );
      d->SmoothingScaleSliderWidget->setMaximum( 50 );
      d->SmoothingScaleSliderWidget->setValue( 0.5 );
      d->SmoothingScaleSliderWidget->setToolTip( "Standard deviation" );
      d->SmoothingScaleLabel->setToolTip( "Standard deviation" );
      break;
    default:
      d->SmoothingScaleSliderWidget->setDecimals( 0 );
      d->SmoothingScaleSliderWidget->setSingleStep( 0 );
      d->SmoothingScaleSliderWidget->setMinimum( 0 );
      d->SmoothingScaleSliderWidget->setMaximum( 0 );
      d->SmoothingScaleSliderWidget->setValue( 0 );
      d->SmoothingScaleSliderWidget->setToolTip( "Unknown Smoothing Method" );
      d->SmoothingScaleLabel->setToolTip( "Unknown Smoothing Method" );
      break;
    }
}
