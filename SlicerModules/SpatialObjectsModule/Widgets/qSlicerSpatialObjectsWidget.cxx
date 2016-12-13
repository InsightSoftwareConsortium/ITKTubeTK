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

// qMRML includes
#include "qSlicerSpatialObjectsWidget.h"
#include "ui_qSlicerSpatialObjectsWidget.h"

// CTK
#include <ctkDoubleRangeSlider.h>

// MRML includes
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLSpatialObjectsDisplayNode.h>
#include <vtkMRMLSpatialObjectsLineDisplayNode.h>
#include <vtkMRMLSpatialObjectsDisplayPropertiesNode.h>

// VTK includes
#include <vtkDataArray.h>
#include <vtkNew.h>

//------------------------------------------------------------------------------
class qSlicerSpatialObjectsWidgetPrivate : public Ui_qSlicerSpatialObjectsWidget
{
  Q_DECLARE_PUBLIC( qSlicerSpatialObjectsWidget );

protected:
  qSlicerSpatialObjectsWidget* const q_ptr;

public:
  qSlicerSpatialObjectsWidgetPrivate( qSlicerSpatialObjectsWidget& object );
  void init();
  bool centeredOrigin( double* origin )const;

  vtkMRMLSpatialObjectsNode* SpatialObjectsNode;
  vtkMRMLSpatialObjectsDisplayNode* SpatialObjectsDisplayNode;
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    SpatialObjectsDisplayPropertiesNode;
};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsWidgetPrivate::qSlicerSpatialObjectsWidgetPrivate
                                      ( qSlicerSpatialObjectsWidget& object )
  : q_ptr( &object )
{
  this->SpatialObjectsDisplayNode = 0;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidgetPrivate::init()
{
  Q_Q( qSlicerSpatialObjectsWidget );
  this->setupUi( q );

  this->ColorBySolidColorPicker->setDialogOptions(
    ctkColorPickerButton::UseCTKColorDialog );

  QObject::connect( this->VisibilityCheckBox,
                   SIGNAL( clicked( bool ) ), q,
                   SLOT( setVisibility( bool ) ) );
  QObject::connect( this->ColorBySolidColorCheckBox,
                   SIGNAL( clicked() ), q,
                   SLOT( setColorBySolid() ) );

  QObject::connect( this->SliceIntersectionCheckBox,
                   SIGNAL( stateChanged( int ) ), q,
                   SLOT( setSliceIntersection( int ) ) );
  QObject::connect( this->SliceIntersectionThicknessSpinBox,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setSliceIntersectionThickness( double ) ) );

  QObject::connect( this->ColorBySolidColorPicker,
                   SIGNAL( colorChanged( QColor ) ), q,
                   SLOT( onColorBySolidChanged( QColor ) ) );
  QObject::connect( this->ColorByScalarsColorTableComboBox,
                   SIGNAL( currentNodeChanged( vtkMRMLNode* ) ), q,
                   SLOT( setColorByCellScalarsColorTable( vtkMRMLNode* ) ) );

  QObject::connect( this->ColorByScalarRadioButton,
                   SIGNAL( clicked() ), q,
                   SLOT( setColorByScalar() ) );
  QObject::connect( this->ColorByScalarComboBox,
                   SIGNAL( currentIndexChanged( int ) ), q,
                   SLOT( onColorByScalarChanged( int ) ) );

  QObject::connect( this->ScalarRangeWidget,
                   SIGNAL( valuesChanged( double,double ) ), q,
                   SLOT( onColorByScalarValuesChanged( double,double ) ) );

  QObject::connect( this->OpacitySlider,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setOpacity( double ) ) );

  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( colorChanged( QColor ) ), q,
                   SLOT( setColor( QColor ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( opacityChanged( double ) ), q,
                   SLOT( setOpacity( double ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( ambientChanged( double ) ), q,
                   SLOT( setAmbient( double ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( diffuseChanged( double ) ), q,
                   SLOT( setDiffuse( double ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( specularChanged( double ) ), q,
                   SLOT( setSpecular( double ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( specularPowerChanged( double ) ), q,
                   SLOT( setSpecularPower( double ) ) );
  QObject::connect( this->MaterialPropertyWidget,
                   SIGNAL( backfaceCullingChanged( bool ) ), q,
                   SLOT( setBackfaceCulling( bool ) ) );

  this->MaterialPropertyWidget->setHidden( true );
  this->MaterialPropertyGroupBox->setHidden( true );

  this->ScalarRangeWidget->setRange( 0., 0. );
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsWidget::qSlicerSpatialObjectsWidget( QWidget *_parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerSpatialObjectsWidgetPrivate( *this ) )
{
  Q_D( qSlicerSpatialObjectsWidget );

  d->init();
  this->inProcess = 0;
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsWidget::~qSlicerSpatialObjectsWidget()
{}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode* qSlicerSpatialObjectsWidget::
SpatialObjectsNode() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->SpatialObjectsNode;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* qSlicerSpatialObjectsWidget::
SpatialObjectsDisplayNode() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->SpatialObjectsDisplayNode;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayPropertiesNode* qSlicerSpatialObjectsWidget::
SpatialObjectsDisplayPropertiesNode() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->SpatialObjectsDisplayNode ?
    d->SpatialObjectsDisplayNode->GetSpatialObjectsDisplayPropertiesNode() : 0;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setSpatialObjectsNode( vtkMRMLNode* node, int DisplayNodeIndex )
{
  this->setSpatialObjectsNode( vtkMRMLSpatialObjectsNode::SafeDownCast( node ), DisplayNodeIndex );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::
setSpatialObjectsNode( vtkMRMLSpatialObjectsNode* SpatialObjectsNode, int DisplayNodeIndex )
{
  Q_D( qSlicerSpatialObjectsWidget );

  vtkMRMLSpatialObjectsNode *oldNode = this->SpatialObjectsNode();
  d->SpatialObjectsNode = SpatialObjectsNode;

  this->setSpatialObjectsDisplayNode(
    d->SpatialObjectsNode ? d->SpatialObjectsNode->GetNthDisplayNode( DisplayNodeIndex ) : NULL );

  qvtkReconnect( oldNode, this->SpatialObjectsNode(),
                vtkCommand::ModifiedEvent, this,
                SLOT( updateWidgetFromMRML() ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::
setSpatialObjectsDisplayNode( vtkMRMLNode* node )
{
  this->setSpatialObjectsDisplayNode(
    vtkMRMLSpatialObjectsDisplayNode::SafeDownCast( node ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::
setSpatialObjectsDisplayNode( vtkMRMLSpatialObjectsDisplayNode*
                             SpatialObjectsDisplayNode )
{
  Q_D( qSlicerSpatialObjectsWidget );

  vtkMRMLSpatialObjectsDisplayNode *oldDisplayNode =
    this->SpatialObjectsDisplayNode();
  vtkMRMLSpatialObjectsDisplayPropertiesNode *oldDisplayPropertiesNode =
    this->SpatialObjectsDisplayPropertiesNode();

  d->SpatialObjectsDisplayNode = SpatialObjectsDisplayNode;

  if( vtkMRMLSpatialObjectsLineDisplayNode::
        SafeDownCast( d->SpatialObjectsDisplayNode ) )
    {
    d->MaterialPropertyWidget->setHidden( true );
    d->MaterialPropertyGroupBox->setHidden( true );
    }
  else
    {
    d->MaterialPropertyWidget->setHidden( false );
    d->MaterialPropertyGroupBox->setHidden( false );
    }

  qvtkReconnect( oldDisplayNode,
                this->SpatialObjectsDisplayNode(),
                vtkCommand::ModifiedEvent, this,
                SLOT( updateWidgetFromMRML() ) );
  qvtkReconnect( oldDisplayPropertiesNode,
                this->SpatialObjectsDisplayPropertiesNode(),
                vtkCommand::ModifiedEvent, this,
                SLOT( updateWidgetFromMRML() ) );

  this->updateWidgetFromMRML();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setVisibility( bool state )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetVisibility( state );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::clickColorBySolid( bool checked )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  if( checked )
    {
    d->SpatialObjectsDisplayNode->SetColorModeToSolid();
    d->SpatialObjectsDisplayNode->SetScalarVisibility( 0 );
    }
  else
    {
    d->SpatialObjectsDisplayNode->SetColorModeToScalarData();
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::onColorBySolidChanged( const QColor &color )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetColor( color.redF(),
                                         color.greenF(),
                                         color.blueF() );
}
//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setColorBySolid()
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetColorModeToSolid();
  d->SpatialObjectsDisplayNode->SetScalarVisibility( 0 );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setSliceIntersection( int intersection )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetSliceIntersectionVisibility( intersection );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setSliceIntersectionThickness( double value )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetSliceIntersectionThickness( value );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setColorByScalar()
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetScalarRangeFlag( vtkMRMLDisplayNode::UseDisplayNodeScalarRange );
  d->SpatialObjectsDisplayNode->SetColorModeToScalarData();
  d->SpatialObjectsDisplayNode->SetScalarVisibility( 1 );
  this->onColorByScalarChanged( d->ColorByScalarComboBox->currentIndex() );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::onColorByScalarChanged( int scalarIndex )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode ||
      !d->SpatialObjectsNode ||
      !d->SpatialObjectsNode->GetPolyData() ||
      !d->SpatialObjectsNode->GetPolyData()->GetPointData() )
    {
    return;
    }

  QString activeScalarName = d->ColorByScalarComboBox->itemText( scalarIndex );
  d->SpatialObjectsDisplayNode->
    SetActiveScalarName( activeScalarName.toLatin1() );

  double range[2] = {0.,0.};
  if( scalarIndex < 0 )
    {
    d->ScalarRangeWidget->setRange( range[0], range[1] );
    }
  else
    {
    // Avoid having to process to an update request
    vtkDataArray* currentScalar = d->SpatialObjectsNode->GetPolyData()->
      GetPointData()->GetScalars( activeScalarName.toLatin1() );
    currentScalar->GetRange( range, -1 );

    double singleStep;
    // If range is 0, set singlestep to 1 to avoid having a singleStep of 0
    if( range[0] == range[1] )
      {
      singleStep = 1;
      }
    else
      {
      singleStep = ( range[1]-range[0] )/100;
      }

    // Workaround a bug of ctkRangeWidget
    // TO FIX in CTK
    if( range[1]-range[0] < d->ScalarRangeWidget->singleStep() )
      {
      d->ScalarRangeWidget->setSingleStep( singleStep );
      d->ScalarRangeWidget->setRange( range[0], range[1] );
      }
    else
      {
      d->ScalarRangeWidget->setRange( range[0], range[1] );
      d->ScalarRangeWidget->setSingleStep( singleStep );
      }

    d->ScalarRangeWidget->setValues( range[0], range[1] );
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::onColorByScalarValuesChanged( double minValue,
                                                              double maxValue )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode ||
      !d->SpatialObjectsNode ||
      !d->SpatialObjectsNode->GetPolyData() ||
      !d->SpatialObjectsNode->GetPolyData()->GetPointData() ||
      !d->SpatialObjectsDisplayNode->GetScalarVisibility() )
    {
    return;
    }

    // Set the values for coloring by scalar
    double values[2];
    values[0] = minValue;
    values[1] = maxValue;
    d->SpatialObjectsDisplayNode->SetScalarRange( values );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setColorByCellScalars()
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetColorModeToUseCellScalars();
  d->SpatialObjectsDisplayNode->SetScalarVisibility( 1 );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::
setColorByCellScalarsColorTable( vtkMRMLNode* colortableNode )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode || !colortableNode )
    {
    return;
    }

  Q_ASSERT( vtkMRMLColorNode::SafeDownCast( colortableNode ) );
  d->SpatialObjectsDisplayNode->
    SetAndObserveColorNodeID( colortableNode->GetID() );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setOpacity( double opacity )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetOpacity( opacity );
}

//------------------------------------------------------------------------------
QColor qSlicerSpatialObjectsWidget::color() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->color();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setColor( const QColor& color )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  // QColors loose precision in the numbers, don't reset the color if it didn't
  // "really" change.
  double* oldColor = d->SpatialObjectsDisplayNode->GetColor();
  if( QColor::fromRgbF( oldColor[0], oldColor[1], oldColor[2] ) != color )
    {
    d->SpatialObjectsDisplayNode->SetColor( color.redF(),
                                           color.greenF(),
                                           color.blueF() );
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setAmbient( double ambient )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetAmbient( ambient );
}

//------------------------------------------------------------------------------
double qSlicerSpatialObjectsWidget::ambient() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->ambient();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setDiffuse( double diffuse )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetDiffuse( diffuse );
}

//------------------------------------------------------------------------------
double qSlicerSpatialObjectsWidget::diffuse() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->diffuse();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setSpecular( double specular )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetSpecular( specular );
}

//------------------------------------------------------------------------------
double qSlicerSpatialObjectsWidget::specular()const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->specular();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setSpecularPower( double specularPower )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetPower( specularPower );
}

//------------------------------------------------------------------------------
double qSlicerSpatialObjectsWidget::specularPower() const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->specularPower();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::setBackfaceCulling( bool backfaceCulling )
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  d->SpatialObjectsDisplayNode->SetBackfaceCulling( backfaceCulling );
}

//------------------------------------------------------------------------------
bool qSlicerSpatialObjectsWidget::backfaceCulling()const
{
  Q_D( const qSlicerSpatialObjectsWidget );
  return d->MaterialPropertyWidget->backfaceCulling();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsWidget::updateWidgetFromMRML()
{
  Q_D( qSlicerSpatialObjectsWidget );

  if( !d->SpatialObjectsNode ||
      !d->SpatialObjectsDisplayNode ||
      this->inProcess )
    {
    return;
    }

  this->inProcess = 1;

  d->VisibilityCheckBox->setChecked(
    d->SpatialObjectsDisplayNode->GetVisibility() );
  d->OpacitySlider->setValue( d->SpatialObjectsDisplayNode->GetOpacity() );

  d->ColorByScalarsColorTableComboBox->setCurrentNodeID
    ( d->SpatialObjectsDisplayNode->GetColorNodeID() );

  // We have to make that check, otherwise, if the datasets are the same,
  // when calling d->ColorByScalarComboBox->setCurrentArray(), it will always actually
  // change the current array ( when it's the same name, it doesn't change it ) and fire
  // signals ( because the currentArray was previously "" )
  // If the signals are fired, it calls back onColorByScalarChanged, which resets the
  // range to the min and max of the dataArray. We have an override situation, and it
  // makes it impossible to move the slider.
  if( d->ColorByScalarComboBox->dataSet() != d->SpatialObjectsNode->GetPolyData() )
    {
    // Block the signals, because setDatSet modifies the comboBox which triggers
    // updateWidgetFromMRML, and onColorByScalarChanged which also triggers
    // updateWidgetFromMRML.
    // In total, updateWidgetFromMRML will be called 5 times...
    // the currentIndexChanged( int ) of the combobox is fired twice when doing
    // setDataset(). This has to be fixed in CTK.
    bool wasBlocking = d->ColorByScalarComboBox->blockSignals( true );

    // This function has to be fixed in CTK, it shouldn't fire 2 signals ( see above )
    d->ColorByScalarComboBox->setDataSet(
      vtkDataSet::SafeDownCast( d->SpatialObjectsNode->GetPolyData() ) );
    // Reset the current Array to blank, because setDataSet sets
    // the currentIndex to "TubeRadius" by default, and if we leave it like this,
    // the signal currentIndexChanged won't be fired when we switch to a node which
    // activeScalarName is TubeRadius ( and the rangeWidget won't be updated )
    d->ColorByScalarComboBox->setCurrentArray( "" );

    d->ColorByScalarComboBox->blockSignals( wasBlocking );
    }

  d->ColorByScalarComboBox->setCurrentArray(
    d->SpatialObjectsDisplayNode->GetActiveScalarName() );


  switch( d->SpatialObjectsDisplayNode->GetColorMode() )
    {
      case vtkMRMLSpatialObjectsDisplayNode::colorModeScalar:
         d->SpatialObjectsDisplayNode->SetColorMode(
           vtkMRMLSpatialObjectsDisplayNode::colorModeSolid );
         break;
      case vtkMRMLSpatialObjectsDisplayNode::colorModeSolid:
        {
        double color[3];

        d->SpatialObjectsDisplayNode->GetColor( color );
        d->ColorBySolidColorPicker->setColor(
          QColor::fromRgbF( color[0],color[1],color[2] ) );
        d->ColorBySolidColorCheckBox->setChecked( 1 );
        }
        break;

      case vtkMRMLSpatialObjectsDisplayNode::colorModeScalarData:
        {
        d->ColorByScalarRadioButton->setChecked( 1 );
        }
        break;
   }

  d->MaterialPropertyWidget->setColor(
    QColor::fromRgbF( d->SpatialObjectsDisplayNode->GetColor()[0],
                     d->SpatialObjectsDisplayNode->GetColor()[1],
                     d->SpatialObjectsDisplayNode->GetColor()[2] ) );

  d->MaterialPropertyWidget->
    setOpacity( d->SpatialObjectsDisplayNode->GetOpacity() );
  d->MaterialPropertyWidget->
    setAmbient( d->SpatialObjectsDisplayNode->GetAmbient() );
  d->MaterialPropertyWidget->
    setDiffuse( d->SpatialObjectsDisplayNode->GetDiffuse() );
  d->MaterialPropertyWidget->
    setSpecular( d->SpatialObjectsDisplayNode->GetSpecular() );
  d->MaterialPropertyWidget->
    setSpecularPower( d->SpatialObjectsDisplayNode->GetPower() );
  d->MaterialPropertyWidget->
    setBackfaceCulling( d->SpatialObjectsDisplayNode->GetBackfaceCulling() );

  this->inProcess = 0;
}
