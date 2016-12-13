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
#include "qSlicerSpatialObjectsGlyphWidget.h"
#include "ui_qSlicerSpatialObjectsGlyphWidget.h"

// MRML includes
#include <vtkMRMLSpatialObjectsDisplayNode.h>
#include <vtkMRMLSpatialObjectsDisplayPropertiesNode.h>

//------------------------------------------------------------------------------
class qSlicerSpatialObjectsGlyphWidgetPrivate
  : public Ui_qSlicerSpatialObjectsGlyphWidget
{
Q_DECLARE_PUBLIC( qSlicerSpatialObjectsGlyphWidget );

protected:
  qSlicerSpatialObjectsGlyphWidget* const q_ptr;

public:
  qSlicerSpatialObjectsGlyphWidgetPrivate(
    qSlicerSpatialObjectsGlyphWidget& object );

  void init();
  bool centeredOrigin( double* origin ) const;

  vtkMRMLSpatialObjectsDisplayNode* SpatialObjectsDisplayNode;
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    SpatialObjectsDisplayPropertiesNode;
};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsGlyphWidgetPrivate::
qSlicerSpatialObjectsGlyphWidgetPrivate(
  qSlicerSpatialObjectsGlyphWidget& object )
  : q_ptr( &object )
{
  this->SpatialObjectsDisplayNode = 0;
  this->SpatialObjectsDisplayPropertiesNode = 0;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidgetPrivate::init()
{
  Q_Q( qSlicerSpatialObjectsGlyphWidget );
  this->setupUi( q );

  QObject::connect( this->GlyphTypeSelector,
                   SIGNAL( currentIndexChanged( int ) ), q,
                   SLOT( setGlyphType( int ) ) );
  QObject::connect( this->ScaleFactorSlider,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setGlyphScaleFactor( double ) ) );
  QObject::connect( this->SpacingSlider,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setGlyphSpacing( double ) ) );
  QObject::connect( this->GlyphSidesSlider,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setTubeGlyphNumberOfSides( double ) ) );
  QObject::connect( this->GlyphRadiusSlider,
                   SIGNAL( valueChanged( double ) ), q,
                   SLOT( setTubeGlyphRadius( double ) ) );
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsGlyphWidget::
qSlicerSpatialObjectsGlyphWidget( QWidget *_parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerSpatialObjectsGlyphWidgetPrivate( *this ) )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );
  d->init();
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsGlyphWidget::~qSlicerSpatialObjectsGlyphWidget()
{}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* qSlicerSpatialObjectsGlyphWidget::
spatialObjectsDisplayNode() const
{
  Q_D( const qSlicerSpatialObjectsGlyphWidget );

  return d->SpatialObjectsDisplayNode;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayPropertiesNode* qSlicerSpatialObjectsGlyphWidget::
spatialObjectsDisplayPropertiesNode() const
{
  Q_D( const qSlicerSpatialObjectsGlyphWidget );

  return d->SpatialObjectsDisplayNode ?
    d->SpatialObjectsDisplayNode->GetSpatialObjectsDisplayPropertiesNode() : 0;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::
setSpatialObjectsDisplayNode( vtkMRMLNode* node )
{
  this->setSpatialObjectsDisplayNode(
    vtkMRMLSpatialObjectsDisplayNode::SafeDownCast( node ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::
setSpatialObjectsDisplayNode(
  vtkMRMLSpatialObjectsDisplayNode* SpatialObjectsDisplayNode )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  vtkMRMLSpatialObjectsDisplayNode *oldDisplayNode =
    d->SpatialObjectsDisplayNode;

  d->SpatialObjectsDisplayNode = SpatialObjectsDisplayNode;

  qvtkReconnect( oldDisplayNode,
                d->SpatialObjectsDisplayNode,
                vtkCommand::ModifiedEvent,
                this,
                SLOT( updateWidgetFromMRMLDisplayNode() ) );

  this->updateWidgetFromMRMLDisplayNode();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::
setSpatialObjectsDisplayPropertiesNode( vtkMRMLNode* node )
{
  this->setSpatialObjectsDisplayPropertiesNode
    ( vtkMRMLSpatialObjectsDisplayPropertiesNode::SafeDownCast( node ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::
setSpatialObjectsDisplayPropertiesNode(
  vtkMRMLSpatialObjectsDisplayPropertiesNode* node )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  d->SpatialObjectsDisplayPropertiesNode = node;

  // Selection of the Glyph type
  d->GlyphTypeSelector->blockSignals( true );
  d->GlyphTypeSelector->clear();

  int i = d->SpatialObjectsDisplayPropertiesNode->GetFirstGlyphGeometry();
  for( ; i <= d->SpatialObjectsDisplayPropertiesNode->GetLastGlyphGeometry();
       ++i )
    {
      std::cout << "Glyph" << i << ": "
                << d->SpatialObjectsDisplayPropertiesNode->GetGlyphGeometryAsString( i )
                << std::endl;
    d->GlyphTypeSelector->addItem(
      d->SpatialObjectsDisplayPropertiesNode->GetGlyphGeometryAsString( i ), i );
    }

  d->GlyphTypeSelector->blockSignals( false );
}
//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::setGlyphType( int type )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->SpatialObjectsDisplayPropertiesNode->SetGlyphGeometry( type );
  QWidget* widget = d->GlyphSubPropertiesWidget->findChild<QWidget*>(
    d->SpatialObjectsDisplayPropertiesNode->GetGlyphGeometryAsString() );

  if( widget )
    {
    d->GlyphSubPropertiesWidget->setCurrentWidget( widget );
    d->GlyphSubPropertiesWidget->setEnabled( true );
    }
  else
    {
    d->GlyphSubPropertiesWidget->setEnabled( false );
    }
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::setGlyphScaleFactor( double scale )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->SpatialObjectsDisplayPropertiesNode->SetGlyphScaleFactor( scale );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::setGlyphSpacing( double spacing )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->SpatialObjectsDisplayPropertiesNode->SetLineGlyphResolution(
    static_cast<int>( spacing ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::setTubeGlyphNumberOfSides( double sides )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->SpatialObjectsDisplayPropertiesNode->SetTubeGlyphNumberOfSides(
    static_cast<int>( sides ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::setTubeGlyphRadius( double radius )
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->SpatialObjectsDisplayPropertiesNode->SetTubeGlyphRadius( radius );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::updateWidgetFromMRMLDisplayNode()
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayNode )
    {
    return;
    }

  if( d->SpatialObjectsDisplayPropertiesNode !=
      d->SpatialObjectsDisplayNode->GetSpatialObjectsDisplayPropertiesNode() )
    {
    this->setSpatialObjectsDisplayPropertiesNode(
      d->SpatialObjectsDisplayNode->GetSpatialObjectsDisplayPropertiesNode() );
    }

  this->updateWidgetFromMRMLDisplayPropertiesNode();
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsGlyphWidget::
updateWidgetFromMRMLDisplayPropertiesNode()
{
  Q_D( qSlicerSpatialObjectsGlyphWidget );

  if( !d->SpatialObjectsDisplayNode ||
      !d->SpatialObjectsDisplayPropertiesNode )
    {
    return;
    }

  d->GlyphTypeSelector->setCurrentIndex(
    d->SpatialObjectsDisplayPropertiesNode->GetGlyphGeometry() );

  d->ScaleFactorSlider->setValue(
    d->SpatialObjectsDisplayPropertiesNode->GetGlyphScaleFactor() );
  d->SpacingSlider->setValue(
    d->SpatialObjectsDisplayPropertiesNode->GetLineGlyphResolution() );
  d->GlyphSidesSlider->setValue(
    d->SpatialObjectsDisplayPropertiesNode->GetTubeGlyphNumberOfSides() );
  d->GlyphRadiusSlider->setValue(
    d->SpatialObjectsDisplayPropertiesNode->GetTubeGlyphRadius() );
}
