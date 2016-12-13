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

// Qt includes
#include "qSlicerSpatialObjectsModuleWidget.h"
#include "ui_qSlicerSpatialObjectsModule.h"
#include "qMRMLSceneSpatialObjectsModel.h"

// MRML includes
#include <vtkMRMLNode.h>
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsDisplayNode.h"
#include "vtkMRMLSpatialObjectsStorageNode.h"
#include "vtkMRMLSpatialObjectsTubeDisplayNode.h"
#include <vtkMRMLScene.h>

//------------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SpatialObjects
class qSlicerSpatialObjectsModuleWidgetPrivate :
 public Ui_qSlicerSpatialObjectsModule
{
  Q_DECLARE_PUBLIC( qSlicerSpatialObjectsModuleWidget );

protected:
  qSlicerSpatialObjectsModuleWidget* const q_ptr;

public:
  qSlicerSpatialObjectsModuleWidgetPrivate(
    qSlicerSpatialObjectsModuleWidget& object );
  void init();

  vtkMRMLSpatialObjectsNode* spatialObjectsNode;
};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsModuleWidgetPrivate
::qSlicerSpatialObjectsModuleWidgetPrivate(
  qSlicerSpatialObjectsModuleWidget& object )
  : q_ptr( &object )
{
  this->spatialObjectsNode = NULL;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModuleWidgetPrivate::init()
{
  Q_Q( qSlicerSpatialObjectsModuleWidget );

  this->setupUi( q );

  QObject::connect( this->SpatialObjectsTreeView,
                   SIGNAL( visibilityChanged( int ) ),
                   this->TractDisplayModesTabWidget,
                   SLOT( setCurrentIndex ( int ) ) );

  // Hide the GlyphTab since, the glyph representations
  // are not build/implemented yet.
  this->TractDisplayModesTabWidget->removeTab( 2 );
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsModuleWidget::
qSlicerSpatialObjectsModuleWidget( QWidget* _parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerSpatialObjectsModuleWidgetPrivate( *this ) )
{
  Q_D( qSlicerSpatialObjectsModuleWidget );
  d->init();
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsModuleWidget::~qSlicerSpatialObjectsModuleWidget()
{}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModuleWidget::setup()
{}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModuleWidget::
setSpatialObjectsNode( vtkMRMLNode* inputNode )
{
  this->setSpatialObjectsNode(
    vtkMRMLSpatialObjectsNode::SafeDownCast( inputNode ) );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModuleWidget::
setSpatialObjectsNode( vtkMRMLSpatialObjectsNode* spatialObjectsNode )
{
  Q_D( qSlicerSpatialObjectsModuleWidget );

  if( d->spatialObjectsNode == spatialObjectsNode )
    return;

  d->spatialObjectsNode = spatialObjectsNode;

  d->LineDisplayWidget->setSpatialObjectsNode( spatialObjectsNode, 0 );
  d->TubeDisplayWidget->setSpatialObjectsNode( spatialObjectsNode, 1 );

  emit currentNodeChanged( d->spatialObjectsNode );
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModuleWidget::setSolidTubeColor( bool solid )
{
  std::vector<vtkMRMLNode *> nodes;
  this->mrmlScene()->
    GetNodesByClass( "vtkMRMLSpatialObjectsTubeDisplayNode", nodes );

  vtkMRMLSpatialObjectsTubeDisplayNode* node = 0;
  for( unsigned int i = 0; i < nodes.size(); ++i )
    {
    node = vtkMRMLSpatialObjectsTubeDisplayNode::SafeDownCast( nodes[i] );
    if( solid )
      {
      node->SetColorMode( vtkMRMLSpatialObjectsDisplayNode::colorModeSolid );
      }
    else
      {
      node->SetColorMode( vtkMRMLSpatialObjectsDisplayNode::colorModeScalar );
      }
    }
}
