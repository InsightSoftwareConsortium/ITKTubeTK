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
#include <QDebug>
#include <QHeaderView>
#include <QMessageBox>
#include <QMouseEvent>

// SpatialObjectss includes
#include "qMRMLSpatialObjectsTreeView.h"
#include "qMRMLSceneSpatialObjectsModel.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsDisplayNode.h"
#include "vtkMRMLSpatialObjectsLineDisplayNode.h"
#include "vtkMRMLSpatialObjectsGlyphDisplayNode.h"
#include "vtkMRMLSpatialObjectsTubeDisplayNode.h"

//------------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SpatialObjects
class qMRMLSpatialObjectsTreeViewPrivate
{
  Q_DECLARE_PUBLIC( qMRMLSpatialObjectsTreeView );
protected:
  qMRMLSpatialObjectsTreeView* const q_ptr;

public:
  qMRMLSpatialObjectsTreeViewPrivate( qMRMLSpatialObjectsTreeView& object );
  void init();

  qMRMLSceneSpatialObjectsModel* SceneModel;
  qMRMLSortFilterProxyModel* SortFilterModel;
};

//------------------------------------------------------------------------------
qMRMLSpatialObjectsTreeViewPrivate::
qMRMLSpatialObjectsTreeViewPrivate( qMRMLSpatialObjectsTreeView& object )
  : q_ptr( &object )
{
  this->SceneModel = 0;
  this->SortFilterModel = 0;
}

//------------------------------------------------------------------------------
void qMRMLSpatialObjectsTreeViewPrivate::init()
{
  Q_Q( qMRMLSpatialObjectsTreeView );

  this->SceneModel = new qMRMLSceneSpatialObjectsModel( q );
  q->setSceneModel( this->SceneModel, "SpatialObjects" );

  // We only want to show vtkMRMLSpatialObjectsNodes.
  QStringList nodeTypes = QStringList();
  nodeTypes.append( "vtkMRMLSpatialObjectsNode" );

  q->setNodeTypes( nodeTypes );
  this->SortFilterModel = q->sortFilterProxyModel();

  q->header()->setStretchLastSection( false );
  q->header()->setResizeMode( QHeaderView::ResizeToContents );
  q->header()->setResizeMode( 0, QHeaderView::Stretch );

  q->setUniformRowHeights( true );
}

//------------------------------------------------------------------------------
qMRMLSpatialObjectsTreeView::qMRMLSpatialObjectsTreeView( QWidget *_parent )
  :qMRMLTreeView( _parent )
  , d_ptr( new qMRMLSpatialObjectsTreeViewPrivate( *this ) )
{
  Q_D( qMRMLSpatialObjectsTreeView );
  d->init();

  // We need to enable mouse tracking to set the appropriate
  // cursor while mouseMove occurs.
  this->setMouseTracking( true );
}

//------------------------------------------------------------------------------
qMRMLSpatialObjectsTreeView::~qMRMLSpatialObjectsTreeView()
{}

//------------------------------------------------------------------------------
#ifndef QT_NO_CURSOR
bool qMRMLSpatialObjectsTreeView::viewportEvent( QEvent* e )
{
  // reset the cursor if we leave the viewport
  if( e->type() == QEvent::Leave )
    {
    setCursor( QCursor() );
    }

  return QTreeView::viewportEvent( e );
}
#endif

//------------------------------------------------------------------------------
void qMRMLSpatialObjectsTreeView::onVisibilityColumnClicked( vtkMRMLNode* node )
{
  if( !node )
    {
    return;
    }
}

//-----------------------------------------------------------------------------
/// Set and observe the logic
//-----------------------------------------------------------------------------
void qMRMLSpatialObjectsTreeView::setLogic( vtkSlicerSpatialObjectsLogic* logic )
{
  if( !logic )
    {
    return;
    }

  this->Logic = logic;
}

//------------------------------------------------------------------------------
void qMRMLSpatialObjectsTreeView::setMRMLScene( vtkMRMLScene* scene )
{
  this->Superclass::setMRMLScene( scene );
  this->setRootIndex( this->sortFilterProxyModel()->mrmlSceneIndex() );
}

//------------------------------------------------------------------------------
bool qMRMLSpatialObjectsTreeView::clickDecoration( const QModelIndex& index )
{
  bool res = false;
  int type = -1;
  QModelIndex sourceIndex = this->sortFilterProxyModel()->mapToSource( index );

  vtkMRMLSpatialObjectsNode* soNode = vtkMRMLSpatialObjectsNode::SafeDownCast( 
    this->sortFilterProxyModel()->mrmlNodeFromIndex( index ) );

  vtkMRMLSpatialObjectsDisplayNode* lineDisplayNode =
    soNode->GetLineDisplayNode();
  vtkMRMLSpatialObjectsDisplayNode* tubeDisplayNode =
    soNode->GetTubeDisplayNode();
  vtkMRMLSpatialObjectsDisplayNode* glyphDisplayNode =
    soNode->GetGlyphDisplayNode();

  qMRMLSceneSpatialObjectsModel* model =
    dynamic_cast<qMRMLSceneSpatialObjectsModel*>( this->sceneModel() );

  if( !( sourceIndex.flags() & Qt::ItemIsEnabled ) )
    {
    res = false;
    }
  else if( sourceIndex.column() == model->lineVisibilityColumn() )
    {
    type = 0;
    if( lineDisplayNode )
      {
      lineDisplayNode->SetVisibility( lineDisplayNode->GetVisibility() ? 0 : 1 );
      res = true;
      }
    }
  else if( sourceIndex.column() == model->tubeVisibilityColumn() )
    {
    type = 1;
    if( tubeDisplayNode )
      {
      tubeDisplayNode->SetVisibility( tubeDisplayNode->GetVisibility() ? 0 : 1 );
      res = true;
      }
    }
  else if( sourceIndex.column() == model->glyphVisibilityColumn() )
    {
    type = 2;
    if( glyphDisplayNode )
      {
      glyphDisplayNode->
        SetVisibility( glyphDisplayNode->GetVisibility() ? 0 : 1 );
      res = true;
      }
    }

  if( glyphDisplayNode->GetVisibility() == 0 &&
      tubeDisplayNode->GetVisibility() == 0  &&
      lineDisplayNode->GetVisibility() == 1 )
    {
    type = 0;
    }
  if( glyphDisplayNode->GetVisibility() == 0 &&
      tubeDisplayNode->GetVisibility() == 1  &&
      lineDisplayNode->GetVisibility() == 0 )
    {
    type = 1;
    }
  if( glyphDisplayNode->GetVisibility() == 1 &&
      tubeDisplayNode->GetVisibility() == 0  &&
      lineDisplayNode->GetVisibility() == 0 )
    {
    type = 2;
    }

  if( res )
    {
    emit decorationClicked( index );
    }
  if( type > -1 )
    {
    emit visibilityChanged( type );
    }
  return res;
}
