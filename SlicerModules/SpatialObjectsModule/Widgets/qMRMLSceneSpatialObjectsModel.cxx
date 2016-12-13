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
#include "qMRMLSceneSpatialObjectsModel.h"
#include <qMRMLSceneDisplayableModel_p.h>

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLSpatialObjectsDisplayNode.h>
#include <vtkMRMLSpatialObjectsLineDisplayNode.h>
#include <vtkMRMLSpatialObjectsTubeDisplayNode.h>
#include <vtkMRMLSpatialObjectsGlyphDisplayNode.h>


class qMRMLSceneSpatialObjectsModelPrivate
  : public qMRMLSceneDisplayableModelPrivate
{
protected:
  Q_DECLARE_PUBLIC( qMRMLSceneSpatialObjectsModel );

public:
  typedef qMRMLSceneDisplayableModelPrivate Superclass;
  qMRMLSceneSpatialObjectsModelPrivate( qMRMLSceneSpatialObjectsModel& object );
  virtual void init();

  int ColorColumn;
  int LineVisibilityColumn;
  int TubeVisibilityColumn;
  int GlyphVisibilityColumn;

};

//------------------------------------------------------------------------------
qMRMLSceneSpatialObjectsModelPrivate
::qMRMLSceneSpatialObjectsModelPrivate( qMRMLSceneSpatialObjectsModel& object )
  : Superclass( object )
{
  this->LineVisibilityColumn = -1;
  this->TubeVisibilityColumn = -1;
  this->GlyphVisibilityColumn = -1;
  this->ColorColumn = -1;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModelPrivate::init()
{}

//------------------------------------------------------------------------------
qMRMLSceneSpatialObjectsModel::qMRMLSceneSpatialObjectsModel( QObject *vparent )
  : Superclass( new qMRMLSceneSpatialObjectsModelPrivate( *this ), vparent )
{
  Q_D( qMRMLSceneSpatialObjectsModel );
  d->init();

  this->setListenNodeModifiedEvent( qMRMLSceneModel::AllNodes );
  this->setVisibilityColumn( -1 );
  this->setOpacityColumn( -1 );
  this->setCheckableColumn( -1 );

  this->setNameColumn( 0 );
  this->setLineVisibilityColumn( 1 );
  this->setTubeVisibilityColumn( 2 );
  //this->setGlyphVisibilityColumn( 3 );
  //this->setColorColumn( 4 );
  this->setColorColumn( 3 );

  //this->setColumnCount( 5 );
  this->setColumnCount( 4 );

  this->setHorizontalHeaderLabels( 
    QStringList() << "Name" << "Lines" << "Tubes" << "Tubes" );
  // TODO Reactivate Glyphs controls since the representations are availabe
  //this->setHorizontalHeaderLabels( 
  //  QStringList() << "Name" << "Lines" << "Tubes" << "Glyphs" << "Tubes" );
}

//------------------------------------------------------------------------------
qMRMLSceneSpatialObjectsModel::qMRMLSceneSpatialObjectsModel( 
  qMRMLSceneSpatialObjectsModelPrivate* pimpl, QObject *vparent )
  : Superclass( pimpl, vparent )
{}

//------------------------------------------------------------------------------
qMRMLSceneSpatialObjectsModel::~qMRMLSceneSpatialObjectsModel()
{}

//------------------------------------------------------------------------------
QFlags<Qt::ItemFlag> qMRMLSceneSpatialObjectsModel::nodeFlags( vtkMRMLNode* node,
                                                              int column ) const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );

  QFlags<Qt::ItemFlag> flags = this->Superclass::nodeFlags( node, column );
  if( column == this->lineVisibilityColumn() &&
      d->displayNode( node ) != 0 )
    {
    flags |= Qt::ItemIsEditable;
    }
  if( column == this->tubeVisibilityColumn() &&
      d->displayNode( node ) != 0 )
    {
    flags |= Qt::ItemIsEditable;
    }
  if( column == this->glyphVisibilityColumn() &&
      d->displayNode( node ) != 0 )
    {
    flags |= Qt::ItemIsEditable;
    }
  if( column == this->colorColumn() &&
      d->displayNode( node ) != 0 )
    {
    flags |= Qt::ItemIsEditable;
    }

  return flags;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel
::updateItemDataFromNode( QStandardItem* item, vtkMRMLNode* node, int column )
{
  vtkMRMLSpatialObjectsNode *soNode =
    vtkMRMLSpatialObjectsNode::SafeDownCast( node );
  if( !soNode )
    {
    return;
    }

  if( column == this->colorColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
        soNode->GetTubeDisplayNode();

    if( displayNode )
      {
      double* rgbF = displayNode->GetColor();
      QColor color =
        QColor::fromRgbF( rgbF[0], rgbF[1], rgbF[2], displayNode->GetOpacity() );
      item->setData( color, Qt::DecorationRole );
      item->setToolTip( "Color" );

      if( displayNode->GetColorMode() ==
            vtkMRMLSpatialObjectsDisplayNode::colorModeSolid )
        {
        item->setEnabled( true );
        }
      else
        {
        item->setEnabled( false );
        }
      }
    }
  else if( column == this->lineVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetLineDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromNode( item, displayNode );
      }
  }
  else if( column == this->tubeVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetTubeDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromNode( item, displayNode );
      }
    }
  else if( column == this->glyphVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetGlyphDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromNode( item, displayNode );
      }
    }
  else
    {
    this->Superclass::updateItemDataFromNode( item, node, column );
    }
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::
updateVilibilityFromNode( QStandardItem* item,
                         vtkMRMLNode* node,
                         bool vtkNotUsed( slice ) )
{
  Q_D( qMRMLSceneSpatialObjectsModel );

  vtkMRMLDisplayNode* displayNode = vtkMRMLDisplayNode::SafeDownCast( node );

  int visible = -1;
  if( displayNode )
    {
    visible = displayNode->GetVisibility();
    }
  switch( visible )
    {
    case 0:
      // It should be fine to set the icon even if it is the same, but due
      // to a bug in Qt ( http://bugreports.qt.nokia.com/browse/QTBUG-20248 ),
      // it would fire a superflous itemChanged() signal.
      if( item->icon().cacheKey() != d->HiddenIcon.cacheKey() )
        {
        item->setIcon( d->HiddenIcon );
        }
      break;
    case 1:
      // It should be fine to set the icon even if it is the same, but due
      // to a bug in Qt ( http://bugreports.qt.nokia.com/browse/QTBUG-20248 ),
      // it would fire a superflous itemChanged() signal.
      if( item->icon().cacheKey() != d->VisibleIcon.cacheKey() )
        {
        item->setIcon( d->VisibleIcon );
        }
      break;
    case 2:
      // It should be fine to set the icon even if it is the same, but due
      // to a bug in Qt ( http://bugreports.qt.nokia.com/browse/QTBUG-20248 ),
      // it would fire a superflous itemChanged() signal.
      if( item->icon().cacheKey() != d->PartiallyVisibleIcon.cacheKey() )
        {
        item->setIcon( d->PartiallyVisibleIcon );
        }
      break;
    default:
      break;
    }
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::
updateNodeFromItemData( vtkMRMLNode* node, QStandardItem* item )
{
  vtkMRMLSpatialObjectsNode *soNode =
    vtkMRMLSpatialObjectsNode::SafeDownCast( node );
  if( !soNode )
    {
    return;
    }

  if( item->column() == this->colorColumn() )
    {
    QColor color = item->data( Qt::DecorationRole ).value<QColor>();
    // Invalid color can happen when the item hasn't been initialized yet
    if( color.isValid() )
      {
      vtkMRMLSpatialObjectsDisplayNode* displayNode =
        soNode->GetTubeDisplayNode();
      if( displayNode )
        {
        int wasModifying = displayNode->StartModify();

        // QColor looses precision, don't change color/opacity if not "really"
        // changed.
        QColor oldColor = QColor::fromRgbF( displayNode->GetColor()[0],
                                           displayNode->GetColor()[1],
                                           displayNode->GetColor()[2],
                                           displayNode->GetOpacity() );
        if( oldColor != color )
          {
          displayNode->SetColor( color.redF(), color.greenF(), color.blueF() );
          displayNode->SetOpacity( color.alphaF() );
          }

        displayNode->EndModify( wasModifying );
        }
      }
    }

  else if( item->column() == this->lineVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetLineDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromItem( item, displayNode );
      }
  }
  else if( item->column() == this->tubeVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetTubeDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromItem( item, displayNode );
      }
    }
  else if( item->column() == this->glyphVisibilityColumn() )
    {
    vtkMRMLSpatialObjectsDisplayNode* displayNode =
      soNode->GetGlyphDisplayNode();
    if( displayNode )
      {
      this->updateVilibilityFromItem( item, displayNode );
      }
    }

  return this->Superclass::updateNodeFromItemData( node, item );
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::
updateVilibilityFromItem( QStandardItem* item,
                         vtkMRMLNode* node,
                         bool vtkNotUsed( slice ) )
{
  Q_D( qMRMLSceneSpatialObjectsModel );

  vtkMRMLDisplayNode* displayNode = vtkMRMLDisplayNode::SafeDownCast( node );

  int visible = -1;
  if( item->icon().cacheKey() == d->HiddenIcon.cacheKey() )
    {
    visible = 0;
    }
  else if( item->icon().cacheKey() == d->VisibleIcon.cacheKey() )
    {
    visible = 1;
    }
  else if( item->icon().cacheKey() == d->PartiallyVisibleIcon.cacheKey() )
    {
    visible = 2;
    }
  if( displayNode )
    {
    displayNode->SetVisibility( visible );
    }
  }

//------------------------------------------------------------------------------
int qMRMLSceneSpatialObjectsModel::colorColumn()const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );
  return d->ColorColumn;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::setColorColumn( int column )
{
  Q_D( qMRMLSceneSpatialObjectsModel );
  d->ColorColumn = column;
  this->updateColumnCount();
}

//------------------------------------------------------------------------------
int qMRMLSceneSpatialObjectsModel::lineVisibilityColumn()const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );
  return d->LineVisibilityColumn;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::setLineVisibilityColumn( int column )
{
  Q_D( qMRMLSceneSpatialObjectsModel );
  d->LineVisibilityColumn = column;
  this->updateColumnCount();
}
//------------------------------------------------------------------------------
int qMRMLSceneSpatialObjectsModel::tubeVisibilityColumn()const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );
  return d->TubeVisibilityColumn;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::setTubeVisibilityColumn( int column )
{
  Q_D( qMRMLSceneSpatialObjectsModel );
  d->TubeVisibilityColumn = column;
  this->updateColumnCount();
}

//------------------------------------------------------------------------------
int qMRMLSceneSpatialObjectsModel::glyphVisibilityColumn()const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );
  return d->GlyphVisibilityColumn;
}

//------------------------------------------------------------------------------
void qMRMLSceneSpatialObjectsModel::setGlyphVisibilityColumn( int column )
{
  Q_D( qMRMLSceneSpatialObjectsModel );
  d->GlyphVisibilityColumn = column;
  this->updateColumnCount();
}

//------------------------------------------------------------------------------
int qMRMLSceneSpatialObjectsModel::maxColumnId()const
{
  Q_D( const qMRMLSceneSpatialObjectsModel );

  int maxId = this->Superclass::maxColumnId();
  maxId = qMax( maxId, d->LineVisibilityColumn );
  maxId = qMax( maxId, d->TubeVisibilityColumn );
  maxId = qMax( maxId, d->GlyphVisibilityColumn );
  maxId = qMax( maxId, d->ColorColumn );

  return maxId;
}
