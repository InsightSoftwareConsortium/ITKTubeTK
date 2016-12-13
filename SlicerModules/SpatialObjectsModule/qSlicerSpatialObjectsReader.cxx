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
#include <QFileInfo>
#include <QDir>

// SlicerQt includes
#include "qSlicerSpatialObjectsReader.h"

// TODO Create the Reader OptionWidget ?
//#include "qSlicerSpatialObjectsReaderOptionsWidget.h"

// Logic includes
#include <vtkSlicerApplicationLogic.h>
#include "vtkSlicerSpatialObjectsLogic.h"

// MRML includes
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

//------------------------------------------------------------------------------
class qSlicerSpatialObjectsReaderPrivate
{
public:
  vtkSmartPointer<vtkSlicerSpatialObjectsLogic> Logic;
};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsReader::qSlicerSpatialObjectsReader( QObject* _parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerSpatialObjectsReaderPrivate )
{}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsReader::
qSlicerSpatialObjectsReader( vtkSlicerSpatialObjectsLogic* logic,
                            QObject* _parent )
 : Superclass( _parent )
 , d_ptr( new qSlicerSpatialObjectsReaderPrivate )
{
  this->setLogic( logic );
}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsReader::~qSlicerSpatialObjectsReader()
{}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsReader::setLogic( vtkSlicerSpatialObjectsLogic* logic )
{
  Q_D( qSlicerSpatialObjectsReader );
  d->Logic = logic;
}

//------------------------------------------------------------------------------
vtkSlicerSpatialObjectsLogic* qSlicerSpatialObjectsReader::logic() const
{
  Q_D( const qSlicerSpatialObjectsReader );
  return d->Logic.GetPointer();
}

//-----------------------------------------------------------------------------
QString qSlicerSpatialObjectsReader::description() const
{
  return "SpatialObjects";
}

//-----------------------------------------------------------------------------
qSlicerSpatialObjectsReader::
IOFileType qSlicerSpatialObjectsReader::fileType() const
{
  return QString( "SpatialObjectFile" );
}

//-----------------------------------------------------------------------------
QStringList qSlicerSpatialObjectsReader::extensions() const
{
  return QStringList()
    << "VesselTubeSpatialObject ( *.tre )"
    << "DTITubeSpatialObject ( *.tre )"
    << "TubeSpatialObject ( *.tre )"
    << "All Files ( * )";
}

//-----------------------------------------------------------------------------
qSlicerIOOptions* qSlicerSpatialObjectsReader::options() const
{
  return 0;
  // TODO add SPATIALOBjects Options
  //return new qSlicerSpatialObjectsReaderOptionsWidget;
}

//-----------------------------------------------------------------------------
bool qSlicerSpatialObjectsReader::load( const IOProperties& properties )
{
  Q_D( qSlicerSpatialObjectsReader );
  Q_ASSERT( d->Logic );
  Q_ASSERT( properties.contains( "fileName" ) );
  QString fileName = properties["fileName"].toString();

  QStringList fileNames;
  if( properties.contains( "suffix" ) )
    {
    QStringList suffixList = properties["suffix"].toStringList();
    suffixList.removeDuplicates();

    // here filename describes a directory
    Q_ASSERT( QFileInfo( fileName ).isDir() );
    QDir dir( fileName );

    // suffix should be of style: *.png
    fileNames = dir.entryList( suffixList );
    }
  else
    {
    fileNames << fileName;
    }

  QStringList nodes;
  foreach( QString file, fileNames )
    {
    vtkMRMLSpatialObjectsNode* node =
      d->Logic->AddSpatialObject( file.toLatin1() );

    if( node )
      {
      if( properties.contains( "name" ) )
        {
        std::string uname = this->mrmlScene()->GetUniqueNameByString( 
          properties["name"].toString().toLatin1() );
        node->SetName( uname.c_str() );
        }
      nodes << node->GetID();
      }
    }
  this->setLoadedNodes( nodes );

  return nodes.size() > 0;
}
