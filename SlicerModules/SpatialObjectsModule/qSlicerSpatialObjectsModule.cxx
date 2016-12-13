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
#include <QtPlugin>

// SlicerQt includes
#include <qSlicerCoreApplication.h>
#include <qSlicerIOManager.h>
#include <qSlicerModuleManager.h>
#include <qSlicerNodeWriter.h>

// SpatialObjects Logic includes
#include "vtkSlicerSpatialObjectsLogic.h"

// SpatialObjects QTModule includes
#include "qSlicerSpatialObjectsModule.h"
#include "qSlicerSpatialObjectsModuleWidget.h"
#include "qSlicerSpatialObjectsReader.h"

// MRML Logic includes
#include <vtkMRMLColorLogic.h>

//------------------------------------------------------------------------------
Q_EXPORT_PLUGIN2( qSlicerSpatialObjectsModule, qSlicerSpatialObjectsModule );

//------------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SpatialObjects
class qSlicerSpatialObjectsModulePrivate
{};

//------------------------------------------------------------------------------
qSlicerSpatialObjectsModule::qSlicerSpatialObjectsModule( QObject* _parent )
  : Superclass( _parent )
  , d_ptr( new qSlicerSpatialObjectsModulePrivate )
{}

//------------------------------------------------------------------------------
qSlicerSpatialObjectsModule::~qSlicerSpatialObjectsModule()
{}

//------------------------------------------------------------------------------
QString qSlicerSpatialObjectsModule::helpText() const
{
  QString help = QString(
    "The SpatialObjects Module loads and adjusts display parameters of "
    "spatial object data.<br>"
    "<a href=\"%1/Documentation/%2.%3/Modules/SpatialObjects\">"
    "%1/Documentation/%2.%3/Modules/SpatialObjects</a><br>"
    "The SpatialObjects Editor allows modifying "
    "parameters ( gradients, bValues, measurement frame ) of the data and "
    "provides a quick way to interpret them. "
    "For that it shows glyphs, line and tubes "
    "for visual exploration.<br><br>"
    "Help for the SpatialObjects Editor: "
    "<a href=\"%1/Modules:SpatialObjects:SpatialObjects_Editor-Documentation\">"
    "%1/Modules:SpatialObjects:SpatialObjects_Editor-Documentation</a>" );

  return help.arg( this->slicerWikiUrl() ).
                    arg( Slicer_VERSION_MAJOR ).arg( Slicer_VERSION_MINOR );
}

//------------------------------------------------------------------------------
QString qSlicerSpatialObjectsModule::acknowledgementText() const
{
  QString acknowledgement = QString(
    "<center><table border=\"0\"><tr>"
    "<td><img src=\":Logos/NAMIC.png\" alt\"NA-MIC\"></td>"
    "</tr></table></center>"
    "This work was supported by NA-MIC and the Slicer "
    "Community. See <a href=\"http://www.slicer.org\">http://www.slicer.org"
    "</a> for details.<br>"
    "The SpatialObjects module was contributed by "
    "Michael Jeulin-L, Kitware Inc. " );

  return acknowledgement;
}

//------------------------------------------------------------------------------
QStringList qSlicerSpatialObjectsModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString( "Michael Jeulin-L ( Kitware )" );

  return moduleContributors;
}

//------------------------------------------------------------------------------
QIcon qSlicerSpatialObjectsModule::icon()const
{
  return QIcon( ":/Icons/SpatialObjects.png" );
}

//------------------------------------------------------------------------------
QStringList qSlicerSpatialObjectsModule::categories() const
{
  return QStringList() << "Surface Models";
}

//------------------------------------------------------------------------------
QStringList qSlicerSpatialObjectsModule::dependencies() const
{
  QStringList moduleDependencies;
  moduleDependencies << "Colors";

  return moduleDependencies;
}

//------------------------------------------------------------------------------
void qSlicerSpatialObjectsModule::setup()
{
  this->Superclass::setup();

  vtkSlicerSpatialObjectsLogic* spatialObjectsLogic =
    vtkSlicerSpatialObjectsLogic::SafeDownCast( this->logic() );

  qSlicerCoreIOManager* coreIOManager =
    qSlicerCoreApplication::application()->coreIOManager();
  coreIOManager->registerIO(
    new qSlicerSpatialObjectsReader( spatialObjectsLogic, this ) );
  coreIOManager->registerIO( new qSlicerNodeWriter(
    "SpatialObjects", QString( "SpatialObjectFile" ),
    QStringList() << "vtkMRMLSpatialObjectsNode", true, this ) );
}

//------------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation*
qSlicerSpatialObjectsModule::createWidgetRepresentation()
{
  return new qSlicerSpatialObjectsModuleWidget;
}

//------------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSpatialObjectsModule::createLogic()
{
  return vtkSlicerSpatialObjectsLogic::New();
}
