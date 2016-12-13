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

#ifndef __qSlicerSpatialObjectsReader_h
#define __qSlicerSpatialObjectsReader_h

// SlicerQt includes
#include <qSlicerFileReader.h>

class qSlicerSpatialObjectsReaderPrivate;
class vtkSlicerSpatialObjectsLogic;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SpatialObjects
class qSlicerSpatialObjectsReader : public qSlicerFileReader
{
  Q_OBJECT

public:
  typedef qSlicerFileReader Superclass;

  qSlicerSpatialObjectsReader( QObject* parent = 0 );
  qSlicerSpatialObjectsReader( vtkSlicerSpatialObjectsLogic* logic,
                              QObject* parent = 0 );

  virtual ~qSlicerSpatialObjectsReader();

  vtkSlicerSpatialObjectsLogic* logic()const;
  void setLogic( vtkSlicerSpatialObjectsLogic* logic );

  virtual QString description()const;
  virtual IOFileType fileType()const;
  virtual QStringList extensions()const;
  virtual qSlicerIOOptions* options()const;

  virtual bool load( const IOProperties& properties );

protected:
  QScopedPointer<qSlicerSpatialObjectsReaderPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE( qSlicerSpatialObjectsReader );
  Q_DISABLE_COPY( qSlicerSpatialObjectsReader );
};

#endif
