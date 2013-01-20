/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Michael Jeulin-L, Kitware Inc.

==============================================================================*/

#ifndef __qSlicerSpatialObjectsReader
#define __qSlicerSpatialObjectsReader

// SlicerQt includes
#include "qSlicerFileReader.h"

class qSlicerSpatialObjectsReaderPrivate;
class vtkSlicerSpatialObjectsLogic;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SpatialObjects
class qSlicerSpatialObjectsReader : public qSlicerFileReader
{
Q_OBJECT

public:
  typedef qSlicerFileReader Superclass;

  qSlicerSpatialObjectsReader(QObject* parent = 0);
  qSlicerSpatialObjectsReader(vtkSlicerSpatialObjectsLogic* logic,
                              QObject* parent = 0);

  virtual ~qSlicerSpatialObjectsReader();

  vtkSlicerSpatialObjectsLogic* logic()const;
  void setLogic(vtkSlicerSpatialObjectsLogic* logic);

  virtual QString description()const;
  virtual IOFileType fileType()const;
  virtual QStringList extensions()const;
  virtual qSlicerIOOptions* options()const;

  virtual bool load(const IOProperties& properties);

protected:
  QScopedPointer<qSlicerSpatialObjectsReaderPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSpatialObjectsReader);
  Q_DISABLE_COPY(qSlicerSpatialObjectsReader);
};

#endif
