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


#ifndef __tubeMetaObjectDocument_h
#define __tubeMetaObjectDocument_h

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "tubeMetaDocument.h"
#include "itkImageDocument.h"
#include "itkBlobSpatialObjectDocument.h"

namespace tube
{

class MetaObjectDocument : public MetaDocument
{
public:

  typedef itk::tube::ObjectDocument                   ObjectDocumentType;
  typedef itk::tube::BlobSpatialObjectDocument        BlobSpatialObjectDocumentType;
  typedef itk::tube::SpatialObjectDocument            SpatialObjectDocumentType;
  typedef itk::tube::ImageDocument                    ImageDocumentType;
  typedef std::vector<ObjectDocumentType::Pointer>    ObjectListType;

  MetaObjectDocument();
  ~MetaObjectDocument();

  void  PrintInfo() const;

  bool  Read(const std::string & _fileName = std::string());

  bool  ReadStream(int _nDims, std::ifstream & _stream);

  bool  Write(const std::string & _fileName = std::string());

  /** Clear tube information */
  void Clear();

  void AddObject( ObjectDocumentType::Pointer );

  /** Overrides any previously added objects */
  void SetObjectList( ObjectListType& list );

  ObjectListType * GetObjectList(void);


protected:
  bool M_Write(void);
  bool M_Read(void);

  void M_SetupReadFields(void);
  void M_SetupObjectReadFields(void);

  void M_SetupWriteFields(void);
  void M_SetupObjectWriteFields(unsigned int);

  static const std::string LABEL_NOBJECTS;
  static const std::string LABEL_TYPE;
  static const std::string LABEL_NAME;
  static const std::string LABEL_NUM_TRANS;
  static const std::string LABEL_TRANSFORM;

  /** Label ID names of possible object types */
  static const std::string ID_LABEL_BLOBTYPE;
  static const std::string ID_LABEL_IMAGETYPE;
  static const std::string ID_LABEL_SPATIALOBJTYPE;

  int                      m_NObjects;

  ObjectListType           m_objects;
  const unsigned int       m_MaxNumTransforms;  //Maximum number of transforms

private:
};

} // End namespace tube

#endif
