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


#ifndef __MetaObjectDocument_h
#define __MetaObjectDocument_h

#include "MetaDocument.h"
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

    bool  Read(const char * _fileName=NULL);

    bool  ReadStream(int _nDims, std::ifstream * _stream);

    bool  Write(const char * _fileName=NULL);

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

    const static char * LABEL_NOBJECTS;
    const static char * LABEL_TYPE;
    const static char * LABEL_NAME;
    const static char * LABEL_NUM_TRANS;
    const static char * LABEL_TRANSFORM;

    /** Label ID names of possible object types */
    const static char * ID_LABEL_BLOBTYPE;
    const static char * ID_LABEL_IMAGETYPE;
    const static char * ID_LABEL_SPATIALOBJTYPE;

    int m_NObjects;

    ObjectListType  m_objects;
    const unsigned int m_MaxNumTransforms;  //Maximum number of transforms

  private:
};

} // End namespace tube

#endif
