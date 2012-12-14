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

#ifndef __MetaDocument_H
#define __MetaDocument_H

#include <iostream>
#include <fstream>
#include <time.h>

#include <metaUtils.h>
#include <metaTypes.h>

extern int META_DEBUG;

namespace tube
{

class MetaDocument
{
  public:

    /** CTOR, DTOR */
    MetaDocument(void);
    MetaDocument(const char * _fileName);

    virtual ~MetaDocument(void);

    void  FileName(const char *_fileName);
    const char  * FileName(void) const;

    void  CopyInfo(const MetaDocument * _object);

    virtual bool  Read(const char * _fileName=NULL);

    virtual bool  Write(const char * _fileName=NULL);

    /** Writes image parameters to stdout */
    virtual void  PrintInfo(void) const;

    const char* DateLastModified(void) const;
    void DateLastModified(const char* _dateModified);

    /** Comment(...), Optional Field, Arbitrary String */
    const char  * Comment(void) const;
    void Comment(const char * _comment);

    /** Name(...), Optional Field, Name of the current MetaDocument */
    virtual void  Name(const char *_Name);
    virtual const char * Name(void) const;

    virtual void Clear(void);

    void ClearFields(void);

  protected:

    std::ifstream* m_ReadStream;
    std::ofstream* m_WriteStream;

    char m_Comment[255];
    char m_DateLastModified[255];
    char m_Name[255];
    char m_FileName[255];

    virtual void M_SetupReadFields(void);
    void M_PrepareNewReadStream(void);

    virtual void M_SetupWriteFields(void);

    virtual bool M_Read(void);

    virtual bool M_Write(void);

    std::vector<MET_FieldRecordType *> m_Fields;

  private:

};

} // Namespace tube end

#endif
