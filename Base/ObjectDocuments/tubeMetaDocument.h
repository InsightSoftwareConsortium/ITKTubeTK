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

#ifndef __tubeMetaDocument_h
#define __tubeMetaDocument_h

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
    MetaDocument(const std::string & _fileName);

    virtual ~MetaDocument(void);

    void  FileName(const std::string & _fileName);
    std::string FileName(void) const;

    void  CopyInfo(const MetaDocument * _object);

    virtual bool  Read(const std::string & _fileName = std::string());

    virtual bool  Write(const std::string & _fileName = std::string());

    /** Writes image parameters to stdout */
    virtual void  PrintInfo(void) const;

    std::string DateLastModified(void) const;
    void DateLastModified(const std::string & _dateModified);

    /** Comment(...), Optional Field, Arbitrary String */
    std::string Comment(void) const;
    void Comment(const std::string & _comment);

    /** Name(...), Optional Field, Name of the current MetaDocument */
    virtual void  Name(const std::string & _Name);
    virtual std::string Name(void) const;

    virtual void Clear(void);

    void ClearFields(void);

  protected:

    std::ifstream m_ReadStream;
    std::ofstream m_WriteStream;

    std::string m_Comment;
    std::string m_DateLastModified;
    std::string m_Name;
    std::string m_FileName;

    virtual void M_SetupReadFields(void);
    void M_PrepareNewReadStream(void);

    virtual void M_SetupWriteFields(void);
    void M_PrepareNewWriteStream(void);

    virtual bool M_Read(void);

    virtual bool M_Write(void);

    std::vector<MET_FieldRecordType *> m_Fields;

  private:

};

} // Namespace tube end

#endif
