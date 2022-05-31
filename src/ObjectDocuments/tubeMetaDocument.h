/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeMetaDocument_h
#define __tubeMetaDocument_h

#include "tubeObject.h"

#include <metaUtils.h>

#include <fstream>

namespace tube
{

/**
 * Encodes a meta document file name and its information.
 *
 * \ingroup  ObjectDocuments
 */
class MetaDocument : public Object
{
public:

  typedef MetaDocument                Self;
  typedef Object                      Superclass;
  typedef Self *                      Pointer;
  typedef const Self *                ConstPointer;

  typedef MET_FieldRecordType         FieldType;
  typedef std::vector< FieldType * >  FieldListType;

  /** Constructor. */
  MetaDocument( void );

  /**  Destructor. */
  virtual ~MetaDocument( void );

  /** Return the name of this class. */
  tubeTypeMacro( MetaDocument );

  /** Return the comment. */
  tubeGetStringMacro( Comment );

  /** Return the date modified. */
  tubeGetStringMacro( DateModified );

  /** Return the file name. */
  tubeGetStringMacro( FileName );

  /** Return the name. */
  tubeGetStringMacro( Name );

  /** Set the comment. */
  tubeSetStringMacro( Comment );

  /** Set the date modified. */
  tubeSetStringMacro( DateModified );

  /** Set the file name. */
  tubeSetStringMacro( FileName );

  /** Set the name. */
  tubeSetStringMacro( Name );

  /** Clear all the information. */
  virtual void Clear( void );

  /** Copy information from the specified meta document. */
  virtual void CopyInformation( const Self * self );

  /** Read the information from the specified file. */
  virtual bool Read( const std::string & fileName = "" );

  /** Write the information to the specified file. */
  virtual bool Write( const std::string & fileName = "" );

protected:

  /** Clear all the fields. */
  void ClearFields( void );

  /** Read the fields. */
  virtual bool ReadFields( void );

  /** Initialize the read fields. */
  virtual void SetupReadFields( void );

  /** Initialize the write fields. */
  virtual void SetupWriteFields( void );

  /** Write the fields. */
  virtual bool WriteFields( void );

  /** Print information about this object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

  std::ifstream  m_ReadStream;
  std::ofstream  m_WriteStream;
  FieldListType  m_FieldList;

private:

  // Copy constructor not implemented.
  MetaDocument( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  std::string  m_Comment;
  std::string  m_DateModified;
  std::string  m_FileName;
  std::string  m_Name;

}; // End class MetaDocument

} // End namespace tube

#endif // End !defined( __tubeMetaDocument_h )
