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

#include "itktubeBlobSpatialObjectDocument.h"
#include "itktubeImageDocument.h"
#include "tubeMetaDocument.h"

namespace tube
{

/**
 * \ingroup  ObjectDocuments
 */
class MetaObjectDocument : public MetaDocument
{
public:

  typedef MetaObjectDocument Self;
  typedef MetaDocument       Superclass;
  typedef Self *             Pointer;
  typedef const Self *       ConstPointer;

  typedef Superclass::FieldType                FieldType;
  typedef Superclass::FieldListType            FieldListType;

  typedef itk::tube::BlobSpatialObjectDocument BlobSpatialObjectDocumentType;
  typedef itk::tube::ImageDocument             ImageDocumentType;
  typedef itk::tube::ObjectDocument            ObjectDocumentType;
  typedef itk::tube::SpatialObjectDocument     SpatialObjectDocumentType;
  typedef std::vector< ObjectDocumentType::Pointer >
                                               ObjectDocumentListType;

  /** Constructor. */
  MetaObjectDocument( void );

  /** Destructor. */
  virtual ~MetaObjectDocument( void );

  /** Return the name of this class. */
  tubeTypeMacro( MetaObjectDocument );

  /** Return the maximum number of transforms. **/
  tubeGetMacro( MaximumNumberOfTransforms, unsigned int );

  /** Return the number of object documents. **/
  tubeGetMacro( NumberOfObjectDocuments, int );

  /** Return a reference to the list of object documents. **/
  virtual ObjectDocumentListType & GetObjectDocumentList( void );

  /** Set the list of object documents. **/
  virtual void SetObjectDocumentList( ObjectDocumentListType &
    objectDocumentList );

  /** Add the specified object document to the back of the list. */
  virtual void AddObjectDocument( ObjectDocumentType::Pointer
    objectDocument );

  /** Clear all the information. */
  virtual void Clear( void );

protected:

  /** Set the maximum number of transforms. **/
  tubeSetMacro( MaximumNumberOfTransforms, unsigned int );

  /** Set the number of object documents. **/
  tubeSetMacro( NumberOfObjectDocuments, int );

  /** Read the fields. */
  virtual bool ReadFields( void );

  /** Initialize the read fields for objects. */
  virtual void SetupObjectReadFields( void );

  /** Initialize the write fields for objects. */
  virtual void SetupObjectWriteFields( unsigned int index );

  /** Initialize the read fields. */
  virtual void SetupReadFields( void );

  /** Initialize the write fields. */
  virtual void SetupWriteFields( void );

  /** Write the fields. */
  virtual bool WriteFields( void );

  /** Print information about this object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  // Copy constructor not implemented.
  MetaObjectDocument( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  unsigned int            m_MaximumNumberOfTransforms;
  int                     m_NumberOfObjectDocuments;
  ObjectDocumentListType  m_ObjectDocumentList;

}; // End class MetaObjectDocument

} // End namespace tube

#endif // End !defined( __tubeMetaObjectDocument_h )
