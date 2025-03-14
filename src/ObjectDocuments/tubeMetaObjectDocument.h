/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef tubeMetaObjectDocument_h
#define tubeMetaObjectDocument_h

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
  using Self = MetaObjectDocument;
  using Superclass = MetaDocument;
  using Pointer = Self *;
  using ConstPointer = const Self *;

  using FieldType = Superclass::FieldType;
  using FieldListType = Superclass::FieldListType;

  using BlobSpatialObjectDocumentType = itk::tube::BlobSpatialObjectDocument;
  using ImageDocumentType = itk::tube::ImageDocument;
  using ObjectDocumentType = itk::tube::ObjectDocument;
  using SpatialObjectDocumentType = itk::tube::SpatialObjectDocument;
  using ObjectDocumentListType = std::vector<ObjectDocumentType::Pointer>;

  /** Constructor. */
  MetaObjectDocument(void);

  /** Destructor. */
  ~MetaObjectDocument(void) override;

  /** Return the name of this class. */
  tubeTypeMacro(MetaObjectDocument);

  /** Return the maximum number of transforms. **/
  tubeGetMacro(MaximumNumberOfTransforms, unsigned int);

  /** Return the number of object documents. **/
  tubeGetMacro(NumberOfObjectDocuments, int);

  /** Return a reference to the list of object documents. **/
  virtual ObjectDocumentListType &
  GetObjectDocumentList(void);

  /** Set the list of object documents. **/
  virtual void
  SetObjectDocumentList(ObjectDocumentListType & objectDocumentList);

  /** Add the specified object document to the back of the list. */
  virtual void
  AddObjectDocument(ObjectDocumentType::Pointer objectDocument);

  /** Clear all the information. */
  void
  Clear(void) override;

protected:
  /** Set the maximum number of transforms. **/
  tubeSetMacro(MaximumNumberOfTransforms, unsigned int);

  /** Set the number of object documents. **/
  tubeSetMacro(NumberOfObjectDocuments, int);

  /** Read the fields. */
  bool
  ReadFields(void) override;

  /** Initialize the read fields for objects. */
  virtual void
  SetupObjectReadFields(void);

  /** Initialize the write fields for objects. */
  virtual void
  SetupObjectWriteFields(unsigned int index);

  /** Initialize the read fields. */
  void
  SetupReadFields(void) override;

  /** Initialize the write fields. */
  void
  SetupWriteFields(void) override;

  /** Write the fields. */
  bool
  WriteFields(void) override;

  /** Print information about this object. */
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // Copy constructor not implemented.
  MetaObjectDocument(const Self & self);

  // Copy assignment operator not implemented.
  void
  operator=(const Self & self);

  unsigned int           m_MaximumNumberOfTransforms;
  int                    m_NumberOfObjectDocuments;
  ObjectDocumentListType m_ObjectDocumentList;

}; // End class MetaObjectDocument

} // End namespace tube

#endif // End !defined( __tubeMetaObjectDocument_h )
