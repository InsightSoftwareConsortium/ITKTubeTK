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

#ifndef itktubeSpatialObjectDocument_h
#define itktubeSpatialObjectDocument_h

#include "itktubeObjectDocument.h"

namespace itk
{

namespace tube
{

/**
 * Encodes a spatial object file name and its ordered transform file names.
 * Spatial object documents store the file name of a spatial object and the file
 * names of the transforms that are to be applied consecutively to the spatial
 * object.
 *
 * \ingroup  ObjectDocuments
 */
class SpatialObjectDocument : public ObjectDocument
{
public:
  using Self = SpatialObjectDocument;
  using Superclass = ObjectDocument;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using TransformNameListType = Superclass::TransformNameListType;

  itkNewMacro(Self);
  itkOverrideGetNameOfClassMacro(SpatialObjectDocument);

protected:
  /** Constructor. */
  SpatialObjectDocument(void) { this->SetObjectType("SpatialObject"); }

  /** Destructor. */
  ~SpatialObjectDocument(void) override {}

private:
  // Copy constructor not implemented.
  SpatialObjectDocument(const Self & self);

  // Copy assignment operator not implemented.
  void
  operator=(const Self & self);

}; // End class SpatialObjectDocument

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSpatialObjectDocument_h )
