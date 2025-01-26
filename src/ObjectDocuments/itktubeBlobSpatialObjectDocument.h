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

#ifndef __itktubeBlobSpatialObjectDocument_h
#define __itktubeBlobSpatialObjectDocument_h

#include "itktubeSpatialObjectDocument.h"

namespace itk
{

namespace tube
{

/**
 * Encodes a blob spatial object file name and its ordered transform file names.
 * Blob spatial object documents store the file name of a blob spatial object
 * and the file names of the transforms that are to be applied consecutively to
 * the blob spatial object.
 *
 * \ingroup  ObjectDocuments
 */
class BlobSpatialObjectDocument : public SpatialObjectDocument
{
public:

  typedef BlobSpatialObjectDocument          Self;
  typedef SpatialObjectDocument              Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  typedef Superclass::TransformNameListType  TransformNameListType;

  itkNewMacro( Self );
  itkOverrideGetNameOfClassMacro( BlobSpatialObjectDocument);

protected:

  /** Constructor. */
  BlobSpatialObjectDocument( void )
    {
    this->SetObjectType( "Blob" );
    }

  /** Destructor. */
  virtual ~BlobSpatialObjectDocument( void )
    {
    }

private:

  // Copy constructor not implemented.
  BlobSpatialObjectDocument( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class BlobSpatialObjectDocument

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeBlobSpatialObjectDocument_h )
