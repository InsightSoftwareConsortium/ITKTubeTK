/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeImageDocument_h
#define __itktubeImageDocument_h

#include "itktubeObjectDocument.h"

namespace itk
{

namespace tube
{

/** \class ImageDocument
 * \brief Encodes an image file name and its ordered transform file names
 *
 *  Image Documents will store the file name of an image
 *    and set file names for the transforms that are to be applied consecutively to the image.
 *
 *  IO is done through MetaObjectDocument.h
 *
 *  \ingroup Document
 */
class ImageDocument : public ObjectDocument
{
public:

  typedef ImageDocument   Self;
  typedef ObjectDocument  Superclass;

  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  typedef Superclass::DateType      DateType;
  typedef Superclass::CommentsType  CommentsType;

  //*** Not Implemented, but would allow for Document objects to be held by other documents
  typedef Superclass::ChildrenListType    ChildrenListType;
  typedef Superclass::ChildrenListPointer ChildrenListPointer;

  /** list that holds the ordered transform Names */
  typedef Superclass::TransformNameListType TransformNameListType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( Self, Superclass );

  /** Return the type of the object within the Document (ie. "Image") */
  std::string GetObjectType( void ) const { return LABEL_IMAGETYPE; }

  ~ImageDocument( void ) {}

protected:

  ImageDocument( void ) : LABEL_IMAGETYPE("Image") {}

  const std::string LABEL_IMAGETYPE;

}; // End class ImageDocument

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeImageDocument_h)
