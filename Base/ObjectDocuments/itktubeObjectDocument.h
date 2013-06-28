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

#ifndef __itktubeObjectDocument_h
#define __itktubeObjectDocument_h

#include "itktubeDocument.h"

namespace itk
{

namespace tube
{


/** \class ObjectDocument
 * \brief Encodes an object file name and its ordered transform file names
 *
 *  Object Documents will store the file name of an object type (eg. image, spatial object, etc)
 *     and a set file names for the transforms that are to be applied consecutively to the object.
 *
 *  IO is done through MetaObjectDocument.h
 *
 *  \ingroup Document
 */
class ObjectDocument : public Document
{
public:

  typedef ObjectDocument          Self;
  typedef Document                Superclass;

  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  typedef Superclass::DateType                  DateType;
  typedef Superclass::CommentsType              CommentsType;
  typedef std::string                           ObjectNameType;
  typedef std::string                           TransformNameType;

  /** Not Implemented, but would allow for Document objects to be held by other documents */
  typedef Superclass::ChildrenListType          ChildrenListType;
  typedef Superclass::ChildrenListPointer       ChildrenListPointer;

  /** list that holds the ordered transform Names */
  typedef std::vector< TransformNameType >      TransformNameListType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ObjectDocument, Document );

  /** To be implemented by the object type that inherits this class */
  virtual std::string  GetObjectType( void ) const{ return "Object";}

  /** Get the Object file name -- default is undefined */
  itkGetConstReferenceMacro( ObjectName, ObjectNameType );

  /** Set the Object file name */
  itkSetMacro( ObjectName, ObjectNameType );

  /** Get the number of transforms associated with the object */
  unsigned int GetNumberOfTransforms( void ) const
    {
    return static_cast<unsigned int>(m_TransformList.size());
    }

  /** Get a std::vector of all the transform file names in order */
  TransformNameListType GetTransformNames( void ) const { return m_TransformList; }

  /** Add a transform name to the end of the transform list */
  void AddTransformNameToBack( const std::string & trans ) { m_TransformList.push_back( trans ); }

  /** Remove last transform from the list -- Does nothing if there are no transforms */
  void RemoveTransformNameFromBack( void )
    {
    if( !m_TransformList.empty() )
      {
      m_TransformList.pop_back();
      }
    }

  void Print( std::ostream& os ) const
    {
    os << "ObjectName: " << m_ObjectName << std::endl;

    os<< "Transform List: " << std::endl;
    TransformNameListType::const_iterator it = m_TransformList.begin();
    while( it != m_TransformList.end() )
      {
      std::cout << *it << std::endl;
      ++it;
      }
    }

  ~ObjectDocument( void ) {}

protected:

  ObjectDocument( void ) {}

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    Superclass::PrintSelf(os,indent);

    os << indent << "ObjectName: " << m_ObjectName << std::endl;

    os << indent << "Transform List: " << std::endl;
    TransformNameListType::const_iterator it = m_TransformList.begin();
    while( it != m_TransformList.end() )
      {
      std::cout << *it << std::endl;
      ++it;
      }
    }

private:

  ObjectNameType                    m_ObjectName;
  TransformNameListType             m_TransformList;

}; // End class ObjectDocument

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeObjectDocument_h)
