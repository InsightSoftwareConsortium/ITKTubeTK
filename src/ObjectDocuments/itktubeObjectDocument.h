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

#ifndef __itktubeObjectDocument_h
#define __itktubeObjectDocument_h

#include "itktubeDocument.h"

namespace itk
{

namespace tube
{

/**
 * Encodes an object file name and its ordered transform file names.
 *
 * \ingroup  ObjectDocuments
 */
class ObjectDocument : public Document
{
public:

  typedef ObjectDocument                   Self;
  typedef Document                         Superclass;
  typedef SmartPointer< Self >             Pointer;
  typedef SmartPointer< const Self >       ConstPointer;

  typedef std::vector< std::string >       TransformNameListType;

  itkNewMacro( Self );
  itkOverrideGetNameOfClassMacro( ObjectDocument);

  /** Return the object type. */
  itkGetStringMacro( ObjectType );

  /** Return the number of transform names. */
  virtual unsigned int GetNumberOfTransformNames( void ) const
    {
    return static_cast< unsigned int >( m_TransformNameList.size() );
    }

  /** Return the list of transform names. */
  virtual TransformNameListType GetTransformNames( void ) const
    {
    return m_TransformNameList;
    }

  /** Add the specified transform name to the back of the list. */
  virtual void AddTransformNameToBack( const std::string & transformName )
    {
    m_TransformNameList.push_back( transformName );
    }

  /** Remove the transform name from the back of the list. */
  virtual void RemoveTransformNameFromBack( void )
    {
    if( !m_TransformNameList.empty() )
      {
      m_TransformNameList.pop_back();
      }
    }

protected:

  /** Constructor. */
  ObjectDocument( void )
    {
    this->SetObjectType( "Object" );
    }

  /** Destructor. */
  virtual ~ObjectDocument( void )
    {
    }

  /** Set the object type. */
  itkSetMacro( ObjectType, std::string );

  /** Print information about the object. */
  virtual void PrintSelf( std::ostream& os, Indent indent ) const override
    {
    Superclass::PrintSelf( os, indent );

    os << indent << "ObjectType:        " << m_ObjectType << std::endl;
    os << indent << "TransformNameList:" << std::endl;

    for( TransformNameListType::const_iterator it = m_TransformNameList.begin();
         it != m_TransformNameList.end(); ++it )
      {
      std::cout << indent << *it << std::endl;
      }
    }

private:

  // Copy constructor not implemented.
  ObjectDocument( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  std::string            m_ObjectType;
  TransformNameListType  m_TransformNameList;

}; // End class ObjectDocument

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeObjectDocument_h )
