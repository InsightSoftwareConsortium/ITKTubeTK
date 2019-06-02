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

#ifndef __itktubeDocument_h
#define __itktubeDocument_h

#include <itkDataObject.h>
#include <itkObjectFactory.h>

namespace itk
{

namespace tube
{

/**
 * Allows text based documentation of meta information.
 *
 * \ingroup  ObjectDocuments
 */
class Document : public DataObject
{
public:

  typedef Document                    Self;
  typedef DataObject                  Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( Document, DataObject );

  /** Return the date modified. */
  itkGetStringMacro( DateModified );

  /** Set the date modified. */
  itkSetStringMacro( DateModified );

  /** Return the comment. */
  itkGetStringMacro( Comment );

  /** Set the comment. */
  itkSetStringMacro( Comment );

  /** Copy the information from the specified data object. */
  virtual void CopyInformation( const DataObject * itkNotUsed( data ) )
    {
    }

  /** Update the output information. */
  virtual void UpdateOutputInformation( void )
    {
    }

  /** Verify that the requested region is within the largest possible region,
      but note that this object does not use region information. */
  virtual bool VerifyRequestedRegion( void )
    {
    return true;
    }

  /** Determine whether the requested region is outside of the buffered region,
      but note that this object does not use region information. */
  virtual bool RequestedRegionIsOutsideOfTheBufferedRegion( void )
    {
    return false;
    }

  /** Set the requested region to match the requested region of the specified
      data object, but note that this object does not use region information. */
  virtual void SetRequestedRegion( const DataObject * itkNotUsed( data ) )
    {
    }

  /** Set the requested region to the largest possible region, but note that
      this object does not use region information. */
  virtual void SetRequestedRegionToLargestPossibleRegion( void )
    {
    }

protected:

  /** Constructor. */
  Document( void )
    {
    }

  /** Destructor. */
  virtual ~Document( void )
    {
    }

  /** Print information about the object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const override
    {
    this->Superclass::PrintSelf( os, indent );

    os << indent << "Comment:      " << m_Comment << std::endl;
    os << indent << "DateModified: " << m_DateModified << std::endl;
    }

private:

  // Copy constructor not implemented.
  Document( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  std::string  m_Comment;
  std::string  m_DateModified;

}; // End class Document

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeDocument_h )
