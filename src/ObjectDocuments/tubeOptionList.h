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

/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: tubeOptionList.h,v $
Language:  C++
Date:      $Date: 2004/05/06 13:45:52 $
Version:   $Revision: 1.1 $

Copyright ( c ) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __tubeOptionList_h
#define __tubeOptionList_h

#include "tubeObject.h"

#include <itkArray.h>

#include <map>
#include <vector>

namespace tube
{

/**
 * Maps arguments to options consisting of tags and values.
 *
 * \ingroup  ObjectDocuments
 */
class OptionList : public Object
{
public:

  typedef OptionList                                 Self;
  typedef Object                                     Superclass;
  typedef Self *                                     Pointer;
  typedef const Self *                               ConstPointer;

  typedef std::multimap< std::string, std::string >  OptionMapType;

  /** Encodes a required option that is missing.  */
  class RequiredOptionMissing
    {
  public:

    /** Constructor. */
    RequiredOptionMissing( const std::string & tag ) : m_Tag( tag )
      {
      }

    /** Return the name of this class. */
    tubeTypeMacroNoOverride( RequiredOptionMissing );

    /** Return the tag of the missing option. */
    tubeGetStringMacro( Tag );

  private:

    std::string m_Tag;

    }; // End class RequiredOptionMissing

  /** Constructor. */
  OptionList( void );

  /** Destructor. */
  virtual ~OptionList( void );

  /** Return the name of this class. */
  tubeTypeMacro( OptionList );

  /** Create the option map from the specified arguments. */
  virtual void CreateOptionMap( int argc, char * argv[] );

  /** Return the values of the specified tag. */
  virtual int GetOptions( const std::string & tag,
    std::vector< std::string > & values ) const;

  /** Dump the values of the specified tag to the standard output. */
  virtual int DumpOptions( const std::string & tag, bool withTag = true,
    bool withNewLine = false ) const;

  /** Return the value of the specified tag as a boolean. */
  virtual bool GetBooleanOption( const std::string & tag, bool defaultValue,
    bool required ) const;

  /** Return the value of the specified tag as a double-precision float. */
  virtual double GetDoubleOption( const std::string & tag,
    double defaultValue, bool required ) const;

  /** Return the value of the specified tag as a signed integer. */
  virtual int GetIntegerOption( const std::string & tag, int defaultValue,
    bool required ) const;

  /** Return the values of the specified tag as a list of unsigned
   * characters. */
  virtual int GetMultipleCharactersOption( const std::string & tag,
    std::vector< unsigned char > & values, bool required ) const;

  /** Return the values of the specified tag as an ITK array of
   * double-precision floats. */
  virtual int GetMultipleDoublesOption( const std::string & tag,
    itk::Array< double > & values, bool required ) const;

  /** Return the values of the specified tag as a list of double-precision
   * floats. */
  virtual int GetMultipleDoublesOption( const std::string & tag,
    std::vector< double > & values, bool required ) const;

  /** Return the values of the specified tag as a list of signed
   * integers. */
  virtual int GetMultipleIntegersOption( const std::string & tag,
    std::vector< int > & values, bool required ) const;

  /** Return the values of the specified tag as a list of strings. */
  virtual int GetMultipleStringsOption( const std::string & tag,
    std::vector< std::string > & values, bool required ) const;

  /** Return the values of the specified tag as a list of unsigned
   * integers. */
  virtual int GetMultipleUnsignedIntegersOption( const std::string & tag,
    std::vector< unsigned int > & values, bool required ) const;

  /** Return the value of the specified tag as a string. */
  virtual int GetStringOption( const std::string & tag, std::string & value,
    bool required ) const;

  /** Return the value of the specified tag as a unsigned character. */
  virtual unsigned char GetCharacterOption( const std::string & tag,
    unsigned char defaultValue, bool required ) const;

  /** Return the value of the specified tag as a unsigned integer. */
  virtual unsigned int GetUnsignedIntegerOption( const std::string & tag,
    unsigned int defaultValue, bool required ) const;

protected:

  /** Return a reference to the option map. */
  virtual OptionMapType & GetOptionMap( void );

  /** Print information about this object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Copy constructor not implemented.
  OptionList( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  OptionMapType m_OptionMap;

}; // End class OptionList

} // End namespace tube

#endif // End !defined( __tubeOptionList_h )
