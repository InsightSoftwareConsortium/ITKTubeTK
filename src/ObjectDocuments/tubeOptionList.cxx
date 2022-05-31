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

/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: tubeOptionList.cxx,v $
Language:  C++
Date:      $Date: 2004/05/06 13:45:52 $
Version:   $Revision: 1.1 $

Copyright ( c ) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or https://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "tubeOptionList.h"

namespace tube
{

OptionList
::OptionList( void )
{
}


OptionList
::~OptionList( void )
{
}


OptionList::OptionMapType &
OptionList
::GetOptionMap( void )
{
  return m_OptionMap;
}


void
OptionList
::CreateOptionMap( int argc, char * argv[] )
{
  std::string tag;
  int index = 1;

  while( index < argc )
    {
    if( argv[index][0] == '-' && argv[index][1] == '-' )
      {
      tag = argv[index];
      tag = tag.erase( 0, 2 ); // Remove '--'.
      }
    else
      {
      m_OptionMap.insert( std::make_pair( tag, argv[index] ) );
      }

    ++index;
    }
}


int
OptionList
::GetOptions( const std::string & tag,
                            std::vector< std::string > & values ) const
{
  values.clear();

  typedef OptionMapType::const_iterator OptionIterator;

  std::pair< OptionIterator, OptionIterator > bound
    = m_OptionMap.equal_range( tag );
  int numberOfValues = 0;

  for( OptionIterator iter = bound.first; iter != bound.second; ++iter )
    {
    values.push_back( iter->second );
    ++numberOfValues;
    }

  return numberOfValues;
}


int
OptionList
::DumpOptions( const std::string & tag, bool withTag, bool withNewLine ) const
{
  typedef OptionMapType::const_iterator OptionIterator;

  std::pair< OptionIterator, OptionIterator > bound
    = m_OptionMap.equal_range( tag );

  if( bound.first != bound.second )
    {
    if( withTag )
      {
      std::cout << "--" << tag << " ";
      }

    int numberOfValues = 0;

    for( OptionIterator iter = bound.first; iter != bound.second; ++iter )
      {
      std::cout << iter->second << " ";
      ++numberOfValues;
      }

    if( withNewLine )
      {
      std::cout << std::endl;
      }

    return numberOfValues++;
    }

  return 0;
}


bool
OptionList
::GetBooleanOption( const std::string & tag,
  bool defaultValue,
  bool required ) const
{
  std::vector< std::string > values;
  const int numberOfValues = this->GetOptions( tag, values );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return defaultValue;
    }

  if( values[0] == "yes" )
    {
      return true;
    }

  return false;
}


unsigned char
OptionList
::GetCharacterOption( const std::string & tag,
  unsigned char defaultValue,
  bool required ) const
{
  std::vector< std::string > values;
  const int numberOfValues = this->GetOptions( tag, values );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return defaultValue;
    }

  return ( unsigned char )std::atoi( values[0].c_str() );
}


double
OptionList
::GetDoubleOption( const std::string & tag,
  double defaultValue,
  bool required ) const
{
  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return defaultValue;
    }

  return std::strtod( tempValues[0].c_str(), NULL );
}


int
OptionList
::GetIntegerOption( const std::string & tag,
  int defaultValue,
  bool required ) const
{
  std::vector< std::string > values;
  const int numberOfValues = this->GetOptions( tag, values );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return defaultValue;
    }

  return ( int )std::strtol( values[0].c_str(), NULL, 0 );
}


int
OptionList
::GetMultipleCharactersOption( const std::string & tag,
  std::vector< unsigned char > & values,
  bool required ) const
{
  values.clear();

  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  if( tempValues[0] == "-" )
    {
    return -2;
    }

  for( int i = 0; i < numberOfValues; ++i )
    {
    values.push_back( ( unsigned char )std::atoi( tempValues[i].c_str() ) );
    }

  return numberOfValues;
}


int
OptionList
::GetMultipleDoublesOption( const std::string & tag,
  itk::Array< double > & values,
  bool required ) const
{
  std::vector< double > tempValues;
  const int numberOfValues = this->GetMultipleDoublesOption( tag, tempValues,
                                                             required );

  if( numberOfValues <= 0 )
    {
    return numberOfValues;
    }

  itk::Array< double > array( numberOfValues );

  for( int i = 0; i < numberOfValues; ++i )
    {
    array[i] = tempValues[i];
    }

  values = array;

  return numberOfValues;
}


int
OptionList
::GetMultipleDoublesOption( const std::string & tag,
  std::vector< double > & values,
  bool required ) const
{
  values.clear();

  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  if( tempValues[0] == "-" )
    {
    return -2;
    }

  for( int i = 0; i < numberOfValues; ++i )
    {
    values.push_back( std::strtod( tempValues[i].c_str(), NULL ) );
    }

  return numberOfValues;
}


int
OptionList
::GetMultipleIntegersOption( const std::string & tag,
  std::vector< int > & values,
  bool required ) const
{
  values.clear();

  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  if( tempValues[0] == "-" )
    {
    return -2;
    }

  for( int i = 0; i < numberOfValues; ++i )
    {
    values.push_back( ( int )std::strtol( tempValues[i].c_str(), NULL, 0 ) );
    }

  return numberOfValues;
}


int
OptionList
::GetMultipleStringsOption( const std::string & tag,
  std::vector< std::string > & values,
  bool required ) const
{
  values.clear();

  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  for( int i = 0; i < numberOfValues; ++i )
    {
    values.push_back( tempValues[i] );
    }

  return numberOfValues;
}


int
OptionList
::GetMultipleUnsignedIntegersOption( const std::string & tag,
  std::vector< unsigned int > & values,
  bool required ) const
{
  values.clear();

  std::vector< std::string > tempValues;
  const int numberOfValues = this->GetOptions( tag, tempValues );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  if( tempValues[0] == "-" )
    {
    return -2;
    }

  for( int i = 0; i < numberOfValues; ++i )
    {
    values.push_back( ( unsigned int )std::strtoul( tempValues[i].c_str(), NULL,
                                                  0 ) );
    }

  return numberOfValues;
}


int
OptionList
::GetStringOption( const std::string & tag,
  std::string & value,
  bool required ) const
{
  std::vector< std::string > values;
  const int numberOfValues = this->GetOptions( tag, values );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return -1;
    }

  value = values[0];

  return numberOfValues;
}


unsigned int
OptionList
::GetUnsignedIntegerOption( const std::string & tag,
  unsigned int defaultValue,
  bool required ) const
{
  std::vector< std::string > values;
  const int numberOfValues = this->GetOptions( tag, values );

  if( required && numberOfValues == 0 )
    {
    throw RequiredOptionMissing( tag );
    }

  if( numberOfValues == 0 )
    {
    return defaultValue;
    }

  return ( unsigned int )std::strtoul( values[0].c_str(), NULL, 0 );
}


void
OptionList
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "OptionMap:" << std::endl;

  for( OptionMapType::const_iterator it = m_OptionMap.begin();
       it != m_OptionMap.end(); ++it )
    {
    os << indent << it->first << " -> " << it->second << std::endl;
    }
}

} // End namespace tube
