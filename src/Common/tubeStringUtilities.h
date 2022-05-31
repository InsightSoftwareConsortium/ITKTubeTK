/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#ifndef __tubeStringUtilities_h
#define __tubeStringUtilities_h

#include "tubeMacro.h"

#include <sstream>
#include <vector>

namespace tube
{

template< class T >
void StringToVector( const std::string & s, std::vector< T > & vec,
  std::string sep = "," )
{
  vec.clear();

  std::string::size_type prevPos = 0;
  std::string::size_type pos = s.find( sep, prevPos );
  T tVal;
  while( pos != std::string::npos )
    {
    // stringstream use whitespace as a delimiter, we don't want that
    // behavior.
    std::string substr = s.substr( prevPos, pos-prevPos );
    while( isspace( substr[ 0 ] ) )
      {
      substr = substr.substr( 1, substr.size()-1 );
      }
    while( isspace( substr[ substr.size()-1 ] ) )
      {
      substr = substr.substr( 0, substr.size()-1 );
      }
    for( std::string::size_type checkspace = 0; checkspace < substr.size();
      ++checkspace )
      {
      if( isspace( substr[ checkspace ] ) )
        {
        substr[ checkspace ] = '_';
        }
      }
    std::stringstream ss;
    ss << substr;
    ss >> tVal;
    vec.push_back( tVal );

    prevPos = pos+1;
    while( isspace( s[ prevPos ] ) && prevPos < s.size() )
      {
      ++prevPos;
      }
    pos = s.find( sep, prevPos );
    }

  std::string substr = s.substr( prevPos, s.size()-prevPos );
  while( isspace( substr[ 0 ] ) )
    {
    substr = substr.substr( 1, substr.size()-1 );
    }
  while( isspace( substr[ substr.size()-1 ] ) )
    {
    substr = substr.substr( 0, substr.size()-1 );
    }
  for( std::string::size_type checkspace = 0; checkspace < substr.size();
    ++checkspace )
    {
    if( isspace( substr[ checkspace ] ) )
      {
      substr[ checkspace ] = '_';
      }
    }
  std::stringstream ss;
  ss << substr;
  ss >> tVal;
  vec.push_back( tVal );
}

} // End namespace tube

#endif // End !defined( __tubeStringUtilities_h )
