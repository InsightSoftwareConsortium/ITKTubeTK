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

#ifndef __tubeStringUtilities_h
#define __tubeStringUtilities_h

#include <sstream>
#include <string>
#include <vector>

namespace tube
{

template <class T>
void StringToVector( const std::string & s, std::vector< T > & vec )
{
  vec.clear();

  std::string::size_type prevPos = 0;
  std::string::size_type pos = s.find( ",", prevPos );
  T tVal;
  while( pos != std::string::npos )
    {
    std::stringstream ss;
    // stringstream use whitespaces as delimiters, we don't want that
    // behavior.
    std::string substr = s.substr( prevPos, pos-prevPos );
    while( substr[0] == ' ' )
      {
      substr = substr.substr( 1, substr.size()-1 );
      }
    while( substr[substr.size()-1] == ' ' )
      {
      substr = substr.substr( 0, substr.size()-1 );
      }
    std::string::size_type whitespacePos = substr.find( ' ' );
    while( whitespacePos != std::string::npos )
      {
      substr[whitespacePos] = '_';
      whitespacePos = substr.find( ' ' );
      }
    ss << substr;
    ss >> tVal;
    vec.push_back( tVal );

    prevPos = pos+1;
    pos = s.find( ",", prevPos );
    }
  std::string substr = s.substr( prevPos, s.size()-prevPos );
  while( substr[0] == ' ' )
    {
    substr = substr.substr( 1, substr.size()-1 );
    }
  while( substr[substr.size()-1] == ' ' )
    {
    substr = substr.substr( 0, substr.size()-1 );
    }
  std::string::size_type whitespacePos = substr.find( ' ' );
  while( whitespacePos != std::string::npos )
    {
    substr[whitespacePos] = '_';
    whitespacePos = substr.find( ' ' );
    }
  std::stringstream ss;
  ss << substr;
  ss >> tVal;
  vec.push_back( tVal );
}

} // End namespace tube

#endif // End !defined(__tubeStringUtilities_h)
