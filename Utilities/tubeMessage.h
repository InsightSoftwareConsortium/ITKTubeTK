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

#ifndef __tubeMessage_h
#define __tubeMessage_h

#include <iostream>

namespace tube {

namespace MessageLevel {
  enum { Information, Warning, Error, Debug };
}

template <class T>
void Message( const T& str, int level = 0 )
{
  #ifndef NDEBUG
  if( level == MessageLevel::Debug )
    {
    std::cout << "<debug>" << str << "</debug>" << std::endl;
    return;
    }
  #endif
  switch( level )
    {
    case MessageLevel::Information:
      std::cout << "<info>" << str << "</info>" << std::endl;
      break;
    case MessageLevel::Warning:
      std::cout << "<warning>" << str << "</warning>" << std::endl;
      break;
    case MessageLevel::Error:
      std::cout << "<error>" << str << "</error>" << std::endl;
      break;
    default:
      break;
    }
}

template <class T>
void InfoMessage( const T& str )
{
  Message( str, MessageLevel::Information );
}

template <class T>
void InformationMessage( const T& str )
{
  Message( str, MessageLevel::Information );
}

template <class T>
void WarningMessage( const T& str )
{
  Message( str, MessageLevel::Warning );
}

template <class T>
void ErrorMessage( const T& str )
{
  Message( str, MessageLevel::Error );
}

template <class T>
void DebugMessage( const T& str )
{
  Message( str, MessageLevel::Debug );
}

};

#endif
