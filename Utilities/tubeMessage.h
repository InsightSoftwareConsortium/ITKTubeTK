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
#include <ostream>
#include <string>

#include <cstdarg>
#include <cstring>
#include <cstdio>


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

/** Formatted logging, as with printf(...) */
inline void FmtMessage(int level, const char* fmt, va_list args )
{
	std::ostream & os = level >= MessageLevel::Information ? std::cerr : std::cout;

	char text[1024];
	memset(text, 0, sizeof(char)*1023);
	vsnprintf(text, 1023, fmt, args);

  #ifndef NDEBUG
  if( level == MessageLevel::Debug )
    {
    os << "<debug>" << text << "</debug>" << std::endl;
    return;
    }
  #endif

  switch( level )
    {
    case MessageLevel::Information:
      os << "<info>" << text << "</info>" << std::endl;
      break;
    case MessageLevel::Warning:
      os << "<warning>" << text << "</warning>" << std::endl;
      break;
    case MessageLevel::Error:
      os << "<error>" << text << "</error>" << std::endl;
      break;
    default:
      break;
    }
  os.flush();
}

inline void FmtInfoMessage(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    FmtMessage( MessageLevel::Information, fmt, args );
    va_end(args);
}

inline void FmtInformationMessage(const char* fmt, ...)
{
  va_list args;
  va_start( args, fmt );
  FmtMessage( MessageLevel::Information, fmt, args );
  va_end( args );

}

inline void FmtWarningMessage(const char* fmt, ...)
{
  va_list args;
  va_start( args, fmt );
  FmtMessage( MessageLevel::Warning, fmt, args );
  va_end( args );
}

inline void FmtErrorMessage(const char* fmt, ...)
{
  va_list args;
  va_start( args, fmt );
  FmtMessage( MessageLevel::Error, fmt, args );
  va_end( args );
}

inline void FmtDebugMessage(const char* fmt, ...)
{
  va_list args;
  va_start( args, fmt );
  FmtMessage( MessageLevel::Debug, fmt, args );
  va_end( args );
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
