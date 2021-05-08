/*=========================================================================

  Copyright 2014 Kitware, Inc. All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   JSON Serialization Utilities
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __JsonSerializationUtilities_h
#define __JsonSerializationUtilities_h

#include <json/json.h>

namespace {

//----------------------------------------------------------------------------
template <typename T>
Json::Value JsonSerialize( const T & value )
{ return value; }

//----------------------------------------------------------------------------
template <typename T>
Json::Value JsonSerialize( const std::vector<T> & value )
{
  Json::Value array( Json::arrayValue );
  const size_t k = value.size();
  for( size_t i = 0; i < k; ++i )
    {
    array.append( JsonSerialize( value[i] ) );
    }
  return array;
}

}

#endif
