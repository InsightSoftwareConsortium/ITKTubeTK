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

#include "tubeARFFParser.h"

#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <limits>
#include <iostream>

namespace tube
{

ARFFParser
::ARFFParser()
: m_Filename(""),
  m_XLabel("x"),
  m_YLabel("y"),
  m_ClassLabel("class"),
  m_XIndex(0),
  m_YIndex(0),
  m_ClassIndex(0),
  m_ClassNames(),
  m_MinX(std::numeric_limits<float>::max()),
  m_MinY(std::numeric_limits<float>::max()),
  m_MaxX(std::numeric_limits<float>::min()),
  m_MaxY(std::numeric_limits<float>::min()),
  m_ARFFData()
{
}

ARFFParser
::~ARFFParser()
{
  std::list<float*>::iterator itr;
  for( itr = m_ARFFData.begin(); itr != m_ARFFData.end(); ++itr )
    {
    delete [] *itr;
    }
}

void
ARFFParser
::SetFilename( const std::string& filename )
{
  m_Filename = filename;
}

void
ARFFParser
::SetXLabel( const std::string& xLabel )
{
  m_XLabel = xLabel;
}

void
ARFFParser
::SetYLabel( const std::string& yLabel )
{
  m_YLabel = yLabel;
}

void
ARFFParser
::SetClassLabel( const std::string& classLabel )
{
  m_ClassLabel = classLabel;
}

void
ARFFParser
::Parse()
{
  std::ifstream file( m_Filename.c_str() );
  std::string line;
  this->determineHeaderParameters( file );
  while( std::getline( file, line ) )
    {
    float* values = new float[3];
    if( line != "" && line.find( "%" ) == std::string::npos )
      {
      this->getValuesFromDataLine( line, values );
      this->adjustMinAndMaxBasedOnNewData( values );
      m_ARFFData.push_back( values );
      }
    else
      {
      delete [] values;
      }
    }
  
}

void 
ARFFParser
::determineHeaderParameters( std::ifstream& file )
{
  std::string line;
  unsigned int index = 0;
  while( std::getline( file, line ) )
    {
    if( line.find( "@DATA" ) != std::string::npos ||
        line.find( "@data" ) != std::string::npos )
      {
      break;
      }
    else if( line.find( "@ATTRIBUTE" ) != std::string::npos || 
             line.find( "@attribute" ) != std::string::npos )
      {
      std::string name;
      this->getAttributeName( line, name );
      if( name == this->m_XLabel )
        {
        m_XIndex = index;
        }
      else if( name == this->m_YLabel )
        {
        m_YIndex = index;
        }
      else if( name == this->m_ClassLabel )
        {
        this->determineClassificationsFromAttributeLine( line );
        m_ClassIndex = index;
        }
      else
        { // do nothing
        }
      ++index;
      }
    else
      { // do nothing
      }   
    }
}

void 
ARFFParser
::getAttributeName( const std::string& line, std::string& name ) const
{
  std::stringstream stringStream( line );
  std::getline( stringStream, name, ' ' );
  std::getline( stringStream, name, ' ' );
}

void 
ARFFParser
::getValuesFromDataLine( const std::string& line, float* values ) const
{
  std::stringstream stringStream( line );
  std::string cell;
  for( unsigned int i = 0; 
       i <= m_XIndex || i <= m_YIndex || i <= m_ClassIndex; 
       ++i )
    {
    std::getline( stringStream, cell, ',' );
    std::stringstream cellStream( cell );
    if( i == m_XIndex )
      {
      cellStream >> values[0];
      }
    else if( i == m_YIndex )
      {
      cellStream >> values[1];
      }
    else if( i == m_ClassIndex )
      {
      std::string name;
      cellStream >> name;
      values[2] = m_ClassNames.find(name)->second;
      }
    else
      { // do nothing
      }
    }
}

void
ARFFParser
::determineClassificationsFromAttributeLine( const std::string& line )
{
  size_t opening = line.find_first_of( '{' );
  ++opening;
  size_t closing = line.find_first_of( '}' );
  std::stringstream stringStream( line.substr( opening, closing-opening ) );
  std::string name;
  float counter = 0;
  while( std::getline( stringStream, name, ',' ) )
    {
    m_ClassNames[name] = counter;
    counter += 1;
    }
}

void 
ARFFParser
::adjustMinAndMaxBasedOnNewData( float* values )
{
  if( values[0] > m_MaxX )
    {
    m_MaxX = values[0];
    }

  if( values[0] < m_MinX )
    {
    m_MinX = values[0];
    }

  if( values[1] > m_MaxY )
    {
    m_MaxY = values[1];
    }

  if( values[1] < m_MinY )
    {
    m_MinY = values[1];
    }
}

const std::list<float*>& 
ARFFParser
::GetARFFData() const
{
  return m_ARFFData;
}

float 
ARFFParser
::GetMinX() const
{
  return m_MinX;
}

float 
ARFFParser
::GetMinY() const
{
  return m_MinY;
}

float 
ARFFParser
::GetMaxX() const
{
  return m_MaxX;
}

float 
ARFFParser
::GetMaxY() const
{
  return m_MaxY;
}

} // end namespace tube
