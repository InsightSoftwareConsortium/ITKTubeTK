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
  m_MinX(0),
  m_MinY(0),
  m_MaxX(0),
  m_MaxY(0),
  m_ARFFData()
{
}

ARFFParser
::~ARFFParser()
{
  std::list<float*>::iterator itr;
  for( itr = m_ARFFData.begin(); itr != m_ARFFData.end(); ++itr )
    {
    delete *itr;
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
  while( std::getline( file, line ) )
    {
    
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
    if( line.find( "@DATA" ) != std::string::npos )
      {
      break;
      }
    if( line.find( "@ATTRIBUTE" ) != std::string::npos )
      {
      ++index;
      }
   
   }
}

const std::list<float*>& 
ARFFParser
::GetARFFData() const
{
  return m_ARFFData;
}

} // end namespace tube
