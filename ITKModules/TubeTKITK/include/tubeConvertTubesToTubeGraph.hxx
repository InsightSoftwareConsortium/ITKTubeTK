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
#ifndef __tubeConvertTubesToTubeGraph_hxx
#define __tubeConvertTubesToTubeGraph_hxx

#include "tubeConvertTubesToTubeGraph.h"

namespace tube
{

template< class TPixel, unsigned int Dimension >
ConvertTubesToTubeGraph< TPixel, Dimension >
::ConvertTubesToTubeGraph( void )
{
  m_Filter = FilterType::New();
}

template< class TPixel, unsigned int Dimension >
vnl_matrix< double >
ConvertTubesToTubeGraph< TPixel, Dimension >
::GetAdjacencyMatrix( void )
{
  return m_Filter->GetAdjacencyMatrix();
}

template< class TPixel, unsigned int Dimension >
vnl_vector< int >
ConvertTubesToTubeGraph< TPixel, Dimension >
::GetRootNodes( void )
{
  return m_Filter->GetRootNodes();
}

template< class TPixel, unsigned int Dimension >
vnl_vector< double >
ConvertTubesToTubeGraph< TPixel, Dimension >
::GetBranchNodes( void )
{
  return m_Filter->GetBranchNodes();
}

template< class TPixel, unsigned int Dimension >
void
ConvertTubesToTubeGraph< TPixel, Dimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  m_filter->PrintSelf();
}

}

#endif
