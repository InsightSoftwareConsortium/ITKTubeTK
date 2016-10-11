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
#ifndef __tubeComputeImageStatistics_hxx
#define __tubeComputeImageStatistics_hxx

#include "tubeComputeImageStatistics.h"

namespace tube
{

template< class TPixel, unsigned int VDimension >
ComputeImageStatistics< TPixel, VDimension >
::ComputeImageStatistics( void )
{
  m_Filter = FilterType::New();
}


/** Write statistics to a CSV formatted file */
template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::WriteCSVStatistics( std::string csvStatisticsFile ) const
{
  m_Filter->WriteCSVStatistics( csvStatisticsFile );
}


template< class TPixel, unsigned int VDimension >
void
ComputeImageStatistics< TPixel, VDimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  
  os << indent << m_Filter << std::endl;
}

} // End namespace tubetk

#endif
