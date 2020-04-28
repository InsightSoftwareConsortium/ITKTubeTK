/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeShortestPathKernel_h
#define __tubeShortestPathKernel_h

#include "GraphKernel.h"

#include <algorithm>

namespace tube
{

/**
 * \brief Implementation of Borgwardt's Shortest-Path Kernel
 *
 * This class implements the shortest-path kernel, proposed in
 *
 * [1] K.M. Borgwardt and H.P. Kriegel, "Shortest-Path Kernels on
 *     Graphs", In: IEEE Int. Conf. on Data Mining, 2005
 */
class ShortestPathKernel : public GraphKernel
{
public:

  /** Edge kernel types */
  static const int EDGE_KERNEL_DEL = 0;

  /** CTOR - Consumer sets graphs */
  ShortestPathKernel( const GraphType & g0, const GraphType & g1 )
    : GraphKernel( g0, g1 ), m_EdgeKernelType( EDGE_KERNEL_DEL )
    {
    }

  /** Sets edge-kernel type */
  void SetEdgeKernel( int edgeKernelType )
    {
    m_EdgeKernelType = edgeKernelType;
    }

  /** Computes the SP kernel value, see [1], Section 4.2 */
  double Compute( void );

private:

  /** Computes a Floyd-transformed graph, see [1], Section 4.1 */
  GraphType FloydTransform( const GraphType & in );

  void ensureOrder( int & first, int & second )
    {
    if( first > second )
      {
      std::swap( first, second );
      }
    }

  /** Floyd-transformed graphs */
  GraphType  m_FG0;
  GraphType  m_FG1;
  int        m_EdgeKernelType;

}; // End class ShortestPathKernel

} // End namespace tube

#endif // End !defined( __tubeShortestPathKernel_h )
