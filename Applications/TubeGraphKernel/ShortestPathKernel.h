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


#ifndef __ShortestPathKernel_h
#define __ShortestPathKernel_h

#include <iostream>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <exception>

#include "GraphKernel.h"

#include "tubeMessage.h"
#include "itkTimeProbesCollectorBase.h"


namespace tube
{


/**
 *
 * \class ShortestPathKernel
 * \brief Impementation of Borgwardt's Shortest-Path Kernel
 *
 * This class implements the shortest-path kernel, proposed in
 *
 * [1] K.M. Borgwardt and H.P. Kriegel, "Shortest-Path Kernels on
 *     Graphs", In: IEEE Int. Conf. on Data Mining, 2005
 *
 */
class ShortestPathKernel : public GraphKernel
  {
  public:

    /** Edge kernel types */
    static const int EDGE_KERNEL_DEL = 0;

    /** CTOR - Consumer sets graphs */
    ShortestPathKernel(const GraphType &G0, const GraphType &G1) :
      GraphKernel(G0, G1), m_edgeKernelType(EDGE_KERNEL_DEL) {}

    /** Sets edge-kernel type */
    // cppcheck-suppress unusedFunction
    void SetEdgeKernel(int type)
      { m_edgeKernelType = type; }

    /** Computes the SP kernel value, see [1], Section 4.2 */
    double Compute(void);


  private:

    /** Floyd-transformed graphs */
    GraphType m_FG0, m_FG1;

    int m_edgeKernelType;

    /** Computes a Floyd-transformed graph, see [1], Section 4.1 */
    GraphType FloydTransform(const GraphType &in);

    template <typename T>
    void ensureOrder(T& first, T& second)
    {
      if ( first > second )
        {
        std::swap( first, second );
        }
    }
  };


} // End of namespace tube
#endif
