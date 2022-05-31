/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeWLSubtreeKernel_h
#define __tubeWLSubtreeKernel_h

#include "GraphKernel.h"
#include "tubeMessage.h"

#include <itkTimeProbesCollectorBase.h>

#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace tube
{

/** \class WLSubtreeKernel
 * \brief WLSubtreekernel implements the Weisfeiler-Lehman
 * subtree kernel proposed in
 *
 * [1]  N. Shervashidze, P. Schweitzer, E.-J. Van Leeuwen and
 *      K. M. Borgwardt, "Weisfeiler-Lehman graph kernels", In:
 *      JMLR 12( Sep ), pp. 2539-2561, 2011.
 *
 * Please read this article for any further details on this
 * kind of graph kernel. Naming of variables in this class
 * is close to the original publication.
 */
class WLSubtreeKernel : public GraphKernel
{
public:

  typedef std::map<std::string, int>  LabelMapType;
  typedef std::vector<LabelMapType>   LabelMapVectorType;

  /** CTOR - Variant with no vertex label information */
  WLSubtreeKernel( const GraphType &G0,
                   const GraphType &G1,
                   const LabelMapVectorType & labelMap,
                   const int & labelCount,
                   int subtreeHeight )
                     : GraphKernel( G0, G1 ),
                       m_SubtreeHeight( subtreeHeight ),
                       m_LabelMap( labelMap ),
                       m_LabelCount( labelCount )
    {
    }

  /** Compute the WLSubtree kernel */
  double Compute( void );

  /**
   * Take graph information and update
   *
   *  1 ) 'labelMap' for each subtree height with compressed label mapping
   *  2 ) the number of compressed labels per subtree level ( in 'cLabCounter' )
   */
  static void UpdateLabelCompression( GraphType &G,
                             std::vector<LabelMapType> & labelMap,
                             int & cLabCounter,
                             int subtreeHeight );

private:
  /**
   * Take a graph 'G' and use the label map information and the number of
   * compressed labels per subtree level to compute a feature mapping phi
   * for the graph, see [1]
   */
  std::vector<int> BuildPhi( GraphType &G );

  /** Our initial set of vertex labels */
  std::set<int>              m_InitialLabelSet;

  /** Subtree height */
  int                        m_SubtreeHeight;

  /** Label map + Count */
  const LabelMapVectorType & m_LabelMap;
  const int                  m_LabelCount;

}; // End class WLSubtreeKernel

} // End namespace tube

#endif // End !defined( __tubeWLSubtreeKernel_h )
