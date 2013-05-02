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


#ifndef __tubeWLSubtreeKernel_h
#define __tubeWLSubtreeKernel_h

#include <iostream>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <map>
#include <set>
#include <exception>

#include "GraphKernel.h"

#include "tubeMessage.h"
#include "itkTimeProbesCollectorBase.h"


namespace tube
{

/**
 *
 * \class WLSubtreeKernel
 * \brief WLSubtreekernel implements the Weisfeiler-Lehman
 * subtree kernel proposed in
 *
 * [1]  N. Shervashidze, P. Schweitzer, E.-J. Van Leeuwen and
 *      K. M. Borgwardt, "Weisfeiler-Lehman graph kernels", In:
 *      JMLR 12(Sep), pp. 2539−2561, 2011.
 *
 * Please read this article for any further details on this
 * kind of graph kernel. Naming of variables in this class
 * is close to the original publication.
 *
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
                     int subtreeHeight) :
                        GraphKernel(G0, G1),
                        m_subtreeHeight(subtreeHeight),
                        m_labelMap(labelMap),
                        m_labelCount(labelCount)
    {
    }

    /** Compute the WLSubtree kernel */
    double Compute(void);

    /*
     * Take graph information and update
     *
     *  1) 'labelMap' for each subtree height with compressed label mapping
     *  2) the number of compressed labels per subtree level (in 'cLabCounter')
     *
     */
    static void UpdateLabelCompression( GraphType &G,
                               std::vector<LabelMapType> & labelMap,
                               int & cLabCounter,
                               int subtreeHeight);


  private:

    /** Our initial set of vertex labels */
    std::set<int> m_initialLabelSet;

    /** Subtree height */
    int m_subtreeHeight;

    /** Label map + Count */
    const LabelMapVectorType & m_labelMap;
    const int m_labelCount;

    /*
     * Take a graph 'G' and use the label map information and the number of
     * compressed labels per subtree level to compute a feature mapping phi
     * for the graph, see [1]
     */
    std::vector<int> BuildPhi( GraphType &G );
  };


} // End of namespace tube
#endif
