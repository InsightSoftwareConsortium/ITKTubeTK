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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include "tubeMessage.h"
#include "itkTimeProbesCollectorBase.h"


namespace tube
{

/** \class ShortestPathKernel
 * \brief Impementation of a Shortest-Path Kernel
 *
 * This class implements the shortest-path kernel as proposed in
 *
 * [1] K.M. Borgwardt and H.P. Kriegel, "Shortest-Path Kernels on
 *     Graphs", In: IEEE Int. Conf. on Data Mining, 2005
 */
class ShortestPathKernel
  {
  public:

    /** Edge kernel types */
    static const int EDGE_KERNEL_DEL = 0;   // k(x,x') = Delta kernel

    /** Node meta-information */
    struct nodeInfoType
    {
      int type;
    };

    /** Some typedefs for BGL */
    typedef boost::adjacency_list<
      boost::listS,     // Use a list for vertices
      boost::vecS,      // Use a vector for edges
      boost::directedS, // We have a directed graph
      nodeInfoType,     // Own vertex information
      boost::property<boost::edge_weight_t, double> > GraphType; // Edge-weights

    typedef boost::exterior_vertex_property<GraphType, double> DistancePropertyType;
    typedef DistancePropertyType::matrix_type DistanceMatrixType;
    typedef DistancePropertyType::matrix_map_type DistanceMatrixMapType;

    /** Vertex descriptor and vertex iterator */
    typedef boost::graph_traits<GraphType>::vertex_descriptor VertexType;
    typedef boost::graph_traits<GraphType>::vertex_iterator VertexIteratorType;

    /** Edge descriptor, edge iterator and edge weight map */
    typedef boost::graph_traits<GraphType>::edge_descriptor EdgeDescriptorType;
    typedef boost::graph_traits<GraphType>::edge_iterator EdgeIteratorType;
    typedef boost::property_map<GraphType, boost::edge_weight_t>::type EdgeWeightMapType;

    /** Types used to hold all vetex (+ property) information  */
    typedef boost::property_map<GraphType, boost::vertex_all_t>::const_type ConstVertexAllMapType;
    typedef boost::property_map<GraphType, boost::vertex_all_t>::type VertexAllMapType;


  private:

    /*
     * g0, g1   ... original graphs
     * fg0, fg1 ... Floyd-transformed graphs
     */
    GraphType m_g0, m_g1, m_fg0, m_fg1;

    /** Graph information */
    int m_nVerticesG0, m_nVerticesG1;

    /** Did we already Floyd-transform the graphs */
    bool m_isFloyd;


  public:

    /** CTOR - Consumer sets graphs */
    ShortestPathKernel(const GraphType &g0, const GraphType &g1) :
      m_isFloyd(false)
    {
      m_g0 = g0;
      m_g1 = g1;
    }

    /** Runs Floyd-transformation on both graphs */
    void ComputeFTGraphs(void);

    /** Reads graph information from JSON file */
    static GraphType GraphFromJSONFile(const char *fileName);

    /** Reads graph information from adjacency matrix */
    static GraphType GraphFromAdjFile(const char *fileName);

    /** Computes the SP kernel value, see [1], Section 4.2 */
    double Compute(int edgeKernelType);


  private:

    /** Computes a Floyd-transformed graph, see [1], Section 4.1 */
    GraphType FloydTransform(const GraphType &in);
  };
} // End of namespace tube
#endif
