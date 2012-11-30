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


#ifndef SPKernel_h
#define SPKernel_h

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


class SPKernel
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
      boost::listS, // Use a list for vertices
      boost::vecS,  // Use a vector for edges
      boost::directedS, // We have a directed graph
      nodeInfoType, // Own vertex information
      boost::property<boost::edge_weight_t, double> > graphType; // Edge-weights

    typedef boost::exterior_vertex_property<graphType, double> distancePropertyType;
    typedef distancePropertyType::matrix_type distanceMatrixType;
    typedef distancePropertyType::matrix_map_type distanceMatrixMapType;

    /** Vertex descriptor and vertex iterator */
    typedef boost::graph_traits<graphType>::vertex_descriptor vertexType;
    typedef boost::graph_traits<graphType>::vertex_iterator vertexIterType;

    /** Edge descriptor, edge iterator and edge weight map */
    typedef boost::graph_traits<graphType>::edge_descriptor edgeDescriptorType;
    typedef boost::graph_traits<graphType>::edge_iterator edgeIteratorType;
    typedef boost::property_map<graphType, boost::edge_weight_t>::type edgeWeightMapType;

    /** Types used to hold all vetex (+ property) information  */
    typedef boost::property_map<graphType, boost::vertex_all_t>::const_type constVertexAllMapType;
    typedef boost::property_map<graphType, boost::vertex_all_t>::type vertexAllMapType;


  private:

    /** Graphs */
    graphType g0, g1, fg0, fg1;

    /** Graph information */
    int nVerticesG0, nVerticesG1;

    /** Did we already Floyd-transform the graphs */
    bool isFloyd;


  private:

    /** Computes a Floyd-transformed graph, see [1], Section 4.1 */
    graphType floydTransform(const graphType &in);


  public:

    /** Avoid problems */
    explicit SPKernel();

    /** CTOR - Consumer sets graphs */
    SPKernel(const graphType &_g0, const graphType &_g1) :
      isFloyd(false)
    {
      g0 = _g0;
      g1 = _g1;
    }

    /** Runs Floyd-transformation on both graphs */
    void computeFTGraphs(void);

    /** Reads graph information from JSON file */
    static graphType graphFromJSONFile(const char *fileName);

    /** Reads graph information from adjacency matrix */
    static graphType graphFromAdjFile(const char *fileName);

    /** Computes the SP kernel value, see [1], Section 4.2 */
    double compute(int edgeKernelType);
};


} // End of namespace tube
#endif
