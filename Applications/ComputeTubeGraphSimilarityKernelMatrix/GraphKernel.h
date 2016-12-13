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

#ifndef __GraphKernel_h
#define __GraphKernel_h

#include "tubeMessage.h"

#include <boost/foreach.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

namespace tube
{

/**
 * \brief Base class for all graph kernels
 */
class GraphKernel
{

public:

  /** Default strategy for node labeling without label file */
  typedef enum
    {
    LABEL_BY_NUM = 0,
    LABEL_BY_DEG = 1
    } DefaultNodeLabelingType;

  /** Node meta-information */
  struct nodeInfoType
    {
    int type;
    }; // End struct nodeInfoType

  /** The graph type */
  typedef boost::adjacency_list< boost::listS,
                                 boost::vecS,
                                 boost::undirectedS,
                                 nodeInfoType,
                                 boost::property< boost::edge_weight_t,
                                                  double> > GraphType;
  /** Access vertex index information */
  typedef boost::property_map<
    GraphType, boost::vertex_index_t>::type     IndexMapType;

  /** Access vertex neighborhoods */
  typedef boost::graph_traits<GraphType>
    ::adjacency_iterator                        AdjacencyIteratorType;

  typedef std::pair< AdjacencyIteratorType,
                     AdjacencyIteratorType>     VertexNeighborType;

  /** Access to shortest-path information */
  typedef boost::exterior_vertex_property<
    GraphType, double>                          DistancePropertyType;
  typedef DistancePropertyType::matrix_type     DistanceMatrixType;
  typedef DistancePropertyType::matrix_map_type DistanceMatrixMapType;

  /** Vertex descriptor and vertex iterator */
  typedef boost::graph_traits<
    GraphType>::vertex_descriptor               VertexType;
  typedef boost::graph_traits<
    GraphType>::vertex_iterator                 VertexIteratorType;

  /** Edge descriptor, edge iterator and edge-weight information */
  typedef boost::graph_traits<
    GraphType>::edge_descriptor                 EdgeDescriptorType;
  typedef boost::graph_traits<
    GraphType>::edge_iterator                   EdgeIteratorType;
  typedef boost::property_map<
    GraphType, boost::edge_weight_t>::type      EdgeWeightMapType;

  /** Types used to access all vertex ( +property ) information  */
  typedef boost::property_map<
    GraphType,
    boost::vertex_all_t>::const_type            ConstVertexAllMapType;
  typedef boost::property_map<
    GraphType,
    boost::vertex_all_t>::type                  VertexAllMapType;

  /** CTOR */
  GraphKernel( const GraphType &G0, const GraphType &G1 )
    {
    m_G0 = G0;
    m_G1 = G1;
    }

  virtual ~GraphKernel( void )
    {
    }

  /** Check if the desired default node labeling is supported */
  static bool IsValidDefaultNodeLabeling( int desiredType );

  /** Read graph from adj file */
  static GraphType GraphFromAdjFile( const char *graphFile,
                                     const char *labelFile,
                                     DefaultNodeLabelingType defNodeLabel );
  /** Read graph from JSON file */
  static GraphType GraphFromJSONFile( const char *graphFile );

  /** Compute kernel value among graphs G0,G1 */
  virtual double Compute( void ) { return 0.0; }


protected:

  /** Build a string representation of the neighbors of v-th vertex */
  static std::string BuildNeighborStr( const GraphType &G, int v );

  /** Builds and returns a prefix "<i>," string, see [1] for details */
  static std::string BuildPrefixFromVertexID( int v );

  /** Format vector of int as string */
  static std::string LabelVectorToString( const std::vector<int> & );

  /** In-place counting sort on int vector */
  static void CountingSort( std::vector<int> & vec );

  /** Two input graphs */
  GraphType m_G0;
  GraphType m_G1;

}; // End class GraphKernel

} // End namespace tube

#endif // End !defined( __GraphKernel_h )
