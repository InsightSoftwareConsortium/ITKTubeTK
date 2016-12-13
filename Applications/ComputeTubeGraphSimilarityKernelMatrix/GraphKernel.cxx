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

#include "GraphKernel.h"

namespace tube
{


std::string GraphKernel::LabelVectorToString( const std::vector<int> &labVec )
{
  std::stringstream res;
  std::copy( 
    labVec.begin(),
    labVec.end(),
    std::ostream_iterator<int>( res, "" ) );
  return res.str();
}


std::string GraphKernel::BuildPrefixFromVertexID( int v )
{
  std::stringstream res;
  res << v << ",";
  return res.str();
}


void GraphKernel::CountingSort( std::vector<int>& vec )
{
  int maxVal = *std::max_element( vec.begin(), vec.end() );

  std::vector<int> C( maxVal + 1 );
  BOOST_FOREACH( int a, vec )
    {
    ++C[a];
    }
  int current = 0;
  for( int i = 0; i <= maxVal; ++i )
    {
    for( int j =0; j < C[i]; ++j )
      {
      vec[current++] = i;
      }
    }
}


std::string GraphKernel::BuildNeighborStr( const GraphType &G, int v )
{
  // Algorithm:
  //
  // 1 ) Get index map for the current graph 'G'
  // 2 ) Compute the #neighbors for vertex 'v'
  // 3 ) Add all neighbors of v to a 'nbVec'
  // 4 ) Sort the neighbor vector ( CountingSort )
  // 5 ) Create string from vector
  IndexMapType index = boost::get( boost::vertex_index, G );

  VertexNeighborType nb = adjacent_vertices( vertex( v, G ), G );
  int nNeighbors = std::distance( nb.first, nb.second );

  // No neighbors
  if( !nNeighbors )
    {
    return BuildPrefixFromVertexID( G[vertex( v, G )].type );
    }

  std::vector<int> nbVec( nNeighbors );
  for( int cnt=0; nb.first != nb.second; ++nb.first, ++cnt )
    {
    int vertexIndex = index[*nb.first];
    int tp = G[vertex( vertexIndex, G )].type;
    nbVec[cnt] = tp;
    }
  CountingSort( nbVec );

  return BuildPrefixFromVertexID( G[vertex( v, G )].type ) +
       LabelVectorToString( nbVec );
}


bool GraphKernel::IsValidDefaultNodeLabeling( int desiredType )
{
  switch( desiredType )
    {
    case LABEL_BY_NUM:
    case LABEL_BY_DEG:
      return true;
    default:
      return false;
    }
}


GraphKernel::GraphType
GraphKernel::GraphFromAdjFile( const char *graphFile,
                               const char *labelFile,
                               DefaultNodeLabelingType defNodeLabel )
{
  std::ifstream reader;

  int nVertices = 0;
  reader.open( graphFile, std::ios::binary | std::ios::in );
  reader >> nVertices;
  reader.get();

  tube::FmtInfoMessage( "Reading graph with %d vertices",
    nVertices );

  GraphKernel::GraphType g( nVertices );
  for( int i=0; i<nVertices; ++i )
      {
      // Defaults to ID
      g[vertex( i, g )].type = i;
      }
  for( int i=0; i<nVertices; ++i )
    {
    for( int j=0; j<nVertices; ++j )
      {
      double tf = 0.0;
      reader >> tf;
      reader.get();
      if( tf > 0 )
        {
        add_edge( i, j, 1, g );
        }
      }
    }
  reader.close();


  // We support either the node ID or the node degree as a
  // label ( degree is suggested in [Shervashidze11a]
  switch( defNodeLabel )
    {
    case LABEL_BY_NUM:
      tube::InfoMessage( "Using node number ( ID ) as label!" );
      break;
    case LABEL_BY_DEG:
      tube::InfoMessage( " Using node degree as label!" );
      break;
    }

  for( int i=0; i<nVertices; ++i )
    {
    switch( defNodeLabel )
      {
      // ... use node number
      case LABEL_BY_NUM:
        g[vertex( i, g )].type = i;
        break;
      // ... use node degree
      case LABEL_BY_DEG:
        g[vertex( i, g )].type = out_degree( i, g );
        break;
      }
    }

  // In case no label file is given
  if( !labelFile )
    {
    return g;
    }

  tube::FmtInfoMessage( "Label file given - reading %s",
    labelFile );

  int nLabels = 0;
  reader.open( labelFile, std::ios::binary | std::ios::in );
  reader >> nLabels;
  reader.get();

  assert( nLabels == nVertices );
  for( int i=0; i<nLabels; ++i )
    {
    int lab;
    reader >> lab;
    reader.get();
    g[vertex( i, g )].type = lab;

    }
  reader.close();
  return g;
}


GraphKernel::GraphType
GraphKernel::GraphFromJSONFile( const char *graphFile )
{
  try
    {
    boost::property_tree::ptree pt;
    read_json( graphFile, pt );

    // Gets #vertices
    int nVertices = boost::lexical_cast<int>( 
      pt.get<std::string>( "nVertices" ) );


    // Parses linkage information
    int cnt = 0;
    std::vector<std::vector<int> > linkInfo( nVertices );
    BOOST_FOREACH( boost::property_tree::ptree::value_type &v,
      pt.get_child( "adjM" ) )
      {
      boost::property_tree::ptree adjInfoTree =
        ( boost::property_tree::ptree )v.second;

      BOOST_FOREACH( boost::property_tree::ptree::value_type &w,
        adjInfoTree.get_child( "" ) )
        {
        int vertexId = boost::lexical_cast<int>( w.second.data() );
        linkInfo[cnt].push_back( vertexId );
        }
      ++cnt;
      }


    // Parses distance information
    cnt = 0;
    std::vector<std::vector<double> > distInfo( nVertices );
    BOOST_FOREACH( boost::property_tree::ptree::value_type &v,
      pt.get_child( "dist" ) )
      {
      boost::property_tree::ptree distInfoTree =
        ( boost::property_tree::ptree )v.second;

      BOOST_FOREACH( boost::property_tree::ptree::value_type &w,
        distInfoTree.get_child( "" ) )
        {
        double distToNeighbor = boost::lexical_cast<double>( w.second.data() );
        distInfo[cnt].push_back( distToNeighbor );
        }
      ++cnt;
      }


    // Parses vertex type information
    std::vector<int> typeInfo;
    BOOST_FOREACH( boost::property_tree::ptree::value_type &v,
      pt.get_child( "type" ) )
      {
      int type = boost::lexical_cast<int>( v.second.data() );
      typeInfo.push_back( type );
      }

    assert( linkInfo.size() == distInfo.size() &&
        static_cast< int >( linkInfo.size() ) <= nVertices &&
        static_cast< int >( distInfo.size() ) <= nVertices &&
        static_cast< int >( typeInfo.size() ) <= nVertices );

    // Now build the graph
    GraphKernel::GraphType g( nVertices );
    for( int i=0; i<nVertices;++i )
      {
      g[vertex( i, g )].type = typeInfo[i];
      }

    for( unsigned int i=0; i<linkInfo.size(); ++i )
      {
      assert( linkInfo[i].size() == distInfo[i].size() );
      for( unsigned int j=1; j<linkInfo[i].size(); ++j )
        {
        tube::FmtDebugMessage( "Adding edge ( %d,%d ) with weight %.5f",
          linkInfo[i][0], linkInfo[i][j], distInfo[i][j] );

        add_edge( linkInfo[i][0],
                  linkInfo[i][j],
                  distInfo[i][j], g );
        }
      }
    return g;
    }
  catch( const std::exception &e )
    {
    tube::FmtErrorMessage( "Error reading JSON graph file %s ( Msg: %s )",
      graphFile, e.what() );
    throw;
    }
}

} // End namespace tube
