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

#include "tubeShortestPathKernel.h"

#include <vnl/vnl_math.h>

namespace tube
{


//-----------------------------------------------------------------------------
ShortestPathKernel::GraphType
ShortestPathKernel::FloydTransform( const GraphType &in )
{
  int nVertices = num_vertices( in );

  DistanceMatrixType distances( nVertices );
  DistanceMatrixMapType dm( distances, in );

  floyd_warshall_all_pairs_shortest_paths( in, dm );

  GraphType out( nVertices );
  assert( nVertices == static_cast< int >( num_vertices( out ) ) );

  ConstVertexAllMapType mapIn = boost::get( boost::vertex_all, in );
  VertexAllMapType mapOut     = boost::get( boost::vertex_all, out );

  for( int i=0; i<nVertices; ++i )
    {
    // Vertex type of i-th node in input graph
    VertexType v0 = vertex( i, in );

    // Vertex type of i-th node in output graph
    VertexType v1 = vertex( i, out );

    put( mapOut, v1, get( mapIn, v0 ) );
    }

  for( int i=0; i<nVertices; ++i )
    {
    for( int j=0; j<=i; ++j )
      {
      // As long as we do not INF distance between ( i,j ), and ...
      if( dm[i][j] != std::numeric_limits<double>::max() )
        {
        // the edge exists ...
        if( !edge( i, j, out ).second )
          {
          // add an edge with the shortest path length
          add_edge( i, j, dm[i][j], out );
          }
        }
      else
        {
        //tube::FmtWarningMessage( "Numeric limit found at ( %d,%d )!",
        //  i,j );
        }
      }
    }

  return out;
}


//-----------------------------------------------------------------------------
double ShortestPathKernel::Compute( void )
{
  double kernelValue = 0.0;
  long int cntEdgeEvaluations = 0;
  EdgeIteratorType aIt, aEnd;

  tube::FmtDebugMessage( "Computing Floyd transform." );
  m_FG0 = FloydTransform( m_G0 );
  m_FG1 = FloydTransform( m_G1 );


  // Get the weight maps for both Floyd-transformed graphs
  EdgeWeightMapType wmFG0 = boost::get( boost::edge_weight, m_FG0 );
  EdgeWeightMapType wmFG1 = boost::get( boost::edge_weight, m_FG1 );

  // Iterate over all the edges of Floyd-transformed graph fg0
  for( tie( aIt, aEnd ) = edges( m_FG0 ); aIt != aEnd; ++aIt )
    {
    const EdgeDescriptorType &e0 = *aIt;

    int srcLabel = m_FG0[source( e0, m_FG0 )].type; // Type of start vertex
    int dstLabel = m_FG0[target( e0, m_FG0 )].type; // Type of end vertex
    ensureOrder( srcLabel, dstLabel );

    // Iterate over all the edges of Floyd-transformed graph fg1
    EdgeIteratorType bIt, bEnd;
    for( tie( bIt, bEnd ) = edges( m_FG1 ); bIt != bEnd; ++bIt )
      {
      cntEdgeEvaluations++;
      const EdgeDescriptorType &e1 = *bIt;

      // We only consider walks of equal length --- At this point we
      // only support weights of 1, since this gives integer lengths
      // of the shortest paths and makes it easy to check for equality.
      //
      // We could also use a Brownian bridge kernel to bound the max.
      // shortest-path length, e.g., max( 0,c - |len( e )-len( e' )| )

      double weightE0 = wmFG0[*aIt];
      double weightE1 = wmFG1[*bIt];

      // Skip edges, unless the weights are equal
      if( !( vnl_math_abs( weightE0 - weightE1 ) <
          std::numeric_limits<double>::epsilon() ) )
        {
        continue;
        }

      double edgeKernelValue = 0.0;
      switch( m_EdgeKernelType )
        {
        case EDGE_KERNEL_DEL:
          int cmpSrcLabel = m_FG1[source( e1, m_FG1 )].type;
          int cmpDstLabel = m_FG1[target( e1, m_FG1 )].type;
          ensureOrder( cmpSrcLabel, cmpDstLabel );

          bool vertexTypeCheck =
            ( srcLabel == cmpSrcLabel ) &&
            ( dstLabel == cmpDstLabel );
          if( !vertexTypeCheck )
            {
            continue;
            }
          else
            {
            edgeKernelValue = 1.0;
            }
          break;
        }
      kernelValue += edgeKernelValue;
      }
    }

  tube::FmtInfoMessage( "Performed %ld edge evaluations",
    cntEdgeEvaluations );
  return kernelValue;
}


} // End namespace tube
