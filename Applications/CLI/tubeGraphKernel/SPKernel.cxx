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

#include "SPKernel.h"


using namespace boost;
using namespace std;

using boost::property_tree::ptree;
using boost::lexical_cast;


namespace tube
{


//-----------------------------------------------------------------------------
SPKernel::graphType SPKernel::graphFromAdjFile(const char *fileName)
{
  ifstream reader;

  int nVertices = 0;
  reader.open(fileName, std::ios::binary | std::ios::in);
  reader >> nVertices;
  reader.get(); // Move forward

  tube::FmtInfoMessage("Reading graph with %d vertices", nVertices);

  SPKernel::graphType g(nVertices);
  for ( int i=0; i<nVertices; ++i )
      {
      // Type is ID for now
      g[vertex(i, g)].type = i;
      }
  for ( int i=0; i<nVertices; ++i )
    {
    for ( int j=0; j<nVertices; ++j )
      {
      double tf;
      reader >> tf;
      reader.get(); // Move forward
      if (tf > 0)
        {
        // Edge between i and j
        add_edge( i, j, 1, g );
        }
      }
    }
    reader.close();
    return g;
}


//-----------------------------------------------------------------------------
SPKernel::graphType SPKernel::graphFromJSONFile(const char *fileName)
{
  try
    {
    ptree pt;
    read_json(fileName, pt);

    // Gets #vertices
    int nVertices = boost::lexical_cast<int>(pt.get<string>("nVertices"));

    // Parses linkage information
    int cnt = 0;
    vector<vector<int> > linkInfo(nVertices);
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("adjM"))
    {
      ptree adjInfoTree = (ptree)v.second;
      BOOST_FOREACH(ptree::value_type &w, adjInfoTree.get_child(""))
        {
        int vertexId = lexical_cast<int>(w.second.data());
        linkInfo[cnt].push_back(vertexId);
        }
      ++cnt;
    }

    // Parses distance information
    cnt = 0;
    vector<vector<double> > distInfo(nVertices);
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("dist"))
      {
      ptree distInfoTree = (ptree)v.second;
      BOOST_FOREACH(ptree::value_type &w, distInfoTree.get_child(""))
        {
        double distToNeighbor = lexical_cast<double>(w.second.data());
        distInfo[cnt].push_back(distToNeighbor);
        }
      ++cnt;
      }

    // Parses vertex type information
    vector<int> typeInfo;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("type"))
      {
      int type = lexical_cast<int>(v.second.data());
      typeInfo.push_back(type);
      }

    assert(linkInfo.size() == distInfo.size() &&
        linkInfo.size() <= nVertices &&
        distInfo.size() <= nVertices &&
        typeInfo.size() <= nVertices);

    // Now build the graph
    SPKernel::graphType g(nVertices);
    for (int i=0; i<nVertices;++i)
      {
      g[vertex(i, g)].type = typeInfo[i];
      }

    for (int i=0; i<linkInfo.size(); ++i)
      {
      assert(linkInfo[i].size() == distInfo[i].size());
      for (int j=1; j<linkInfo[i].size(); ++j)
        {
        tube::FmtDebugMessage("Adding edge (%d,%d) with weight %.5f",
          linkInfo[i][0], linkInfo[i][j], distInfo[i][j]);
        add_edge(linkInfo[i][0], linkInfo[i][j], distInfo[i][j], g);
      }
    }
    return g;
  }
  catch (std::exception &e)
  {
    tube::FmtErrorMessage("Error reading JSON graph file %s (Msg: %s)",
      fileName, e.what());
    throw std::exception();
  }
}


//-----------------------------------------------------------------------------
SPKernel::graphType SPKernel::floydTransform(const graphType &in)
{
  int nVertices = num_vertices(in);

  distanceMatrixType distances(nVertices);
  distanceMatrixMapType dm(distances, in);

  floyd_warshall_all_pairs_shortest_paths(in, dm);

  graphType out(nVertices);
  assert(nVertices == num_vertices(out));

  constVertexAllMapType mapIn = get(vertex_all, in);
  vertexAllMapType mapOut = get(vertex_all, out);

  for (int i=0; i<nVertices; ++i)
    {
    vertexType v0 = vertex(i,in);   // Vertex type of i-th node in input graph
    vertexType v1 = vertex(i,out);  // Vertex type of i-th node in output graph
    put(mapOut, v1, get(mapIn, v0));
    }

  for(int i=0; i<nVertices; ++i)
    {
    for (int j=0; j<nVertices; ++j)
      {
      // As long as we do not INF distance between (i,j), and ...
      if (dm[i][j] != numeric_limits<double>::max())
        {
        // the edge exists ...
        if (!edge(i, j, out).second)
          {
          // add an edge with the shortest path length
          add_edge(i, j, dm[i][j], out);
          }
        }
      }
    }
    return out;
}


//-----------------------------------------------------------------------------
void SPKernel::computeFTGraphs(void)
{
  tube::FmtDebugMessage("Computing Floyd transform.");
  fg0 = floydTransform(g0);
  fg1 = floydTransform(g1);
  isFloyd = true;
}


//-----------------------------------------------------------------------------
double SPKernel::compute(int edgeKernelType)
{
  double kernelValue = 0.0;
  edgeIteratorType aIt, aEnd, bIt, bEnd;

  if (!isFloyd)
    {
    computeFTGraphs();
    }

  long int cnt = 0;
  for (tie(aIt,aEnd) = edges(fg0); aIt != aEnd; ++aIt)
    {
    const edgeDescriptorType &e0 = *aIt;
    int src_type = fg0[source(e0, fg0)].type;
    int dst_type = fg0[target(e0, fg0)].type;
    edgeWeightMapType wmFG0 = get(edge_weight, fg0);

    for (tie(bIt, bEnd) = edges(fg1); bIt != bEnd; ++bIt)
      {
      cnt++;
      const edgeDescriptorType &e1 = *bIt;
      edgeWeightMapType wmFG1 = get(edge_weight, fg1);

      /*
       * We only consider walks of equal length --- At this point we
       * only support weights of 1, since this gives integer lengths
       * of the shortest paths and makes it easy to check for equality.
       *
       * We could also use a Brownian bridge kernel to bound the max.
       * shortest-path length, e.g., max(0,c - |len(e)-len(e')|)
       */
      double weightE0 = wmFG0[*aIt];
      double weightE1 = wmFG1[*bIt];
      if (!fabs(weightE0 - weightE1) < numeric_limits<double>::epsilon())
        {
        continue;
        }

      double edgeKernelValue = 0.0;
      switch (edgeKernelType)
        {
        // Our only supported kernel so far
        case EDGE_KERNEL_DEL:
          /*
           * Delta kernel on edge start/end types, i.e.,
           * 1 if types are equal, 0 otherwise
           */
          bool vertexTypeCheck =
            (src_type == fg1[source(e1, fg1)].type) &&
            (dst_type == fg1[target(e1, fg1)].type);

          if (!vertexTypeCheck)
            {
            continue;
            }
          edgeKernelValue = 1.0;
          break;
        }
      kernelValue += edgeKernelValue;
      ++cnt;
      }
    }
  tube::FmtInfoMessage("Performed %ld edge evaluations", cnt);
  return kernelValue;
}


} // End of namespace tube
