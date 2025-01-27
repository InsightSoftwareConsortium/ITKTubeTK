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

#include "tubeWLSubtreeKernel.h"

namespace tube
{

void
WLSubtreeKernel::UpdateLabelCompression(GraphType &                 G,
                                        std::vector<LabelMapType> & labelMap,
                                        int &                       cLabCounter,
                                        int                         subtreeHeight)
{
  const int N = num_vertices(G);
  for (int i = 0; i < N; ++i)
  {
    /* At height = 0, we relabel the vertex types to start with 0 ( since the
       type can be an arbitrary number ). This will allow convenient
       indexing. */
    const int                    height = 0;
    const int                    type = G[vertex(i, G)].type;
    const std::string            vertexStr = boost::lexical_cast<std::string>(type);
    LabelMapType::const_iterator it = labelMap[height].find(vertexStr);
    if (it == labelMap[height].end())
    {
      labelMap[height][vertexStr] = cLabCounter;
      G[vertex(i, G)].type = cLabCounter;
      ++cLabCounter;
    }
    else
    {
      G[vertex(i, G)].type = it->second;
    }
  }

  /* For heights > 0, we have to be careful with relabeling immediately,
   * since we need neighbor information. We record the relabeling result
   * while fetching neighbor information and building the compressed label.
   * Once we are done with all vertices, we relabel. */
  for (int height = 1; height < subtreeHeight; ++height)
  {
    std::vector<int> relabel(N, -1);
    for (int i = 0; i < N; ++i)
    {
      const std::string            nbStr = BuildNeighborStr(G, i);
      LabelMapType::const_iterator it = labelMap[height].find(nbStr);
      if (it == labelMap[height].end())
      {
        labelMap[height][nbStr] = cLabCounter;
        relabel[i] = cLabCounter;
        ++cLabCounter;
      }
      else
      {
        relabel[i] = labelMap[height][nbStr];
      }
    }

    // Relabel now.
    for (int i = 0; i < N; ++i)
    {
      if (relabel[i] > 0)
      {
        G[vertex(i, G)].type = relabel[i];
      }
    }
  }
}

std::vector<int>
WLSubtreeKernel::BuildPhi(GraphType & G)
{
  std::vector<int> phi(m_LabelCount, 0);
  const int        N = num_vertices(G);

  for (int i = 0; i < N; ++i)
  {
    const int                    height = 0;
    const int                    type = G[vertex(i, G)].type;
    LabelMapType::const_iterator it = m_LabelMap[height].find(boost::lexical_cast<std::string>(type));
    if (it != m_LabelMap[height].end())
    {
      const int cLab = it->second;
      G[vertex(i, G)].type = cLab;
      ++phi[cLab];
    }
  }

  for (int height = 1; height < m_SubtreeHeight; ++height)
  {
    std::vector<int> relabel(N, -1);
    for (int i = 0; i < N; ++i)
    {
      const std::string            nbStr = BuildNeighborStr(G, i);
      LabelMapType::const_iterator it = m_LabelMap[height].find(nbStr);
      if (it != m_LabelMap[height].end())
      {
        const int cLab = it->second;
        relabel[i] = cLab;
        ++phi[cLab];
      }
    }
    for (int i = 0; i < N; ++i)
    {
      if (relabel[i] > 0)
      {
        G[vertex(i, G)].type = relabel[i];
      }
    }
  }
  return phi;
}

double
WLSubtreeKernel::Compute(void)
{
  const std::vector<int> phiG0 = this->BuildPhi(m_G0);
  const std::vector<int> phiG1 = this->BuildPhi(m_G1);

  if (phiG0.size() != phiG1.size())
  {
    tube::ErrorMessage("Mismatch in feature mapping.");
    throw std::exception();
  }
  return static_cast<double>(inner_product(phiG0.begin(), phiG0.end(), phiG1.begin(), 0.0));
}

} // End namespace tube
