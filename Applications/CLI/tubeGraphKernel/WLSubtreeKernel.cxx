#include "WLSubtreeKernel.h"


namespace tube
{


//-----------------------------------------------------------------------------
void WLSubtreeKernel::UpdateLabelCompression( GraphType &G,
                                              std::vector<LabelMapType> & labelMap,
                                              int &cLabCounter,
                                              int subtreeHeight)
{
  int h = 0;
  int N = num_vertices(G);
  LabelMapType::const_iterator it;
  for(int i=0; i<N; ++i)
    {
    /*
     * On level h = 0, we relabel the vertex types to start
     * with 0 (since the type can be an arbitrary number).
     * This will allow convenient indexing!
     */
    int type = G[vertex(i,G)].type;
    std::string vertexStr = boost::lexical_cast<std::string>(type);
    it = labelMap[h].find(vertexStr);
    if(it == labelMap[h].end())
      {
      tube::FmtDebugMessage("Relabel %d (%s) -> %d",
        i, vertexStr.c_str(), cLabCounter);
      labelMap[h][vertexStr] = cLabCounter;
      G[vertex(i,G)].type = cLabCounter;
      ++cLabCounter;
      }
    else
      {
      tube::FmtDebugMessage("Relabel %d (%s) -> %d",
        i, vertexStr.c_str(), (*it).second);
      int cLab =(*it).second;
      G[vertex(i,G)].type = cLab;
      }
    }

  /*
   * For levels > 0, we have to be careful with relabeling immediately, since
   * we need neighbor information. Thus, we recored the relabeling result and
   * while fetching neighbor information and building the compressed label.
   * Once we are done with all vertices, we relabel!
   */
  for( h=1; h<subtreeHeight; ++h )
    {
    std::vector<int> relabel(N,-1);
    LabelMapType::iterator it;
    for( int i=0; i<N; ++i )
      {
      std::string nbStr = BuildNeighborStr(G, i);
      it = labelMap[h].find(nbStr);
      if( it == labelMap[h].end() )
        {
        tube::FmtDebugMessage("(N) Relabel %d (%s) -> %d",
          i, nbStr.c_str(), cLabCounter);
        labelMap[h][nbStr] = cLabCounter;
        relabel[i] = cLabCounter;
        ++cLabCounter;
        }
      else
        {
        tube::FmtDebugMessage("(F) Relabel %d (%s) -> %d",
          i, nbStr.c_str(), labelMap[h][nbStr]);
        relabel[i] = labelMap[h][nbStr];
        }
      }
      // RELABEL now!
      for( int i=0; i<N; ++i )
        {
        if( relabel[i] > 0 )
          {
          G[vertex(i,G)].type = relabel[i];
          }
        }
    }
}


//-----------------------------------------------------------------------------
std::vector<int> WLSubtreeKernel::BuildPhi( GraphType &G )
{
  int dim = m_labelCount;
  std::vector<int> phi(dim, 0);
  int N = num_vertices(G);

  int h = 0;
  LabelMapType::const_iterator it;
  for( int i=0; i<N; ++i )
    {
    int tp = G[vertex(i,G)].type;
    it = m_labelMap[h].find(boost::lexical_cast<std::string>(tp));
    if( it != m_labelMap[h].end() )
      {
      int cLab = (*it).second;
      G[vertex(i,G)].type = cLab;
      phi[cLab]++;
      }
    }

  for( h=1; h<m_subtreeHeight; ++h )
    {
    std::vector<int> relabel(N,-1);
    for( int i=0; i<N; ++i )
      {
      std::string nbStr = BuildNeighborStr(G, i);
      it = m_labelMap[h].find(nbStr);
      if( it != m_labelMap[h].end() )
        {
        int cLab = (*it).second;
        relabel[i] = cLab;
        phi[cLab]++;
        }
      }
      for( int i=0; i<N; ++i )
        {
        if (relabel[i]>0)
          {
          G[vertex(i,G)].type = relabel[i];
          }
        }
    }
  return phi;
}


//-----------------------------------------------------------------------------
double WLSubtreeKernel::Compute()
{
  std::vector<int> phiG0 = BuildPhi(m_G0);
  std::vector<int> phiG1 = BuildPhi(m_G1);

  if( phiG0.size() != phiG1.size() )
    {
    tube::ErrorMessage("Mismatch in feature mapping!");
    throw std::exception();
    }
  int kVal = inner_product( phiG0.begin(),
                            phiG0.end(),
                            phiG1.begin(),
                            0.0);
  return (double)kVal;
}


} // End of namespace tube
