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

#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkMatrix.h"
#include "ShortestPathKernel.h"
#include "tubeGraphKernelCLP.h"

using namespace boost;
using namespace tube;
using namespace std;

using boost::property_tree::ptree;
using boost::lexical_cast;


//-----------------------------------------------------------------------------
/** Reads a list of input graphs from JSON file */
void readGraphList(const string &fileName,
                   vector<string> &list,
                   vector<int> &labels)
{
  list.clear();
  try
    {
    ptree pt;
    read_json(fileName, pt);

    int nGraphs = boost::lexical_cast<int>(pt.get<string>("nGraphs"));

    int graphCount = 0;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("graphList"))
      {
      const string &file = lexical_cast<string>(v.second.data());
      list.push_back( file );
      ++graphCount;
      }

    int labelCount = 0;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("labels"))
      {
      int label = lexical_cast<int>(v.second.data());
      labels.push_back( label );
      ++labelCount;
      }
    assert(graphCount == nGraphs);
    assert(labelCount == nGraphs);
    }
  catch (std::exception &e)
    {
    tube::FmtErrorMessage("Error reading JSON graph file %s (Msg: %s)",
      fileName.c_str(), e.what());
    throw std::exception();
    }
}


//-----------------------------------------------------------------------------
/** Writes kernel to binary file */
void writeKernel(const string &baseFileName, const vnl_matrix<double> &K)
{
  string outFileName = baseFileName + ".bin";
  ofstream ofs;
  ofs.open( outFileName.c_str(), ios::binary | ios::out );
  if (!ofs)
    {
    tube::FmtErrorMessage("Could not open kernel matrix %s for writing!",
      outFileName.c_str());
    throw std::exception();
    }
  ofs.write( (const char *)K.data_block(), K.size()*sizeof(double) );
  if (ofs.bad())
    {
    tube::FmtErrorMessage("Could not write kernel matrix %s!",
      outFileName.c_str());
    throw std::exception();
    }
  ofs.close();
  if (ofs.fail())
    {
    tube::FmtErrorMessage("Could not close kernel matrix file %s!",
      outFileName.c_str());
    throw std::exception();
    }
}


//-----------------------------------------------------------------------------
/** Writes kernel file in LIBSVM compatible format */
void writeKernelLibSVM(const string &baseFileName,
                       const vnl_matrix<double> &K,
                       const vector<int> &labels)
{
  string outFileName = baseFileName + ".libsvm";
  ofstream ofs;
  ofs.open( outFileName.c_str(), ios::binary | ios::out );
  if (!ofs)
    {
    tube::FmtErrorMessage("Could not open kernel matrix %s for writing!",
      outFileName.c_str());
    throw std::exception();
    }

  assert(K.rows() == labels.size());

  for (unsigned int r = 0; r < K.rows(); ++r)
    {
    ofs << labels[r] << " " << "0:" << r+1 << " ";
    for (unsigned int c = 0; c < K.cols()-1; ++c)
      {
      ofs << c+1 << ":" << K[r][c] << " ";
      }
    ofs << K.cols() << ":" << K[r][K.cols()-1] << endl;
    }
  ofs.close();
  if (ofs.fail())
    {
    tube::FmtErrorMessage("Could not close kernel matrix file %s!",
      outFileName.c_str());
    throw std::exception();
    }
}


//-----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  PARSE_ARGS;
  try
    {
    tube::FmtInfoMessage("Input graph list %s.", graphListA.c_str());

    vector<string> listA, listB;
    vector<int> notUsed, labelsA;
    readGraphList( graphListA, listA, labelsA ); // First list contributes labels
    readGraphList( graphListB, listB, notUsed ); // Labels are not used!!

    int N = listA.size();
    int M = listB.size();

    assert( N > 0 && M > 0);

    tube::FmtDebugMessage("Read N=%d entries from %s.", N, graphListA.c_str());
    tube::FmtDebugMessage("Read M=%d entries from %s.", M, graphListB.c_str());

    vnl_matrix<double> K(N,M); // Kernel matrix
    K.fill(0);

    for ( int i = 0; i < N; ++i)
      {
      // Reads i-th graph from file
      ShortestPathKernel::GraphType f =
        ShortestPathKernel::GraphFromAdjFile( listA[i].c_str() );


      // Now, iterate over the 'compare-to' graphs
      for ( int j = 0; j < M; ++j)
        {
        tube::FmtInfoMessage("Computing (%d,%d)-th kernel entry",
          i, j);

        // Reads j-th 'compare-to' graphs
        ShortestPathKernel::GraphType g =
          ShortestPathKernel::GraphFromAdjFile( listB[j].c_str() );

        // Creates the Shortest-Path Kernel and sets K[i][j]-th kernel entry
        ShortestPathKernel spk(f,g);
        K[i][j] = spk.Compute(ShortestPathKernel::EDGE_KERNEL_DEL);
        }
      }

    writeKernel( outputKernel, K ); // in binary format
    writeKernelLibSVM( outputKernel, K, labelsA ); // in LibSVM format
    }
  catch (std::exception &e)
    {
    tube::ErrorMessage("Exiting ...");
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
