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

#include "tubeMessage.h"
#include "tubeMetaObjectDocument.h"

#include <vnl/algo/vnl_matrix_inverse.h>

#include "MergeTubeGraphsCLP.h"

int
DoIt(int argc, char * argv[]);

using namespace tube;

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  return DoIt(argc, argv);
}

int
DoIt(int argc, char * argv[])
{
  PARSE_ARGS;

  // Generic stream for log messages
  std::stringstream logMsg;

  typedef MetaObjectDocument                         DocumentReaderType;
  typedef DocumentReaderType::ObjectDocumentListType DocumentListType;

  int numberOfCentroids = nCentroids;
  assert(numberOfCentroids > 0);

  logMsg.str("");
  logMsg << "Number of centroids: " << numberOfCentroids;
  tube::InfoMessage(logMsg.str());

  DocumentReaderType * reader = new DocumentReaderType();
  reader->SetFileName(listFile.c_str());
  if (!reader->Read(listFile.c_str()))
  {
    tube::ErrorMessage("Could not read ObjectDocument file!");
    delete reader;
    return EXIT_FAILURE;
  }

  DocumentListType graphObjects = reader->GetObjectDocumentList();

  int numberOfGraphs = graphObjects.size();
  logMsg.str("");
  logMsg << "Number of graphs " << numberOfGraphs;
  tube::InfoMessage(logMsg.str());

  vnl_matrix<double> aMat(numberOfCentroids, numberOfCentroids);
  aMat.fill(0);
  vnl_vector<double> bVect(numberOfCentroids);
  bVect.fill(0);
  vnl_vector<double> rVect(numberOfCentroids);
  rVect.fill(0);

  double      tf;
  int         numberOfCentroids2;
  std::string filename;

  DocumentListType::const_iterator graphIt = graphObjects.begin();
  while (graphIt != graphObjects.end())
  {
    filename = (*graphIt)->GetObjectName();

    std::string matrixFilename = filename + ".mat";
    tube::InfoMessage("Reading file " + matrixFilename);
    std::ifstream readMatrixStream;
    readMatrixStream.open(matrixFilename.c_str(), std::ios::in);
    readMatrixStream >> numberOfCentroids2;
    readMatrixStream.get();
    if (numberOfCentroids != numberOfCentroids2)
    {
      std::cerr << "Error: fileList's #Centroids != matrix #Centroids" << std::endl;
      std::cerr << numberOfCentroids << " != " << numberOfCentroids2 << std::endl;
      delete reader;
      return 0;
    }
    for (int i = 0; i < numberOfCentroids; i++)
    {
      for (int j = 0; j < numberOfCentroids; j++)
      {
        readMatrixStream >> tf;
        readMatrixStream.get();
        aMat[i][j] += tf;
      }
    }
    readMatrixStream.close();

    std::string   branchFilename = filename + ".brc";
    std::ifstream readBranchStream;
    readBranchStream.open(branchFilename.c_str(), std::ios::in);
    readBranchStream >> numberOfCentroids2;
    readBranchStream.get();
    if (numberOfCentroids != numberOfCentroids2)
    {
      std::cerr << "Error: fileList's #Centroids != branch #Centroids" << std::endl;
      delete reader;
      return 0;
    }
    for (int i = 0; i < numberOfCentroids; i++)
    {
      readBranchStream >> tf;
      readBranchStream.get();
      bVect[i] += tf;
    }
    readBranchStream.close();

    std::string   rootFilename = filename + ".rot";
    std::ifstream readRootStream;
    readRootStream.open(rootFilename.c_str(), std::ios::in);
    readRootStream >> numberOfCentroids2;
    readRootStream.get();
    if (numberOfCentroids != numberOfCentroids2)
    {
      std::cerr << "Error: fileList's #Centroids != root #Centroids" << std::endl;
      delete reader;
      return 0;
    }
    for (int i = 0; i < numberOfCentroids; i++)
    {
      readRootStream >> tf;
      readRootStream.get();
      rVect[i] += tf;
    }
    readRootStream.close();

    ++graphIt;
  }

  std::string   matrixFile = graphFile + ".mat";
  std::ofstream writeStream;
  writeStream.open(matrixFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for (int i = 0; i < numberOfCentroids; i++)
  {
    for (int j = 0; j < numberOfCentroids; j++)
    {
      aMat[i][j] = aMat[i][j] / numberOfGraphs;
      writeStream << aMat[i][j];
      if (j < numberOfCentroids - 1)
      {
        writeStream << " ";
      }
    }
    writeStream << std::endl;
  }
  writeStream.close();

  vnl_matrix<double> iMat(numberOfCentroids, numberOfCentroids);
  iMat.set_identity();

  vnl_matrix<double> cntMat(numberOfCentroids, numberOfCentroids);
  cntMat = iMat - 0.1 * aMat;
  vnl_vector<double> e(numberOfCentroids);
  e.fill(1);
  vnl_vector<double> cnt(numberOfCentroids);
  vnl_matrix<double> cntMatI(numberOfCentroids, numberOfCentroids);
  cntMatI = vnl_matrix_inverse<double>(cntMat).inverse();
  cnt = cntMatI * e;
  std::string cntFile = graphFile + ".cnt";
  writeStream.open(cntFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for (int i = 0; i < numberOfCentroids; i++)
  {
    writeStream << cnt[i] << std::endl;
  }
  writeStream.close();

  std::string branchFile = graphFile + ".brc";
  writeStream.open(branchFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for (int i = 0; i < numberOfCentroids; i++)
  {
    writeStream << bVect[i] / numberOfGraphs << std::endl;
  }
  writeStream.close();

  std::string rootFile = graphFile + ".rot";
  writeStream.open(rootFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for (int i = 0; i < numberOfCentroids; i++)
  {
    writeStream << rVect[i] / numberOfGraphs << std::endl;
  }
  writeStream.close();

  delete reader;
  return 1;
}
