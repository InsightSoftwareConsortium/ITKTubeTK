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

#include "tubeMessage.h"

#include <itkVesselTubeSpatialObject.h>

#include <fstream>

#include "ComputeTubeGraphProbabilityCLP.h"

int DoIt( int argc, char * argv[] );

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}


int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::stringstream logMsg;
  std::ifstream readStream;

  double tf;
  unsigned int numberOfCentroids;

  // Graph connectivity information file
  std::string filename = inGraphFile + ".mat";

  logMsg.str( "" );
  logMsg << "Reading file: " << filename;
  tube::InfoMessage( logMsg.str() );

  readStream.open( filename.c_str(), std::ios::binary | std::ios::in );
  readStream >> numberOfCentroids;
  readStream.get();

  vnl_matrix<double> cMat( numberOfCentroids, numberOfCentroids );
  cMat.fill( 0 );
  vnl_vector<double> bVect( numberOfCentroids );
  bVect.fill( 0 );

  vnl_matrix<double> meanCMat( numberOfCentroids, numberOfCentroids );
  meanCMat.fill( 0 );
  vnl_vector<double> meanBVect( numberOfCentroids );
  meanBVect.fill( 0 );
  vnl_vector<double> meanCVect( numberOfCentroids );
  meanCVect.fill( 0 );

  // Fill-up connectivity matrix
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    for( unsigned int j=0; j<numberOfCentroids; j++ )
      {
      readStream >> tf;
      readStream.get();
      cMat[i][j] = tf;
      }
    }
  readStream.close();

  unsigned int numberOfCentroids2;

  // Branch information
  filename = inGraphFile + ".brc";

  logMsg.str( "" );
  logMsg << "Reading file: " << filename;
  tube::InfoMessage( logMsg.str() );

  readStream.open( filename.c_str(), std::ios::binary | std::ios::in );
  readStream >> numberOfCentroids2;
  readStream.get();
  if( numberOfCentroids != numberOfCentroids2 )
    {
    tube::ErrorMessage(
      "Error: fileList's #Centroids != branch #Centroids" );
    return EXIT_FAILURE;
    }
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    readStream >> tf;
    readStream.get();
    bVect[i] = tf;
    }
  readStream.close();

  // MEAN graph connectivity file
  filename = meanGraphFile + ".mat";

  logMsg.str( "" );
  logMsg << "Reading file: " << filename;
  tube::InfoMessage( logMsg.str() );

  readStream.open( filename.c_str(), std::ios::binary | std::ios::in );
  readStream >> numberOfCentroids2;
  readStream.get();
  if( numberOfCentroids != numberOfCentroids2 )
    {
    tube::ErrorMessage(
      "Error: fileList's #Centroids != mean matrix #Centroids" );
    return EXIT_FAILURE;
    }
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    for( unsigned int j=0; j<numberOfCentroids; j++ )
      {
      readStream >> tf;
      readStream.get();
      meanCMat[i][j] = tf;
      }
    }
  readStream.close();

  // MEAN branch file
  filename = meanGraphFile + ".brc";

  logMsg.str( "" );
  logMsg << "Reading file: " << filename;
  tube::InfoMessage( logMsg.str() );

  readStream.open( filename.c_str(), std::ios::binary | std::ios::in );
  readStream >> numberOfCentroids2;
  readStream.get();
  if( numberOfCentroids != numberOfCentroids2 )
    {
    tube::ErrorMessage(
      "Error: fileList's #Centroids != mean branch #Centroids" );
    return EXIT_FAILURE;
    }
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    readStream >> tf;
    readStream.get();
    meanBVect[i] = tf;
    }
  readStream.close();

  // MEAN centrality information
  filename = meanGraphFile + ".cnt";

  logMsg.str( "" );
  logMsg << "Reading file: " << filename;
  tube::InfoMessage( logMsg.str() );

  readStream.open( filename.c_str(), std::ios::binary | std::ios::in );
  readStream >> numberOfCentroids2;
  readStream.get();
  if( numberOfCentroids != numberOfCentroids2 )
    {
    tube::ErrorMessage(
      "Error: fileList's #Centroids != mean centrality #Centroids" );
    return EXIT_FAILURE;
    }
  for( unsigned int i=0; i<numberOfCentroids; i++ )
    {
    readStream >> tf;
    readStream.get();
    meanCVect[i] = tf;
    }
  readStream.close();

  std::ofstream writeMatStream;
  std::ofstream writeCntStream;
  std::ofstream writeBrcStream;
  filename = outGraphFile + ".matProb.txt";
  writeMatStream.open( filename.c_str(), std::ios::binary | std::ios::out );
  filename = outGraphFile + ".cntProb.txt";
  writeCntStream.open( filename.c_str(), std::ios::binary | std::ios::out );
  filename = outGraphFile + ".brcProb.txt";
  writeBrcStream.open( filename.c_str(), std::ios::binary | std::ios::out );

  for( unsigned int i=0; i < numberOfCentroids; i++ )
    {
    bool used = false;
    for( unsigned int j = 0; j < numberOfCentroids; j++ )
      {
      // ( i,j ) connected
      if( cMat[i][j] > 0 )
        {
        writeMatStream << meanCMat[i][j] << std::endl;
        used = true;
        }
      }
    if( used )
      {
      writeCntStream << meanCVect[i] << std::endl;
      if( bVect[i] > 0 )
        {
        writeBrcStream << meanBVect[i] << std::endl;
        }
      }
    }
  writeMatStream.close();
  writeBrcStream.close();
  writeCntStream.close();

  return EXIT_SUCCESS;
}
