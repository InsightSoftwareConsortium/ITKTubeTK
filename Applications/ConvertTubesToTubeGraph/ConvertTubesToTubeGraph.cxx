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

#include <itkImageFileReader.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>
#include <tubeConvertTubesToTubeGraph.h>

#include "ConvertTubesToTubeGraphCLP.h"

int DoIt( int argc, char * argv[] );

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}

int DoIt( int argc, char * argv[] )
{
  enum { Dimension = 3 };

  typedef short                                      PixelType;
  typedef itk::Image< PixelType, Dimension >         ImageType;
  typedef itk::ImageFileReader< ImageType >          ImageReaderType;
  typedef itk::SpatialObjectReader< >                SpatialObjectReaderType;

  typedef tube::ConvertTubesToTubeGraph< PixelType, Dimension >
    FilterType;
  FilterType::Pointer filter = FilterType::New();
  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;
  std::stringstream logMsg;

  tube::InfoMessage( "Reading spatial objects..." );
  timeCollector.Start( "Load tubes" );
  SpatialObjectReaderType::Pointer soReader = SpatialObjectReaderType::New();
  soReader->SetFileName( tubeFile.c_str() );
  try
    {
    soReader->Update();
    filter->SetInputTubeGroup( soReader->GetGroup() );
    }
  catch(...)
    {
    tube::ErrorMessage( "ERROR: Cannot read spatial object file" );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load tubes" );

  tube::InfoMessage( "Reading image..." );
  timeCollector.Start( "Load CVT" );
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( voronoiFile.c_str() );
  try
    {
    imReader->Update();
    filter->SetCVTImage( imReader->GetOutput() );
    }
  catch( itk::ExceptionObject & ex )
    {
    tube::FmtErrorMessage( "Cannot read image file: %s",
        ex.what());
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load CVT" );
  timeCollector.Start( "Processing" );

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    tube::FmtErrorMessage( "TubeSpatialObjectToTubeGraphFilter Update error: %s",
        ex.what());
    return EXIT_FAILURE;
    }

  int numberOfCentroids = filter->GetNumberOfCenteroids();

  logMsg.str( "" );
  logMsg << "Number of Centroids = " << numberOfCentroids;
  tube::InfoMessage( logMsg.str() );

  vnl_matrix< double > aMat( numberOfCentroids, numberOfCentroids );
  aMat = filter->GetAdjacencyMatrix();

  vnl_vector< int > rootNodes( numberOfCentroids );
  rootNodes = filter->GetRootNodes();
  vnl_vector< double > branchNodes( numberOfCentroids );
  branchNodes = filter->GetBranchNodes();

  timeCollector.Stop( "Processing" );

  timeCollector.Start( "Save data" );

  std::string matrixFile = graphFile + ".mat";
  std::ofstream writeStream;
  writeStream.open( matrixFile.c_str(), std::ios::binary | std::ios::out );
  writeStream << numberOfCentroids << std::endl;
  for( int i = 0; i < numberOfCentroids; i++ )
    {
    for( int j = 0; j < numberOfCentroids; j++ )
      {
      writeStream << aMat[i][j];
      if( j < numberOfCentroids - 1 )
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  std::string branchFile = graphFile + ".brc";
  writeStream.open( branchFile.c_str(), std::ios::binary | std::ios::out );
  writeStream << numberOfCentroids << std::endl;
  for( int i = 0; i < numberOfCentroids; i++ )
    {
    writeStream << branchNodes[i] << std::endl;
    }
  writeStream.close();

  std::string rootFile = graphFile + ".rot";
  writeStream.open( rootFile.c_str(), std::ios::binary | std::ios::out );
  writeStream << numberOfCentroids << std::endl;
  for( int i = 0; i < numberOfCentroids; i++ )
    {
    writeStream << rootNodes[i] << std::endl;
    }
  writeStream.close();
  timeCollector.Stop( "Save data" );

  timeCollector.Report();

  return EXIT_SUCCESS;
}
