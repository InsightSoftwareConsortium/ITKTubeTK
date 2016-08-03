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
#include <itkMinimumMaximumImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

#include <metaTubeGraph.h>

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
  typedef itk::GroupSpatialObject< Dimension >       GroupType;
  typedef itk::ImageFileReader< ImageType >          ImageReaderType;
  typedef itk::SpatialObjectReader< >                SpatialObjectReaderType;
  typedef itk::VesselTubeSpatialObject< Dimension >  TubeSpatialObjectType;
  typedef TubeSpatialObjectType::TubePointType       TubePointType;
  typedef TubeSpatialObjectType::TransformType       TubeTransformType;

  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;
  std::stringstream logMsg;

  tube::InfoMessage( "Reading spatial objects..." );
  timeCollector.Start( "Load tubes" );
  SpatialObjectReaderType::Pointer soReader = SpatialObjectReaderType::New();
  soReader->SetFileName(tubeFile.c_str());
  GroupType::Pointer group;
  try
    {
    soReader->Update();
    group = soReader->GetGroup();
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
  imReader->SetFileName(voronoiFile.c_str());
  ImageType::Pointer image;
  try
    {
    imReader->Update();
    image = imReader->GetOutput();
    }
  catch( itk::ExceptionObject & ex )
    {
    tube::FmtErrorMessage( "Cannot read image file: %s",
        ex.what());
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load CVT" );
  timeCollector.Start( "Processing" );

  typedef itk::MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
  MinMaxFilterType::Pointer mmFilter = MinMaxFilterType::New();
  mmFilter->SetInput(image);
  mmFilter->Update();

  int numberOfCentroids = mmFilter->GetMaximum();

  logMsg.str( "" );
  logMsg << "Number of Centroids = " << numberOfCentroids;
  tube::InfoMessage( logMsg.str() );

  vnl_matrix<int> aMat(numberOfCentroids, numberOfCentroids);
  aMat.fill(0);

  vnl_matrix<double> cMat(3, 3);
  vnl_vector<double> cVect(3);

  vnl_vector<int> rootNodes(numberOfCentroids);
  rootNodes.fill(0);
  vnl_vector<double> branchNodes(numberOfCentroids);
  branchNodes.fill(0);

  int count = 0;
  char tubeName[10];
  std::sprintf( tubeName, "Tube" );
  TubeSpatialObjectType::ChildrenListType *
    tubeList = group->GetChildren( group->GetMaximumDepth(), tubeName );
  TubeSpatialObjectType::ChildrenListType::const_iterator
           tubeIt = tubeList->begin();
  int numTubes = tubeList->size();
  TubePointType tubePoint;
  MetaScene scene(3);
  MetaTubeGraph * graph;
  TubeTransformType::Pointer tubeTransform;
  while(tubeIt != tubeList->end())
    {
    TubeSpatialObjectType::Pointer tube =
          dynamic_cast<TubeSpatialObjectType *>((*tubeIt).GetPointer());

    tube->RemoveDuplicatePoints();
    tube->ComputeTangentAndNormals();

    int numberOfPoints = tube->GetNumberOfPoints();

    graph = new MetaTubeGraph(3);

    itk::Point<double, 3> pnt;
    itk::Index< 3 > indx;
    tubePoint = static_cast<TubePointType>(tube->GetPoints()[0]);
    pnt = tubePoint.GetPosition();
    tube->ComputeObjectToWorldTransform();
    tubeTransform = tube->GetIndexToWorldTransform();
    pnt = tubeTransform->TransformPoint(pnt);
    image->TransformPhysicalPointToIndex(pnt, indx);
    double cCount = 1;
    int cNode = image->GetPixel(indx);
    double cRadius = tubePoint.GetRadius();
    for(int i=0; i<3; i++)
      {
      cVect[i] = tubePoint.GetTangent()[i];
      }
    cMat = outer_product(cVect, cVect);
    if(tube->GetRoot())
      {
      rootNodes[cNode-1] = rootNodes[cNode-1]+1;
      }
    branchNodes[cNode-1] = branchNodes[cNode-1]+1.0/numTubes;
    int numberOfNodesCrossed = 0;
    for(int p=1; p<numberOfPoints; p++)
      {
      tubePoint = static_cast<TubePointType>(tube->GetPoints()[p]);
      pnt = tubePoint.GetPosition();
      pnt = tubeTransform->TransformPoint(pnt);
      image->TransformPhysicalPointToIndex(pnt, indx);
      int tNode = image->GetPixel(indx);
      if(tNode == cNode)
        {
        cCount++;
        cRadius += tubePoint.GetRadius();
        for(int i=0; i<3; i++)
          {
          cVect[i] = tubePoint.GetTangent()[i];
          }
        cMat = cMat + outer_product(cVect, cVect);
        }
      else
        {
        int len = graph->GetPoints().size();
        if(graph->GetPoints().size()>3
          && graph->GetPoints().at(len-1)->m_GraphNode == tNode
          && graph->GetPoints().at(len-2)->m_GraphNode == cNode)
          {
          logMsg.str( "" );
          logMsg  << "Oscillation detected"
                  << " : tube = " << cNode
                  << " : seq = " << graph->GetPoints().at(len-3)->m_GraphNode
                  << " " << graph->GetPoints().at(len-2)->m_GraphNode
                  << " " << graph->GetPoints().at(len-1)->m_GraphNode
                  << " " << cNode << " " << tNode;
          tube::WarningMessage( logMsg.str() );

          TubeGraphPnt * tgP = graph->GetPoints().back();
          cNode = tNode;
          cRadius = tgP->m_R;
          for(int i=0; i<3; i++)
            {
            for(int j=0; j<3; j++)
              {
              cMat[i][j] = tgP->m_T[i*3+j];
              }
            }
          cCount = tgP->m_P;
          graph->GetPoints().pop_back();
          /* Memory allocated for each element of list returned by
          graph->GetPoints() usually released when destructor of graph called,
          but since tgP is popped off back of list, memory would not be
          released without explicit delete. */
          delete tgP;
          }
        else
          {
          numberOfNodesCrossed++;
          aMat[cNode-1][tNode-1] = aMat[cNode-1][tNode-1]+1;
          TubeGraphPnt * tgP = new TubeGraphPnt(3);
          tgP->m_GraphNode = cNode;
          tgP->m_R = cRadius/cCount;
          tgP->m_P = cCount;
          for(int i=0; i<3; i++)
            {
            for(int j=0; j<3; j++)
              {
              tgP->m_T[i*3+j] = cMat[i][j] / cCount;
              }
            }
          graph->GetPoints().push_back(tgP);
          cNode = tNode;
          cRadius = tubePoint.GetRadius();
          for(int i=0; i<3; i++)
            {
            cVect[i] = tubePoint.GetTangent()[i];
            }
          cMat = outer_product(cVect, cVect);
          cCount = 1;
          }
        }
      }
    if(numberOfNodesCrossed>0)
      {
      TubeGraphPnt * tgP = new TubeGraphPnt(3);
      tgP->m_GraphNode = cNode;
      tgP->m_R = cRadius/cCount;
      for(int i=0; i<3; i++)
        {
        for(int j=0; j<3; j++)
          {
          tgP->m_T[i*3+j] = cMat[i][j] / cCount;
          }
        }
      graph->GetPoints().push_back(tgP);
      scene.AddObject(graph);
      }
    else
      {
      delete graph;
      }
    ++tubeIt;
    ++count;
    }

  timeCollector.Stop( "Processing" );

  timeCollector.Start( "Save data" );
  scene.Write(graphFile.c_str());

  std::string matrixFile = graphFile + ".mat";
  std::ofstream writeStream;
  writeStream.open(matrixFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for(int i=0; i<numberOfCentroids; i++)
    {
    for(int j=0; j<numberOfCentroids; j++)
      {
      writeStream << aMat[i][j];
      if(j<numberOfCentroids-1)
        {
        writeStream << " ";
        }
      }
    writeStream << std::endl;
    }
  writeStream.close();

  std::string branchFile = graphFile + ".brc";
  writeStream.open(branchFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for(int i=0; i<numberOfCentroids; i++)
    {
    writeStream << branchNodes[i] << std::endl;
    }
  writeStream.close();

  std::string rootFile = graphFile + ".rot";
  writeStream.open(rootFile.c_str(), std::ios::binary | std::ios::out);
  writeStream << numberOfCentroids << std::endl;
  for(int i=0; i<numberOfCentroids; i++)
    {
    writeStream << rootNodes[i] << std::endl;
    }
  writeStream.close();
  timeCollector.Stop( "Save data" );

  delete tubeList;
  timeCollector.Report();

  return EXIT_SUCCESS;
}
