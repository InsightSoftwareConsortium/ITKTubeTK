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


#include <iostream>
#include <fstream>
#include <list>
#include <vector>

#include "itkSpatialObjectReader.h"
#include "itkImageFileReader.h"
#include "itkVesselTubeSpatialObject.h"

#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

#include "tubeDensityProbabilityCLP.h"


/** Forward decl. */
int DoIt( int, char **argv );


int main(int argc, char **argv)
{
  PARSE_ARGS;

  return DoIt( argc, argv );
}


int DoIt(int argc, char **argv)
{
  PARSE_ARGS;

  typedef itk::Image< short, 3 >             ImageType;
  typedef itk::GroupSpatialObject<3>         GroupType;
  typedef itk::ImageFileReader< ImageType >  ImageReaderType;
  typedef itk::SpatialObjectReader<3>        SOReaderType;
  typedef itk::VesselTubeSpatialObject< 3 >  TubeType;
  typedef TubeType::TubePointType            TubePointType;
  typedef TubeType::TransformType            TubeTransformType;


  tube::CLIProgressReporter progressReporter(
    "tubeDensityProbability",
    CLPProcessInformation );

  progressReporter.Start();

  /*
   * Read in spatial object file (tubes)
   */
  SOReaderType::Pointer soReader = SOReaderType::New();
  soReader->SetFileName( inTubeFile.c_str() );
  soReader->Update();
  GroupType::Pointer group = soReader->GetGroup();

  progressReporter.Report(0.1);

  /*
   * Read in ATLAS EMD image
   */
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( inMeanImageFile.c_str() );
  imReader->Update();
  ImageType::Pointer meanImage = imReader->GetOutput();

  progressReporter.Report(0.2);


  TubeType::ChildrenListType * tubeList = group->GetChildren(99999, "Tube");
  TubeType::ChildrenListType::const_iterator tubeIt = tubeList->begin();
  TubePointType tubePoint;
  TubeTransformType::Pointer tubeTransform;
  std::ofstream writeStream;
  writeStream.open(outFile.c_str(), std::ios::binary | std::ios::out);
  while(tubeIt != tubeList->end()) // Iterate over tubes
    {
    TubeType::Pointer tube = dynamic_cast<TubeType *>((*tubeIt).GetPointer());

    tube->RemoveDuplicatePoints();
    tube->ComputeTangentAndNormals();

    itk::Point<double, 3> pnt;
    itk::Index< 3 > indx;
    tube->ComputeObjectToWorldTransform();
    tubeTransform = tube->GetIndexToWorldTransform();
    for(int i=0; i<tube->GetNumberOfPoints() ; i++)
      {
      tubePoint = static_cast<TubePointType>(tube->GetPoints()[i]); // Get point
      pnt = tubePoint.GetPosition(); // Get point's position
      pnt = tubeTransform->TransformPoint(pnt); // Point coords to physical coords
      meanImage->TransformPhysicalPointToIndex(pnt, indx); // Get closest voxel
      writeStream << meanImage->GetPixel(indx) << std::endl; // Write value of ATLAS EMD file at voxel
      }

    ++tubeIt;
    }
  writeStream.close();

  progressReporter.Report(1.0);
  progressReporter.End();

  return 1;
  }
