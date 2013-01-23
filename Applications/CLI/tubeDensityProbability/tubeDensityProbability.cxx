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


template< unsigned int VDimension >
int DoIt( int argc, char **argv );

// This needs to be declared for tubeCLIHelperFunctions.
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char **argv ){ return 0; }
#include "tubeCLIHelperFunctions.h"

int main(int argc, char **argv)
{
  PARSE_ARGS;

  itk::ImageIOBase::IOComponentType componentType;
  unsigned int dimension;

  try
    {
    tube::GetImageInformation( inMeanImageFile, componentType, dimension );
    switch( dimension )
      {
      case 2:
        return DoIt< 2 >( argc, argv );
      case 3:
        return DoIt< 3 >( argc, argv );
      default:
        return EXIT_FAILURE;
      }
    }
  catch( const std::exception & exception )
    {
    tube::ErrorMessage( exception.what() );
    return EXIT_FAILURE;
    }
  return EXIT_FAILURE;
}


template< unsigned int VDimension >
int DoIt(int argc, char **argv)
{
  PARSE_ARGS;

  const unsigned int Dimension = VDimension;

  typedef itk::Image< short, Dimension >             ImageType;
  typedef itk::GroupSpatialObject< Dimension >       GroupType;
  typedef itk::ImageFileReader< ImageType >          ImageReaderType;
  typedef itk::SpatialObjectReader< Dimension >      SOReaderType;
  typedef itk::VesselTubeSpatialObject< Dimension >  TubeType;
  typedef typename TubeType::TubePointType           TubePointType;
  typedef typename TubeType::TransformType           TubeTransformType;


  tube::CLIProgressReporter progressReporter(
    "tubeDensityProbability",
    CLPProcessInformation );

  progressReporter.Start();

  /*
   * Read in spatial object file (tubes)
   */
  typename SOReaderType::Pointer soReader = SOReaderType::New();
  soReader->SetFileName( inTubeFile.c_str() );
  soReader->Update();
  typename GroupType::Pointer group = soReader->GetGroup();

  progressReporter.Report(0.1);

  /*
   * Read in ATLAS EMD image
   */
  typename ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( inMeanImageFile.c_str() );
  imReader->Update();
  typename ImageType::Pointer meanImage = imReader->GetOutput();

  progressReporter.Report(0.2);

  char childName[] = "Tube";
  typename TubeType::ChildrenListType * tubeList = group->GetChildren(99999, childName);
  typename TubeType::ChildrenListType::const_iterator tubeIt = tubeList->begin();
  TubePointType tubePoint;
  typename TubeTransformType::Pointer tubeTransform;
  std::ofstream writeStream;
  writeStream.open(outFile.c_str(), std::ios::binary | std::ios::out);
  while(tubeIt != tubeList->end()) // Iterate over tubes
    {
    typename TubeType::Pointer tube = dynamic_cast<TubeType *>((*tubeIt).GetPointer());

    tube->RemoveDuplicatePoints();
    tube->ComputeTangentAndNormals();

    itk::Point<double, Dimension> pnt;
    itk::Index< Dimension > indx;
    tube->ComputeObjectToWorldTransform();
    tubeTransform = tube->GetIndexToWorldTransform();
    for( unsigned int i=0; i<tube->GetNumberOfPoints() ; i++)
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
