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

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"
#include "itktubeTubeAngleOfIncidenceWeightFunction.h"
#include "itktubeTubePointWeightsCalculator.h"

#include <itkOptimizerParameters.h>
#include <itkSpatialObjectReader.h>
#include <itkTubeSpatialObjectPoint.h>
#include <itksys/SystemTools.hxx>

#include <fstream>

int itktubeTubeAngleOfIncidenceWeightFunctionTest( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: "
      << argv[0]
      << " inputTubeNetwork.tre"
      << " outputDirectory"
      << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputTubeNetwork = argv[1];
  const char * outputDirectory = argv[2];

  enum { Dimension = 3 };
  typedef itk::VesselTubeSpatialObject< Dimension > TubeSpatialObjectType;
  typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;
  typedef itk::TubeSpatialObjectPoint< Dimension >  TubePointType;
  typedef float                                     WeightType;

  // Read input tube tree.
  typedef itk::SpatialObjectReader< Dimension >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputTubeNetwork );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();
  std::cout << "Number of children = "
    << groupSpatialObject->GetNumberOfChildren()
    << std::endl;

  // Sub-sample the tube tree.
  typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< GroupSpatialObjectType,
    TubeSpatialObjectType >
      SubSampleTubeTreeFilterType;
  SubSampleTubeTreeFilterType::Pointer subSampleTubeTreeFilter =
    SubSampleTubeTreeFilterType::New();
  subSampleTubeTreeFilter->SetInput( reader->GetGroup() );
  subSampleTubeTreeFilter->SetSampling( 10 );
  try
    {
    subSampleTubeTreeFilter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::tube::Function::TubeAngleOfIncidenceWeightFunction<
    TubePointType, WeightType > WeightFunctionType;
  WeightFunctionType::Pointer weightFunction = WeightFunctionType::New();
  weightFunction->SetFractionalImportance( 0.8 );
  weightFunction->SetAngleDependence( 2.5 );
  WeightFunctionType::PointType probeOrigin;
  probeOrigin[0] = 354.88;
  probeOrigin[1] = -214.83;
  probeOrigin[2] = 185.739;
  weightFunction->SetUltrasoundProbeOrigin( probeOrigin );

  typedef itk::OptimizerParameters< WeightType > PointWeightsType;
  typedef itk::tube::TubePointWeightsCalculator< Dimension,
    TubeSpatialObjectType, WeightFunctionType,
    PointWeightsType > PointWeightsCalculatorType;
  PointWeightsCalculatorType::Pointer weightsCalculator
    = PointWeightsCalculatorType::New();
  weightsCalculator->SetTubeTreeSpatialObject(
    subSampleTubeTreeFilter->GetOutput() );
  weightsCalculator->SetPointWeightFunction( weightFunction );
  weightsCalculator->Compute();
  const PointWeightsType & weights
    = weightsCalculator->GetPointWeights();

  itksys::SystemTools::MakeDirectory( outputDirectory );
  std::string outputPrefix( outputDirectory );
  outputPrefix += "/";

  std::ostringstream ostrm;
  char childName[] = "Tube";
  GroupSpatialObjectType::Pointer tubeTree
    = subSampleTubeTreeFilter->GetOutput();
  GroupSpatialObjectType::ChildrenListType * tubeList =
    tubeTree->GetChildren( tubeTree->GetMaximumDepth(), childName );
  itk::SizeValueType tubeIndex = 0;
  typedef GroupSpatialObjectType::ChildrenListType::const_iterator
    TubesIteratorType;
  for( TubesIteratorType tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeSpatialObjectType * currentTube =
      dynamic_cast< TubeSpatialObjectType * >( ( *tubeIterator ).GetPointer() );
    if( currentTube != NULL )
      {
      const TubeSpatialObjectType::PointListType & currentTubePoints
        = currentTube->GetPoints();
      typedef TubeSpatialObjectType::PointListType::const_iterator
        TubePointIteratorType;
      for( TubePointIteratorType tubePointIterator = currentTubePoints.begin();
            tubePointIterator != currentTubePoints.end();
            ++tubePointIterator )
        {
        ostrm.str( "" );
        ostrm << outputPrefix;
        ostrm << "Weight" << tubeIndex << ".acsv";

        std::ofstream outputFile( ostrm.str().c_str() );
        if( !outputFile.is_open() )
          {
          std::cerr << "Could not open output file: " << ostrm.str() << std::endl;
          return EXIT_FAILURE;
          }

        outputFile << "# Name = W" << tubeIndex << "\n";
        outputFile << "# pointColor = "
                   << weights[tubeIndex] << ","
                   << weights[tubeIndex] << ","
                   << 0.3 << "\n";
        outputFile << "# pointColumns = type|x|y|z|sel|vis\n";
        outputFile << "point|";
        const TubePointType::PointType & position
          = tubePointIterator->GetPosition();
        outputFile << position[0] << "|";
        outputFile << position[1] << "|";
        outputFile << position[2] << "|";
        outputFile << "1|1\n";
        std::cout << position << ": " << weights[tubeIndex] << '\n';

        ++tubeIndex;
        }
      }
    }
  delete tubeList;

  // Regression test on a few points
  WeightType tolerance = 0.01;
  if( std::fabs( weights[3] - 0.818436 ) > tolerance ||
      std::fabs( weights[4] - 0.909169 ) > tolerance ||
      std::fabs( weights[5] - 0.951959 ) > tolerance )
    {
    std::cerr << "Did not generate the expected weights." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
