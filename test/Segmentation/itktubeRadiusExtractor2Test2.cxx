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

#include "itktubeRadiusExtractor2.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkSpatialObjectReader.h>

int
itktubeRadiusExtractor2Test2(int argc, char * argv[])
{
  if (argc != 3)
  {
    std::cout << "itktubeRadiusExtractor2Test2 <inputImage> <vessel.tre>" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Image<float, 3> ImageType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer                imReader = ImageReaderType::New();
  imReader->SetFileName(argv[1]);
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::tube::RadiusExtractor2<ImageType> RadiusOpType;
  RadiusOpType::Pointer                          radiusOp = RadiusOpType::New();

  radiusOp->SetInputImage(im);

  bool returnStatus = EXIT_SUCCESS;

  radiusOp->SetMinMedialness(0.005);
  if (radiusOp->GetMinMedialness() != 0.005)
  {
    tube::ErrorMessage("MinMedialness != 0.005");
    returnStatus = EXIT_FAILURE;
  }

  radiusOp->SetMinMedialnessStart(0.002);
  if (radiusOp->GetMinMedialnessStart() != 0.002)
  {
    tube::ErrorMessage("MinMedialnessStart != 0.002");
    returnStatus = EXIT_FAILURE;
  }

  typedef itk::SpatialObjectReader<>             ReaderType;
  typedef itk::SpatialObject<>::ChildrenListType ObjectListType;
  typedef itk::GroupSpatialObject<>              GroupType;
  typedef itk::TubeSpatialObject<>               TubeType;
  typedef TubeType::TubePointListType            TubePointListType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  reader->Update();
  GroupType::Pointer group = reader->GetGroup();

  std::cout << "Number of children = " << group->GetNumberOfChildren() << std::endl;

  char tubeName[17];
  std::strcpy(tubeName, "Tube");
  ObjectListType * tubeList = group->GetChildren(-1, tubeName);

  unsigned int numTubes = tubeList->size();
  std::cout << "Number of tubes = " << numTubes << std::endl;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator RandGenType;
  RandGenType::Pointer                                           rndGen = RandGenType::New();
  rndGen->Initialize(); // set seed here

  double       avgFailures = 0;
  double       avgAvgDiff = 0;
  double       avgMaxDiff = 0;
  unsigned int numMCRuns = 5;
  for (unsigned int mcRun = 0; mcRun < numMCRuns; mcRun++)
  {
    std::cout << std::endl;
    std::cout << "*** RUN = " << mcRun << std::endl;
    unsigned int rndTubeNum = 1;
    // Use the shorter tube
    // rndGen->GetUniformVariate( 0, 1 ) * numTubes;
    if (rndTubeNum >= numTubes)
    {
      rndTubeNum = numTubes - 1;
    }
    ObjectListType::iterator tubeIter = tubeList->begin();
    for (unsigned int i = 0; i < rndTubeNum; i++)
    {
      ++tubeIter;
    }

    std::cout << "Test tube = " << rndTubeNum << std::endl;

    TubeType *        tubep = static_cast<TubeType *>(tubeIter->GetPointer());
    TubeType::Pointer tube = static_cast<TubeType *>(tubeIter->GetPointer());

    tube->ComputeTangentsAndNormals();

    unsigned int numPoints = tube->GetPoints().size();

    unsigned int startPoint = rndGen->GetUniformVariate(0, 1) * (numPoints - 10) + 5;
    double       radiusStart = rndGen->GetUniformVariate(0, 1) * 2 + 1;

    double              failures = 0;
    std::vector<double> idealR;
    idealR.clear();
    TubePointListType::iterator pntIter = tube->GetPoints().begin();
    for (unsigned int i = 0; i < numPoints; i++)
    {
      if (i == startPoint)
      {
        pntIter->SetId(0);
        radiusStart = (radiusStart + pntIter->GetRadiusInObjectSpace()) / 2;
      }
      else
      {
        pntIter->SetId(i + 1);
      }

      if (std::fabs(pntIter->GetTangentInObjectSpace().GetVnlVector().magnitude() - 1) > 0.01)
      {
        std::cout << "Point: " << i << ": Tangent not of unit length." << std::endl;
        ++failures;
      }
      if (std::fabs(pntIter->GetNormal1InObjectSpace().GetVnlVector().magnitude() - 1) > 0.01)
      {
        std::cout << "Point: " << i << ": Normal1 not of unit length." << std::endl;
        ++failures;
      }
      if (std::fabs(pntIter->GetNormal2InObjectSpace().GetVnlVector().magnitude() - 1) > 0.01)
      {
        std::cout << "Point: " << i << ": Normal2 not of unit length." << std::endl;
        ++failures;
      }
      if (std::fabs(dot_product(pntIter->GetTangentInObjectSpace().GetVnlVector(),
                                pntIter->GetNormal1InObjectSpace().GetVnlVector())) > 0.001)
      {
        std::cout << "Point: " << i << ": dot_product( Tangent, Normal1 ) != 0." << std::endl;
        ++failures;
      }
      if (std::fabs(dot_product(pntIter->GetTangentInObjectSpace().GetVnlVector(),
                                pntIter->GetNormal2InObjectSpace().GetVnlVector())) > 0.001)
      {
        std::cout << "Point: " << i << ": dot_product( Tangent, Normal2 ) != 0." << std::endl;
        ++failures;
      }
      if (std::fabs(dot_product(pntIter->GetNormal1InObjectSpace().GetVnlVector(),
                                pntIter->GetNormal2InObjectSpace().GetVnlVector())) > 0.001)
      {
        std::cout << "Point: " << i << ": dot_product( Normal1, Normal2 ) != 0." << std::endl;
        ++failures;
      }

      idealR.push_back(pntIter->GetRadiusInObjectSpace());
      if (pntIter->GetRadiusInObjectSpace() < 0)
      {
        std::cout << "Point: " << i << ": radius < 0. ( " << pntIter->GetRadiusInObjectSpace() << " )" << std::endl;
        ++failures;
      }

      ++pntIter;
    }

    radiusOp->SetRadiusStart(radiusStart);
    radiusOp->SetRadiusMin(0.33);
    radiusOp->SetRadiusMax(15.0);
    radiusOp->SetRadiusStep(0.25);
    radiusOp->SetRadiusTolerance(0.125);

    radiusOp->ExtractRadii(tubep);

    double avgDiff = 0;
    double maxDiff = 0;
    pntIter = tube->GetPoints().begin();
    for (unsigned int i = 0; i < numPoints; i++)
    {
      double diff = std::fabs(pntIter->GetRadiusInObjectSpace() - idealR[i]);
      avgDiff += diff;
      if (diff > maxDiff)
      {
        maxDiff = diff;
      }
      if (diff > 1)
      {
        std::cout << "Point: " << i << "  idealR = " << idealR[i]
                  << "  estimatedR = " << pntIter->GetRadiusInObjectSpace() << " : FAIL " << std::endl;
        ++failures;
      }
      else
      {
        std::cout << "Point: " << i << "  idealR = " << idealR[i]
                  << "  estimatedR = " << pntIter->GetRadiusInObjectSpace() << std::endl;
      }

      // reset radius for re-testing
      pntIter->SetRadiusInObjectSpace(idealR[i]);
      ++pntIter;
      if (i < numPoints - 1)
      {
        pntIter->SetRadiusInObjectSpace(idealR[i]);
        ++pntIter;
        ++i;
      }
    }
    avgDiff /= numPoints;
    failures = failures / numPoints;
    std::cout << "Failures = " << failures << std::endl;
    std::cout << "Average diff = " << avgDiff << std::endl;
    std::cout << "Max diff = " << maxDiff << std::endl;
    avgAvgDiff += avgDiff;
    avgMaxDiff += maxDiff;
    avgFailures += failures;
  }

  delete tubeList;

  avgFailures /= numMCRuns;
  avgAvgDiff /= numMCRuns;
  avgMaxDiff /= numMCRuns;
  std::cout << "Average Failures = " << avgFailures << std::endl;
  std::cout << "Average Average diff = " << avgAvgDiff << std::endl;
  std::cout << "Average Max diff = " << avgMaxDiff << std::endl;

  if (avgFailures > 0.15)
  {
    return EXIT_FAILURE;
  }

  return returnStatus;
}
