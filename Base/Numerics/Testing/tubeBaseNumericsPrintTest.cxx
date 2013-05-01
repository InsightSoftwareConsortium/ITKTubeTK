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
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#define ITK_LEAN_AND_MEAN

#include "itkImage.h"

#include "itkImageRegionMomentsCalculator.h"
#include "itkJointHistogramImageFunction.h"
#include "itkTubeBlurImageFunction.h"
#include "itkTubeLDAGenerator.h"
#include "itkTubeNJetLDAGenerator.h"
#include "itkTubeNJetImageFunction.h"
#include "itkVectorImageToListGenerator.h"
#include "itkVotingResampleImageFunction.h"

int tubeBaseNumericsPrintTest( int, char*[] )
{
  typedef itk::Image< float, 2 >                 ImageType;
  typedef itk::Image< itk::Vector<float, 2>, 2 > VectorImageType;

  itk::ImageRegionMomentsCalculator< ImageType >::Pointer
    regionMomentsObject =
    itk::ImageRegionMomentsCalculator< ImageType >::New();
  std::cout << "-------------itkImageRegionMomentsCalculator"
            << regionMomentsObject
            << std::endl;

  itk::JointHistogramImageFunction< ImageType >::Pointer
    jointHistoObject =
    itk::JointHistogramImageFunction< ImageType >::New();
  std::cout << "-------------itkJointHistogramImageFunction"
            << jointHistoObject
            << std::endl;

  itk::tube::BlurImageFunction< ImageType >::Pointer tbif =
    itk::tube::BlurImageFunction< ImageType > ::New();
  std::cout << "-------------tbif" << tbif << std::endl;

  itk::tube::NJetImageFunction< ImageType >::Pointer
    njetObject =
    itk::tube::NJetImageFunction< ImageType >::New();
  std::cout << "-------------itkTubeNJetImageFunction"
            << njetObject
            << std::endl;

  itk::tube::LDAGenerator< ImageType, ImageType >::Pointer
    ldaGenerator =
    itk::tube::LDAGenerator< ImageType, ImageType >::New();
  std::cout << "-------------LDAGenerator"
    << ldaGenerator << std::endl;

  itk::tube::NJetLDAGenerator< ImageType, ImageType >::Pointer
    njetLDAGenerator =
    itk::tube::NJetLDAGenerator< ImageType, ImageType >::New();
  std::cout << "-------------LDANJetGenerator"
    << njetLDAGenerator << std::endl;

  itk::Statistics::VectorImageToListGenerator< VectorImageType,
      ImageType >::Pointer
    vectorImageToListObject =
    itk::Statistics::VectorImageToListGenerator< VectorImageType,
      ImageType >::New();
  std::cout << "-------------itkVectorImageToListGenerator"
            << vectorImageToListObject
            << std::endl;

  itk::VotingResampleImageFunction< ImageType >::Pointer
    votingResampleObject =
    itk::VotingResampleImageFunction< ImageType >::New();
  std::cout << "-------------itkVotingResampleImageFunction"
            << votingResampleObject
            << std::endl;

  return EXIT_SUCCESS;
}

