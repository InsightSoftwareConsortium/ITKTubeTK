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

#include "itkOrientedImage.h"
#include "itkMatrix.h"

#include "itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h"
#include "itkAnisotropicDiffusionTensorFunction.h"
#include "itkAnisotropicDiffusionTensorImageFilter.h"
#include "itkAnisotropicEdgeEnhancementDiffusionImageFilter.h"
#include "itkAnisotropicHybridDiffusionImageFilter.h"
#include "itkStructureTensorRecursiveGaussianImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"
#include "itkTubeBlurImageFunction.h"
#include "itkTubeEnhancingDiffusion2DImageFilter.h"
#include "itkTubeLDAGenerator.h"
#include "itkTubeNJetLDAGenerator.h"


int tubeBaseFilteringPrintTest( int, char* [] )
{
  typedef itk::OrientedImage< float, 3 > ImageType;
  typedef itk::Matrix< float, 3 >        MatrixType;

  itk::AnisotropicCoherenceEnhancingDiffusionImageFilter< ImageType, ImageType >
    ::Pointer acedif =
    itk::AnisotropicCoherenceEnhancingDiffusionImageFilter< ImageType,
    ImageType >::New();
  std::cout << "-------------acedif" << acedif << std::endl;

  itk::AnisotropicDiffusionTensorFunction< ImageType >
    ::Pointer adtf =
    itk::AnisotropicDiffusionTensorFunction< ImageType >::New();
  std::cout << "-------------adtf" << adtf << std::endl;

  //itk::AnisotropicDiffusionTensorImageFilter< ImageType, ImageType >
    //::Pointer adtif =
    //itk::AnisotropicDiffusionTensorImageFilter< ImageType, ImageType >
    //::New();
  //std::cout << "-------------adtif" << adtif << std::endl;

  itk::AnisotropicEdgeEnhancementDiffusionImageFilter< ImageType, ImageType >
    ::Pointer aeedif =
    itk::AnisotropicEdgeEnhancementDiffusionImageFilter< ImageType, ImageType >
    ::New();
  std::cout << "-------------aeedif" << aeedif << std::endl;

  itk::AnisotropicHybridDiffusionImageFilter< ImageType, ImageType >
    ::Pointer ahdif =
    itk::AnisotropicHybridDiffusionImageFilter< ImageType, ImageType >
    ::New();
  std::cout << "-------------ahdif" << ahdif << std::endl;

  //itk::StructureTensorRecursiveGaussianImageFilter< ImageType, ImageType >
    //::Pointer strgif =
    //itk::StructureTensorRecursiveGaussianImageFilter< ImageType, ImageType >
    //::New();
  //std::cout << "-------------strgif" << strgif << std::endl;

  //itk::SymmetricEigenVectorAnalysisImageFilter< ImageType, ImageType,
    //ImageType >::Pointer sevaif =
    //itk::SymmetricEigenVectorAnalysisImageFilter< ImageType, ImageType,
    //ImageType >::New();
  //std::cout << "-------------sevaif" << sevaif << std::endl;

  itk::tube::BlurImageFunction< ImageType >::Pointer tbif =
    itk::tube::BlurImageFunction< ImageType > ::New();
  std::cout << "-------------tbif" << tbif << std::endl;

  itk::TubeEnhancingDiffusion2DImageFilter< float, 2 >::Pointer
    vesselEnahncingObj =
    itk::TubeEnhancingDiffusion2DImageFilter< float, 2 >::New();
  std::cout << "-------------TubeEnhancingDiffusion2DImageFilter"
    << vesselEnahncingObj << std::endl;

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

  return EXIT_SUCCESS;
}
