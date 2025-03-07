/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#include "itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.h"
#include "itktubeAnisotropicDiffusionTensorFunction.h"
#include "itktubeAnisotropicDiffusionTensorImageFilter.h"
#include "itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h"
#include "itktubeAnisotropicHybridDiffusionImageFilter.h"
#include "itktubeCropImageFilter.h"
#include "itktubeCVTImageFilter.h"
#include "itktubeExtractTubePointsSpatialObjectFilter.h"
#include "itktubeFFTGaussianDerivativeIFFTFilter.h"
#include "itktubeMinimumSpanningTreeVesselConnectivityFilter.h"
#include "itktubeRidgeFFTFilter.h"
#include "itktubeReResampleImageFilter.h"
#include "itktubeSheetnessMeasureImageFilter.h"
#include "itktubeShrinkWithBlendingImageFilter.h"
#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"
#include "itktubeSymmetricEigenVectorAnalysisImageFilter.h"
#include "itktubeTubeEnhancingDiffusion2DImageFilter.h"

#include <itkImage.h>
#include <itkMatrix.h>
#include <itkTubeSpatialObject.h>

int tubeFilteringPrintTest( int itkNotUsed( argc ), char * itkNotUsed(
  argv )[] )
{
  const unsigned int Dimension = 3;
  using ImageType = itk::Image< float, Dimension >;
  using TubeSpatialObjectType = itk::TubeSpatialObject< Dimension >;

  itk::tube::AnisotropicCoherenceEnhancingDiffusionImageFilter< ImageType,
    ImageType > ::Pointer acedif =
    itk::tube::AnisotropicCoherenceEnhancingDiffusionImageFilter< ImageType,
    ImageType >::New();
  std::cout << "-------------acedif" << acedif << std::endl;

  itk::tube::AnisotropicDiffusionTensorFunction< ImageType >
    ::Pointer adtf = itk::tube::AnisotropicDiffusionTensorFunction<
    ImageType >::New();
  std::cout << "-------------adtf" << adtf << std::endl;

  itk::tube::AnisotropicEdgeEnhancementDiffusionImageFilter< ImageType,
    ImageType >::Pointer aeedif =
    itk::tube::AnisotropicEdgeEnhancementDiffusionImageFilter< ImageType,
    ImageType >::New();
  std::cout << "-------------aeedif" << aeedif << std::endl;

   itk::tube::AnisotropicHybridDiffusionImageFilter< ImageType, ImageType >
    ::Pointer ahdif =
    itk::tube::AnisotropicHybridDiffusionImageFilter< ImageType, ImageType >
    ::New();
  std::cout << "-------------ahdif" << ahdif << std::endl;

  using CropImageFilter = itk::tube::CropImageFilter< ImageType,
    ImageType >;
  CropImageFilter::Pointer cropImage = CropImageFilter::New();
  std::cout << "-------------cropImage" << cropImage << std::endl;

  using CVTImageFilter = itk::tube::CVTImageFilter< ImageType >;
  CVTImageFilter::Pointer cvtImage = CVTImageFilter::New();
  std::cout << "-------------cvtImage" << cvtImage << std::endl;

  using ExtractTubePointsSpatialObjectFilterType = itk::tube::ExtractTubePointsSpatialObjectFilter<
    TubeSpatialObjectType >;
  ExtractTubePointsSpatialObjectFilterType::Pointer etpsof =
    ExtractTubePointsSpatialObjectFilterType::New();
  std::cout << "-------------etpsof" << etpsof << std::endl;

  using FFTGaussianDerivativeIFFTFilter = itk::tube::FFTGaussianDerivativeIFFTFilter< ImageType >;
  FFTGaussianDerivativeIFFTFilter::Pointer fgdif =
    FFTGaussianDerivativeIFFTFilter::New();
  std::cout << "-------------fgdif " << fgdif << std::endl;

  using VesselConnectivityFilterType = itk::tube::MinimumSpanningTreeVesselConnectivityFilter<
    Dimension >;
  VesselConnectivityFilterType::Pointer mstvcf =
    VesselConnectivityFilterType::New();
  std::cout << "-------------mstvcf " << mstvcf << std::endl;

  using RidgeFFTFilterType = itk::tube::RidgeFFTFilter< ImageType >;
  RidgeFFTFilterType::Pointer rfif = RidgeFFTFilterType::New();
  std::cout << "-------------rfif " << rfif << std::endl;

  using ReResampleImageFilterType = itk::tube::ReResampleImageFilter< float, 3 >;
  ReResampleImageFilterType::Pointer rrif = ReResampleImageFilterType::New();
  std::cout << "-------------rrif " << rrif << std::endl;

  itk::tube::StructureTensorRecursiveGaussianImageFilter< ImageType >
    ::Pointer strgif =
    itk::tube::StructureTensorRecursiveGaussianImageFilter< ImageType >
    ::New();
  std::cout << "-------------strgif" << strgif << std::endl;

  using EigenValueArrayType = itk::FixedArray< float, 3 >;
  using EigenVectorMatrixType = itk::Matrix< float, 3, 3 >;
  using EigenValueImageType = itk::Image< EigenValueArrayType, 3>;
  using EigenVectorImageType = itk::Image< EigenVectorMatrixType, 3>;

  typedef itk::tube::StructureTensorRecursiveGaussianImageFilter<
    ImageType >::OutputImageType SymmetricSecondRankTensorImageType;
  itk::tube::SymmetricEigenVectorAnalysisImageFilter<
    SymmetricSecondRankTensorImageType,
    EigenValueImageType,
    EigenVectorImageType >::Pointer sevaif =
    itk::tube::SymmetricEigenVectorAnalysisImageFilter<
    SymmetricSecondRankTensorImageType, EigenValueImageType,
    EigenVectorImageType >::New();
  std::cout << "-------------sevaif" << sevaif << std::endl;

  itk::tube::TubeEnhancingDiffusion2DImageFilter< float, 2 >::Pointer
    vesselEnahncingObj = itk::tube::TubeEnhancingDiffusion2DImageFilter<
    float, 2 >::New();
  std::cout << "-------------TubeEnhancingDiffusion2DImageFilter"
    << vesselEnahncingObj << std::endl;

  itk::tube::SheetnessMeasureImageFilter< float >::Pointer
    sheetnessMeasureImageFilterObj =
    itk::tube::SheetnessMeasureImageFilter< float >::New();
  std::cout << "-------------SheetnessMeasureImageFilter"
    << sheetnessMeasureImageFilterObj << std::endl;

  itk::tube::ShrinkWithBlendingImageFilter< ImageType, ImageType >::Pointer
    shrinkUsingMaxImageFilterObj =
    itk::tube::ShrinkWithBlendingImageFilter< ImageType, ImageType >::New();
  std::cout << "-------------ShrinkWithBlendingImageFilter"
    << shrinkUsingMaxImageFilterObj << std::endl;

  return EXIT_SUCCESS;
}
