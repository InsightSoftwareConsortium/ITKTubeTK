##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTK_Filtering_H_Files
  Filtering/itkGeneralizedDistanceTransformImageFilter.h
  Filtering/itkImageRegionSplitter.h
  Filtering/itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.h
  Filtering/itktubeAnisotropicDiffusionTensorFunction.h
  Filtering/itktubeAnisotropicDiffusionTensorImageFilter.h
  Filtering/itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h
  Filtering/itktubeAnisotropicHybridDiffusionImageFilter.h
  Filtering/itktubeBinaryThinningImageFilter3D.h
  Filtering/itktubeComputeTubeFlyThroughImageFilter.h
  Filtering/itktubeComputeTubeMeasuresFilter.h
  Filtering/itktubeContrastCostFunction.h
  Filtering/itktubeConvertSpatialGraphToImageFilter.h
  Filtering/itktubeCropImageFilter.h
  Filtering/itktubeCropTubesFilter.h
  Filtering/itktubeCVTImageFilter.h
  Filtering/itktubeDifferenceImageFilter.h
  Filtering/itktubeEnhanceContrastUsingPriorImageFilter.h
  Filtering/itktubeExtractTubePointsSpatialObjectFilter.h
  Filtering/itktubeFFTGaussianDerivativeIFFTFilter.h
  Filtering/itktubeGaussianDerivativeFilter.h
  Filtering/itktubeGaussianDerivativeImageSource.h
  Filtering/itktubeInverseIntensityImageFilter.h
  Filtering/itktubeLimitedMinimumMaximumImageFilter.h
  Filtering/itktubeMinimumSpanningTreeVesselConnectivityFilter.h
  Filtering/itktubePadImageFilter.h
  Filtering/itktubeRegionFromReferenceImageFilter.h
  Filtering/itktubeReResampleImageFilter.h
  Filtering/itktubeSheetnessMeasureImageFilter.h
  Filtering/itktubeShrinkWithBlendingImageFilter.h
  Filtering/itktubeSpatialObjectSource.h
  Filtering/itktubeSpatialObjectToSpatialObjectFilter.h
  Filtering/itktubeStructureTensorRecursiveGaussianImageFilter.h
  Filtering/itktubeSubSampleTubeSpatialObjectFilter.h
  Filtering/itktubeSubSampleSpatialObjectFilter.h
  Filtering/itktubeSymmetricEigenVectorAnalysisImageFilter.h
  Filtering/itktubeTortuositySpatialObjectFilter.h
  Filtering/itktubeTubeEnhancingDiffusion2DImageFilter.h
  Filtering/itktubeTubeSpatialObjectToDensityImageFilter.h
  Filtering/itktubeTubeSpatialObjectToImageFilter.h
  Filtering/itktubeTubeSpatialObjectToTubeGraphFilter.h
  Filtering/tubeImageMathFilters.h
  Filtering/tubeTubeMathFilters.h )

set( TubeTK_Filtering_HXX_Files
  Filtering/itkGeneralizedDistanceTransformImageFilter.txx
  Filtering/itkImageRegionSplitter.hxx
  Filtering/itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.hxx
  Filtering/itktubeAnisotropicDiffusionTensorFunction.hxx
  Filtering/itktubeAnisotropicDiffusionTensorImageFilter.hxx
  Filtering/itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.hxx
  Filtering/itktubeAnisotropicHybridDiffusionImageFilter.hxx
  Filtering/itktubeBinaryThinningImageFilter3D.hxx
  Filtering/itktubeComputeTubeFlyThroughImageFilter.hxx
  Filtering/itktubeComputeTubeMeasuresFilter.hxx
  Filtering/itktubeContrastCostFunction.h
  Filtering/itktubeConvertSpatialGraphToImageFilter.hxx
  Filtering/itktubeCropImageFilter.hxx
  Filtering/itktubeCropTubesFilter.hxx
  Filtering/itktubeCVTImageFilter.hxx
  Filtering/itktubeDifferenceImageFilter.hxx
  Filtering/itktubeEnhanceContrastUsingPriorImageFilter.hxx
  Filtering/itktubeExtractTubePointsSpatialObjectFilter.hxx
  Filtering/itktubeFFTGaussianDerivativeIFFTFilter.hxx
  Filtering/itktubeGaussianDerivativeFilter.hxx
  Filtering/itktubeGaussianDerivativeImageSource.hxx
  Filtering/itktubeInverseIntensityImageFilter.hxx
  Filtering/itktubeLimitedMinimumMaximumImageFilter.hxx
  Filtering/itktubeMinimumSpanningTreeVesselConnectivityFilter.hxx
  Filtering/itktubePadImageFilter.hxx
  Filtering/itktubeRegionFromReferenceImageFilter.hxx
  Filtering/itktubeReResampleImageFilter.hxx
  Filtering/itktubeSheetnessMeasureImageFilter.hxx
  Filtering/itktubeShrinkWithBlendingImageFilter.hxx
  Filtering/itktubeSpatialObjectSource.hxx
  Filtering/itktubeSpatialObjectToSpatialObjectFilter.hxx
  Filtering/itktubeStructureTensorRecursiveGaussianImageFilter.hxx
  Filtering/itktubeSubSampleTubeSpatialObjectFilter.hxx
  Filtering/itktubeSubSampleSpatialObjectFilter.hxx
  Filtering/itktubeTortuositySpatialObjectFilter.h
  Filtering/itktubeTubeEnhancingDiffusion2DImageFilter.hxx
  Filtering/itktubeTubeSpatialObjectToDensityImageFilter.hxx
  Filtering/itktubeTubeSpatialObjectToImageFilter.hxx
  Filtering/itktubeTubeSpatialObjectToTubeGraphFilter.hxx
  Filtering/tubeImageMathFilters.hxx
  Filtering/tubeTubeMathFilters.hxx )

set( TubeTK_Filtering_CXX_Files )

list( APPEND TubeTK_SRCS
  ${TubeTK_Filtering_H_Files}
  ${TubeTK_Filtering_HXX_Files}
  ${TubeTK_Filtering_CXX_Files} )
