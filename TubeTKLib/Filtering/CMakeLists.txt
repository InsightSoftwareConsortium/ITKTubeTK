##############################################################################
#
# Library:   TubeTKLib
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTKLib_Base_Filtering_H_Files
  itkGeneralizedDistanceTransformImageFilter.h
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.h
  itktubeAnisotropicDiffusionTensorFunction.h
  itktubeAnisotropicDiffusionTensorImageFilter.h
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h
  itktubeAnisotropicHybridDiffusionImageFilter.h
  itktubeBinaryThinningImageFilter3D.h
  itktubeComputeTubeFlyThroughImageFilter.h
  itktubeComputeTubeMeasuresFilter.h
  itktubeContrastCostFunction.h
  itktubeConvertSpatialGraphToImageFilter.h
  itktubeCropImageFilter.h
  itktubeCropTubesFilter.h
  itktubeCVTImageFilter.h
  itktubeEnhanceContrastUsingPriorImageFilter.h
  itktubeExtractTubePointsSpatialObjectFilter.h
  itktubeFFTGaussianDerivativeIFFTFilter.h
  itktubeGaussianDerivativeFilter.h
  itktubeGaussianDerivativeImageSource.h
  itktubeInverseIntensityImageFilter.h
  itktubeMinimumSpanningTreeVesselConnectivityFilter.h
  itktubePadImageFilter.h
  itktubeRegionFromReferenceImageFilter.h
  itktubeResampleImageFilter.h
  itktubeSheetnessMeasureImageFilter.h
  itktubeShrinkWithBlendingImageFilter.h
  itktubeSpatialObjectSource.h
  itktubeSpatialObjectToSpatialObjectFilter.h
  itktubeStructureTensorRecursiveGaussianImageFilter.h
  itktubeSubSampleTubeSpatialObjectFilter.h
  itktubeSubSampleTubeTreeSpatialObjectFilter.h
  itktubeSymmetricEigenVectorAnalysisImageFilter.h
  itktubeTortuositySpatialObjectFilter.h
  itktubeTubeEnhancingDiffusion2DImageFilter.h
  itktubeTubeSpatialObjectToDensityImageFilter.h
  itktubeTubeSpatialObjectToImageFilter.h
  itktubeTubeSpatialObjectToTubeGraphFilter.h
  tubeImageFilters.h
  tubeTreeFilters.h )

set( TubeTKLib_Base_Filtering_HXX_Files
  itkGeneralizedDistanceTransformImageFilter.hxx
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.hxx
  itktubeAnisotropicDiffusionTensorFunction.hxx
  itktubeAnisotropicDiffusionTensorImageFilter.hxx
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.hxx
  itktubeAnisotropicHybridDiffusionImageFilter.hxx
  itktubeBinaryThinningImageFilter3D.hxx
  itktubeComputeTubeFlyThroughImageFilter.hxx
  itktubeComputeTubeMeasuresFilter.hxx
  itktubeContrastCostFunction.h
  itktubeConvertSpatialGraphToImageFilter.hxx
  itktubeCropImageFilter.hxx
  itktubeCropTubesFilter.hxx
  itktubeCVTImageFilter.hxx
  itktubeEnhanceContrastUsingPriorImageFilter.hxx
  itktubeExtractTubePointsSpatialObjectFilter.hxx
  itktubeFFTGaussianDerivativeIFFTFilter.hxx
  itktubeGaussianDerivativeFilter.hxx
  itktubeGaussianDerivativeImageSource.hxx
  itktubeInverseIntensityImageFilter.hxx
  itktubeMinimumSpanningTreeVesselConnectivityFilter.hxx
  itktubePadImageFilter.hxx
  itktubeRegionFromReferenceImageFilter.hxx
  itktubeResampleImageFilter.hxx
  itktubeSheetnessMeasureImageFilter.hxx
  itktubeShrinkWithBlendingImageFilter.hxx
  itktubeSpatialObjectSource.hxx
  itktubeSpatialObjectToSpatialObjectFilter.hxx
  itktubeStructureTensorRecursiveGaussianImageFilter.hxx
  itktubeSubSampleTubeSpatialObjectFilter.hxx
  itktubeSubSampleTubeTreeSpatialObjectFilter.hxx
  itktubeTortuositySpatialObjectFilter.h
  itktubeTubeEnhancingDiffusion2DImageFilter.hxx
  itktubeTubeSpatialObjectToDensityImageFilter.hxx
  itktubeTubeSpatialObjectToImageFilter.hxx
  itktubeTubeSpatialObjectToTubeGraphFilter.hxx
  tubeImageFilters.hxx
  tubeTreeFilters.hxx )

# Add GPU Array Fire based filter source files here
if( TubeTK_USE_ARRAYFIRE )
  list( APPEND TubeTKLib_Base_Filtering_H_Files
    itktubeArrayFireGlueUtilities.h
    itktubeGPUArrayFireGaussianDerivativeFilter.h )

  list( APPEND TubeTKLib_Base_Filtering_HXX_Files
    itktubeGPUArrayFireGaussianDerivativeFilter.hxx )
endif()

add_library( TubeTKLibFiltering INTERFACE )
target_include_directories( TubeTKLibFiltering INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  )
target_link_libraries( TubeTKLibFiltering INTERFACE
  TubeTKLibNumerics )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_Filtering_H_Files}
  ${TubeTKLib_Base_Filtering_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS TubeTKLibFiltering )
