##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
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

project( TubeTKFiltering )

find_package( ITK REQUIRED )
set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_Filtering_H_Files
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.h
  itktubeAnisotropicDiffusionTensorFunction.h
  itktubeAnisotropicDiffusionTensorImageFilter.h
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h
  itktubeAnisotropicHybridDiffusionImageFilter.h
  itktubeExtractTubePointsSpatialObjectFilter.h
  itktubeFFTGaussianDerivativeIFFTFilter.h
  itktubeGaussianDerivativeImageSource.h
  itktubePadImageFilter.h
  itktubeRegionFromReferenceImageFilter.h
  itktubeSheetnessMeasureImageFilter.h
  itktubeSpatialObjectSource.h
  itktubeSpatialObjectToSpatialObjectFilter.h
  itktubeStructureTensorRecursiveGaussianImageFilter.h
  itktubeSubSampleTubeSpatialObjectFilter.h
  itktubeSubSampleTubeTreeSpatialObjectFilter.h
  itktubeSymmetricEigenVectorAnalysisImageFilter.h
  itktubeTortuositySpatialObjectFilter.h
  itktubeTubeEnhancingDiffusion2DImageFilter.h )

set( TubeTK_Base_Filtering_HXX_Files
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.hxx
  itktubeAnisotropicDiffusionTensorFunction.hxx
  itktubeAnisotropicDiffusionTensorImageFilter.hxx
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.hxx
  itktubeAnisotropicHybridDiffusionImageFilter.hxx
  itktubeExtractTubePointsSpatialObjectFilter.hxx
  itktubeFFTGaussianDerivativeIFFTFilter.hxx
  itktubeGaussianDerivativeImageSource.hxx
  itktubePadImageFilter.hxx
  itktubeRegionFromReferenceImageFilter.hxx
  itktubeSheetnessMeasureImageFilter.hxx
  itktubeSpatialObjectSource.hxx
  itktubeSpatialObjectToSpatialObjectFilter.hxx
  itktubeStructureTensorRecursiveGaussianImageFilter.hxx
  itktubeSubSampleTubeSpatialObjectFilter.hxx
  itktubeSubSampleTubeTreeSpatialObjectFilter.hxx
  itktubeTortuositySpatialObjectFilter.h
  itktubeTubeEnhancingDiffusion2DImageFilter.hxx )

add_custom_target( TubeTKFiltering SOURCES
  ${TubeTK_Base_Filtering_H_Files}
  ${TubeTK_Base_Filtering_HXX_Files} )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Filtering_H_Files}
  ${TubeTK_Base_Filtering_HXX_Files}
  DESTINATION include )
