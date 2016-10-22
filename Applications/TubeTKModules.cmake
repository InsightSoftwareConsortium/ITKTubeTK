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

set( proj Applications )

set( TubeTK_${proj}_MODULES
  AtlasBuilderUsingIntensity
  ComputeBinaryImageSimilarityMetrics
  ComputeImageSimilarityMetrics
  ComputeImageStatistics
  ComputeImageToTubeRigidMetricImage
  ComputeSegmentTubesParameters
  ComputeTrainingMask
  ComputeTubeFlyThroughImage
  ComputeTubeGraphProbability
  ComputeTubeMeasures
  ComputeTubeProbability
  ConvertCSVToImages
  ConvertImagesToCSV
  ConvertShrunkenSeedImageToList
  ConvertToMetaImage
  ConvertTRE
  ConvertSpatialGraphToImage
  ConvertTubesToDensityImage
  ConvertTubesToImage
  ConvertTubesToTubeTree
  ConvertTubesToTubeGraph
  CropImage
  CropTubes
  DeblendTomosynthesisSlicesUsingPrior
  EnhanceCoherenceAndEdgesUsingDiffusion
  EnhanceCoherenceUsingDiffusion
  EnhanceContrastUsingPrior
  EnhanceEdgesUsingDiffusion
  EnhanceTubesUsingDiffusion
  EnhanceUsingDiscriminantAnalysis
  EnhanceUsingNJetDiscriminantAnalysis
  ExtractMetricImageSlice
  ImageMath
  MergeAdjacentImages
  MergeTubeGraphs
  ResampleImage
  ResampleTubes
  RegisterImages
  RegisterImageToTubesUsingRigidTransform
  RegisterUsingImageCenters
  SampleCLIApplication
  SegmentBinaryImageSkeleton
  SegmentConnectedComponents
  SegmentConnectedComponentsUsingParzenPDFs
  SegmentTubes
  SegmentTubeUsingMinimalPath
  SegmentUsingOtsuThreshold
  ShrinkImage
  SimulateAcquisitionArtifactsUsingPrior
  TubeMath
  TreeMath )

set( TubeTK_${proj}_Boost_MODULES )
if( TubeTK_USE_BOOST )
  set( TubeTK_${proj}_Boost_MODULES
    ComputeRegionSignatures
    SegmentUsingQuantileThreshold
    TransferLabelsToRegions )
  if( TubeTK_USE_JsonCpp )
    list(TubeTK_${proj}_Boost_MODULES
      ComputeTubeGraphSimilarityKernelMatrix )
  endif( TubeTK_USE_JsonCpp )
  list( APPEND TubeTK_${proj}_MODULES
    ${TubeTK_${proj}_Boost_MODULES} )
endif( TubeTK_USE_BOOST )

set( TubeTK_${proj}_LIBSVM_MODULES )
if( TubeTK_USE_LIBSVM )
  set( TubeTK_${proj}_LIBSVM_MODULES
    EnhanceTubesUsingDiscriminantAnalysis )
  list( APPEND TubeTK_${proj}_MODULES
    ${TubeTK_${proj}_LIBSVM_MODULES} )
endif( TubeTK_USE_LIBSVM )

set( TubeTK_${proj}_VTK_MODULES )
if( TubeTK_USE_VTK )
  set( TubeTK_${proj}_VTK_MODULES
    ConvertTubesToSurface
    ComputeTubeTortuosityMeasures
    RegisterUsingSlidingGeometries )
  list( APPEND TubeTK_${proj}_MODULES
    ${TubeTK_${proj}_VTK_MODULES} )
endif( TubeTK_USE_VTK )

if( NOT TubeTK_SOURCE_DIR )
  find_package( TubeTK REQUIRED )
  include( ${TubeTK_USE_FILE} )
endif( NOT TubeTK_SOURCE_DIR )

if( ${CMAKE_PROJECT_NAME} STREQUAL "Slicer" )
  include( ${CMAKE_CURRENT_SOURCE_DIR}/../CMake/TubeTKMacroAddModules.cmake )
else()
  include(TubeTKMacroAddModules)
endif( ${CMAKE_PROJECT_NAME} STREQUAL "Slicer" )
TubeTKMacroAddModules( MODULES ${TubeTK_${proj}_MODULES} )
