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
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

project( TubeTKApplications )

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/test/Baseline )
if( NOT ITK_TEST_OUTPUT_DIR )
  make_directory( ${TubeTK_BINARY_DIR}/Testing/Temporary )
  set( ITK_TEST_OUTPUT_DIR ${TubeTK_BINARY_DIR}/Testing/Temporary )
endif()

set( TubeTK_CompareTextFiles_EXE 
  "tubeCompareTextFilesTestDriver")
set( TubeTK_CompareImages_EXE 
  "tubeCompareImagesTestDriver" )

set( SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/bin" )
set( SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/lib" )
set( SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/lib" )

find_package( SlicerExecutionModel REQUIRED
  HINTS ./SlicerExecutionModel-build ../SlicerExecutionModel-build
        ../../SlicerExecutionModel-build )
include( ${SlicerExecutionModel_USE_FILE} )
mark_as_advanced( GENERATECLP_EXE )

add_subdirectory( AtlasBuilderUsingIntensity )
add_subdirectory( ComputeBinaryImageSimilarityMetrics )
add_subdirectory( ComputeImageSimilarityMetrics )
add_subdirectory( ComputeImageStatistics )
add_subdirectory( ComputeTrainingMask )
add_subdirectory( ComputeTubeFlyThroughImage )
add_subdirectory( ComputeTubeGraphProbability )
add_subdirectory( ComputeTubeMeasures )
add_subdirectory( ComputeTubeProbability )
add_subdirectory( ConvertCSVToImages )
add_subdirectory( Convert3DImagesTo4DImage )
add_subdirectory( Convert4DImageTo3DImages )
add_subdirectory( ConvertImagesToCSV )
add_subdirectory( ConvertShrunkenSeedImageToList )
add_subdirectory( ConvertTRE )
add_subdirectory( ConvertSpatialGraphToImage )
add_subdirectory( ConvertTubesToDensityImage )
add_subdirectory( ConvertTubesToImage )
add_subdirectory( ConvertTubesToTubeTree )
add_subdirectory( ConvertTubesToTubeGraph )
add_subdirectory( CropImage )
add_subdirectory( CropTubes )
add_subdirectory( DeblendTomosynthesisSlicesUsingPrior )
#add_subdirectory( DICOMSeriesConverter )
add_subdirectory( EnhanceCoherenceAndEdgesUsingDiffusion )
add_subdirectory( EnhanceCoherenceUsingDiffusion )
add_subdirectory( EnhanceContrastUsingAHE )
add_subdirectory( EnhanceContrastUsingPrior )
add_subdirectory( EnhanceEdgesUsingDiffusion )
add_subdirectory( EnhanceTubesUsingDiffusion )
add_subdirectory( EnhanceTubesUsingDiscriminantAnalysis )
add_subdirectory( EnhanceUsingDiscriminantAnalysis )
add_subdirectory( EnhanceUsingNJetDiscriminantAnalysis )
add_subdirectory( ImageMath )
add_subdirectory( MergeAdjacentImages )
add_subdirectory( MergeTubeGraphs )
add_subdirectory( ResampleImage )
add_subdirectory( ResampleTubes )
add_subdirectory( RegisterImages )
add_subdirectory( RegisterTubesToImage )
add_subdirectory( SampleCLIApplication )
add_subdirectory( SegmentBinaryImageSkeleton3D )
add_subdirectory( SegmentConnectedComponents )
add_subdirectory( SegmentConnectedComponentsUsingParzenPDFs )
add_subdirectory( SegmentTubes )
add_subdirectory( SegmentTubeUsingMinimalPath )
add_subdirectory( SegmentUsingOtsuThreshold )
add_subdirectory( ShrinkImageWithBlending )
add_subdirectory( SimulateAcquisitionArtifactsUsingPrior )
add_subdirectory( TubeMath )

if( TubeTK_USE_VTK )
  add_subdirectory( ConvertTubesToSurface )
  add_subdirectory( ComputeTubeTortuosityMeasures )
  add_subdirectory( RegisterUsingSlidingGeometries )
endif()


####
# NOT CONVERTED YET
####
#add_subdirectory( ComputeRegionSignatures )
#add_subdirectory( SegmentUsingQuantileThreshold )
#add_subdirectory( TransferLabelsToRegions )
#add_subdirectory( ComputeTubeGraphSimilarityKernelMatrix )
