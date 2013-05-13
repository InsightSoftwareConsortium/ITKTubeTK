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

add_subdirectory( AtlasBuilderUsingIntensity )
add_subdirectory( ComputeBinaryImageSimilarityMetrics )
add_subdirectory( ComputeImageSimilarityMetrics )
add_subdirectory( ComputeImageStatisticsUsingMask )
if( TubeTK_USE_Boost )
  add_subdirectory( ComputeImageQuantiles )
endif( TubeTK_USE_Boost )
add_subdirectory( ComputeTubeGraphProbability )
add_subdirectory( ComputeTubeProbability )
# TODO: Has not been updated for GDCM 2.0.
if( 0 )
  add_subdirectory( ConvertDICOMSeriesToImages )
endif()
add_subdirectory( ConvertToMetaImage )
add_subdirectory( CropImage )
add_subdirectory( DeblendTomosynthesisSlicesUsingPrior )
add_subdirectory( EnhanceCoherenceAndEdgesUsingDiffusion )
add_subdirectory( EnhanceCoherenceUsingDiffusion )
add_subdirectory( EnhanceContrastUsingPrior )
add_subdirectory( EnhanceEdgesUsingDiffusion )
add_subdirectory( EnhanceTubesUsingDiffusion )
add_subdirectory( EnhanceUsingDiscriminantAnalysis )
add_subdirectory( EnhanceUsingNJetDiscriminantAnalysis )
add_subdirectory( ImageMath )
add_subdirectory( MergeAdjacentImages )
add_subdirectory( MergeTubeGraphs )
add_subdirectory( ResampleImage )
add_subdirectory( RegisterImageToTubesUsingRigidTransform )
add_subdirectory( RegisterUsingImageCenters )
if( TubeTK_USE_VTK )
  add_subdirectory( RegisterUsingSlidingGeometries )
endif( TubeTK_USE_VTK )
add_subdirectory( SampleCLIApplication )
add_subdirectory( SegmentBinaryImageSkeleton )
add_subdirectory( SegmentConnectedComponents )
add_subdirectory( SegmentConnectedComponentsUsingParzenPDFs )
#add_subdirectory( SegmentTubeSeeds )
#add_subdirectory( SegmentTubes )
add_subdirectory( SegmentUsingOtsuThreshold )
add_subdirectory( SimulateAcquisitionArtifactsUsingPrior )
if( TubeTK_USE_Boost )
  add_subdirectory( TransferLabelsToRegions )
  add_subdirectory( TubeGraphKernel )
endif( TubeTK_USE_Boost )
add_subdirectory( TubeGraphToImage )
add_subdirectory( TubesToDensityImage )
add_subdirectory( TubeToTubeGraph )
add_subdirectory( TubeTransform )
