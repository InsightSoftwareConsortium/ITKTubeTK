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
if( TubeTK_USE_VTK )
  add_subdirectory( tubeAnisotropicDiffusiveDeformableRegistration )
endif( TubeTK_USE_VTK )
add_subdirectory( tubeCoherenceEnhancingAnisotropicDiffusion )
add_subdirectory( tubeComputeImageSimilarity )
add_subdirectory( tubeComputeBinaryImageSimilarity )
add_subdirectory( tubeConnectedComponents )
add_subdirectory( tubeContrastImageUsingPrior )
add_subdirectory( tubeCrop )
add_subdirectory( tubeDeblendImagesUsingPrior )
add_subdirectory( tubeEdgeEnhancingAnisotropicDiffusion )
add_subdirectory( tubeExtractCurves2D )
#add_subdirectory( tubeExtractVessels )
#add_subdirectory( tubeExtractVesselSeeds )
add_subdirectory( tubeHessianTubeness2D )
add_subdirectory( tubeHybridEnhancingAnisotropicDiffusion )
add_subdirectory( tubeLDAGenerator )
add_subdirectory( tubeImageToTubeRigidRegistration )
add_subdirectory( tubeMaskToStats )
add_subdirectory( tubeMatchImageWithPrior )
add_subdirectory( tubeMerge )
add_subdirectory( tubeNJetLDAGenerator )
add_subdirectory( tubeOtsuThreshold )
add_subdirectory( tubePDFSegmenter )
add_subdirectory( tubeResample )
add_subdirectory( tubeSampleCLIApplication )
add_subdirectory( tubeSkeletonize )
add_subdirectory( tubeToGraph )
add_subdirectory( tubeTransform )
add_subdirectory( tubeGraphsToGraph )
add_subdirectory( tubeAtlasSummation )
add_subdirectory( tubeDensityImageRadiusBuilder)
add_subdirectory( tubeDensityProbability )
add_subdirectory( tubeGraphToImage )
add_subdirectory( tubeGraphProbability )
add_subdirectory( tubeCenterPhantom )
if ( TubeTK_USE_Boost )
  add_subdirectory( tubeGraphKernel )
  add_subdirectory( tubeImageQuantiles )
  add_subdirectory( tubeAtlasLabelTransfer )
endif()
