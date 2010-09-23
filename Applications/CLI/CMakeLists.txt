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
add_subdirectory( CoherenceEnhancingAnisotropicDiffusion )
add_subdirectory( ComputeImageSimilarity )
add_subdirectory( Crop )
add_subdirectory( DeblendImagesUsingPrior )
add_subdirectory( EdgeEnhancingAnisotropicDiffusion )
add_subdirectory( ExtractCurves2D )
add_subdirectory( HessianTubeness2D )
add_subdirectory( HybridEnhancingAnisotropicDiffusion )
add_subdirectory( MatchImageWithPrior )
add_subdirectory( Merge )
add_subdirectory( OtsuThreshold )
add_subdirectory( OrImages )
add_subdirectory( PdPfaScoring )
add_subdirectory( RandomlyModifyBinaryMask )
add_subdirectory( Resample )
add_subdirectory( SampleCLIApplication )
add_subdirectory( SimulateTomosynthesis )
if( TubeTK_USE_OpenIGTLink )
  add_subdirectory( SendPositionMessageUsingOpenIGTLink )
  add_subdirectory( ReceivePositionMessageUsingOpenIGTLink )
endif( TubeTK_USE_OpenIGTLink )
add_subdirectory( GenerateFeaturesForWeka )
add_subdirectory( Skeletonize )
add_subdirectory( WekaToImage )
add_subdirectory( ComputeZMap )
add_subdirectory( PDFSegmenter )
add_subdirectory( ComputeContrastImage )
if( TubeTK_USE_FANN )
  add_subdirectory( FANNClassify )
endif()


