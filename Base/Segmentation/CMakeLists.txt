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

set( TubeTK_Base_Segmentation_H_Files
  itkTubeCVTImageFilter.h
  itkTubeGenerateLDA.h
  itkTubePDFSegmenter.h
  itkTubeRadiusExtractor.h
  itkTubeRidgeExtractor.h
  itkTubeTubeExtractor.h
  itkTubeTubeNetExtractor.h
  )

set(TubeTK_Base_Segmentation_TXX_Files
  itkTubeCVTImageFilter.txx
  itkTubeGenerateLDA.txx
  itkTubePDFSegmenter.txx
  itkTubeRadiusExtractor.txx
  itkTubeRidgeExtractor.txx
  itkTubeTubeExtractor.txx
  itkTubeTubeNetExtractor.txx
  )

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES
  #${TubeTK_Base_Segmentation_H_Files}
  #${TubeTK_Base_Segmentation_TXX_Files}
  #DESTINATION include
  #)

