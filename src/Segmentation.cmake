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

set( TubeTK_Segmentation_H_Files
  Segmentation/itktubePDFSegmenterBase.h
  Segmentation/itktubePDFSegmenterParzen.h
  Segmentation/itktubeRadiusExtractor2.h
  Segmentation/itktubeRadiusExtractor3.h
  Segmentation/itktubeRidgeExtractor.h
  Segmentation/itktubeSegmentTubeUsingMinimalPathFilter.h
  Segmentation/itktubeTubeExtractor.h
  Segmentation/itktubeRidgeSeedFilter.h
  Segmentation/itktubeComputeTrainingMaskFilter.h)

set( TubeTK_Segmentation_HXX_Files
  Segmentation/itktubePDFSegmenterBase.hxx
  Segmentation/itktubePDFSegmenterParzen.hxx
  Segmentation/itktubeRadiusExtractor2.hxx
  Segmentation/itktubeRadiusExtractor3.hxx
  Segmentation/itktubeRidgeExtractor.hxx
  Segmentation/itktubeSegmentTubeUsingMinimalPathFilter.hxx
  Segmentation/itktubeTubeExtractor.hxx
  Segmentation/itktubeRidgeSeedFilter.hxx
  Segmentation/itktubeComputeTrainingMaskFilter.hxx)

set( TubeTK_Segmentation_CXX_Files )

list( APPEND TubeTK_SRCS
  ${TubeTK_Segmentation_H_Files}
  ${TubeTK_Segmentation_HXX_Files}
  ${TubeTK_Segmentation_CXX_Files} )
