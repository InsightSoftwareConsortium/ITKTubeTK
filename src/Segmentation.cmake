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
  Segmentation/itktubeComputeSegmentTubesParameters.h
  Segmentation/itktubePDFSegmenterBase.h
  Segmentation/itktubePDFSegmenterParzen.h
  Segmentation/itktubeRadiusExtractor2.h
  Segmentation/itktubeRidgeExtractor.h
  Segmentation/itktubeSegmentBinaryImageSkeleton.h
  Segmentation/itktubeSegmentBinaryImageSkeleton3D.h
  Segmentation/itktubeSegmentTubeUsingMinimalPathFilter.h
  Segmentation/itktubeTubeExtractor.h
  Segmentation/itktubeRidgeSeedFilter.h
  Segmentation/itktubeComputeTrainingMask.h)

set( TubeTK_Segmentation_HXX_Files
  Segmentation/itktubeComputeSegmentTubesParameters.hxx
  Segmentation/itktubePDFSegmenterBase.hxx
  Segmentation/itktubePDFSegmenterParzen.hxx
  Segmentation/itktubePDFSegmenterRandomForest.hxx
  Segmentation/itktubeRadiusExtractor2.hxx
  Segmentation/itktubeRidgeExtractor.hxx
  Segmentation/itktubeSegmentBinaryImageSkeleton.hxx
  Segmentation/itktubeSegmentBinaryImageSkeleton3D.hxx
  Segmentation/itktubeSegmentTubeUsingMinimalPathFilter.hxx
  Segmentation/itktubeTubeExtractor.hxx
  Segmentation/itktubeRidgeSeedFilter.hxx
  Segmentation/itktubeComputeTrainingMask.hxx)

set( TubeTK_Segmentation_CXX_Files )

list( APPEND TubeTK_SRCS
  ${TubeTK_Segmentation_Files}
  ${TubeTK_Segmentation_X_Files}
  ${TubeTK_Segmentation_X_Files} )
