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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( MODULE_NAME AtlasBuilderUsingIntensity )

set( AtlasBuilderUsingIntensity_H_Files
  itktubeCompleteImageResampleFilter.h
  itktubeMeanAndSigmaImageBuilder.h
  itktubeMinimizeImageSizeFilter.h
  itktubeRobustMeanAndSigmaImageBuilder.h
  tubeAtlasSummation.h )

set( AtlasBuilderUsingIntensity_HXX_Files
  itktubeCompleteImageResampleFilter.hxx
  itktubeMeanAndSigmaImageBuilder.hxx
  itktubeMinimizeImageSizeFilter.hxx
  itktubeRobustMeanAndSigmaImageBuilder.hxx )

add_custom_target( AtlasBuilderUsingIntensityInclude SOURCES
  ${AtlasBuilderUsingIntensity_H_Files}
  ${AtlasBuilderUsingIntensity_HXX_Files} )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  ADDITIONAL_SRCS tubeAtlasSummation.cxx
  LOGO_HEADER ${CMAKE_CURRENT_LIST_DIR}/../../../docs/TubeTKLogo.h
  TARGET_LIBRARIES
    TubeTK ${ITK_LIBRARIES}
  )

if( BUILD_TESTING )
  if( PYTHON_EXECUTABLE )
    add_subdirectory( Testing )
  endif( PYTHON_EXECUTABLE )
endif( BUILD_TESTING )
