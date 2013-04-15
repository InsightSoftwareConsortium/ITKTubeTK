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

#-----------------------------------------------------------------------------
set( MODULE_NAME AtlasBuilderUsingIntensity )
project( ${MODULE_NAME} )

set( AtlasBuilderUsingIntensity_H_Files
  AtlasSummation.h
  itkCompleteImageResampleFilter.h
  itkMeanAndSigmaImageBuilder.h
  itkMinimizeImageSizeFilter.h
  itkRobustMeanAndSigmaImageBuilder.h
  )

set( AtlasBuilderUsingIntensity_TXX_Files
  itkCompleteImageResampleFilter.txx
  itkMeanAndSigmaImageBuilder.txx
  itkMinimizeImageSizeFilter.txx
  itkRobustMeanAndSigmaImageBuilder.txx
  )

add_custom_target( AtlasBuilderUsingIntensityInclude SOURCES
  ${AtlasBuilderUsingIntensity_H_Files}
  ${AtlasBuilderUsingIntensity_TXX_Files}
  )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  ADDITIONAL_SRCS AtlasSummation.cxx
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Utilities/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} TubeObjectDocuments
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Utilities
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Registration
    ${TubeTK_SOURCE_DIR}/Base/ObjectDocuments
  )

#-----------------------------------------------------------------------------
if( BUILD_TESTING )
  find_package( PythonInterp )
  if( PYTHON_EXECUTABLE )
    add_subdirectory( Testing )
  endif( PYTHON_EXECUTABLE )
endif( BUILD_TESTING )
