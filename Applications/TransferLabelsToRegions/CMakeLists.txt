##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 ( the "License" );
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
# Style Guidelines
#   * Place a space after ( and before )
#   * Test input and output images should be in .mha format
#   * Test coverage should be > 85% before committing
#   * Test outputs and baselines should be pre-pended with ${MODULE_NAME}
#   * Use appropriate TubeTK comment blocks at the top every file
#   * Follow our naming convention - if you are uncertain, ask
#   * Use the tube namespace for all classes
#   * Your code is a part of TubeTK - don't be "different" - it makes the code
#        harder to maintain and understand, even if you think your way is
#        "better."

set( MODULE_NAME TransferLabelsToRegions )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find Boost
set(Boost_ADDITIONAL_VERSIONS "1.51")
find_package(Boost 1.49 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#-----------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -msse3" ${CMAKE_FLAGS_RELEASE})

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Utilities/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Utilities
    ${TubeTK_SOURCE_DIR}/Base/Filtering
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
