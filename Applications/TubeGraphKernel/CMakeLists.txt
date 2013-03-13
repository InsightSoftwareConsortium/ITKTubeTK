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

set( MODULE_NAME TubeGraphKernel )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find Boost
set(Boost_ADDITIONAL_VERSIONS "1.51")
find_package(Boost 1.49 COMPONENTS system filesystem REQUIRED)

#-----------------------------------------------------------------------------
if( CMAKE_COMPILER_IS_GNUCXX )
  set(CMAKE_CXX_FLAGS "-O3 -ffast-math -msse3")
endif()

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Utilities/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}
  ADDITIONAL_SRCS ShortestPathKernel.cxx WLSubtreeKernel.cxx GraphKernel.cxx
  INCLUDE_DIRECTORIES
    ${Boost_INCLUDE_DIRS}
    ${TubeTK_SOURCE_DIR}/Utilities
    ${TubeTK_SOURCE_DIR}/Base/Filtering
  )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
