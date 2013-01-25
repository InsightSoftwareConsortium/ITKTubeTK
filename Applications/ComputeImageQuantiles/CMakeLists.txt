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
project( tubeImageQuantiles )

cmake_minimum_required( VERSION 2.8 )

# Help finding BOOST
set(BOOST_ROOT  /software/boost_1_51_0/)
set(Boost_ADDITIONAL_VERSIONS "1.51")

find_package(Boost 1.49 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find GenerateCLP
find_package( GenerateCLP REQUIRED )
include( ${GenerateCLP_USE_FILE} )


set(USER_CMAKE_CXX_FLAGS "-ftemplate-depth-100")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USER_CMAKE_CXX_FLAGS}")

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  )

set( PROJECT_SOURCE ${PROJECT_NAME}.cxx )

generateclp( PROJECT_SOURCE ${PROJECT_NAME}.xml )

add_executable( ${PROJECT_NAME} ${PROJECT_SOURCE} )
target_link_libraries( ${PROJECT_NAME} ${ITK_LIBRARIES} )
