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

include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)

if( Slicer3_INSTALL_BIN_DIR )
  set( BASE_USTK_TESTS
    ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseUSTKTests )
else()
  set( BASE_USTK_TESTS
    ${TubeTK_BINARY_DIR}/bin/tubeBaseUSTKTests )
endif()

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseUSTK_SRCS
  itkLabelMapToAcousticImpedanceImageFilterTest.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/USTK
  )

add_executable( tubeBaseUSTKTests
  tubeBaseUSTKTests.cxx
  ${tubeBaseUSTK_SRCS}
  )
target_link_libraries( tubeBaseUSTKTests
  #TubeUSTK
  ${ITK_LIBRARIES}
  )

MIDAS_ADD_TEST( itkLabelMapToAcousticImpedanceImageFilterTest
  ${BASE_USTK_TESTS}
  itkLabelMapToAcousticImpedanceImageFilterTest
  ${CMAKE_CURRENT_SOURCE_DIR}/LabelMapToAcousticImpedanceTestLookupTable.csv
  )
