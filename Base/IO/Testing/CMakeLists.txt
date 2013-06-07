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

include_regular_expression( "^.*$" )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/IO
    )
endif()

set( BASE_IO_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseIOTests )

set( BASE_IO_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseIOHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseIOTests_SRCS
  tubeBaseIOPrintTest.cxx
  itkTubeMetaLDATest.cxx
  itkTubeMetaNJetLDATest.cxx
  itkTubeMetaRidgeSeedTest.cxx
  itkTubeMetaTubeParamsTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/IO )

add_executable( tubeBaseIOHeaderTest
  tubeBaseIOHeaderTest.cxx )
target_link_libraries( tubeBaseIOHeaderTest
  ${ITK_LIBRARIES} )

add_executable( tubeBaseIOTests
  tubeBaseIOTests.cxx
  ${tubeBaseIOTests_SRCS} )
target_link_libraries( tubeBaseIOTests
  ${ITK_LIBRARIES} TubeIO )

add_test( tubeBaseIOHeaderTest
  ${BASE_IO_HEADER_TEST} )

add_test( tubeBaseIOPrintTest
  ${BASE_IO_TESTS}
  tubeBaseIOPrintTest )

add_test( itkTubeMetaLDATest
  ${BASE_IO_TESTS}
  itkTubeMetaLDATest
    ${TEMP}/itkTubeMetaLDATest.mlda
  )

add_test( itkTubeMetaNJetLDATest
  ${BASE_IO_TESTS}
  itkTubeMetaNJetLDATest
    ${TEMP}/itkTubeMetaNJetLDATest.mnda
  )

add_test( itkTubeMetaRidgeSeedTest
  ${BASE_IO_TESTS}
  itkTubeMetaRidgeSeedTest
    ${TEMP}/itkTubeMetaRidgeSeedTest.mrs
  )

add_test( itkTubeMetaTubeParamsTest
  ${BASE_IO_TESTS}
  itkTubeMetaTubeParamsTest
    ${TEMP}/itkTubeMetaTubeParamsTest.mtp
  )
