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

find_package( Tclsh )
if( TCL_TCLSH )
  add_test( PrintSelf-BaseIO ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl
    ${TubeTK_SOURCE_DIR}/Base/IO )
endif( TCL_TCLSH )

set( BASE_IO_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseIOTests )

set( BASE_IO_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseIOHeaderTest )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/IO/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseIOTests_SRCS
  tubeBaseIOPrintTest.cxx
  itkTubeMetaLDATest.cxx
  itkTubeMetaNJetLDATest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/IO )

add_executable( tubeBaseIOHeaderTest
  tubeBaseIOHeaderTest.cxx )
target_link_libraries( tubeBaseIOHeaderTest
  ITKIO ITKNumerics ITKBasicFilters )

add_executable( tubeBaseIOTests
  tubeBaseIOTests.cxx
  ${tubeBaseIOTests_SRCS} )
target_link_libraries( tubeBaseIOTests
  ITKIO ITKNumerics ITKBasicFilters TubeIO )

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
