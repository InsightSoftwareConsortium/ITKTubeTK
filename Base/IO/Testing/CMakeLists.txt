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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/IO )
endif( PYTHON_EXECUTABLE )

set( BASE_IO_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseIOTests> )

set( BASE_IO_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseIOHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseIOTests_SRCS
  tubeBaseIOPrintTest.cxx
  itktubeMetaLDATest.cxx
  itktubeMetaNJetLDATest.cxx
  itktubeMetaClassPDFTest.cxx
  itktubeMetaRidgeSeedTest.cxx
  itktubeMetaTubeParamsTest.cxx
  itktubePDFSegmenterIOTest.cxx
  itktubeRidgeSeedFilterIOTest.cxx )

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
  ${ITK_LIBRARIES} TubeTKIO )

add_test( NAME tubeBaseIOHeaderTest
  COMMAND ${BASE_IO_HEADER_TEST} )

add_test( NAME tubeBaseIOPrintTest
  COMMAND ${BASE_IO_TESTS}
    tubeBaseIOPrintTest )

add_test( NAME itktubeMetaLDATest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaLDATest
      ${TEMP}/itkTubeMetaLDATest.mlda )

add_test( NAME itktubeMetaNJetLDATest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaNJetLDATest
      ${TEMP}/itkTubeMetaNJetLDATest.mnda )

add_test( NAME itktubeMetaClassPDFTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaClassPDFTest
      ${TEMP}/itkTubeMetaClassPDFTest.mha )

add_test( NAME itktubeMetaRidgeSeedTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaRidgeSeedTest
      ${TEMP}/itkTubeMetaRidgeSeedTest.mrs )

add_test( NAME itktubeMetaTubeParamsTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaTubeParamsTest
      ${TEMP}/itkTubeMetaTubeParamsTest.mtp )

Midas3FunctionAddTest( NAME itktubePDFSegmenterIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubePDFSegmenterIOTest.mha
      ${TEMP}/itktubePDFSegmenterIOTest2.mha
    itktubePDFSegmenterIOTest
      MIDAS{ES0015_Large.mha.md5}
      MIDAS{ES0015_Large.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask.mha.md5}
      ${TEMP}/itktubePDFSegmenterIOTest.mha
      ${TEMP}/itktubePDFSegmenterIOTest.mpd
      ${TEMP}/itktubePDFSegmenterIOTest2.mha
      ${TEMP}/itktubePDFSegmenterIOTest2.mpd )

Midas3FunctionAddTest( NAME itktubeRidgeSeedFilterIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output2.mha
    itktubeRidgeSeedFilterIOTest
      MIDAS{im0001.crop.contrast.mha.md5}
      MIDAS{im0001.vk.maskRidge.crop.mha.md5}
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterIOTest.mrs
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output2.mha )
