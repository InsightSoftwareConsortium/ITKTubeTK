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

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )


find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/IO )
endif( PYTHON_EXECUTABLE )

set( IMAGECOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:ImageCompareCommand> )

set( TEXTCOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:TextCompareCommand> )

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
  itktubeMetaTubeExtractorTest.cxx
  itktubePDFSegmenterIOTest.cxx
  itktubeRidgeSeedFilterIOTest.cxx
  itktubeTubeExtractorIOTest.cxx
  itktubeTubeXIOTest.cxx )

add_executable( tubeBaseIOHeaderTest
  tubeBaseIOHeaderTest.cxx )
target_link_libraries( tubeBaseIOHeaderTest
  TubeTKIO ${ITK_LIBRARIES} )

set( tubeBaseIO_External_Libraries )
if( TubeTK_USE_GPU_ARRAYFIRE )
  set( tubeBaseIO_External_Libraries ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_GPU_ARRAYFIRE )

SEMMacroBuildCLI(
  NAME tubeBaseIOTests
  ADDITIONAL_SRCS
    ${tubeBaseIOTests_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${tubeBaseIO_External_Libraries}
    TubeCLI
    TubeTKCommon
    TubeTKNumerics
    TubeTKIO
    ${ITK_LIBRARIES}
    ITKIOMeta
  INSTALL_RUNTIME_DESTINATION bin
  EXECUTABLE_ONLY
  )

add_test( NAME tubeBaseIOHeaderTest
  COMMAND ${BASE_IO_HEADER_TEST} )

add_test( NAME tubeBaseIOPrintTest
  COMMAND ${BASE_IO_TESTS}
    tubeBaseIOPrintTest )

add_test( NAME itktubeMetaLDATest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaLDATest
      ${TEMP}/itktubeMetaLDATest.mlda )

add_test( NAME itktubeMetaNJetLDATest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaNJetLDATest
      ${TEMP}/itktubeMetaNJetLDATest.mnda )

add_test( NAME itktubeMetaClassPDFTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaClassPDFTest
      ${TEMP}/itktubeMetaClassPDFTest.mha )

add_test( NAME itktubeMetaRidgeSeedTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaRidgeSeedTest
      ${TEMP}/itktubeMetaRidgeSeedTest.mrs )

add_test( NAME itktubeMetaTubeExtractorTest
  COMMAND ${BASE_IO_TESTS}
    itktubeMetaTubeExtractorTest
      ${TEMP}/itktubeMetaTubeExtractorTest.mtp )

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

Midas3FunctionAddTest( NAME itktubeTubeExtractorIOTest
  COMMAND ${BASE_IO_TESTS}
    itktubeTubeExtractorIOTest
      MIDAS{TubeExtractorIOTest.mtp.md5}
      ${TEMP}/itktubeTubeExtractorIOTest.mtp
      ${TEMP}/itktubeTubeExtractorIOTest2.mtp )

Midas3FunctionAddTest( NAME itktubeTubeExtractorIOTest-Compare
  COMMAND ${TEXTCOMPARE_EXE}
    -b MIDAS{TubeExtractorIOTest.mtp.md5}
    -t ${TEMP}/itktubeTubeExtractorIOTest.mtp )
set_property( TEST itktubeTubeExtractorIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeExtractorIOTest )

Midas3FunctionAddTest( NAME itktubeTubeExtractorIOTest-Compare2
  COMMAND ${TEXTCOMPARE_EXE}
    -b MIDAS{TubeExtractorIOTest.mtp.md5}
    -t ${TEMP}/itktubeTubeExtractorIOTest2.mtp )
set_property( TEST itktubeTubeExtractorIOTest-Compare2 APPEND PROPERTY DEPENDS
  itktubeTubeExtractorIOTest )

Midas3FunctionAddTest( NAME itktubeTubeXIOTest
  COMMAND ${BASE_IO_TESTS}
    itktubeTubeXIOTest
      MIDAS{TubeXIOTest.tre.md5}
      ${TEMP}/itktubeTubeXIOTest.tre
      512 512 393
      )

Midas3FunctionAddTest( NAME itktubeTubeXIOTest-Compare
  COMMAND ${TEXTCOMPARE_EXE}
    -b MIDAS{TubeXIOTest.tre.md5}
    -t ${TEMP}/itktubeTubeXIOTest.tre )
set_property( TEST itktubeTubeXIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeXIOTest )
