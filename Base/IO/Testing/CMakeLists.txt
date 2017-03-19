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

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RANDOMFOREST )
  find_package( RandomForest REQUIRED )
  include_directories( ${RandomForest_INCLUDE_DIR} )
endif( TubeTK_USE_RANDOMFOREST )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/IO )
endif( PYTHON_EXECUTABLE )

set( CompareImages_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

set( BASE_IO_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseIOTests> )

set( BASE_IO_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseIOHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseIOTests_SRCS
  tubeBaseIOPrintTest.cxx
  itktubePDFSegmenterParzenIOTest.cxx
  itktubeTubeExtractorIOTest.cxx
  itktubeRidgeSeedFilterIOTest.cxx
  itktubeTubeXIOTest.cxx )
if( TubeTK_USE_LIBSVM )
  list( APPEND tubeBaseIOTests_SRCS
    itktubePDFSegmenterSVMIOTest.cxx )
endif( TubeTK_USE_LIBSVM )
if( TubeTK_USE_RANDOMFOREST )
  list( APPEND tubeBaseIOTests_SRCS
    itktubePDFSegmenterRandomForestIOTest.cxx )
endif( TubeTK_USE_RANDOMFOREST )


add_executable( tubeBaseIOHeaderTest
  tubeBaseIOHeaderTest.cxx )
target_link_libraries( tubeBaseIOHeaderTest
  TubeTKIO ${ITK_LIBRARIES} ${LIBSVM_LIBRARIES} )
if( TubeTK_USE_LIBSVM )
  target_link_libraries( tubeBaseIOHeaderTest ${LIBSVM_LIBRARIES} )
endif( TubeTK_USE_LIBSVM )

set( tubeBaseIO_External_Libraries )
if( TubeTK_USE_ARRAYFIRE )
  set( tubeBaseIO_External_Libraries ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_ARRAYFIRE )
if( TubeTK_USE_LIBSVM )
  set( tubeBaseIO_External_Libraries ${LIBSVM_LIBRARIES} )
endif( TubeTK_USE_LIBSVM )

set( no_install_option )
if( NOT TubeTK_INSTALL_DEVELOPMENT )
  set( no_install_option NO_INSTALL )
endif()

SEMMacroBuildCLI(
  NAME tubeBaseIOTests
  ADDITIONAL_SRCS
    ${tubeBaseIOTests_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${tubeBaseIO_External_Libraries}
    TubeCLI
    TubeTKMetaIO
    TubeTKCommon
    TubeTKNumerics
    TubeTKIO
    ${ITK_LIBRARIES}
    ITKIOMeta
  EXECUTABLE_ONLY
  ${no_install_option}
  )

add_test( NAME tubeBaseIOHeaderTest
  COMMAND ${BASE_IO_HEADER_TEST} )

add_test( NAME tubeBaseIOPrintTest
  COMMAND ${BASE_IO_TESTS}
    tubeBaseIOPrintTest )

ExternalData_Add_Test( TubeTKData
  NAME itktubePDFSegmenterParzenIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubePDFSegmenterParzenIOTest.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mha
    itktubePDFSegmenterParzenIOTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
      ${TEMP}/itktubePDFSegmenterParzenIOTest.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest.mpd
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mpd )

if( TubeTK_USE_LIBSVM )
  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterSVMIOTest
    COMMAND ${BASE_IO_TESTS}
      --compareNumberOfPixelsTolerance 20
      --compare ${TEMP}/itktubePDFSegmenterSVMIOTest.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.mha
      itktubePDFSegmenterSVMIOTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${TEMP}/itktubePDFSegmenterSVMIOTest.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest.msvm
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.msvm )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RANDOMFOREST )
  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterRandomForestIOTest
    COMMAND ${BASE_IO_TESTS}
      --compareNumberOfPixelsTolerance 20
      --compare ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mha
      itktubePDFSegmenterRandomForestIOTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mrf
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mrf )
    endif( TubeTK_USE_RANDOMFOREST )

ExternalData_Add_Test( TubeTKData
  NAME itktubeRidgeSeedFilterIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output2.mha
    --compareNumberOfPixelsTolerance 50
    itktubeRidgeSeedFilterIOTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterIOTest.mrs
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output2.mha )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeExtractorIOTest
  COMMAND ${BASE_IO_TESTS}
    itktubeTubeExtractorIOTest
      DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
      ${TEMP}/itktubeTubeExtractorIOTest.mtp
      ${TEMP}/itktubeTubeExtractorIOTest2.mtp )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeExtractorIOTest-Compare
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
    -t ${TEMP}/itktubeTubeExtractorIOTest.mtp )
set_tests_properties( itktubeTubeExtractorIOTest-Compare PROPERTIES DEPENDS
  itktubeTubeExtractorIOTest )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeExtractorIOTest-Compare2
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
    -t ${TEMP}/itktubeTubeExtractorIOTest2.mtp )
set_tests_properties( itktubeTubeExtractorIOTest-Compare2 PROPERTIES DEPENDS
  itktubeTubeExtractorIOTest )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeXIOTest
  COMMAND ${BASE_IO_TESTS}
    itktubeTubeXIOTest
      DATA{${TubeTK_DATA_ROOT}/TubeXIOTest.tre}
      ${TEMP}/itktubeTubeXIOTest.tre
      512 512 393
      )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeXIOTest-Compare
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeXIOTest.tre}
    -t ${TEMP}/itktubeTubeXIOTest.tre )
set_tests_properties( itktubeTubeXIOTest-Compare PROPERTIES DEPENDS
  itktubeTubeXIOTest )
