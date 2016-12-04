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

if( TubeTK_USE_LIBSVM )
  find_package( LIBSVM REQUIRED )
  include_directories( ${LIBSVM_INCLUDE_DIRS} )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RandomForest )
  find_package( RandomForest REQUIRED )
  include_directories( ${RandomForest_INCLUDE_DIR} )
endif( TubeTK_USE_RandomForest )

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
if( TubeTK_USE_RandomForest )
  list( APPEND tubeBaseIOTests_SRCS
    itktubePDFSegmenterRandomForestIOTest.cxx )
endif( TubeTK_USE_RandomForest )


add_executable( tubeBaseIOHeaderTest
  tubeBaseIOHeaderTest.cxx )
target_link_libraries( tubeBaseIOHeaderTest
  TubeTKIO ${ITK_LIBRARIES} ${LIBSVM_LIBRARIES} )
if( TubeTK_USE_LIBSVM )
  target_link_libraries( tubeBaseIOHeaderTest ${LIBSVM_LIBRARIES} )
endif( TubeTK_USE_LIBSVM )

set( tubeBaseIO_External_Libraries )
if( TubeTK_USE_GPU_ARRAYFIRE )
  set( tubeBaseIO_External_Libraries ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_GPU_ARRAYFIRE )
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

Midas3FunctionAddTest( NAME itktubePDFSegmenterParzenIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubePDFSegmenterParzenIOTest.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mha
    itktubePDFSegmenterParzenIOTest
      MIDAS{ES0015_Large.mha.md5}
      MIDAS{ES0015_Large.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask.mha.md5}
      ${TEMP}/itktubePDFSegmenterParzenIOTest.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest.mpd
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mha
      ${TEMP}/itktubePDFSegmenterParzenIOTest2.mpd )

if( TubeTK_USE_LIBSVM )
  Midas3FunctionAddTest( NAME itktubePDFSegmenterSVMIOTest
    COMMAND ${BASE_IO_TESTS}
      --compareNumberOfPixelsTolerance 20
      --compare ${TEMP}/itktubePDFSegmenterSVMIOTest.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.mha
      itktubePDFSegmenterSVMIOTest
        MIDAS{ES0015_Large.mha.md5}
        MIDAS{ES0015_Large.mha.md5}
        MIDAS{GDS0015_Large-TrainingMask.mha.md5}
        ${TEMP}/itktubePDFSegmenterSVMIOTest.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest.msvm
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.mha
        ${TEMP}/itktubePDFSegmenterSVMIOTest2.msvm )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_RandomForest )
  Midas3FunctionAddTest( NAME itktubePDFSegmenterRandomForestIOTest
    COMMAND ${BASE_IO_TESTS}
      --compareNumberOfPixelsTolerance 20
      --compare ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mha
      itktubePDFSegmenterRandomForestIOTest
        MIDAS{ES0015_Large.mha.md5}
        MIDAS{ES0015_Large.mha.md5}
        MIDAS{GDS0015_Large-TrainingMask.mha.md5}
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest.mrf
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mha
        ${TEMP}/itktubePDFSegmenterRandomForestIOTest2.mrf )
endif( TubeTK_USE_RandomForest )

Midas3FunctionAddTest( NAME itktubeRidgeSeedFilterIOTest
  COMMAND ${BASE_IO_TESTS}
    --compare ${TEMP}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterIOTest_Output2.mha
    --compareNumberOfPixelsTolerance 50
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
  COMMAND ${CompareTextFiles_EXE}
    -b MIDAS{TubeExtractorIOTest.mtp.md5}
    -t ${TEMP}/itktubeTubeExtractorIOTest.mtp )
set_property( TEST itktubeTubeExtractorIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeExtractorIOTest )

Midas3FunctionAddTest( NAME itktubeTubeExtractorIOTest-Compare2
  COMMAND ${CompareTextFiles_EXE}
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
  COMMAND ${CompareTextFiles_EXE}
    -b MIDAS{TubeXIOTest.tre.md5}
    -t ${TEMP}/itktubeTubeXIOTest.tre )
set_property( TEST itktubeTubeXIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeXIOTest )
