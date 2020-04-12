##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
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
set( CompareTextFiles_EXE $<TARGET_FILE:tubeCompareTextFiles> )

if( PYTHON_EXECUTABLE )
  itk_add_test( PrintSelf-TubeTKIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/IO )
endif( PYTHON_EXECUTABLE )

CreateTestDriver( tubeIOHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeIOHeaderTest.cxx" )
if( TubeTK_USE_ARRAYFIRE )
  target_include_directories( tubeIOHeaderTestDriver PUBLIC
    ${ArrayFire_INCLUDE_DIRS} )
  target_link_libraries( tubeIOHeaderTestDriver LINK_PUBLIC
    ${ArrayFire_LIBRARIES} )
endif()

itk_add_test( NAME tubeIOHeaderTest
  COMMAND tubeIOHeaderTestDriver
    tubeIOHeaderTest )

set( tubeIOTests_SRCS
  tubeIOPrintTest.cxx
  itktubePDFSegmenterParzenIOTest.cxx
  itktubeTubeExtractorIOTest.cxx
  itktubeRidgeSeedFilterIOTest.cxx
  itktubeTubeXIOTest.cxx )

CreateTestDriver( tubeIO
  "${TubeTK-Test_LIBRARIES}"
  "${tubeIOTests_SRCS}" )
if( TubeTK_USE_ARRAYFIRE )
  target_include_directories( tubeIOTestDriver PUBLIC
    ${ArrayFire_INCLUDE_DIRS} )
  target_link_libraries( tubeIOTestDriver LINK_PUBLIC
    ${ArrayFire_LIBRARIES} )
endif()

itk_add_test( NAME tubeIOPrintTest
  COMMAND tubeIOTestDriver
    tubeIOPrintTest )

itk_add_test( 
  NAME itktubePDFSegmenterParzenIOTest
  COMMAND tubeIOTestDriver
    --compare ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest2.mha
    itktubePDFSegmenterParzenIOTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest.mpd
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest2.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenIOTest2.mpd )

itk_add_test( 
  NAME itktubeRidgeSeedFilterIOTest
  COMMAND tubeIOTestDriver
    --compare ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterIOTest_Output2.mha
    --compareNumberOfPixelsTolerance 50
    itktubeRidgeSeedFilterIOTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterIOTest_Output.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterIOTest.mrs
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterIOTest_Output2.mha )

itk_add_test( 
  NAME itktubeTubeExtractorIOTest
  COMMAND tubeIOTestDriver
    itktubeTubeExtractorIOTest
      DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeExtractorIOTest.mtp
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeExtractorIOTest2.mtp )

itk_add_test( 
  NAME itktubeTubeExtractorIOTest-Compare
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
    -t ${ITK_TEST_OUTPUT_DIR}/itktubeTubeExtractorIOTest.mtp )
set_property( TEST itktubeTubeExtractorIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeExtractorIOTest )

itk_add_test( 
  NAME itktubeTubeExtractorIOTest-Compare2
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeExtractorIOTest.mtp}
    -t ${ITK_TEST_OUTPUT_DIR}/itktubeTubeExtractorIOTest2.mtp )
set_property( TEST itktubeTubeExtractorIOTest-Compare2 APPEND PROPERTY DEPENDS
  itktubeTubeExtractorIOTest )

itk_add_test( 
  NAME itktubeTubeXIOTest
  COMMAND tubeIOTestDriver
    itktubeTubeXIOTest
      DATA{${TubeTK_DATA_ROOT}/TubeXIOTest.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeXIOTest.tre
      512 512 393
      )

itk_add_test( 
  NAME itktubeTubeXIOTest-Compare
  COMMAND ${CompareTextFiles_EXE}
    -b DATA{${TubeTK_DATA_ROOT}/TubeXIOTest.tre}
    -t ${ITK_TEST_OUTPUT_DIR}/itktubeTubeXIOTest.tre )
set_property( TEST itktubeTubeXIOTest-Compare APPEND PROPERTY DEPENDS
  itktubeTubeXIOTest )
