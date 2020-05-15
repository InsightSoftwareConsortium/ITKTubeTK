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

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# expect image comparison perfection on 64-bit linux machines
if( UNIX AND NOT APPLE AND NOT CYGWIN AND CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( CompareImagesTolerance 0 )
# ... and give a bit of tolerance for 32-bit machines or windows/mac
else( UNIX AND NOT APPLE AND NOT CYGWIN AND CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( CompareImagesTolerance 0.000001 )
endif( UNIX AND NOT APPLE AND NOT CYGWIN AND CMAKE_SIZEOF_VOID_P EQUAL 8 )

# Test1
itk_add_test(
            NAME ${MODULE_NAME}-Test1
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1.mha )

# Test1-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Test1-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
               CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test1.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test1 )

# Test2
itk_add_test(
            NAME ${MODULE_NAME}-Test2
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/SyntheticTubePlusSurfaceImage.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test2.mha )

# Test2-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Test2-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
               CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test2.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test2.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test2 )

# Test3
itk_add_test(
            NAME ${MODULE_NAME}-Test3
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/SyntheticTubePlusSurfaceImage_oriented.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test3.mha )

# Test3-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Test3-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
               CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test3.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test3.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test3 )
