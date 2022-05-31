#############################################################################
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
#       https://www.apache.org/licenses/LICENSE-2.0
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

# Test1
itk_add_test(
            NAME ${MODULE_NAME}-Test1
            COMMAND ${PROJ_EXE}
               -a
               --saveTransform ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.tfm
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Left.mha}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Right.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha )

itk_add_test(
            NAME ${MODULE_NAME}-Test1-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha
               -b DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test1 )

itk_add_test(
            NAME ${MODULE_NAME}-Test1-Tfm-Compare
            COMMAND ${TubeTK_CompareTextFiles_EXE}
              CompareTextFiles
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.tfm
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.tfm}
               -d 0.1 )
set_tests_properties( ${MODULE_NAME}-Test1-Tfm-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test1 )

# Test2
itk_add_test(
            NAME ${MODULE_NAME}-Test2
            COMMAND ${PROJ_EXE}
               -i 0
               -b 50,50
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Wo_offset.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha )

itk_add_test(
            NAME ${MODULE_NAME}-Test2-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test2 )

# Test3
itk_add_test(
            NAME ${MODULE_NAME}-Test3
            COMMAND ${PROJ_EXE}
               -i 0
               -k 10
               -b 50,50
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Wo_offset.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha )

itk_add_test(
            NAME ${MODULE_NAME}-Test3-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3.mha} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test3 )

# Test4
itk_add_test(
            NAME ${MODULE_NAME}-Test4
            COMMAND ${PROJ_EXE}
               -b 50,50
               -M
               -i 0
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Wo_offset.mha}
               ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha )

itk_add_test(
            NAME ${MODULE_NAME}-Test4-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-Test4 )
