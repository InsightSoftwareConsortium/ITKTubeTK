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

# Test1
itk_add_test(
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Den.mha
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Rad.mha
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Tan.mha
    --inputTemplateImage DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
    --useSquaredDistance )

# Test1-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare-Den
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -i 0.001
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Den.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Den.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Den PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

itk_add_test( NAME ${MODULE_NAME}-Test1-Compare-Rad
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -i 0.001
    -n 50
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Rad.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Rad.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Rad PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

itk_add_test( NAME ${MODULE_NAME}-Test1-Compare-Tan
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -i 0.00001
    -n 50
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1-Tan.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Tan.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Tan PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )
