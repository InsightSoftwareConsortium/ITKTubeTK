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

include_regular_expression( "^.*$" )

set( PROJ_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1
itk_add_test(
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    -o ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.txt
    DATA{${TubeTK_DATA_ROOT}/GDS0015_Large_Modified.mha}
    DATA{${TubeTK_DATA_ROOT}/GDS0015_Large.mha} )

# Test1-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.txt
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.txt} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )
