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

#
# Test with 2D Cylinder
#

# Run module and generate result on test case
itk_add_test(
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Cylinder_2D.mha}
    DATA{${TubeTK_DATA_ROOT}/Cylinder_2D.tre} 1
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.flyThrough.mha
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.flyThrough.tubeMask.mha )

# compare outputs with groundtruth
itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare-flyThrough
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.flyThrough.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.flyThrough.mha}
    -i 0.001 )

set_tests_properties( ${MODULE_NAME}-Test1-Compare-flyThrough PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare-flyThrough-TubeMask
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.flyThrough.tubeMask.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.flyThrough.tubeMask.mha}
    -i 0.001 )

set_tests_properties( ${MODULE_NAME}-Test1-Compare-flyThrough-TubeMask PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

#
# Test with 3D Cylinder
#

# Run module and generate result on test case
itk_add_test(
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Cylinder_3D.mha}
    DATA{${TubeTK_DATA_ROOT}/Cylinder_3D.tre} 1
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.flyThrough.mha
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.flyThrough.tubeMask.mha )

# compare outputs with groundtruth
itk_add_test(
  NAME ${MODULE_NAME}-Test2-Compare-flyThrough
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.flyThrough.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.flyThrough.mha}
    -i 0.001 )

set_tests_properties( ${MODULE_NAME}-Test2-Compare-flyThrough PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

itk_add_test(
  NAME ${MODULE_NAME}-Test2-Compare-flyThrough-TubeMask
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.flyThrough.tubeMask.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.flyThrough.tubeMask.mha}
    -i 0.001 )

set_tests_properties( ${MODULE_NAME}-Test2-Compare-flyThrough-TubeMask PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )
