############################################################################
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
############################################################################
include_regular_expression( "^.*$" )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1 - NotCropped
itk_add_test(
  NAME ${MODULE_NAME}-Test1_NotCropped
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/CropTubes-TestInput.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1_NotCropped.tre
    --boxCorner 17.5,14,7.25
    --boxSize 5,5,5 )

# Test1 - Compare - NotCropped
itk_add_test(
  NAME ${MODULE_NAME}-Test1_NotCropped-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1_NotCropped.tre
    -b DATA{${TubeTK_DATA_ROOT}/CropTubes-Test1_NotCropped.tre}
    -d 0.01 )
set_tests_properties( ${MODULE_NAME}-Test1_NotCropped-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Test1_NotCropped )

# Test2 - Cropped
itk_add_test(
  NAME ${MODULE_NAME}-Test2_Cropped
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/CropTubes-TestInput.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test2_Cropped.tre
    --boxCorner 17.5,14,7.25
    --boxSize 5,5,5
    --CropTubes )

# Test2 - Compare - Cropped
itk_add_test(
  NAME ${MODULE_NAME}-Test2_Cropped-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test2_Cropped.tre
    -b DATA{${TubeTK_DATA_ROOT}/CropTubes-Test2_Cropped.tre}
    -d 0.01 )

set_tests_properties( ${MODULE_NAME}-Test2_Cropped-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Test2_Cropped )
