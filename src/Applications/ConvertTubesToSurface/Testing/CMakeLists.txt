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
    DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1.vtp
    )

# Test2
itk_add_test(
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test2.stl
    )

# Test3
itk_add_test(
  NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test3.stt
    )

# Test4
itk_add_test(
  NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Input.tre}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test4_Surface.vtp
    --centerline ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test4_Centerline.vtp )

# Test5 - Compare - Surface
#itk_add_test(
  #NAME ${MODULE_NAME}-Test4-Compare-Surface
  #COMMAND ${TubeTK_CompareTextFiles_EXE}
    #CompareTextFiles
    #-t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test4_Surface.vtk
    #-b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test4_Surface.vtk}
    #-d 0.01 )

#set_tests_properties( ${MODULE_NAME}-Test4-Compare-Surface PROPERTIES DEPENDS
  #${MODULE_NAME}-Test4 )
