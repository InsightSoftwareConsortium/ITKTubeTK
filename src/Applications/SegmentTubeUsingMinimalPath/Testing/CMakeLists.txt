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

# Test1 - RegularStepGradientDescent
itk_add_test(
  NAME ${MODULE_NAME}-RegularStepGradientDescent
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Synthetic-04-Speed.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-RegularStepGradientDescent.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Regular_Step_Gradient_Descent
    --terminationValue 2
    --numberOfIterations 1000
    --stepLengthFactor 1
    --stepLengthRelax 0.999 )

# Test1 - Compare - RegularStepGradientDescent
itk_add_test(
  NAME ${MODULE_NAME}-RegularStepGradientDescent-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-RegularStepGradientDescent.tre
    -b DATA{${TubeTK_DATA_ROOT}/SegmentTubeUsingMinimalPath-RegularStepGradientDescent.tre}
    -d 0.01 )
set_tests_properties( ${MODULE_NAME}-RegularStepGradientDescent-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-RegularStepGradientDescent )

# Test2 - GradientDescent
itk_add_test(
  NAME ${MODULE_NAME}-GradientDescent
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Synthetic-04-Speed.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-GradientDescent.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Gradient_Descent
    --terminationValue 2
    --numberOfIterations 1000 )

# Test2 - Compare - GradientDescent
itk_add_test(
  NAME ${MODULE_NAME}-GradientDescent-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-GradientDescent.tre
    -b DATA{${TubeTK_DATA_ROOT}/SegmentTubeUsingMinimalPath-GradientDescent.tre}
    -d 0.01 )
set_tests_properties( ${MODULE_NAME}-GradientDescent-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-GradientDescent )

# Test3 - IterateNeighborhood
itk_add_test(
  NAME ${MODULE_NAME}-IterateNeighborhood
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Synthetic-04-Speed.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-IterateNeighborhood.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Iterate_Neighborhood
    --terminationValue 2
    --stepLengthFactor 1 )

# Test3 - Compare - IterateNeighborhood
itk_add_test(
  NAME ${MODULE_NAME}-IterateNeighborhood-Compare
  COMMAND ${TubeTK_CompareTextFiles_EXE}
    CompareTextFiles
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-IterateNeighborhood.tre
    -b DATA{${TubeTK_DATA_ROOT}/SegmentTubeUsingMinimalPath-IterateNeighborhood.tre}
    -d 0.01 )
set_tests_properties( ${MODULE_NAME}-IterateNeighborhood-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-IterateNeighborhood )
