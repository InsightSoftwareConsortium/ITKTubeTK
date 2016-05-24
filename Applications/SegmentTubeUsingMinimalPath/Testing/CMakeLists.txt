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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( DATA ${CMAKE_CURRENT_SOURCE_DIR}/Data )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

# Test1 - RegularStepGradientDescent
Midas3FunctionAddTest( NAME ${MODULE_NAME}-RegularStepGradientDescent
  COMMAND ${PROJ_EXE}
    MIDAS{Synthetic-04-Speed.mha.md5}
    ${TEMP}/${MODULE_NAME}-RegularStepGradientDescent.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Regular_Step_Gradient_Descent
    --terminationValue 2
    --numberOfIterations 1000
    --stepLengthFactor 1
    --stepLengthRelax 0.999 )

# Test1 - Compare - RegularStepGradientDescent
Midas3FunctionAddTest( NAME ${MODULE_NAME}-RegularStepGradientDescent-Compare
  COMMAND ${CompareTextFiles_EXE}
    -t ${TEMP}/${MODULE_NAME}-RegularStepGradientDescent.tre
    -b MIDAS{MPE-RegularStepGradientDescent.tre.md5}
    -d 0.01 )
set_property( TEST ${MODULE_NAME}-RegularStepGradientDescent-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-RegularStepGradientDescent)

# Test2 - GradientDescent
Midas3FunctionAddTest( NAME ${MODULE_NAME}-GradientDescent
  COMMAND ${PROJ_EXE}
    MIDAS{Synthetic-04-Speed.mha.md5}
    ${TEMP}/${MODULE_NAME}-GradientDescent.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Gradient_Descent
    --terminationValue 2
    --numberOfIterations 1000 )

# Test2 - Compare - GradientDescent
Midas3FunctionAddTest( NAME ${MODULE_NAME}-GradientDescent-Compare
  COMMAND ${CompareTextFiles_EXE}
    -t ${TEMP}/${MODULE_NAME}-GradientDescent.tre
    -b MIDAS{MPE-GradientDescent.tre.md5}
    -d 0.01 )
set_property( TEST ${MODULE_NAME}-GradientDescent-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-GradientDescent)

# Test3 - IterateNeighborhood
Midas3FunctionAddTest( NAME ${MODULE_NAME}-IterateNeighborhood
  COMMAND ${PROJ_EXE}
    MIDAS{Synthetic-04-Speed.mha.md5}
    ${TEMP}/${MODULE_NAME}-IterateNeighborhood.tre
    --startPoint 113.962,58.9872,115
    --targetPoint 66.3718,113.141,76
    --optimizer Iterate_Neighborhood
    --terminationValue 2
    --stepLengthFactor 1 )

# Test3 - Compare - IterateNeighborhood
Midas3FunctionAddTest( NAME ${MODULE_NAME}-IterateNeighborhood-Compare
  COMMAND ${CompareTextFiles_EXE}
    -t ${TEMP}/${MODULE_NAME}-IterateNeighborhood.tre
    -b MIDAS{MPE-IterateNeighborhood.tre.md5}
    -d 0.01 )
set_property( TEST ${MODULE_NAME}-IterateNeighborhood-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-IterateNeighborhood)
