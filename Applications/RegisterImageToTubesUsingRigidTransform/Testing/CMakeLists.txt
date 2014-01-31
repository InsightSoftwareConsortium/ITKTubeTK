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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
  --parameterProgression ${TEMP}/${MODULE_NAME}Test1ParameterProgression.h5
  MIDAS{Branch.n020.mha.md5}
  MIDAS{tube.tre.md5}
  ${TEMP}/${MODULE_NAME}Test1.h5 )

set( METRIC_SAMPLER ComputeImageToTubeRigidMetricImage )
set( METRIC_SAMPLER_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${METRIC_SAMPLER}> )
Midas3FunctionAddTest( NAME ${METRIC_SAMPLER}-Test1
  COMMAND ${METRIC_SAMPLER_EXE}
    MIDAS{Branch.n020.mha.md5}
    MIDAS{tube.tre.md5}
    ${TEMP}/${METRIC_SAMPLER}Test1.mha )

set( METRIC_EXTRACTOR ExtractMetricImageSlice )
set( METRIC_EXTRACTOR_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${METRIC_EXTRACTOR}> )
Midas3FunctionAddTest( NAME ${METRIC_EXTRACTOR}-Test1
  COMMAND ${METRIC_EXTRACTOR_EXE}
    ${TEMP}/${METRIC_SAMPLER}Test1.mha
    ${TEMP}/${METRIC_EXTRACTOR}Test1.mha )
set_property( TEST ${METRIC_EXTRACTOR}-Test1
  APPEND PROPERTY DEPENDS ${METRIC_SAMPLER}-Test1 )

get_target_property( METRIC_SAMPLER_LOC ${METRIC_SAMPLER} LOCATION )
get_target_property( ANALYSIS_EXE ${MODULE_NAME} LOCATION )
configure_file( ${MODULE_NAME}Tuner.json.in
  ${TEMP}/${MODULE_NAME}Tuner.json )

# Needs ParameterSerializer
if( NOT TubeTK_USE_SLICER )
  if( TubeTK_USE_PYQTGRAPH )
    Midas3FunctionAddTest( NAME ${MODULE_NAME}TunerTest
      COMMAND ${PYTHON_TESTING_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../${MODULE_NAME}Tuner.py
        --non-interactive
        ${TEMP}/${MODULE_NAME}Tuner.json
      MIDAS_FETCH_ONLY{Branch.n020.mha.md5}
      MIDAS_FETCH_ONLY{tube.tre.md5} )
  endif( TubeTK_USE_PYQTGRAPH )

  Midas3FunctionAddTest( NAME ${METRIC_SAMPLER}-Test2
    COMMAND ${METRIC_SAMPLER_EXE}
        --parameterstorestore ${TEMP}/${MODULE_NAME}Tuner.json
    MIDAS{Branch.n020.mha.md5}
    MIDAS{tube.tre.md5}
    ${TEMP}/${METRIC_SAMPLER}Test2.mha )
endif( NOT TubeTK_USE_SLICER )
