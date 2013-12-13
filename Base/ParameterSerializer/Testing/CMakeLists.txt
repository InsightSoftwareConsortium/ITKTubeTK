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

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseParameterSerializer ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer )
endif( PYTHON_EXECUTABLE )

set( BASE_PARAMETERSERIALIZER_TESTS
  $<TARGET_FILE:tubeBaseParameterSerializerTests> )

set( BASE_PARAMETERSERIALIZER_HEADER_TEST
  $<TARGET_FILE:tubeBaseParameterSerializerHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseParameterSerializerTests_SRCS
  Filtering/itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest.cxx
  USTK/itkAcousticImpulseResponseImageFilterSerializerTest.cxx
  USTK/itkUltrasoundProbeGeometryCalculatorSerializerTest.cxx
  Registration/itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/USTK
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Registration
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Filtering
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/USTK
  ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Registration )

SEMMacroBuildCLI(
  NAME tubeBaseParameterSerializerHeaderTest
  ADDITIONAL_SRCS
    Filtering/itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest.cxx
    USTK/itkAcousticImpulseResponseImageFilterSerializerTest.cxx
    USTK/itkUltrasoundProbeGeometryCalculatorSerializerTest.cxx
    Registration/itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest.cxx
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/USTK
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Registration
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Filtering
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/USTK
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Registration
  INSTALL_RUNTIME_DESTINATION lib
  INSTALL_LIBRARY_DESTINATION lib
  INSTALL_ARCHIVE_DESTINATION bin )

SEMMacroBuildCLI(
  NAME tubeBaseParameterSerializerTests
  ADDITIONAL_SRCS
    Filtering/itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest.cxx
    USTK/itkAcousticImpulseResponseImageFilterSerializerTest.cxx
    USTK/itkUltrasoundProbeGeometryCalculatorSerializerTest.cxx
    Registration/itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest.cxx
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/USTK
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Registration
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Filtering
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/USTK
    ${TubeTK_SOURCE_DIR}/Base/ParameterSerializer/Registration
  INSTALL_RUNTIME_DESTINATION lib
  INSTALL_LIBRARY_DESTINATION lib
  INSTALL_ARCHIVE_DESTINATION bin )

add_test( NAME tubeBaseParameterSerializerHeaderTest
  COMMAND ${BASE_PARAMETERSERIALIZER_HEADER_TEST} )

add_test( NAME itkAcousticImpulseResponseImageFilterSerializerTest
  COMMAND ${BASE_PARAMETERSERIALIZER_TESTS}
  itkAcousticImpulseResponseImageFilterSerializerTest
  ${TEMP}/itkAcousticImpulseResponseImageFilterSerializerTest.json )

Midas3FunctionAddTest( NAME itkAcousticImpulseResponseImageFilterSerializerTestBaseline
  COMMAND ${TubeTK_BINARY_DIR}/bin/JsonCompareCommand
  -t ${TEMP}/itkAcousticImpulseResponseImageFilterSerializerTest.json
  -b MIDAS{itkAcousticImpulseResponseImageFilterSerializerTest.json.md5} )

add_test( NAME itkUltrasoundProbeGeometryCalculatorSerializerTest
  COMMAND ${BASE_PARAMETERSERIALIZER_TESTS}
  itkUltrasoundProbeGeometryCalculatorSerializerTest
  ${TEMP}/itkUltrasoundProbeGeometryCalculatorSerializerTest.json )

add_test( NAME itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest
  COMMAND ${BASE_PARAMETERSERIALIZER_TESTS}
  itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest
  ${TEMP}/itktubeSubSampleTubeTreeSpatialObjectFilterSerializerTest.json )

add_test( NAME itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest
  COMMAND ${BASE_PARAMETERSERIALIZER_TESTS}
  itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest
  ${TEMP}/itktubeTubeAngleOfIncidenceWeightFunctionSerializerTest.json )
