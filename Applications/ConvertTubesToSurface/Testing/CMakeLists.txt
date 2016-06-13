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

set( CompareTextFiles_EXE
     ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

set( CompareImages_EXE
     ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

# Test1
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    MIDAS{Branch-truth-new.tre.md5}
    ${TEMP}/${MODULE_NAME}-Test1.vtk
    )

# Test2
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    MIDAS{Branch-truth-new.tre.md5}
    ${TEMP}/${MODULE_NAME}-Test2.stl
    )

# Test3
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    MIDAS{Branch-truth-new.tre.md5}
    ${TEMP}/${MODULE_NAME}-Test3.stt
    )
set_tests_properties( ${MODULE_NAME}-Test3
  PROPERTIES WILL_FAIL true )

# Test4
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    MIDAS{ConvertTubesToSurfaceTestFile.tre.md5}
    ${TEMP}/${MODULE_NAME}-Test4_Surface.vtk
    --centerline ${TEMP}/${MODULE_NAME}-Test4_Centerline.vtp)

# Test5 - Compare - Surface
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4-Compare-Surface
  COMMAND ${CompareTextFiles_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test4_Surface.vtk
    -b MIDAS{ConvertTubesToSurfaceTestFile_Surface.vtk.md5}
    -d 0.01 )

set_property( TEST ${MODULE_NAME}-Test4-Compare-Surface
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )
