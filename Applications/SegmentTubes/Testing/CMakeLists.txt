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

set( IMAGECOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:ImageCompare> )

# Test1
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
            COMMAND ${PROJ_EXE}
               -b 0
               -i 30,50,30
               -o ${TEMP}/${MODULE_NAME}Test1.mha
               MIDAS{Branch.n010.mha.md5}
               ${TEMP}/${MODULE_NAME}Test1.tre )

# Test1-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare
            COMMAND ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test1.mha
               -b MIDAS{${MODULE_NAME}Test1.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare
              APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )


# Test2
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2
            COMMAND ${PROJ_EXE}
               -b 0
               --seedP 32.0,50.0,32
               -o ${TEMP}/${MODULE_NAME}Test2.mha
               MIDAS{Branch.n010.mha.md5}
               ${TEMP}/${MODULE_NAME}Test2.tre )

# Test2-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-Compare
            COMMAND ${IMAGECOMPARE_EXE}
               -r 1
               -t ${TEMP}/${MODULE_NAME}Test2.mha
               -b MIDAS{${MODULE_NAME}Test1.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test2-Compare
              APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

# Test 3
# Test duplicate extraction - using prior extraction mask
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3
            COMMAND ${PROJ_EXE}
               -b 0
               --seedP 32.0,50.0,32
               -e ${TEMP}/${MODULE_NAME}Test2.mha
               -o ${TEMP}/${MODULE_NAME}Test3.mha
               MIDAS{Branch.n010.mha.md5}
               ${TEMP}/${MODULE_NAME}Test3.tre )
set_tests_properties( ${MODULE_NAME}-Test3 PROPERTIES WILL_FAIL true )
set_property( TEST ${MODULE_NAME}-Test3
              APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

# Test 4
# Test duplicate extraction - using prior extraction tubes
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4
            COMMAND ${PROJ_EXE}
               -b 0
               --seedP 32.0,50.0,32
               -v ${TEMP}/${MODULE_NAME}Test2.tre
               -o ${TEMP}/${MODULE_NAME}Test4.mha
               MIDAS{Branch.n010.mha.md5}
               ${TEMP}/${MODULE_NAME}Test4.tre )
set_tests_properties( ${MODULE_NAME}-Test4 PROPERTIES WILL_FAIL true )
set_property( TEST ${MODULE_NAME}-Test4
              APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )
