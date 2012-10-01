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

Include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)


set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
 ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${MODULE_NAME} )

set( IMAGECOMPARE_EXE
 ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )


#######
# Test1
#######
MIDAS_ADD_Test( ${MODULE_NAME}-Test1
               ${PROJ_EXE}
               MIDAS{ES0015_Large.mha.md5}
               ${TEMP}/${MODULE_NAME}Test1.mha
               -m 10,10 -M 100,100 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test1-Compare
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test1.mha
               -b MIDAS{${MODULE_NAME}Test1.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test1-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test2
#######
MIDAS_ADD_Test( ${MODULE_NAME}-Test2
               ${PROJ_EXE}
               MIDAS{ES0015_Large.mha.md5}
               ${TEMP}/${MODULE_NAME}Test2.mha
               -m 10,10 -s 91,91 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test2-Compare
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test2.mha
               -b MIDAS{${MODULE_NAME}Test1.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test2-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

#######
# Test3
#######
MIDAS_ADD_Test( ${MODULE_NAME}-Test3
               ${PROJ_EXE}
               MIDAS{ES0015_Large.mha.md5}
               ${TEMP}/${MODULE_NAME}Test3.mha
               -m 15,15 -s 81,81 -b 5,5 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test3-Compare
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test3.mha
               -b MIDAS{${MODULE_NAME}Test1.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test3-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )

#######
# Test4
#######
MIDAS_ADD_Test( ${MODULE_NAME}-Test4
               ${PROJ_EXE}
               MIDAS{ES0015_Large.mha.md5}
               ${TEMP}/${MODULE_NAME}Test4
               -S 2,2 -b 5,5 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test4-Compare00
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test4_00.mha
               -b MIDAS{${MODULE_NAME}Test4_00.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test4-Compare00
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test4-Compare01
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test4_01.mha
               -b MIDAS{${MODULE_NAME}Test4_01.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test4-Compare01
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test4-Compare10
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test4_10.mha
               -b MIDAS{${MODULE_NAME}Test4_10.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test4-Compare10
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )

MIDAS_ADD_Test( ${MODULE_NAME}-Test4-Compare11
               ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}Test4_11.mha
               -b MIDAS{${MODULE_NAME}Test4_11.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test4-Compare11
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )
