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
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( TEXTCOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/TextCompareCommand )

#######
# Test1
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test1
            ${PROJ_EXE}
               --useNumberOfLDA 3
               --objectId 0,255
               --zeroScales 1
               --firstScales 3
               --secondScales 5
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveFeatureImages ${TEMP}/${PROJECT_NAME}Test1.fi
               MIDAS{ES0015_Large_Subs.mha.md5}
               MIDAS{GDS0015_Large_Subs2.mha.md5}
               ${TEMP}/${PROJECT_NAME}Test1 )

#######
# Test1-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test1-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.01
               -t ${TEMP}/${PROJECT_NAME}Test1.lda00.mha
               -b MIDAS{${PROJECT_NAME}Test1.lda00.mha.md5} )
set_tests_properties( ${PROJECT_NAME}-Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test1 )

#######
# Test2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test2
            ${PROJ_EXE}
               --useNumberOfLDA 2
               --objectId 0,255
               --usePCA
               --zeroScales 1
               --firstScales 3
               --secondScales 5
               MIDAS{ES0015_Large_Subs.mha.md5}
               MIDAS{GDS0015_Large_Subs2.mha.md5}
               ${TEMP}/${PROJECT_NAME}Test2 )

#######
# Test2-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test2-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.01
               -t ${TEMP}/${PROJECT_NAME}Test2.lda01.mha
               -b MIDAS{${PROJECT_NAME}Test2.lda01.mha.md5} )
set_tests_properties( ${PROJECT_NAME}-Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test2 )

#######
# Test3
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test3
            ${PROJ_EXE}
               --useNumberOfLDA 4
               --objectId 0,255
               --loadLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test3.mnda
               MIDAS{ES0015_Large_Subs.mha.md5}
               MIDAS{GDS0015_Large_Subs2.mha.md5}
               ${TEMP}/${PROJECT_NAME}Test3 )
set_tests_properties( ${PROJECT_NAME}-Test3
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test1 )

#######
# Test3-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test3-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.02
               -t ${TEMP}/${PROJECT_NAME}Test3.lda00.mha
               -b MIDAS{${PROJECT_NAME}Test1.lda00.mha.md5} )
set_tests_properties( ${PROJECT_NAME}-Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test3 )

#######
# Test3-Compare2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test3-Compare2
            ${TEXTCOMPARE_EXE}
               -d 0.01
               -t ${TEMP}/${PROJECT_NAME}Test3.mnda
               -b MIDAS{${PROJECT_NAME}Test1.mnda.md5} )
set_tests_properties( ${PROJECT_NAME}-Test3-Compare2
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test3 )

#######
# Test4
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test4
            ${PROJ_EXE}
               --useNumberOfLDA 1
               --objectId 127,255
               --zeroScales 1
               --firstScales 3
               --secondScales 5
               --loadLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test4.mnda
               MIDAS{im0001_n05.crop.mha.md5}
               MIDAS{im0001.vk.crop.mha.md5}
               ${TEMP}/${PROJECT_NAME}Test4 )
set_tests_properties( ${PROJECT_NAME}-Test4
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test1 )

#######
# Test4-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test4-Compare
            ${TEXTCOMPARE_EXE}
               -d 0.01
               -t ${TEMP}/${PROJECT_NAME}Test4.mnda
               -b MIDAS{${PROJECT_NAME}Test1.mnda.md5} )
set_tests_properties( ${PROJECT_NAME}-Test4-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}-Test4 )
