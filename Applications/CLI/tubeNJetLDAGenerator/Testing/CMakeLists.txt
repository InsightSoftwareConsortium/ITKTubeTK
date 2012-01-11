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

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )

set( BASELINE
  ${TubeTK_SOURCE_DIR}/Applications/CLI/${PROJECT_NAME}/Testing/Baselines )

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
add_test( ${PROJECT_NAME}Test1
            ${PROJ_EXE}
               --useNumberOfLDA 3
               --objectId 0,255
               --zeroScales 1
               --firstScales 1,3
               --secondScales 5
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveFeatureImages ${TEMP}/${PROJECT_NAME}Test1.fi
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               ${INPUT_DATA}/GDS0015_Large_Subs2.mha
               ${TEMP}/${PROJECT_NAME}Test1 )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}Test1-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.01
               -t ${TEMP}/${PROJECT_NAME}Test1.lda00.mha
               -b ${BASELINE}/${PROJECT_NAME}Test1.lda00.mha )
set_tests_properties( ${PROJECT_NAME}Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test2
#######
add_test( ${PROJECT_NAME}Test2
            ${PROJ_EXE}
               --useNumberOfLDA 2
               --objectId 0,255
               --usePCA
               --zeroScales 1
               --firstScales 1,3
               --secondScales 5
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               ${INPUT_DATA}/GDS0015_Large_Subs2.mha
               ${TEMP}/${PROJECT_NAME}Test2 )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}Test2-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.01
               -t ${TEMP}/${PROJECT_NAME}Test2.lda01.mha
               -b ${BASELINE}/${PROJECT_NAME}Test2.lda01.mha )
set_tests_properties( ${PROJECT_NAME}Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test2 )

#######
# Test3
#######
add_test( ${PROJECT_NAME}Test3
            ${PROJ_EXE}
               --useNumberOfLDA 4
               --objectId 0,255
               --loadLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test3.mnda
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               ${INPUT_DATA}/GDS0015_Large_Subs2.mha
               ${TEMP}/${PROJECT_NAME}Test3 )
set_tests_properties( ${PROJECT_NAME}Test3
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}Test3-Compare
            ${IMAGECOMPARE_EXE}
               -i 0.02
               -t ${TEMP}/${PROJECT_NAME}Test3.lda00.mha
               -b ${BASELINE}/${PROJECT_NAME}Test1.lda00.mha )
set_tests_properties( ${PROJECT_NAME}Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )

#######
# Test3-Compare2
#######
add_test( ${PROJECT_NAME}Test3-Compare2
            ${TEXTCOMPARE_EXE}
               -d 0.01
               -t ${TEMP}/${PROJECT_NAME}Test3.mnda
               -b ${BASELINE}/${PROJECT_NAME}Test1.mnda )
set_tests_properties( ${PROJECT_NAME}Test3-Compare2
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )

#######
# Test4
#######
add_test( ${PROJECT_NAME}Test4
            ${PROJ_EXE}
               --useNumberOfLDA 1
               --objectId 127,255
               --zeroScales 1
               --firstScales 1,3
               --secondScales 5
               --loadLDAInfo ${TEMP}/${PROJECT_NAME}Test1.mnda
               --saveLDAInfo ${TEMP}/${PROJECT_NAME}Test4.mnda
               ${INPUT_DATA}/im0001_n05.crop.mha
               ${INPUT_DATA}/im0001.vk.crop.mha
               ${TEMP}/${PROJECT_NAME}Test4 )
set_tests_properties( ${PROJECT_NAME}Test4
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test4-Compare
#######
add_test( ${PROJECT_NAME}Test4-Compare
            ${TEXTCOMPARE_EXE}
               -d 0.01
               -t ${TEMP}/${PROJECT_NAME}Test4.mnda
               -b ${BASELINE}/${PROJECT_NAME}Test1.mnda )
set_tests_properties( ${PROJECT_NAME}Test4-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test4 )
