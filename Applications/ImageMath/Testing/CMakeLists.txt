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
  ${TubeTK_SOURCE_DIR}/Applications/${PROJECT_NAME}/Testing/Baselines )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE 
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand )

#######
# Test1
#######
add_test( ${PROJECT_NAME}Test1
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -i -4000 60000 -4 60
               -w ${TEMP}/${PROJECT_NAME}Test1.mha )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test1.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test1.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test2
#######
add_test( ${PROJECT_NAME}Test2
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -i -4000 60000 -4 60
               -W 2 ${TEMP}/${PROJECT_NAME}Test2.mha )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test2.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test2.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test2 )

#######
# Test3
#######
add_test( ${PROJECT_NAME}Test3
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -n -4000 60000 -100 10
               -W 2 ${TEMP}/${PROJECT_NAME}Test3.mha )

#######
# Test4
#######
add_test( ${PROJECT_NAME}Test4
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -N -4000 60000 0 100 )

#######
# Test5
#######
add_test( ${PROJECT_NAME}Test5
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -f ${INPUT_DATA}/ES0015_Large.mha ${BASELINE}/${PROJECT_NAME}Test3.mha 
               -W 2 ${TEMP}/${PROJECT_NAME}Test5.mha )
set_tests_properties( ${PROJECT_NAME}Test5
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )

#######
# Test5-Compare
#######
add_test( ${PROJECT_NAME}Test5-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test5.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test5.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test5-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test5 )

#######
# Test6
#######
add_test( ${PROJECT_NAME}Test6
            ${PROJ_EXE}
               -a 0.5 0.5 ${INPUT_DATA}/ES0015_Large.mha  ${BASELINE}/${PROJECT_NAME}Test1.mha 
               -w ${TEMP}/${PROJECT_NAME}Test6.mha )
set_tests_properties( ${PROJECT_NAME}Test6
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test6-Compare
#######
add_test( ${PROJECT_NAME}Test6-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test6.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test6.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test6-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test6 )

#######
# Test7
#######
add_test( ${PROJECT_NAME}Test7
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -t -1 1 0 100
               -w ${TEMP}/${PROJECT_NAME}Test7.mha )

#######
# Test7-Compare
#######
add_test( ${PROJECT_NAME}Test7-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test7.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test7.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test7-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test7 )

#######
# Test8
#######
add_test( ${PROJECT_NAME}Test8
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -t -1 1 ${TEMP}/${PROJECT_NAME}Test2.mha 100
               -w ${TEMP}/${PROJECT_NAME}Test8.mha )

#######
# Test8-Compare
#######
add_test( ${PROJECT_NAME}Test8-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test8.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test8.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test8-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test8 )

#######
# Test9
#######
add_test( ${PROJECT_NAME}Test9
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -m -1 1 ${TEMP}/${PROJECT_NAME}Test2.mha 100
               -w ${TEMP}/${PROJECT_NAME}Test9.mha )

#######
# Test9-Compare
#######
add_test( ${PROJECT_NAME}Test9-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test9.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test9.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test9-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test9 )

#######
# Test10
#######
add_test( ${PROJECT_NAME}Test10
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -M 1 5 1 1
               -w ${TEMP}/${PROJECT_NAME}Test10.mha )

#######
# Test10-Compare
#######
add_test( ${PROJECT_NAME}Test10-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test10.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test10.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test10-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test10 )

#######
# Test11
#######
add_test( ${PROJECT_NAME}Test11
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -b 10
               -w ${TEMP}/${PROJECT_NAME}Test11.mha )

#######
# Test11-Compare
#######
add_test( ${PROJECT_NAME}Test11-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test11.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test11.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test11-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test11 )

#######
# Test12
#######
add_test( ${PROJECT_NAME}Test12
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -B 5 2 0
               -w ${TEMP}/${PROJECT_NAME}Test12.mha )

#######
# Test12-Compare
#######
add_test( ${PROJECT_NAME}Test12-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test12.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test12.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test12-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test12 )

#######
# Test13
#######
add_test( ${PROJECT_NAME}Test13
            ${PROJ_EXE}
               ${INPUT_DATA}/im0001.crop.mha
               -z 1 2 3
               -w ${TEMP}/${PROJECT_NAME}Test13.mha )

#######
# Test13-Compare
#######
add_test( ${PROJECT_NAME}Test13-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test13.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test13.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test13-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test13 )

#######
# Test14
#######
add_test( ${PROJECT_NAME}Test14
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -l 10 ${TEMP}/${PROJECT_NAME}Test14.txt )

#######
# Test14-Compare
#######
add_test( ${PROJECT_NAME}Test14-Compare
            cmake -E compare_files 
            ${TEMP}/${PROJECT_NAME}Test14.txt
            ${BASELINE}/${PROJECT_NAME}Test14.txt )
set_tests_properties( ${PROJECT_NAME}Test14-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test14 )

#######
# Test15
#######
add_test( ${PROJECT_NAME}Test15
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -L 10 1 10000 ${TEMP}/${PROJECT_NAME}Test15.txt )

#######
# Test15-Compare
#######
add_test( ${PROJECT_NAME}Test15-Compare
            cmake -E compare_files 
            ${TEMP}/${PROJECT_NAME}Test15.txt
            ${BASELINE}/${PROJECT_NAME}Test15.txt )
set_tests_properties( ${PROJECT_NAME}Test15-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test15 )

#######
# Test16
#######
add_test( ${PROJECT_NAME}Test16
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -c 10 1000
               -w ${TEMP}/${PROJECT_NAME}Test16.mha )

#######
# Test16-Compare
#######
add_test( ${PROJECT_NAME}Test16-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test16.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test16.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test16-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test16 )

#######
# Test17
#######
add_test( ${PROJECT_NAME}Test17
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -C 10 1000 ${BASELINE}/${PROJECT_NAME}Test7.mha
               -w ${TEMP}/${PROJECT_NAME}Test17.mha )

#######
# Test17-Compare
#######
add_test( ${PROJECT_NAME}Test17-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test17.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test17.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test17-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test17 )

#######
# Test18
#######
add_test( ${PROJECT_NAME}Test18
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -r 2 
               -w ${TEMP}/${PROJECT_NAME}Test18.mha )

#######
# Test18-Compare
#######
add_test( ${PROJECT_NAME}Test18-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test18.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test18.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test18-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test18 )

#######
# Test19
#######
add_test( ${PROJECT_NAME}Test19
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -R ${BASELINE}/${PROJECT_NAME}Test18.mha
               -w ${TEMP}/${PROJECT_NAME}Test19.mha )

#######
# Test19-Compare
#######
add_test( ${PROJECT_NAME}Test19-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test19.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test19.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test19-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test19 )

#######
# Test20
#######
add_test( ${PROJECT_NAME}Test20
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -s 0 -1000 1000 3 100 100 0
               -w ${TEMP}/${PROJECT_NAME}Test20.mha )

#######
# Test20-Compare
#######
add_test( ${PROJECT_NAME}Test20-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test20.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test20.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test20-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test20 )

#######
# Test21
#######
add_test( ${PROJECT_NAME}Test21
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -O -1000 -1000 -2000
               -w ${TEMP}/${PROJECT_NAME}Test21.mha )

#######
# Test21-Compare
#######
add_test( ${PROJECT_NAME}Test21-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test21.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test21.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test21-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test21 )

#######
# Test22
#######
add_test( ${PROJECT_NAME}Test22
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -S 7278492372
               -w ${TEMP}/${PROJECT_NAME}Test22.mha )

#######
# Test22-Compare
#######
add_test( ${PROJECT_NAME}Test22-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test22.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test22.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test22-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test22 )


#######
# Test23
#######
add_test( ${PROJECT_NAME}Test23
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -P 0 ${INPUT_DATA}/ES0015_Large.mha
               -w ${TEMP}/${PROJECT_NAME}Test23.mha )

#######
# Test23-Compare
#######
add_test( ${PROJECT_NAME}Test23-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test23.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test23.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test23-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test23 )

#######
# Test24
#######
add_test( ${PROJECT_NAME}Test24
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -u ${INPUT_DATA}/ES0015_Large.mha
               -w ${TEMP}/${PROJECT_NAME}Test24.mha )

#######
# Test24-Compare
#######
add_test( ${PROJECT_NAME}Test24-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test24.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test24.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test24-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test24 )

#######
# Test25
#######
add_test( ${PROJECT_NAME}Test25
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large.mha
               -p 0 ${INPUT_DATA}/ES0015_Large.mha
               -w ${TEMP}/${PROJECT_NAME}Test25.mha )

#######
# Test25-Compare
#######
add_test( ${PROJECT_NAME}Test25-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test25.mha 
               -b ${BASELINE}/${PROJECT_NAME}Test25.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test25-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test25 )
