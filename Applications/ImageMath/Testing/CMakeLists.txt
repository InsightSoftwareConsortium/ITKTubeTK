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

set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)


set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )

set( BASELINE
  ${TubeTK_SOURCE_DIR}/Applications/${PROJECT_NAME}/Testing/Baselines )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_BINARY_DIR}/bin/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( TEXTCOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/TextCompareCommand )

#######
# Test1
#######
add_test( ${PROJECT_NAME}Test1
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -i -1 1 -100 100
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -b 5
               -w ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -a 0 1 ${INPUT_DATA}/ES0015_Large_Subs.mha
               -I ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -w ${TEMP}/${PROJECT_NAME}Test2.mha )

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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -S 1234
               -n -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test3.mha )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test3.mha
               -b ${BASELINE}/${PROJECT_NAME}Test3.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )


#######
# Test4
#######
add_test( ${PROJECT_NAME}Test4
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -S 1234
               -N -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test4.mha )

#######
# Test4-Compare
#######
add_test( ${PROJECT_NAME}Test4-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test4.mha
               -b ${BASELINE}/${PROJECT_NAME}Test4.mha
               -i 0.01 )
set_tests_properties( ${PROJECT_NAME}Test4-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test4 )


#######
# Test5
#######
add_test( ${PROJECT_NAME}Test5
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -f 10
                 ${BASELINE}/${PROJECT_NAME}Test2.mha
               -w ${TEMP}/${PROJECT_NAME}Test5.mha )

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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -a 0.95 0.05 ${BASELINE}/${PROJECT_NAME}Test5.mha
               -w ${TEMP}/${PROJECT_NAME}Test6.mha )

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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -t -1 -0.33 -1 1
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -m -1 -0.33 ${INPUT_DATA}/ES0015_Large_Subs.mha 0
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
               ${INPUT_DATA}/GDS0015_1.mha
               -M 0 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test9.mha )

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
               ${INPUT_DATA}/GDS0015_1.mha
               -M 1 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test10.mha )

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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/im0001.crop2.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -l 10 ${TEMP}/${PROJECT_NAME}Test14.txt )

#######
# Test14-Compare
#######
add_test( ${PROJECT_NAME}Test14-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test14.txt
            -b ${BASELINE}/${PROJECT_NAME}Test14.txt )
set_tests_properties( ${PROJECT_NAME}Test14-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test14 )

#######
# Test15
#######
add_test( ${PROJECT_NAME}Test15
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -L 100 -1 0.01 ${TEMP}/${PROJECT_NAME}Test15.txt )

#######
# Test15-Compare
#######
add_test( ${PROJECT_NAME}Test15-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test15.txt
            -b ${BASELINE}/${PROJECT_NAME}Test15.txt )
set_tests_properties( ${PROJECT_NAME}Test15-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test15 )

#######
# Test16
#######
add_test( ${PROJECT_NAME}Test16
            ${PROJ_EXE}
               ${INPUT_DATA}/CroppedWholeLungCTScan.mhd
               -c 100 10
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -C 100 3 ${BASELINE}/${PROJECT_NAME}Test2.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -s 0 -1 -0.33 1 70 50 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test20.mha )

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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -O -10 -10 0
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
# Test23
#######
add_test( ${PROJECT_NAME}Test23
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -P 0 ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -u ${INPUT_DATA}/ES0015_Large_Subs.mha
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
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -p 0 ${INPUT_DATA}/ES0015_Large_Subs.mha
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

#######
# Test26
#######
add_test( ${PROJECT_NAME}Test26
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -d 0
               -w ${TEMP}/${PROJECT_NAME}Test26.mha )

#######
# Test26-Compare
#######
add_test( ${PROJECT_NAME}Test26-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test26.mha
               -b ${BASELINE}/${PROJECT_NAME}Test26.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test26-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test26 )

#######
# Test27
#######
add_test( ${PROJECT_NAME}Test27
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -x 6
               -w ${TEMP}/${PROJECT_NAME}Test27.mha )

#######
# Test27-Compare
#######
add_test( ${PROJECT_NAME}Test27-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test27.mha
               -b ${BASELINE}/${PROJECT_NAME}Test27.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test27-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test27 )

#######
# Test28
#######
add_test( ${PROJECT_NAME}Test28
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -d 1
               -w ${TEMP}/${PROJECT_NAME}Test28.mha )

#######
# Test28-Compare
#######
add_test( ${PROJECT_NAME}Test28-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test28.mha
               -b ${BASELINE}/${PROJECT_NAME}Test28.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test28-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test28 )

#######
# Test29
#######
add_test( ${PROJECT_NAME}Test29
            ${PROJ_EXE}
               ${INPUT_DATA}/ES0015_Large_Subs.mha
               -d 2
               -w ${TEMP}/${PROJECT_NAME}Test29.mha )

#######
# Test29-Compare
#######
add_test( ${PROJECT_NAME}Test29-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test29.mha
               -b ${BASELINE}/${PROJECT_NAME}Test29.mha
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}Test29-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test29 )


######
#MIDAS TESTS
########

#######
# MIDAS-Test1
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test1
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -i -1 1 -100 100
               -w ${TEMP}/${PROJECT_NAME}Test1.mha )

#######
# MIDAS-Test1-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test1.mha
               -b MIDAS{${PROJECT_NAME}Test1.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test1 )

#######
# MIDAS-Test2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test2
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 5
               -w ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -a 0 1 MIDAS{ES0015_Large_Subs.mha.md5}
               -I ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -w ${TEMP}/${PROJECT_NAME}Test2.mha )

#######
# MIDAS-Test2-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test2.mha
               -b MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test2 )

#######
# MIDAS-Test3
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test3
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -n -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test3.mha )

#######
# MIDAS-Test3-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test3.mha
               -b MIDAS{${PROJECT_NAME}Test3.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test3 )


#######
# MIDAS-Test4
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test4
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -N -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test4.mha )

#######
# MIDAS-Test4-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test4-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test4.mha
               -b MIDAS{${PROJECT_NAME}Test4.mha.md5}
               -i 0.01 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test4-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test4 )


#######
# MIDAS-Test5
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test5
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -f 10
                 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test5.mha )

#######
# MIDAS-Test5-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test5-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test5.mha
               -b MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test5-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test5 )

#######
# MIDAS-Test6
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test6
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -a 0.95 0.05 MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test6.mha )

#######
# MIDAS-Test6-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test6-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test6.mha
               -b MIDAS{${PROJECT_NAME}Test6.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test6-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test6 )

#######
# MIDAS-Test7
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test7
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -t -1 -0.33 -1 1
               -w ${TEMP}/${PROJECT_NAME}Test7.mha )

#######
# MIDAS-Test7-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test7-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test7.mha
               -b MIDAS{${PROJECT_NAME}Test7.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test7-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test7 )

#######
# MIDAS-Test8
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test8
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -m -1 -0.33 MIDAS{ES0015_Large_Subs.mha.md5} 0
               -w ${TEMP}/${PROJECT_NAME}Test8.mha )

#######
# MIDAS-Test8-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test8-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test8.mha
               -b MIDAS{${PROJECT_NAME}Test8.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test8-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test8 )

#######
# MIDAS-Test9
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test9
            ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 0 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test9.mha )

#######
# MIDAS-Test9-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test9-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test9.mha
               -b MIDAS{${PROJECT_NAME}Test9.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test9-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test9 )

#######
# MIDAS-Test10
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test10
            ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 1 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test10.mha )

#######
# MIDAS-Test10-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test10-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test10.mha
               -b MIDAS{${PROJECT_NAME}Test10.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test10-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test10 )

#######
# MIDAS-Test11
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test11
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 10
               -w ${TEMP}/${PROJECT_NAME}Test11.mha )

#######
# MIDAS-Test11-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test11-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test11.mha
               -b MIDAS{${PROJECT_NAME}Test11.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test11-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test11 )

#######
# MIDAS-Test12
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test12
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -B 5 2 0
               -w ${TEMP}/${PROJECT_NAME}Test12.mha )

#######
# MIDAS-Test12-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test12-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test12.mha
               -b MIDAS{${PROJECT_NAME}Test12.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test12-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test12 )

#######
# MIDAS-Test13
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test13
            ${PROJ_EXE}
               MIDAS{im0001.crop2.mha.md5}
               -z 1 2 3
               -w ${TEMP}/${PROJECT_NAME}Test13.mha )

#######
# MIDAS-Test13-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test13-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test13.mha
               -b MIDAS{${PROJECT_NAME}Test13.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test13-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test13 )

#######
# MIDAS-Test14
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test14
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -l 10 ${TEMP}/${PROJECT_NAME}Test14.txt )

#######
# MIDAS-Test14-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test14-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test14.txt
            -b MIDAS{${PROJECT_NAME}Test14.txt.md5})
set_tests_properties( ${PROJECT_NAME}MIDAS-Test14-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test14 )

#######
# MIDAS-Test15
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test15
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -L 100 -1 0.01 ${TEMP}/${PROJECT_NAME}Test15.txt )

#######
# MIDAS-Test15-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test15-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test15.txt
            -b MIDAS{${PROJECT_NAME}Test15.txt.md5})
set_tests_properties( ${PROJECT_NAME}MIDAS-Test15-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test15 )

#######
# Test16
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test16
            ${PROJ_EXE}
               MIDAS{CroppedWholeLungCTScan.mhd.md5}
               -c 100 10
               -w ${TEMP}/${PROJECT_NAME}Test16.mha
               MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5}
               )

#######
# MIDAS-Test16-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test16-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test16.mha
               -b MIDAS{${PROJECT_NAME}Test16.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test16-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test16 )

#######
# MIDAS-Test17
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test17
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -C 100 3 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test17.mha )

#######
# MIDAS-Test17-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test17-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test17.mha
               -b MIDAS{${PROJECT_NAME}Test17.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test17-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test17 )

#######
# MIDAS-Test18
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test18
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -r 2
               -w ${TEMP}/${PROJECT_NAME}Test18.mha )

#######
# MIDAS-Test18-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test18-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test18.mha
               -b MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test18-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test18 )

#######
# MIDAS-Test19
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test19
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -R MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test19.mha )

#######
# MIDAS-Test19-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test19-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test19.mha
               -b MIDAS{${PROJECT_NAME}Test19.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test19-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test19 )

#######
# MIDAS-Test20
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test20
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -s 0 -1 -0.33 1 70 50 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test20.mha )

#######
# MIDAS-Test20-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test20-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test20.mha
               -b MIDAS{${PROJECT_NAME}Test20.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test20-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test20 )

#######
# MIDAS-Test21
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test21
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -O -10 -10 0
               -w ${TEMP}/${PROJECT_NAME}Test21.mha )

#######
# MIDAS-Test21-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test21-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test21.mha
               -b MIDAS{${PROJECT_NAME}Test21.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test21-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test21 )

#######
# MIDAS-Test23
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test23
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -P 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test23.mha )

#######
# MIDAS-Test23-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test23-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test23.mha
               -b MIDAS{${PROJECT_NAME}Test23.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test23-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test23 )

#######
# MIDAS-Test24
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test24
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -u MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test24.mha )

#######
# MIDAS-Test24-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test24-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test24.mha
               -b MIDAS{${PROJECT_NAME}Test24.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test24-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test24 )

#######
# MIDAS-Test25
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test25
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -p 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test25.mha )

#######
# MIDAS-Test25-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test25-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test25.mha
               -b MIDAS{${PROJECT_NAME}Test25.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test25-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test25 )

#######
# MIDAS-Test26
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test26
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 0
               -w ${TEMP}/${PROJECT_NAME}Test26.mha )

#######
# MIDAS-Test26-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test26-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test26.mha
               -b MIDAS{${PROJECT_NAME}Test26.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test26-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test26 )

#######
# MIDAS-Test27
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test27
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -x 6
               -w ${TEMP}/${PROJECT_NAME}Test27.mha )

#######
# MIDAS-Test27-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test27-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test27.mha
               -b MIDAS{${PROJECT_NAME}Test27.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test27-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test27 )

#######
# MIDAS-Test28
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test28
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 1
               -w ${TEMP}/${PROJECT_NAME}Test28.mha )

#######
# MIDAS-Test28-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test28-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test28.mha
               -b MIDAS{${PROJECT_NAME}Test28.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test28-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test28 )

#######
# MIDAS-Test29
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test29
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 2
               -w ${TEMP}/${PROJECT_NAME}Test29.mha )

#######
# MIDAS-Test29-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test29-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test29.mha
               -b MIDAS{${PROJECT_NAME}Test29.mha.md5}
               -i 0.001 )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test29-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test29 )
