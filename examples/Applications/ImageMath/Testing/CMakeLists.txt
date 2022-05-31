##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

include_regular_expression( "^.*$" )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1
itk_add_test(
            NAME ${MODULE_NAME}Test1
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -i -1 1 -100 100
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha )

# Test1-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test1-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test1-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test1 )

# Test2
itk_add_test(
            NAME ${MODULE_NAME}Test2
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -b 5
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2_meanfield.mha
               -a 0 1 DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -I ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2_meanfield.mha
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha )

# Test2-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test2-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test2-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test2 )

# Test3
itk_add_test(
            NAME ${MODULE_NAME}Test3
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -S 1234
               -n -1 1 10 0.2
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha )

# Test3-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test3-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test3-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test3 )

# Test4
itk_add_test(
            NAME ${MODULE_NAME}Test4
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -S 1234
               -N -1 1 9.8 10.2
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha )

# Test4-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test4-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4.mha}
               -i 0.01 )
set_tests_properties( ${MODULE_NAME}Test4-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test4 )

# Test5
itk_add_test(
            NAME ${MODULE_NAME}Test5
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -f 10
                 DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test5.mha )

# Test5-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test5-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test5.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test5.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test5-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test5 )

# Test6
itk_add_test(
            NAME ${MODULE_NAME}Test6
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -a 0.95 0.05 DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test5.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test6.mha )

# Test6-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test6-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test6.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test6.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test6-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test6 )

# Test7
itk_add_test(
            NAME ${MODULE_NAME}Test7
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -t -1 -0.33 -1 1
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test7.mha )

# Test7-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test7-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test7.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test7.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test7-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test7 )

# Test8
itk_add_test(
            NAME ${MODULE_NAME}Test8
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -m -1 -0.33 DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha} 0
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test8.mha )

# Test8-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test8-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test8.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test8.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test8-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test8 )

# Test9
itk_add_test(
            NAME ${MODULE_NAME}Test9
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/GDS0015_1.mha}
               -M 0 11 255 0
               -W 2 ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test9.mha )

# Test9-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test9-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test9.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test9.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test9-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test9 )

# Test10
itk_add_test(
            NAME ${MODULE_NAME}Test10
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/GDS0015_1.mha}
               -M 1 5 255 0
               -W 2 ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test10.mha )

# Test10-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test10-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test10.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test10.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test10-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test10 )

# Test11
itk_add_test(
            NAME ${MODULE_NAME}Test11
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -b 10
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test11.mha )

# Test11-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test11-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test11.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test11.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test11-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test11 )

# Test12
itk_add_test(
            NAME ${MODULE_NAME}Test12
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -B 5 2 0
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test12.mha )

# Test12-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test12-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test12.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test12.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test12-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test12 )

# Test13
itk_add_test(
            NAME ${MODULE_NAME}Test13
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/im0001.crop2.mha}
               -z 0.7 1.4 3
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test13.mha )

# Test13-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test13-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test13.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test13.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test13-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test13 )

# Test14
itk_add_test(
            NAME ${MODULE_NAME}Test14
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -l 10 ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test14.txt )

# Test14-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test14-Compare
            COMMAND ${TubeTK_CompareTextFiles_EXE}
              CompareTextFiles
              -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test14.txt
              -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test14.txt} )
set_tests_properties( ${MODULE_NAME}Test14-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test14 )

# Test15
itk_add_test(
            NAME ${MODULE_NAME}Test15
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -L 100 -1 0.01 ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test15.txt )

# Test15-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test15-Compare
            COMMAND ${TubeTK_CompareTextFiles_EXE}
              CompareTextFiles
              -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test15.txt
              -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test15.txt} )
set_tests_properties( ${MODULE_NAME}Test15-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test15 )

# Test16
itk_add_test(
            NAME ${MODULE_NAME}Test16
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
               -c 100 10
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test16.mha )

# Test16-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test16-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test16.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test16.mha}
               -i 0.5
               -n 300 )
set_tests_properties( ${MODULE_NAME}Test16-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test16 )

# Test17
itk_add_test(
            NAME ${MODULE_NAME}Test17
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -C 100 3 DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test17.mha )

# Test17-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test17-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test17.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test17.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test17-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test17 )

# Test18
itk_add_test(
            NAME ${MODULE_NAME}Test18
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -r 2
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test18.mha )

# Test18-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test18-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test18.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test18.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test18-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test18 )

# Test19
itk_add_test(
            NAME ${MODULE_NAME}Test19
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -R DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test18.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test19.mha )

# Test19-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test19-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test19.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test19.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test19-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test19 )

# Test20
itk_add_test(
            NAME ${MODULE_NAME}Test20
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -s 0 -1 -0.33 1 70 50 0
               -W 2 ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test20.mha )

# Test20-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test20-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test20.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test20.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test20-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test20 )

# Test21
itk_add_test(
            NAME ${MODULE_NAME}Test21
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -O -10 -10 0
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test21.mha )

# Test21-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test21-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test21.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test21.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test21-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test21 )

# Test23
itk_add_test(
            NAME ${MODULE_NAME}Test23
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -P 0 DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test23.mha )

# Test23-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test23-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test23.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test23.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test23-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test23 )

# Test24
itk_add_test(
            NAME ${MODULE_NAME}Test24
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -P 0 DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test24.mha )

# Test24-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test24-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test24.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test24.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test24-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test24 )

# Test25
itk_add_test(
            NAME ${MODULE_NAME}Test25
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -p 0 DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test25.mha )

# Test25-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test25-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test25.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test25.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test25-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test25 )

# Test26
itk_add_test(
            NAME ${MODULE_NAME}Test26
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -d 0
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test26.mha )

# Test26-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test26-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test26.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test26.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test26-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test26 )

# Test27
itk_add_test(
            NAME ${MODULE_NAME}Test27
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -x 6
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test27.mha )

# Test27-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test27-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test27.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test27.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test27-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test27 )

# Test28
itk_add_test(
            NAME ${MODULE_NAME}Test28
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -d 1
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test28.mha )

# Test28-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test28-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test28.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test28.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test28-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test28 )

# Test29
itk_add_test(
            NAME ${MODULE_NAME}Test29
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -d 2
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test29.mha )

# Test29-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test29-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test29.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test29.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test29-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test29 )

# Test30
itk_add_test(
            NAME ${MODULE_NAME}Test30
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
               -g 5
               -w ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test30.mha )

# Test30-Compare
itk_add_test(
            NAME ${MODULE_NAME}Test30-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test30.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test30.mha}
               -i 0.001 )
set_tests_properties( ${MODULE_NAME}Test30-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}Test30 )
