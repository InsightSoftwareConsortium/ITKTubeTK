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
  ${TubeTK_BINARY_DIR}/bin/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( TEXTCOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/TextCompareCommand )

#######
# Test1
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test1
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -i -1 1 -100 100
               -w ${TEMP}/${PROJECT_NAME}Test1.mha )

#######
# Test1-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test1.mha
               -b MIDAS{${PROJECT_NAME}Test1.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test1-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test1 )

#######
# Test2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test2
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 5
               -w ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -a 0 1 MIDAS{ES0015_Large_Subs.mha.md5}
               -I ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -w ${TEMP}/${PROJECT_NAME}Test2.mha )

#######
# Test2-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test2.mha
               -b MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test2-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test2 )

#######
# Test3
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test3
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -n -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test3.mha )

#######
# Test3-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test3.mha
               -b MIDAS{${PROJECT_NAME}Test3.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test3-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test3 )


#######
# Test4
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test4
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -N -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test4.mha )

#######
# Test4-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test4-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test4.mha
               -b MIDAS{${PROJECT_NAME}Test4.mha.md5}
               -i 0.01 )
set_property(TEST ${PROJECT_NAME}-Test4-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test4 )


#######
# Test5
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test5
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -f 10
                 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test5.mha )

#######
# Test5-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test5-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test5.mha
               -b MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test5-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test5 )

#######
# Test6
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test6
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -a 0.95 0.05 MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test6.mha )

#######
# Test6-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test6-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test6.mha
               -b MIDAS{${PROJECT_NAME}Test6.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test6-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test6 )

#######
# Test7
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test7
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -t -1 -0.33 -1 1
               -w ${TEMP}/${PROJECT_NAME}Test7.mha )

#######
# Test7-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test7-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test7.mha
               -b MIDAS{${PROJECT_NAME}Test7.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test7-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test7 )

#######
# Test8
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test8
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -m -1 -0.33 MIDAS{ES0015_Large_Subs.mha.md5} 0
               -w ${TEMP}/${PROJECT_NAME}Test8.mha )

#######
# Test8-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test8-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test8.mha
               -b MIDAS{${PROJECT_NAME}Test8.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test8-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test8 )

#######
# Test9
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test9
            ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 0 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test9.mha )

#######
# Test9-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test9-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test9.mha
               -b MIDAS{${PROJECT_NAME}Test9.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test9-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test9 )

#######
# Test10
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test10
            ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 1 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test10.mha )

#######
# Test10-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test10-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test10.mha
               -b MIDAS{${PROJECT_NAME}Test10.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test10-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test10 )

#######
# Test11
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test11
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 10
               -w ${TEMP}/${PROJECT_NAME}Test11.mha )

#######
# Test11-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test11-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test11.mha
               -b MIDAS{${PROJECT_NAME}Test11.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test11-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test11 )

#######
# Test12
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test12
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -B 5 2 0
               -w ${TEMP}/${PROJECT_NAME}Test12.mha )

#######
# Test12-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test12-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test12.mha
               -b MIDAS{${PROJECT_NAME}Test12.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test12-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test12 )

#######
# Test13
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test13
            ${PROJ_EXE}
               MIDAS{im0001.crop2.mha.md5}
               -z 1 2 3
               -w ${TEMP}/${PROJECT_NAME}Test13.mha )

#######
# Test13-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test13-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test13.mha
               -b MIDAS{${PROJECT_NAME}Test13.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test13-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test13 )

#######
# Test14
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test14
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -l 10 ${TEMP}/${PROJECT_NAME}Test14.txt )

#######
# Test14-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test14-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test14.txt
            -b MIDAS{${PROJECT_NAME}Test14.txt.md5})
set_property(TEST ${PROJECT_NAME}-Test14-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test14 )

#######
# Test15
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test15
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -L 100 -1 0.01 ${TEMP}/${PROJECT_NAME}Test15.txt )

#######
# Test15-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test15-Compare
            ${TEXTCOMPARE_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test15.txt
            -b MIDAS{${PROJECT_NAME}Test15.txt.md5})
set_property(TEST ${PROJECT_NAME}-Test15-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test15 )

#######
# Test16
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test16
            ${PROJ_EXE}
               MIDAS{CroppedWholeLungCTScan.mhd.md5}
               -c 100 10
               -w ${TEMP}/${PROJECT_NAME}Test16.mha
               MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5}
               )

#######
# Test16-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test16-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test16.mha
               -b MIDAS{${PROJECT_NAME}Test16.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test16-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test16 )

#######
# Test17
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test17
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -C 100 3 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test17.mha )

#######
# Test17-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test17-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test17.mha
               -b MIDAS{${PROJECT_NAME}Test17.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test17-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test17 )

#######
# Test18
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test18
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -r 2
               -w ${TEMP}/${PROJECT_NAME}Test18.mha )

#######
# Test18-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test18-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test18.mha
               -b MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test18-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test18 )

#######
# Test19
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test19
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -R MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test19.mha )

#######
# Test19-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test19-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test19.mha
               -b MIDAS{${PROJECT_NAME}Test19.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test19-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test19 )

#######
# Test20
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test20
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -s 0 -1 -0.33 1 70 50 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test20.mha )

#######
# Test20-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test20-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test20.mha
               -b MIDAS{${PROJECT_NAME}Test20.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test20-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test20 )

#######
# Test21
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test21
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -O -10 -10 0
               -w ${TEMP}/${PROJECT_NAME}Test21.mha )

#######
# Test21-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test21-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test21.mha
               -b MIDAS{${PROJECT_NAME}Test21.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test21-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test21 )

#######
# Test23
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test23
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -P 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test23.mha )

#######
# Test23-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test23-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test23.mha
               -b MIDAS{${PROJECT_NAME}Test23.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test23-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test23 )

#######
# Test24
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test24
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -u MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test24.mha )

#######
# Test24-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test24-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test24.mha
               -b MIDAS{${PROJECT_NAME}Test24.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test24-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test24 )

#######
# Test25
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test25
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -p 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test25.mha )

#######
# Test25-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test25-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test25.mha
               -b MIDAS{${PROJECT_NAME}Test25.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test25-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test25 )

#######
# Test26
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test26
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 0
               -w ${TEMP}/${PROJECT_NAME}Test26.mha )

#######
# Test26-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test26-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test26.mha
               -b MIDAS{${PROJECT_NAME}Test26.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test26-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test26 )

#######
# Test27
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test27
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -x 6
               -w ${TEMP}/${PROJECT_NAME}Test27.mha )

#######
# Test27-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test27-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test27.mha
               -b MIDAS{${PROJECT_NAME}Test27.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test27-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test27 )

#######
# Test28
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test28
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 1
               -w ${TEMP}/${PROJECT_NAME}Test28.mha )

#######
# Test28-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test28-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test28.mha
               -b MIDAS{${PROJECT_NAME}Test28.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test28-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test28 )

#######
# Test29
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test29
            ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 2
               -w ${TEMP}/${PROJECT_NAME}Test29.mha )

#######
# Test29-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}-Test29-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test29.mha
               -b MIDAS{${PROJECT_NAME}Test29.mha.md5}
               -i 0.001 )
set_property(TEST ${PROJECT_NAME}-Test29-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test29 )
