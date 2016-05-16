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
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${PROJECT_NAME}> )

set( CompareImages_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

# Test1
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test1
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -i -1 1 -100 100
               -w ${TEMP}/${PROJECT_NAME}Test1.mha )

# Test1-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test1-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test1.mha
               -b MIDAS{${PROJECT_NAME}Test1.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test1-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test1 )

# Test2
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test2
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 5
               -w ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -a 0 1 MIDAS{ES0015_Large_Subs.mha.md5}
               -I ${TEMP}/${PROJECT_NAME}Test2_meanfield.mha
               -w ${TEMP}/${PROJECT_NAME}Test2.mha )

# Test2-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test2-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test2.mha
               -b MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test2-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test2 )

# Test3
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test3
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -n -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test3.mha )

# Test3-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test3-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test3.mha
               -b MIDAS{${PROJECT_NAME}Test3.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test3-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test3 )

# Test4
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test4
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -S 1234
               -N -1 1 10 0.2
               -w ${TEMP}/${PROJECT_NAME}Test4.mha )

# Test4-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test4-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test4.mha
               -b MIDAS{${PROJECT_NAME}Test4.mha.md5}
               -i 0.01 )
set_property( TEST ${PROJECT_NAME}-Test4-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test4 )

# Test5
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test5
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -f 10
                 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test5.mha )

# Test5-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test5-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test5.mha
               -b MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test5-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test5 )

# Test6
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test6
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -a 0.95 0.05 MIDAS{${PROJECT_NAME}Test5.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test6.mha )

# Test6-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test6-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test6.mha
               -b MIDAS{${PROJECT_NAME}Test6.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test6-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test6 )

# Test7
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test7
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -t -1 -0.33 -1 1
               -w ${TEMP}/${PROJECT_NAME}Test7.mha )

# Test7-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test7-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test7.mha
               -b MIDAS{${PROJECT_NAME}Test7.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test7-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test7 )

# Test8
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test8
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -m -1 -0.33 MIDAS{ES0015_Large_Subs.mha.md5} 0
               -w ${TEMP}/${PROJECT_NAME}Test8.mha )

# Test8-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test8-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test8.mha
               -b MIDAS{${PROJECT_NAME}Test8.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test8-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test8 )

# Test9
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test9
            COMMAND ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 0 10 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test9.mha )

# Test9-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test9-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test9.mha
               -b MIDAS{${PROJECT_NAME}Test9.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test9-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test9 )

# Test10
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test10
            COMMAND ${PROJ_EXE}
               MIDAS{GDS0015_1.mha.md5}
               -M 1 5 255 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test10.mha )

# Test10-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test10-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test10.mha
               -b MIDAS{${PROJECT_NAME}Test10.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test10-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test10 )

# Test11
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test11
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -b 10
               -w ${TEMP}/${PROJECT_NAME}Test11.mha )

# Test11-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test11-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test11.mha
               -b MIDAS{${PROJECT_NAME}Test11.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test11-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test11 )

# Test12
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test12
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -B 5 2 0
               -w ${TEMP}/${PROJECT_NAME}Test12.mha )

# Test12-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test12-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test12.mha
               -b MIDAS{${PROJECT_NAME}Test12.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test12-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test12 )

# Test13
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test13
            COMMAND ${PROJ_EXE}
               MIDAS{im0001.crop2.mha.md5}
               -z 0.7 1.4 3
               -w ${TEMP}/${PROJECT_NAME}Test13.mha )

# Test13-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test13-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test13.mha
               -b MIDAS{${PROJECT_NAME}Test13.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test13-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test13 )

# Test14
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test14
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -l 10 ${TEMP}/${PROJECT_NAME}Test14.txt )

# Test14-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test14-Compare
            COMMAND ${CompareTextFiles_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test14.txt
            -b MIDAS{${PROJECT_NAME}Test14.txt.md5} )
set_property( TEST ${PROJECT_NAME}-Test14-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test14 )

# Test15
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test15
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -L 100 -1 0.01 ${TEMP}/${PROJECT_NAME}Test15.txt )

# Test15-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test15-Compare
            COMMAND ${CompareTextFiles_EXE}
            -t ${TEMP}/${PROJECT_NAME}Test15.txt
            -b MIDAS{${PROJECT_NAME}Test15.txt.md5} )
set_property( TEST ${PROJECT_NAME}-Test15-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test15 )

# Test16
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test16
            COMMAND ${PROJ_EXE}
               MIDAS{CroppedWholeLungCTScan.mhd.md5}
               -c 100 10
               -w ${TEMP}/${PROJECT_NAME}Test16.mha
               MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )

# Test16-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test16-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test16.mha
               -b MIDAS{${PROJECT_NAME}Test16.mha.md5}
               -i 0.5
               -n 300 )
set_property( TEST ${PROJECT_NAME}-Test16-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test16 )

# Test17
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test17
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -C 100 3 MIDAS{${PROJECT_NAME}Test2.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test17.mha )

# Test17-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test17-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test17.mha
               -b MIDAS{${PROJECT_NAME}Test17.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test17-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test17 )

# Test18
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test18
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -r 2
               -w ${TEMP}/${PROJECT_NAME}Test18.mha )

# Test18-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test18-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test18.mha
               -b MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test18-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test18 )

# Test19
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test19
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -R MIDAS{${PROJECT_NAME}Test18.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test19.mha )

# Test19-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test19-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test19.mha
               -b MIDAS{${PROJECT_NAME}Test19.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test19-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test19 )

# Test20
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test20
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -s 0 -1 -0.33 1 70 50 0
               -W 2 ${TEMP}/${PROJECT_NAME}Test20.mha )

# Test20-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test20-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test20.mha
               -b MIDAS{${PROJECT_NAME}Test20.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test20-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test20 )

# Test21
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test21
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -O -10 -10 0
               -w ${TEMP}/${PROJECT_NAME}Test21.mha )

# Test21-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test21-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test21.mha
               -b MIDAS{${PROJECT_NAME}Test21.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test21-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test21 )

# Test23
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test23
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -P 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test23.mha )

# Test23-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test23-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test23.mha
               -b MIDAS{${PROJECT_NAME}Test23.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test23-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test23 )

# Test24
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test24
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -u MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test24.mha )

# Test24-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test24-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test24.mha
               -b MIDAS{${PROJECT_NAME}Test24.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test24-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test24 )

# Test25
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test25
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -p 0 MIDAS{ES0015_Large_Subs.mha.md5}
               -w ${TEMP}/${PROJECT_NAME}Test25.mha )

# Test25-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test25-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test25.mha
               -b MIDAS{${PROJECT_NAME}Test25.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test25-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test25 )

# Test26
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test26
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 0
               -w ${TEMP}/${PROJECT_NAME}Test26.mha )

# Test26-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test26-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test26.mha
               -b MIDAS{${PROJECT_NAME}Test26.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test26-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test26 )

# Test27
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test27
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -x 6
               -w ${TEMP}/${PROJECT_NAME}Test27.mha )

# Test27-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test27-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test27.mha
               -b MIDAS{${PROJECT_NAME}Test27.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test27-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test27 )

# Test28
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test28
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 1
               -w ${TEMP}/${PROJECT_NAME}Test28.mha )

# Test28-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test28-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test28.mha
               -b MIDAS{${PROJECT_NAME}Test28.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test28-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test28 )

# Test29
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test29
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -d 2
               -w ${TEMP}/${PROJECT_NAME}Test29.mha )

# Test29-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test29-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test29.mha
               -b MIDAS{${PROJECT_NAME}Test29.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test29-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test29 )

# Test30
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test30
            COMMAND ${PROJ_EXE}
               MIDAS{ES0015_Large_Subs.mha.md5}
               -g 5
               -w ${TEMP}/${PROJECT_NAME}Test30.mha )

# Test30-Compare
Midas3FunctionAddTest( NAME ${PROJECT_NAME}-Test30-Compare
            COMMAND ${CompareImages_EXE}
               -t ${TEMP}/${PROJECT_NAME}Test30.mha
               -b MIDAS{${PROJECT_NAME}Test30.mha.md5}
               -i 0.001 )
set_property( TEST ${PROJECT_NAME}-Test30-Compare
                      APPEND PROPERTY DEPENDS ${PROJECT_NAME}-Test30 )
