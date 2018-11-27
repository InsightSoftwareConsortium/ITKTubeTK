##############################################################################
#
# Library:   TubeTKLib
#
# Copyright Kitware Inc.
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

set( TEMP ${TubeTKLib_BINARY_DIR}/Temporary )

set( CompareImages_EXE
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:CompareImages> )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test1
  COMMAND ${CompareImages_EXE}
    -t DATA{${TubeTKLib_DATA_ROOT}/ES0015_Large.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/ES0015_Large.mha} )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test2
  COMMAND ${CompareImages_EXE}
    -n 1000
    -r 1
    -i 50
    -t DATA{${TubeTKLib_DATA_ROOT}/im0001.crop.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop.mha} )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test3
  COMMAND ${CompareImages_EXE}
    -t DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop.mha} )
set_tests_properties( CompareImages-Test3 PROPERTIES WILL_FAIL true )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test4
  COMMAND ${CompareImages_EXE}
    -t DATA{${TubeTKLib_DATA_ROOT}/im0001.crop.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop.mha},DATA{${TubeTKLib_DATA_ROOT}/im0001.crop.mha} )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test5
  COMMAND ${CompareImages_EXE}
    -t DATA{${TubeTKLib_DATA_ROOT}/im0001.crop.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop.mha},DATA{${TubeTKLib_DATA_ROOT}/im0001_n20.crop.mha} )
set_tests_properties( CompareImages-Test5 PROPERTIES WILL_FAIL true )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test6
  COMMAND ${CompareImages_EXE}
    -t DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop.mha}
    -b ${TEMP}/filedoesnotexist.mha )
set_tests_properties( CompareImages-Test6 PROPERTIES WILL_FAIL true )

ExternalData_Add_Test( TubeTKLibData
  NAME CompareImages-Test7
  COMMAND ${CompareImages_EXE}
    -c 1.0e-4
    -t DATA{${TubeTKLib_DATA_ROOT}/ResampleImageTest2.mha}
    -b DATA{${TubeTKLib_DATA_ROOT}/ResampleImageTest2SpacingOff.mha} )
