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
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

if( Python3_EXECUTABLE )
  itk_add_test( PrintSelf-TubeTKCommon ${Python3_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/Common )
endif( Python3_EXECUTABLE )

CreateTestDriver( tubeCommonHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeCommonHeaderTest.cxx" )
itk_add_test( NAME tubeCommonHeaderTest
 COMMAND tubeCommonHeaderTestDriver
   tubeCommonHeaderTest )

set( tubeCommonTestDriver_SRCS
  tubeCommonPrintTest.cxx
  tubeMacroTest.cxx
  tubeMessageTest.cxx
  tubeObjectTest.cxx )

CreateTestDriver( tubeCommon
  "${TubeTK-Test_LIBRARIES}"
  "${tubeCommonTestDriver_SRCS}" )

itk_add_test( NAME tubeCommonPrintTest
  COMMAND tubeCommonTestDriver
    tubeCommonPrintTest )

itk_add_test( NAME tubeMacroTest
  COMMAND tubeCommonTestDriver
    tubeMacroTest )

itk_add_test( NAME tubeMessageTest
  COMMAND tubeCommonTestDriver
    tubeMessageTest )

itk_add_test( NAME tubeObjectTest
  COMMAND tubeCommonTestDriver
    tubeObjectTest )
