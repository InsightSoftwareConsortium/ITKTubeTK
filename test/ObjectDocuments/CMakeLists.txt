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

if( Python3_EXECUTABLE )
  itk_add_test( PrintSelf-ObjectDocuments ${Python3_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/ObjectDocuments )
endif( Python3_EXECUTABLE )

set( tubeObjectDocumentsTests_SRCS
  tubeObjectDocumentsPrintTest.cxx )

CreateTestDriver( tubeObjectDocumentsHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeObjectDocumentsHeaderTest.cxx" )
itk_add_test( NAME tubeObjectDocumentsHeaderTest
  COMMAND tubeObjectDocumentsHeaderTestDriver
    tubeObjectdocumentsHeaderTest )

CreateTestDriver( tubeObjectDocuments
  "${TubeTK-Test_LIBRARIES}"
  "${tubeObjectDocumentsTests_SRCS}" )

itk_add_test( NAME tubeObjectDocumentsPrintTest
  COMMAND tubeObjectDocumentsTestDriver
    tubeObjectDocumentsPrintTest )
