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

if( PYTHON_EXECUTABLE )
  itk_add_test( PrintSelf-MetaIO ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/MetaIO )
endif( PYTHON_EXECUTABLE )

CreateTestDriver( tubeMetaIOHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeMetaIOHeaderTest.cxx" )

itk_add_test( NAME tubeMetaIOHeaderTest
  COMMAND tubeMetaIOHeaderTestDriver
    tubeMetaIOHeaderTest )

set( tubeMetaIOTests_SRCS
  tubeMetaIOPrintTest.cxx
  itktubeMetaLDATest.cxx
  itktubeMetaNJetLDATest.cxx
  itktubeMetaClassPDFTest.cxx
  itktubeMetaRidgeSeedTest.cxx
  itktubeMetaTubeExtractorTest.cxx )

CreateTestDriver( tubeMetaIO
  "${TubeTK-Test_LIBRARIES}"
  "${tubeMetaIOTests_SRCS}" )

itk_add_test( NAME tubeMetaIOPrintTest
  COMMAND tubeMetaIOTestDriver
  tubeMetaIOPrintTest )

itk_add_test( NAME itktubeMetaLDATest
  COMMAND tubeMetaIOTestDriver
    itktubeMetaLDATest
      ${ITK_TEST_OUTPUT_DIR}/itktubeMetaLDATest.mlda )

itk_add_test( NAME itktubeMetaNJetLDATest
  COMMAND tubeMetaIOTestDriver
    itktubeMetaNJetLDATest
      ${ITK_TEST_OUTPUT_DIR}/itktubeMetaNJetLDATest.mnda )

itk_add_test( NAME itktubeMetaClassPDFTest
  COMMAND tubeMetaIOTestDriver
    itktubeMetaClassPDFTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeMetaClassPDFTest.mha )

itk_add_test( NAME itktubeMetaRidgeSeedTest
  COMMAND tubeMetaIOTestDriver
    itktubeMetaRidgeSeedTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeMetaRidgeSeedTest.mrs )

itk_add_test( NAME itktubeMetaTubeExtractorTest
  COMMAND tubeMetaIOTestDriver
    itktubeMetaTubeExtractorTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeMetaTubeExtractorTest.mtp )
