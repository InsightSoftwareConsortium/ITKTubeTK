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

find_package( Tclsh )
if( TCL_TCLSH )
  add_test( PrintSelf-BaseClassification ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl  
    ${TubeTK_SOURCE_DIR}/Base/Classification )
endif( TCL_TCLSH )

set( BASE_CLASSIFICATION_TESTS 
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseClassificationTests 
  )

set( BASE_CLASSIFICATION_HEADER_TEST 
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseClassificationHeaderTest
  )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Classification/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseClassification_SRCS
  tubeBaseClassificationPrintTest.cxx
  itkFeatureGeneratingImageFunctionTest.cxx 
  itkStandardFeatureGeneratingImageFunctionTest.cxx 
  )

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Classification
  )

add_executable( tubeBaseClassificationHeaderTest 
  tubeBaseClassificationHeaderTest.cxx 
  )
target_link_libraries( tubeBaseClassificationHeaderTest
  ITKIO 
  ITKBasicFilters 
  )

add_executable( tubeBaseClassificationTests
  tubeBaseClassificationTests.cxx
  ${tubeBaseClassification_SRCS} 
  )
target_link_libraries( tubeBaseClassificationTests
  ITKIO 
  ITKBasicFilters 
  )

add_test( tubeBaseClassificationHeaderTest 
  ${BASE_CLASSIFICATION_HEADER_TEST} 
  )

add_test( tubeBaseClassificationPrintTest 
  ${BASE_CLASSIFICATION_TESTS}
  tubeBaseClassificationPrintTest 
  )