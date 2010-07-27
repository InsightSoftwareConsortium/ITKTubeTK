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
  add_test( PrintSelf-BaseSegmentationVTree ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl  
    ${TubeTK_SOURCE_DIR}/Base/SegmentationVTree 
    )
endif( TCL_TCLSH )

set( BASE_SEGMENTATION_VTREE_TESTS 
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseSegmentationVTreeTests 
  )

set( BASE_SEGMENTATION_HEADER_VTREE_TEST 
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseSegmentationVTreeHeaderTest
  )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Segmentation/VTree/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseSegmentationVTree_SRCS
  tubeBaseSegmentationVTreePrintTest.cxx
  itkBlur3DImageFunctionTest.cxx
  itkMatrixMathTest.cxx
  itkOptBrentTest.cxx
  )

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Numerics 
  ${TubeTK_SOURCE_DIR}/Base/Segmentation/VTree
  )

add_executable( tubeBaseSegmentationVTreeHeaderTest 
  tubeBaseSegmentationVTreeHeaderTest.cxx 
  )
target_link_libraries( tubeBaseSegmentationVTreeHeaderTest
  ITKIO 
  ITKNumerics 
  ITKBasicFilters 
  )

add_executable( tubeBaseSegmentationVTreeTests
  tubeBaseSegmentationVTreeTests.cxx 
  ${tubeBaseSegmentationVTree_SRCS} 
  )
target_link_libraries( tubeBaseSegmentationVTreeTests
  ITKIO 
  ITKNumerics 
  ITKBasicFilters 
  VTree3D
  )

add_test( tubeBaseSegmentationVTreeHeaderTest 
  ${BASE_SEGMENTATION_HEADER_VTREE_TEST} 
  )

add_test( tubeBaseSegmentationVTreePrintTest 
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeBaseSegmentationVTreePrintTest 
  )

add_test( itkBlur3DImageFunctionTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  --compare ${BASELINE}/itkBlur3DImageFunctionTest.mha
            ${TEMP}/itkBlur3DImageFunctionTest.mha
  itkBlur3DImageFunctionTest
    ${TEMP}/itkBlur3DImageFunctionTest.mha 
  )

add_test( itkMatrixMathTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  itkMatrixMathTest
  )

add_test( itkOptBrentTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  itkOptBrentTest
  )
