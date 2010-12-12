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
  itkTubeBlurImageFunctionTest.cxx
  tubeMatrixMathTest.cxx
  tubeOptBrentTest.cxx
  tubeOptBrentNDTest.cxx
  tubeOptGoldenMeanTest.cxx
  tubeOptParabolicFitTest.cxx
  itkTubeRidgeExtractorTest.cxx
  itkTubeRidgeExtractorTest2.cxx
  itkTubeRadiusExtractorTest.cxx
  itkTubeRadiusExtractorTest2.cxx
  tubeSplineApprox1DTest.cxx
  tubeSplineNDTest.cxx
  tubeUserFuncTest.cxx
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
  VTree
  )

add_test( tubeBaseSegmentationVTreeHeaderTest
  ${BASE_SEGMENTATION_HEADER_VTREE_TEST}
  )

add_test( tubeBaseSegmentationVTreePrintTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeBaseSegmentationVTreePrintTest
  )

add_test( itkTubeBlurImageFunctionTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  --compare ${BASELINE}/itkTubeBlurImageFunctionTest.mha
            ${TEMP}/itkTubeBlurImageFunctionTest.mha
  itkTubeBlurImageFunctionTest
    ${TEMP}/itkTubeBlurImageFunctionTest.mha
  )

add_test( tubeMatrixMathTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeMatrixMathTest
  )

add_test( tubeOptBrentTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeOptBrentTest
  )

add_test( tubeOptGoldenMeanTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeOptGoldenMeanTest
  )

add_test( tubeOptParabolicFitTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeOptParabolicFitTest
  )

add_test( tubeOptBrentNDTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeOptBrentNDTest
  )

add_test( tubeUserFuncTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  tubeUserFuncTest
  )

add_test( tubeSplineApprox1DTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  --compare ${BASELINE}/tubeSplineApprox1DTest.mha
            ${TEMP}/tubeSplineApprox1DTest.mha
  tubeSplineApprox1DTest
    ${TEMP}/tubeSplineApprox1DTest.mha
  )

add_test( tubeSplineNDTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  --compare ${BASELINE}/tubeSplineNDTest.mha
            ${TEMP}/tubeSplineNDTest.mha
  tubeSplineNDTest
    ${TEMP}/tubeSplineNDTest.mha
  )

add_test( itkTubeRidgeExtractorTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  --compare ${BASELINE}/itkTubeRidgeExtractorTest.mha
            ${TEMP}/itkTubeRidgeExtractorTest.mha
  itkTubeRidgeExtractorTest
    ${INPUT_DATA}/Branch.n010.sub.mha
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  )

add_test( itkTubeRidgeExtractorTest2
  ${BASE_SEGMENTATION_VTREE_TESTS}
  itkTubeRidgeExtractorTest2
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )

add_test( itkTubeRadiusExtractorTest
  ${BASE_SEGMENTATION_VTREE_TESTS}
  itkTubeRadiusExtractorTest
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )

add_test( itkTubeRadiusExtractorTest2
  ${BASE_SEGMENTATION_VTREE_TESTS}
  itkTubeRadiusExtractorTest2
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )
