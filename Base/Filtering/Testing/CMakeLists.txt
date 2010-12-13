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
  add_test( PrintSelf-BaseFiltering ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl
    ${TubeTK_SOURCE_DIR}/Base/Filter )
endif( TCL_TCLSH )

set( BASE_FILTERING_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseFilteringTests )

set( BASE_FILTERING_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseFilteringHeaderTest )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Filtering/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseFilteringTests_SRCS
  tubeBaseFilteringPrintTest.cxx
  itkTubeEnhancingDiffusion2DImageFilterTest.cxx
  itkStructureTensorRecursiveGaussianImageFilterTest.cxx
  itkStructureTensorRecursiveGaussianImageFilterTestNew.cxx
  itkTubeBlurImageFunctionTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_executable( tubeBaseFilteringHeaderTest
  tubeBaseFilteringHeaderTest.cxx )
target_link_libraries( tubeBaseFilteringHeaderTest
  ITKIO ITKNumerics ITKBasicFilters )

add_executable( tubeBaseFilteringTests
  tubeBaseFilteringTests.cxx
  ${tubeBaseFilteringTests_SRCS} )
target_link_libraries( tubeBaseFilteringTests
  ITKIO ITKNumerics ITKBasicFilters )

add_test( tubeBaseFilteringHeaderTest
  ${BASE_FILTERING_HEADER_TEST} )

add_test( tubeBaseFilteringPrintTest
  ${BASE_FILTERING_TESTS}
  tubeBaseFilteringPrintTest )

add_test( itkTubeEnhancingDiffusion2DImageFilterTest
  ${BASE_FILTERING_TESTS}
    --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      ${INPUT_DATA}/ES0015_1.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha )

add_test( itkTubeEnhancingDiffusion2DImageFilterRetinaTest
  ${BASE_FILTERING_TESTS}
    --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      ${INPUT_DATA}/im0001.crop.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

add_test( itkTubeEnhancingDiffusion2DImageFilterRetina10Test
  ${BASE_FILTERING_TESTS}
    --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      ${INPUT_DATA}/im0001_n10.crop.mha
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

add_test( itkStructureTensorRecursiveGaussianImageFilterTest
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTest
      ${INPUT_DATA}/PrimitiveObjects.mha
      ${TEMP}/PrimaryEigenVectorImage.mha
      ${TEMP}/PrimaryEigenValueImage.mha )

add_test( itkStructureTensorRecursiveGaussianImageFilterTest2
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTest
      ${INPUT_DATA}/CylinderSynthetic.mha
      ${TEMP}/CylinderPrimaryEigenVectorImage.mha
      ${TEMP}/CylinderPrimaryEigenValueImage.mha )

add_test( itkStructureTensorRecursiveGaussianImageFilterTest3
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTestNew )

add_test( itkTubeBlurImageFunctionTest
  ${BASE_FILTERING_TESTS}
  --compare ${BASELINE}/itkTubeBlurImageFunctionTest.mha
            ${TEMP}/itkTubeBlurImageFunctionTest.mha
  itkTubeBlurImageFunctionTest
    ${TEMP}/itkTubeBlurImageFunctionTest.mha
  )
