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

include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseFiltering ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    )
endif()

set( BASE_FILTERING_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseFilteringTests )

set( BASE_FILTERING_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseFilteringHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseFilteringTests_SRCS
  tubeBaseFilteringPrintTest.cxx
  itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx
  itkAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx
  itkAnisotropicHybridDiffusionImageFilterTest.cxx
  itkStructureTensorRecursiveGaussianImageFilterTest.cxx
  itkStructureTensorRecursiveGaussianImageFilterTestNew.cxx
  itkTubeEnhancingDiffusion2DImageFilterTest.cxx
  itkSheetnessMeasureImageFilterTest.cxx
  itkSheetnessMeasureImageFilterTest2.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_executable( tubeBaseFilteringHeaderTest
  tubeBaseFilteringHeaderTest.cxx )
target_link_libraries( tubeBaseFilteringHeaderTest
  ${ITK_LIBRARIES} )

add_executable( tubeBaseFilteringTests
  tubeBaseFilteringTests.cxx
  ${tubeBaseFilteringTests_SRCS} )
target_link_libraries( tubeBaseFilteringTests
  ${ITK_LIBRARIES} )

add_test( tubeBaseFilteringHeaderTest
  ${BASE_FILTERING_HEADER_TEST} )

add_test( tubeBaseFilteringPrintTest
  ${BASE_FILTERING_TESTS}
  tubeBaseFilteringPrintTest )

add_test( itkSheetnessMeasureImageFilterTest
  ${BASE_FILTERING_TESTS}
  itkSheetnessMeasureImageFilterTest )

MIDAS_ADD_TEST( itkSheetnessMeasureImageFilterTest2
  ${BASE_FILTERING_TESTS}
    itkSheetnessMeasureImageFilterTest2
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimitiveObjectsSheetnessMeasure.mha
      ${TEMP}/PrimitiaveObjectsPrimaryEigenVector.mha )

MIDAS_ADD_TEST( itkSheetnessMeasureImageFilterTest3
  ${BASE_FILTERING_TESTS}
    itkSheetnessMeasureImageFilterTest2
      MIDAS{UltrasoundCroppedResampled.mha.md5}
      ${TEMP}/UltrasoundCroppedResampledSheetnessMeasure.mha
      ${TEMP}/UltrasoundCroppedResampledPrimaryEigenVector.mha 2.0)


MIDAS_ADD_TEST( itkStructureTensorRecursiveGaussianImageFilterTest
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimaryEigenVectorImage.mha
      ${TEMP}/PrimaryEigenValueImage.mha )

MIDAS_ADD_TEST( itkStructureTensorRecursiveGaussianImageFilterTest2
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{CylinderSynthetic.mha.md5}
      ${TEMP}/CylinderPrimaryEigenVectorImage.mha
      ${TEMP}/CylinderPrimaryEigenValueImage.mha )

MIDAS_ADD_TEST( itkStructureTensorRecursiveGaussianImageFilterTest3
  ${BASE_FILTERING_TESTS}
    itkStructureTensorRecursiveGaussianImageFilterTestNew )

MIDAS_ADD_TEST( itkTubeEnhancingDiffusion2DImageFilterTest
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterTest.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{ES0015_1_Subs.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha )

MIDAS_ADD_TEST( itkTubeEnhancingDiffusion2DImageFilterRetinaTest
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001.crop2.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

MIDAS_ADD_TEST( itkTubeEnhancingDiffusion2DImageFilterRetina10Test
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    itkTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001_n10.crop2.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

MIDAS_ADD_TEST( itkAnisotropicHybridDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itkAnisotropicHybridDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTScanHybridDiffused.mha )

MIDAS_ADD_TEST( itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTCoherenceEnhanced.mha )

MIDAS_ADD_TEST( itkAnisotropicEdgeEnhancementDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itkAnisotropicEdgeEnhancementDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTEdgeEnhanced.mha )
