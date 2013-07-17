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
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Filtering )
endif()

set( BASE_FILTERING_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseFilteringTests )

set( BASE_FILTERING_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseFilteringHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseFilteringTests_SRCS
  tubeBaseFilteringPrintTest.cxx
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx
  itktubeAnisotropicHybridDiffusionImageFilterTest.cxx
  itktubeExtractTubePointsSpatialObjectFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTestNew.cxx
  itktubeSubSampleTubeSpatialObjectFilterTest.cxx
  itktubeSubSampleTubeTreeSpatialObjectFilterTest.cxx
  itktubeTubeEnhancingDiffusion2DImageFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest2.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
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

add_test( itktubeSheetnessMeasureImageFilterTest
  ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest )

MIDAS_ADD_TEST( itktubeExtractTubePointsSpatialObjectFilterTest
  ${BASE_FILTERING_TESTS}
    itktubeExtractTubePointsSpatialObjectFilterTest
      MIDAS{VascularNetwork.tre.md5} )

MIDAS_ADD_TEST( itktubeSheetnessMeasureImageFilterTest2
  ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimitiveObjectsSheetnessMeasure.mha
      ${TEMP}/PrimitiaveObjectsPrimaryEigenVector.mha )

MIDAS_ADD_TEST( itktubeSheetnessMeasureImageFilterTest3
  ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      MIDAS{UltrasoundCroppedResampled.mha.md5}
      ${TEMP}/UltrasoundCroppedResampledSheetnessMeasure.mha
      ${TEMP}/UltrasoundCroppedResampledPrimaryEigenVector.mha 2.0)

MIDAS_ADD_TEST( itktubeStructureTensorRecursiveGaussianImageFilterTest
  ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimaryEigenVectorImage.mha
      ${TEMP}/PrimaryEigenValueImage.mha )

MIDAS_ADD_TEST( itktubeStructureTensorRecursiveGaussianImageFilterTest2
  ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{CylinderSynthetic.mha.md5}
      ${TEMP}/CylinderPrimaryEigenVectorImage.mha
      ${TEMP}/CylinderPrimaryEigenValueImage.mha )

MIDAS_ADD_TEST( itktubeStructureTensorRecursiveGaussianImageFilterTest3
  ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTestNew )

MIDAS_ADD_TEST( itktubeSubSampleTubeSpatialObjectFilterTest
  ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeSpatialObjectFilterTest
      MIDAS{Branch-truth-new.tre.md5}
      ${TEMP}/itkSubSampleTubeSpatialObjectFilterTest.tre )

MIDAS_ADD_TEST( itktubeSubSampleTubeTreeSpatialObjectFilterTest
  ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeTreeSpatialObjectFilterTest
      MIDAS{Branch-truth-new.tre.md5}
      ${TEMP}/itkSubSampleTubeTreeSpatialObjectFilterTest.tre )

MIDAS_ADD_TEST( itktubeTubeEnhancingDiffusion2DImageFilterTest
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterTest.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
    --compareIntensityTolerance 0.001
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{ES0015_1_Subs.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha )

MIDAS_ADD_TEST( itktubeTubeEnhancingDiffusion2DImageFilterRetinaTest
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    --compareIntensityTolerance 0.001
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001.crop2.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

MIDAS_ADD_TEST( itktubeTubeEnhancingDiffusion2DImageFilterRetina10Test
  ${BASE_FILTERING_TESTS}
    --compare MIDAS{itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    --compareIntensityTolerance 0.001
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001_n10.crop2.mha.md5}
      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

MIDAS_ADD_TEST( itktubeAnisotropicHybridDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itktubeAnisotropicHybridDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTScanHybridDiffused.mha
 MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )


MIDAS_ADD_TEST( itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTCoherenceEnhanced.mha
 MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )


MIDAS_ADD_TEST( itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
 ${BASE_FILTERING_TESTS}
 itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
 MIDAS{CroppedWholeLungCTScan.mhd.md5}
 ${TEMP}/CroppedWholeLungCTEdgeEnhanced.mha
 MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )
