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

include_regular_expression( "^.*$" )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseFiltering ${PYTHON_EXECUTABLE}
    ${TubeTKLib_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTKLib_SOURCE_DIR}/Base/Filtering )
endif( PYTHON_EXECUTABLE )

set( BASE_FILTERING_TESTS
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseFilteringTests> )

set( BASE_FILTERING_HEADER_TEST
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseFilteringHeaderTest> )

set( TEMP ${TubeTKLib_BINARY_DIR}/Temporary )

set( tubeBaseFilteringTests_SRCS
  tubeBaseFilteringTests.cxx
  tubeBaseFilteringPrintTest.cxx
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx
  itktubeAnisotropicHybridDiffusionImageFilterTest.cxx
  itktubeCVTImageFilterTest.cxx
  itktubeExtractTubePointsSpatialObjectFilterTest.cxx
  itktubeFFTGaussianDerivativeIFFTFilterTest.cxx
  itktubeRidgeFFTFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest2.cxx
  itktubeShrinkWithBlendingImageFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTestNew.cxx
  itktubeSubSampleTubeSpatialObjectFilterTest.cxx
  itktubeSubSampleTubeTreeSpatialObjectFilterTest.cxx
  itktubeTortuositySpatialObjectFilterTest.cxx
  itktubeTubeEnhancingDiffusion2DImageFilterTest.cxx )

# Add tests of filters based on the ArrayFire Library
if( TubeTKLib_USE_ARRAYFIRE )
  list( APPEND tubeBaseFilteringTests_SRCS
    itktubeGPUArrayFireGaussianDerivativeFilterTest.cxx )
endif( TubeTKLib_USE_ARRAYFIRE )

include_directories(
  ${TubeTKLib_SOURCE_DIR}/Common
  ${TubeTKLib_SOURCE_DIR}/Numerics
  ${TubeTKLib_SOURCE_DIR}/Filtering )

add_executable( tubeBaseFilteringHeaderTest
  tubeBaseFilteringHeaderTest.cxx )
target_link_libraries( tubeBaseFilteringHeaderTest
  TubeTKLibFiltering ${ITK_LIBRARIES} )
if( TubeTKLib_USE_ARRAYFIRE )
  target_link_libraries( tubeBaseFilteringHeaderTest
    ${ArrayFire_LIBRARIES} )
endif( TubeTKLib_USE_ARRAYFIRE )

add_executable( tubeBaseFilteringTests
  ${tubeBaseFilteringTests_SRCS} )
target_link_libraries( tubeBaseFilteringTests 
  TubeTKLibFiltering
  ${ITK_LIBRARIES}
  ITKIOMeta
  ITKIOSpatialObjects )
if( TubeTKLib_USE_ARRAYFIRE )
  target_link_libraries( tubeBaseFilteringTests ${ArrayFire_LIBRARIES} )
endif( TubeTKLib_USE_ARRAYFIRE )

add_test( NAME tubeBaseFilteringHeaderTest
  COMMAND ${BASE_FILTERING_HEADER_TEST} )

add_test( NAME tubeBaseFilteringPrintTest
  COMMAND ${BASE_FILTERING_TESTS}
    tubeBaseFilteringPrintTest )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeCVTImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeCVTImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/GDS0015_1.mha}
      ${TEMP}/itktubeCVTImageFilterTest.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeExtractTubePointsSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeExtractTubePointsSpatialObjectFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/VascularNetwork.tre} )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest1
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 0 0.1 0.1
      DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 1 0.1 1
      DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      2 2 2 2
      DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha )

# GPU ArrayFire Gaussian Derivative Tests - begin
if( TubeTKLib_USE_ARRAYFIRE )
  ExternalData_Add_Test( TubeTKLibData
    NAME itktubeGPUArrayFireGaussianDerivativeFilterTest1
    COMMAND ${BASE_FILTERING_TESTS}
      --compareIntensityTolerance 0.01
      --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeGPUArrayFireGaussianDerivativeFilterTest1.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest1.mha
      itktubeGPUArrayFireGaussianDerivativeFilterTest
        0 0 0.1 0.1
        DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest1.mha )

  ExternalData_Add_Test( TubeTKLibData
    NAME itktubeGPUArrayFireGaussianDerivativeFilterTest2
    COMMAND ${BASE_FILTERING_TESTS}
      --compareIntensityTolerance 0.01
      --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeGPUArrayFireGaussianDerivativeFilterTest2.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest2.mha
      itktubeGPUArrayFireGaussianDerivativeFilterTest
        0 1 0.1 1
        DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest2.mha )

  ExternalData_Add_Test( TubeTKLibData
    NAME itktubeGPUArrayFireGaussianDerivativeFilterTest3
    COMMAND ${BASE_FILTERING_TESTS}
      --compareIntensityTolerance 0.01
      --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeGPUArrayFireGaussianDerivativeFilterTest3.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest3.mha
      itktubeGPUArrayFireGaussianDerivativeFilterTest
        2 2 2 2
        DATA{${TubeTKLib_DATA_ROOT}/im0001.mha}
        ${TEMP}/itktubeGPUArrayFireGaussianDerivativeFilterTest3.mha )

endif( TubeTKLib_USE_ARRAYFIRE )
# GPU ArrayFire Gaussian Derivative Tests - end

#--compare DATA{${TubeTKLib_DATA_ROOT}/itktubeRidgeFFTFilterTest1.mha}
#${TEMP}/itktubeRidgeFFTFilterTest1.mha
ExternalData_Add_Test( TubeTKLibData
  NAME itktubeRidgeFFTFilterTest1
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeRidgeFFTFilterTest
      2
      DATA{${TubeTKLib_DATA_ROOT}/Branch.n010.mha}
      ${TEMP}/itktubeRidgeFFTFilterTest1_Ridgeness.mha
      ${TEMP}/itktubeRidgeFFTFilterTest1_Roundness.mha
      ${TEMP}/itktubeRidgeFFTFilterTest1_Curvature.mha
      ${TEMP}/itktubeRidgeFFTFilterTest1_Levelness.mha )

add_test( NAME itktubeSheetnessMeasureImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeSheetnessMeasureImageFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      DATA{${TubeTKLib_DATA_ROOT}/PrimitiveObjects.mha}
      ${TEMP}/PrimitiveObjectsSheetnessMeasure.mha
      ${TEMP}/PrimitiaveObjectsPrimaryEigenVector.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeSheetnessMeasureImageFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      DATA{${TubeTKLib_DATA_ROOT}/UltrasoundCroppedResampled.mha}
      ${TEMP}/UltrasoundCroppedResampledSheetnessMeasure.mha
      ${TEMP}/UltrasoundCroppedResampledPrimaryEigenVector.mha 2.0 )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeShrinkWithBlendingImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
  --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeShrinkWithBlendingImageFilterTest.mha}
  ${TEMP}/itktubeShrinkWithBlendingImageFilterTest.mha
  itktubeShrinkWithBlendingImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/im0001.crop.contrast.mha}
      ${TEMP}/itktubeShrinkWithBlendingImageFilterTest.mha
      ${TEMP}/itktubeShrinkWithBlendingImageFilterTest-IndexImage.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/PrimitiveObjects.mha}
      ${TEMP}/PrimaryEigenVectorImage.mha
      ${TEMP}/PrimaryEigenValueImage.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/CylinderSynthetic.mha}
      ${TEMP}/CylinderPrimaryEigenVectorImage.mha
      ${TEMP}/CylinderPrimaryEigenValueImage.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTestNew )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeSubSampleTubeSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeSpatialObjectFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch-truth-new.tre}
      ${TEMP}/itkSubSampleTubeSpatialObjectFilterTest.tre )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeSubSampleTubeTreeSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeTreeSpatialObjectFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch-truth-new.tre}
      ${TEMP}/itkSubSampleTubeTreeSpatialObjectFilterTest.tre )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeEnhancingDiffusion2DImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterTest.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/ES0015_1_Subs.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterTest.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeEnhancingDiffusion2DImageFilterRetinaTest
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/im0001.crop2.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeEnhancingDiffusion2DImageFilterRetina10Test
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTKLib_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/im0001_n10.crop2.mha}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeAnisotropicHybridDiffusionImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicHybridDiffusionImageFilterTest
     DATA{${TubeTKLib_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${TEMP}/CroppedWholeLungCTScanHybridDiffused.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
     DATA{${TubeTKLib_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${TEMP}/CroppedWholeLungCTCoherenceEnhanced.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
     DATA{${TubeTKLib_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${TEMP}/CroppedWholeLungCTEdgeEnhanced.mha )

add_test( NAME itktubeTortuositySpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeTortuositySpatialObjectFilterTest )
