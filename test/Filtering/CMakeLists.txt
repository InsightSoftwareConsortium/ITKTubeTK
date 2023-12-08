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
  itk_add_test( PrintSelf-TubeTKFiltering ${Python3_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/Filtering )
endif( Python3_EXECUTABLE )

CreateTestDriver( tubeFilteringHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeFilteringHeaderTest.cxx" )
itk_add_test( NAME tubeFilteringHeaderTest
  COMMAND tubeFilteringHeaderTestDriver
    tubeFilteringHeaderTest )

set( tubeFilteringTests_SRCS
  tubeFilteringPrintTest.cxx
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
  itktubeSubSampleSpatialObjectFilterTest.cxx
  itktubeTortuositySpatialObjectFilterTest.cxx
  itktubeTubeEnhancingDiffusion2DImageFilterTest.cxx
  tubeTubeMathFiltersTest.cxx )

CreateTestDriver( tubeFiltering
  "${TubeTK-Test_LIBRARIES}"
  "${tubeFilteringTests_SRCS}" )

itk_add_test( NAME tubeFilteringPrintTest
  COMMAND tubeFilteringTestDriver
    tubeFilteringPrintTest )

itk_add_test(
  NAME itktubeCVTImageFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeCVTImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/GDS0015_1.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeCVTImageFilterTest.mha )

itk_add_test(
  NAME itktubeExtractTubePointsSpatialObjectFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeExtractTubePointsSpatialObjectFilterTest
      DATA{${TubeTK_DATA_ROOT}/VascularNetwork.tre} )

itk_add_test(
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest1
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 0 1 1
      DATA{${TubeTK_DATA_ROOT}/im0001.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha )

itk_add_test(
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest2
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 1 0.1 1
      DATA{${TubeTK_DATA_ROOT}/im0001.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha )

itk_add_test(
  NAME itktubeFFTGaussianDerivativeIFFTFilterTest3
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.01
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      2 2 2 2
      DATA{${TubeTK_DATA_ROOT}/im0001.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha )

#--compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeFFTFilterTest1.mha}
#${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFilterTest1.mha
itk_add_test(
  NAME itktubeRidgeFFTFilterTest1
  COMMAND tubeFilteringTestDriver
    itktubeRidgeFFTFilterTest
      2
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFilterTest1_Ridgeness.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFilterTest1_Roundness.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFilterTest1_Curvature.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFilterTest1_Levelness.mha )
set_tests_properties(itktubeRidgeFFTFilterTest1
  PROPERTIES TIMEOUT 10)

add_test( NAME itktubeSheetnessMeasureImageFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeSheetnessMeasureImageFilterTest )

itk_add_test(
  NAME itktubeSheetnessMeasureImageFilterTest2
  COMMAND tubeFilteringTestDriver
    itktubeSheetnessMeasureImageFilterTest2
      DATA{${TubeTK_DATA_ROOT}/PrimitiveObjects.mha}
      ${ITK_TEST_OUTPUT_DIR}/PrimitiveObjectsSheetnessMeasure.mha
      ${ITK_TEST_OUTPUT_DIR}/PrimitiaveObjectsPrimaryEigenVector.mha )

itk_add_test(
  NAME itktubeSheetnessMeasureImageFilterTest3
  COMMAND tubeFilteringTestDriver
    itktubeSheetnessMeasureImageFilterTest2
      DATA{${TubeTK_DATA_ROOT}/UltrasoundCroppedResampled.mha}
      ${ITK_TEST_OUTPUT_DIR}/UltrasoundCroppedResampledSheetnessMeasure.mha
      ${ITK_TEST_OUTPUT_DIR}/UltrasoundCroppedResampledPrimaryEigenVector.mha 2.0 )

itk_add_test(
  NAME itktubeShrinkWithBlendingImageFilterTest
  COMMAND tubeFilteringTestDriver
  --compare DATA{${TubeTK_DATA_ROOT}/itktubeShrinkWithBlendingImageFilterTest.mha}
    ${ITK_TEST_OUTPUT_DIR}/itktubeShrinkWithBlendingImageFilterTest.mha
  itktubeShrinkWithBlendingImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeShrinkWithBlendingImageFilterTest.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeShrinkWithBlendingImageFilterTest-IndexImage.mha )

itk_add_test(
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/PrimitiveObjects.mha}
      ${ITK_TEST_OUTPUT_DIR}/PrimaryEigenVectorImage.mha
      ${ITK_TEST_OUTPUT_DIR}/PrimaryEigenValueImage.mha )

itk_add_test(
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest2
  COMMAND tubeFilteringTestDriver
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/CylinderSynthetic.mha}
      ${ITK_TEST_OUTPUT_DIR}/CylinderPrimaryEigenVectorImage.mha
      ${ITK_TEST_OUTPUT_DIR}/CylinderPrimaryEigenValueImage.mha )

itk_add_test(
  NAME itktubeStructureTensorRecursiveGaussianImageFilterTest3
  COMMAND tubeFilteringTestDriver
    itktubeStructureTensorRecursiveGaussianImageFilterTestNew )

itk_add_test(
  NAME itktubeSubSampleTubeSpatialObjectFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeSubSampleTubeSpatialObjectFilterTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itkSubSampleTubeSpatialObjectFilterTest.tre )

itk_add_test(
  NAME itktubeSubSampleSpatialObjectFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeSubSampleSpatialObjectFilterTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itkSubSampleSpatialObjectFilterTest.tre )

itk_add_test(
  NAME itktubeTubeEnhancingDiffusion2DImageFilterTest
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_1_Subs.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterTest.mha )

itk_add_test(
  NAME itktubeTubeEnhancingDiffusion2DImageFilterRetinaTest
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop2.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

itk_add_test(
  NAME itktubeTubeEnhancingDiffusion2DImageFilterRetina10Test
  COMMAND tubeFilteringTestDriver
    --compareIntensityTolerance 0.001
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/im0001_n10.crop2.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

itk_add_test(
  NAME itktubeAnisotropicHybridDiffusionImageFilterTest
  COMMAND tubeFilteringTestDriver
   itktubeAnisotropicHybridDiffusionImageFilterTest
     DATA{${TubeTK_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${ITK_TEST_OUTPUT_DIR}/CroppedWholeLungCTScanHybridDiffused.mha )

itk_add_test(
  NAME itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
  COMMAND tubeFilteringTestDriver
   itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
     DATA{${TubeTK_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${ITK_TEST_OUTPUT_DIR}/CroppedWholeLungCTCoherenceEnhanced.mha )

itk_add_test(
  NAME itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
  COMMAND tubeFilteringTestDriver
   itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
     DATA{${TubeTK_DATA_ROOT}/CroppedWholeLungCTScan.mhd,CroppedWholeLungCTScan.raw}
     ${ITK_TEST_OUTPUT_DIR}/CroppedWholeLungCTEdgeEnhanced.mha )

add_test( NAME itktubeTortuositySpatialObjectFilterTest
  COMMAND tubeFilteringTestDriver
    itktubeTortuositySpatialObjectFilterTest )

add_test( NAME tubeTubeMathFiltersTest
  COMMAND tubeFilteringTestDriver
  tubeTubeMathFiltersTest )
