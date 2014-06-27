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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseFiltering ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Filtering )
endif( PYTHON_EXECUTABLE )

set( BASE_FILTERING_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseFilteringTests> )

set( BASE_FILTERING_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseFilteringHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseFilteringTests_SRCS
  tubeBaseFilteringPrintTest.cxx
  itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx
  itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx
  itktubeAnisotropicHybridDiffusionImageFilterTest.cxx
  itktubeExtractTubePointsSpatialObjectFilterTest.cxx
  itktubeFFTGaussianDerivativeIFFTFilterTest.cxx
  itktubeRidgeFFTFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest.cxx
  itktubeSheetnessMeasureImageFilterTest2.cxx
  itktubeShrinkUsingMaxImageFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTest.cxx
  itktubeStructureTensorRecursiveGaussianImageFilterTestNew.cxx
  itktubeSubSampleTubeSpatialObjectFilterTest.cxx
  itktubeSubSampleTubeTreeSpatialObjectFilterTest.cxx
  itktubeTubeEnhancingDiffusion2DImageFilterTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_executable( tubeBaseFilteringHeaderTest
  tubeBaseFilteringHeaderTest.cxx )
target_link_libraries( tubeBaseFilteringHeaderTest
  ${ITK_LIBRARIES} )

SEMMacroBuildCLI(
  NAME tubeBaseFilteringTests
  ADDITIONAL_SRCS
    ${tubeBaseFilteringTests_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Numerics
  INSTALL_RUNTIME_DESTINATION lib
  INSTALL_LIBRARY_DESTINATION lib
  INSTALL_ARCHIVE_DESTINATION bin
  EXECUTABLE_ONLY
  )

add_test( NAME tubeBaseFilteringHeaderTest
  COMMAND ${BASE_FILTERING_HEADER_TEST} )

add_test( NAME tubeBaseFilteringPrintTest
  COMMAND ${BASE_FILTERING_TESTS}
    tubeBaseFilteringPrintTest )

Midas3FunctionAddTest( NAME itktubeExtractTubePointsSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeExtractTubePointsSpatialObjectFilterTest
      MIDAS{VascularNetwork.tre.md5} )

Midas3FunctionAddTest( NAME itktubeFFTGaussianDerivativeIFFTFilterTest1
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare MIDAS{itktubeFFTGaussianDerivativeIFFTFilterTest1.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 0 0.1 0.1
      MIDAS{im0001.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest1.mha )

Midas3FunctionAddTest( NAME itktubeFFTGaussianDerivativeIFFTFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare MIDAS{itktubeFFTGaussianDerivativeIFFTFilterTest2.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      0 1 0.1 1
      MIDAS{im0001.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest2.mha )

Midas3FunctionAddTest( NAME itktubeFFTGaussianDerivativeIFFTFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.01
    --compare MIDAS{itktubeFFTGaussianDerivativeIFFTFilterTest3.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha
    itktubeFFTGaussianDerivativeIFFTFilterTest
      2 2 2 2
      MIDAS{im0001.mha.md5}
      ${TEMP}/itktubeFFTGaussianDerivativeIFFTFilterTest3.mha )

    #--compare MIDAS{itktubeRidgeFFTFilterTest1.mha.md5}
      #${TEMP}/itktubeRidgeFFTFilterTest1.mha
Midas3FunctionAddTest( NAME itktubeRidgeFFTFilterTest1
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeRidgeFFTFilterTest
      2
      MIDAS{Branch.n010.mha.md5}
      ${TEMP}/itktubeRidgeFFTFilterTest1_Ridgeness.mha 
      ${TEMP}/itktubeRidgeFFTFilterTest1_Roundness.mha 
      ${TEMP}/itktubeRidgeFFTFilterTest1_Curvature.mha 
      ${TEMP}/itktubeRidgeFFTFilterTest1_Levelness.mha )

add_test( NAME itktubeSheetnessMeasureImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest )

Midas3FunctionAddTest( NAME itktubeSheetnessMeasureImageFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimitiveObjectsSheetnessMeasure.mha
      ${TEMP}/PrimitiaveObjectsPrimaryEigenVector.mha )

Midas3FunctionAddTest( NAME itktubeSheetnessMeasureImageFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSheetnessMeasureImageFilterTest2
      MIDAS{UltrasoundCroppedResampled.mha.md5}
      ${TEMP}/UltrasoundCroppedResampledSheetnessMeasure.mha
      ${TEMP}/UltrasoundCroppedResampledPrimaryEigenVector.mha 2.0 )

Midas3FunctionAddTest( NAME itktubeShrinkUsingMaxImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    --compare MIDAS{itktubeShrinkUsingMaxImageFilterTest.mha.md5}
      ${TEMP}/itktubeShrinkUsingMaxImageFilterTest.mha
    itktubeShrinkUsingMaxImageFilterTest
      MIDAS{im0001.crop.contrast.mha.md5}
      ${TEMP}/itktubeShrinkUsingMaxImageFilterTest.mha
      ${TEMP}/itktubeShrinkUsingMaxImageFilterTest-IndexImage.mha )

Midas3FunctionAddTest( NAME itktubeStructureTensorRecursiveGaussianImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{PrimitiveObjects.mha.md5}
      ${TEMP}/PrimaryEigenVectorImage.mha
      ${TEMP}/PrimaryEigenValueImage.mha )

Midas3FunctionAddTest( NAME itktubeStructureTensorRecursiveGaussianImageFilterTest2
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTest
      MIDAS{CylinderSynthetic.mha.md5}
      ${TEMP}/CylinderPrimaryEigenVectorImage.mha
      ${TEMP}/CylinderPrimaryEigenValueImage.mha )

Midas3FunctionAddTest( NAME itktubeStructureTensorRecursiveGaussianImageFilterTest3
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeStructureTensorRecursiveGaussianImageFilterTestNew )

Midas3FunctionAddTest( NAME itktubeSubSampleTubeSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeSpatialObjectFilterTest
      MIDAS{Branch-truth-new.tre.md5}
      ${TEMP}/itkSubSampleTubeSpatialObjectFilterTest.tre )

Midas3FunctionAddTest( NAME itktubeSubSampleTubeTreeSpatialObjectFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    itktubeSubSampleTubeTreeSpatialObjectFilterTest
      MIDAS{Branch-truth-new.tre.md5}
      ${TEMP}/itkSubSampleTubeTreeSpatialObjectFilterTest.tre )

Midas3FunctionAddTest( NAME itktubeTubeEnhancingDiffusion2DImageFilterTest
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare MIDAS{itktubeEnhancingDiffusion2DImageFilterTest.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{ES0015_1_Subs.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterTest.mha )

Midas3FunctionAddTest( NAME itktubeTubeEnhancingDiffusion2DImageFilterRetinaTest
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare MIDAS{itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001.crop2.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetinaTest.mha
      true )

Midas3FunctionAddTest( NAME itktubeTubeEnhancingDiffusion2DImageFilterRetina10Test
  COMMAND ${BASE_FILTERING_TESTS}
    --compareIntensityTolerance 0.001
    --compare MIDAS{itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
    itktubeTubeEnhancingDiffusion2DImageFilterTest
      MIDAS{im0001_n10.crop2.mha.md5}
      ${TEMP}/itktubeEnhancingDiffusion2DImageFilterRetina10Test.mha
      true )

Midas3FunctionAddTest( NAME itktubeAnisotropicHybridDiffusionImageFilterTest
 COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicHybridDiffusionImageFilterTest
     MIDAS{CroppedWholeLungCTScan.mhd.md5}
     ${TEMP}/CroppedWholeLungCTScanHybridDiffused.mha
     MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )

Midas3FunctionAddTest( NAME itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
 COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest
     MIDAS{CroppedWholeLungCTScan.mhd.md5}
     ${TEMP}/CroppedWholeLungCTCoherenceEnhanced.mha
     MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )

Midas3FunctionAddTest( NAME itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
 COMMAND ${BASE_FILTERING_TESTS}
   itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest
     MIDAS{CroppedWholeLungCTScan.mhd.md5}
     ${TEMP}/CroppedWholeLungCTEdgeEnhanced.mha
     MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )
