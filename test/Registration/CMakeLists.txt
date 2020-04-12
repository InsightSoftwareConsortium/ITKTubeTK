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
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

if( PYTHON_EXECUTABLE )
  itk_add_test( PrintSelf-Registration ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/Registration )
endif( PYTHON_EXECUTABLE )

CreateTestDriver( tubeRegistrationHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeRegistrationHeaderTest.cxx" )
itk_add_test( NAME tubeRegistrationHeaderTest
  COMMAND tubeRegistrationHeaderTestDriver
    tubeRegistrationHeaderTest )

set( tubeRegistrationTests_SRCS
  tubeRegistrationPrintTest.cxx
  itktubeImageToTubeRigidMetricPerformanceTest.cxx
  itktubeImageToTubeRigidMetricTest.cxx
  itktubeImageToTubeRigidRegistrationPerformanceTest.cxx
  itktubeImageToTubeRigidRegistrationTest.cxx
  itktubePointsToImageTest.cxx
  itktubeSyntheticTubeImageGenerationTest.cxx
  itktubeTubeAngleOfIncidenceWeightFunctionTest.cxx
  itktubeTubeExponentialResolutionWeightFunctionTest.cxx
  itktubeTubeParametricExponentialResolutionWeightFunctionTest.cxx
  itktubeTubeParametricExponentialWithBoundsResolutionWeightFunctionTest.cxx
  itktubeTubePointWeightsCalculatorTest.cxx
  itktubeTubeToTubeTransformFilterTest.cxx )

if( TubeTK_USE_VTK )
  list( APPEND tubeRegistrationTests_SRCS
    itktubeAnisotropicDiffusiveRegistrationGenerateTestingImages.cxx
    itktubeAnisotropicDiffusiveRegistrationRegularizationTest.cxx )
endif( TubeTK_USE_VTK )

CreateTestDriver( tubeRegistration
  "${TubeTK-Test_LIBRARIES}"
  "${tubeRegistrationTests_SRCS}" )
if( TubeTK_USE_VTK )
  target_link_libraries( tubeRegistrationTestDriver LINK_PUBLIC
    ${VTK_LIBRARIES} )
endif( TubeTK_USE_VTK )

itk_add_test( NAME tubeRegistrationPrintTest
  COMMAND tubeRegistrationTestDriver
    tubeRegistrationPrintTest )

itk_add_test(
  NAME itktubeTubeToTubeTransformFilterTest
  COMMAND tubeRegistrationTestDriver
    --compare
      DATA{${TubeTK_DATA_ROOT}/itktubeTubeToTubeTransformFilter.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeToTubeTransformFilter.mha
    itktubeTubeToTubeTransformFilterTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeToTubeTransformFilter.tre
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeTubeToTubeTransformFilter.mha
      0.2 0.1 0.1 5 -5 5
      1 )

itk_add_test(
  NAME itktubeImageToTubeRigidRegistrationTest
  COMMAND tubeRegistrationTestDriver
    itktubeImageToTubeRigidRegistrationTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeImageToTubeRigidRegistrationOutputTube.tre
      ${ITK_TEST_OUTPUT_DIR}/itktubeImageToTubeRigidRegistrationOutputImage.mha )

itk_add_test(
  NAME itktubeImageToTubeRigidRegistrationPerformanceTest
  COMMAND tubeRegistrationTestDriver
    itktubeImageToTubeRigidRegistrationPerformanceTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeImageToTubeRigidRegistrationPerformance )

if( TubeTK_USE_VTK )
  itk_add_test(
    NAME itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise
    COMMAND tubeRegistrationTestDriver
      --compare
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise.mha
        DATA{${TubeTK_DATA_ROOT}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise.mha}
      itktubeAnisotropicDiffusiveRegistrationRegularizationTest
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise.mha
        0 0
        5 0.125 1 )

  itk_add_test(
    NAME
      itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraight
    COMMAND tubeRegistrationTestDriver
      --compare
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraight.mha
        DATA{${TubeTK_DATA_ROOT}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraight.mha}
      itktubeAnisotropicDiffusiveRegistrationRegularizationTest
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestStraight.mha
        0.1 0
        5 0.125 1 )

  itk_add_test(
    NAME
      itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
    COMMAND tubeRegistrationTestDriver
      --compare
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise.mha
        DATA{${TubeTK_DATA_ROOT}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise.mha}
      itktubeAnisotropicDiffusiveRegistrationRegularizationTest
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise.mha
        0 0.5
        5 0.125 1 )
  itk_add_test(
    NAME
      itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngled
    COMMAND tubeRegistrationTestDriver
      --compare
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngled.mha
        DATA{${TubeTK_DATA_ROOT}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngled.mha}
      itktubeAnisotropicDiffusiveRegistrationRegularizationTest
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngled.mha
        0.1 0.5
        5 0.125 1 )

  itk_add_test(
    NAME
      itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian
    COMMAND tubeRegistrationTestDriver
      --compare
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian.mha
        DATA{${TubeTK_DATA_ROOT}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian.mha}
      itktubeAnisotropicDiffusiveRegistrationRegularizationTest
        ${ITK_TEST_OUTPUT_DIR}/itktubeAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian.mha
        0.1 0.5
        5 0.125 0 )
endif()

itk_add_test(
  NAME itktubePointsToImageTest
  COMMAND tubeRegistrationTestDriver
    itktubePointsToImageTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/Branch-truth-new-points.mha )

itk_add_test(
  NAME itktubeImageToTubeRigidMetricTest1
  COMMAND tubeRegistrationTestDriver
    itktubeImageToTubeRigidMetricTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      2.30028 )

itk_add_test(
  NAME itktubeImageToTubeRigidMetricTest2
  COMMAND tubeRegistrationTestDriver
    itktubeImageToTubeRigidMetricTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticTransformedVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      0.0739628 )

itk_add_test(
  NAME itktubeImageToTubeRigidMetricPerformanceTest
  COMMAND tubeRegistrationTestDriver
    itktubeImageToTubeRigidMetricPerformanceTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${ITK_TEST_OUTPUT_DIR}/itkImageToTubeRigidMetricPerformance.txt )

itk_add_test(
  NAME itktubeSyntheticTubeImageGenerationTest
  COMMAND tubeRegistrationTestDriver
    itktubeSyntheticTubeImageGenerationTest
      ${ITK_TEST_OUTPUT_DIR}/SyntheticBlurredTubeImage.mha
      ${ITK_TEST_OUTPUT_DIR}/SyntheticVesselTube.tre
      ${ITK_TEST_OUTPUT_DIR}/SyntheticVesselTubeImage.mha
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${ITK_TEST_OUTPUT_DIR}/SyntheticTransformedVesselTubeImage.mha )

itk_add_test(
  NAME itktubeTubeExponentialResolutionWeightFunctionTest
  COMMAND tubeRegistrationTestDriver
    itktubeTubeExponentialResolutionWeightFunctionTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeExponentialResolutionWeightFunctionTest.csv )

itk_add_test(
  NAME
    itktubeTubeParametricExponentialResolutionWeightFunctionTest
  COMMAND tubeRegistrationTestDriver
    itktubeTubeParametricExponentialResolutionWeightFunctionTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeParametricExponentialResolutionWeightFunctionTest.csv )

itk_add_test(
  NAME
    itktubeTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
  COMMAND tubeRegistrationTestDriver
    itktubeTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeParametricExponentialWithBoundsResolutionWeightFunctionTest.csv )

itk_add_test(
  NAME itktubeTubePointWeightsCalculatorTest
  COMMAND tubeRegistrationTestDriver
    itktubeTubePointWeightsCalculatorTest
      ${ITK_TEST_OUTPUT_DIR}/itktubePointWeightsCalculatorTest.csv )

itk_add_test(
  NAME itktubeTubeAngleOfIncidenceWeightFunctionTest
  COMMAND tubeRegistrationTestDriver
    itktubeTubeAngleOfIncidenceWeightFunctionTest
      DATA{${TubeTK_DATA_ROOT}/AlmostMaxInhale01.vessels.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeAngleOfIncidenceWeightFunctionTest )
