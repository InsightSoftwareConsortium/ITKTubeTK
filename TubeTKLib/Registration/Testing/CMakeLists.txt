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
  add_test( PrintSelf-BaseRegistration ${PYTHON_EXECUTABLE}
    ${TubeTKLib_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTKLib_SOURCE_DIR}/Base/Registration )
endif( PYTHON_EXECUTABLE )

set( BASE_REGISTRATION_TESTS
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseRegistrationTests> )

set( CompareImages_EXE
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:CompareImages> )

set( TEMP ${TubeTKLib_BINARY_DIR}/Temporary )

set( tubeBaseRegistrationTests_SRCS
  tubeBaseRegistrationTests.cxx
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

# give a bit of tolerance
set( CompareImagesTolerance 0.000001 )

set( tubeBaseRegistration_ADDITIONAL_LIBRARIES )
if( TubeTKLib_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )

  set( tubeBaseRegistrationTests_SRCS
    ${tubeBaseRegistrationTests_SRCS}
    itktubeAnisotropicDiffusiveRegistrationGenerateTestingImages.cxx
    itktubeAnisotropicDiffusiveRegistrationRegularizationTest.cxx )

  set( tubeBaseRegistration_ADDITIONAL_LIBRARIES
    ${VTK_LIBRARIES} )
endif( TubeTKLib_USE_VTK )

add_executable( tubeBaseRegistrationTests
  ${tubeBaseRegistrationTests_SRCS} )
target_link_libraries( tubeBaseRegistrationTests
  TubeTKLibFiltering TubeTKLibRegistration
  ${ITK_LIBRARIES} ITKIOMeta ITKIOSpatialObjects ITKOptimizers
  ${tubeBaseRegistration_ADDITIONAL_LIBRARIES} )

if( TubeTKLib_USE_VTK )
  add_definitions( -DTubeTKLib_USE_VTK )
endif( TubeTKLib_USE_VTK )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeToTubeTransformFilterTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    --compare
      DATA{${TubeTKLib_DATA_ROOT}/itktubeTubeToTubeTransformFilter.mha}
      ${TEMP}/itktubeTubeToTubeTransformFilter.mha
    itktubeTubeToTubeTransformFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch-truth-new.tre}
      ${TEMP}/itktubeTubeToTubeTransformFilter.tre
      DATA{${TubeTKLib_DATA_ROOT}/Branch.n020.mha}
      ${TEMP}/itktubeTubeToTubeTransformFilter.mha
      0.2 0.1 0.1 5 -5 5
      1 )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeImageToTubeRigidRegistrationTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeImageToTubeRigidRegistrationTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTKLib_DATA_ROOT}/tube.tre}
      ${TEMP}/itktubeImageToTubeRigidRegistrationOutputTube.tre
      ${TEMP}/itktubeImageToTubeRigidRegistrationOutputImage.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeImageToTubeRigidRegistrationPerformanceTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeImageToTubeRigidRegistrationPerformanceTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTKLib_DATA_ROOT}/Branch-truth-new.tre}
      ${TEMP}/itktubeImageToTubeRigidRegistrationPerformance )

if( TubeTKLib_USE_VTK )
  ExternalData_Add_Test( TubeTKLibData
    NAME itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise
    COMMAND ${BASE_REGISTRATION_TESTS}
      itkAnisotropicDiffusiveRegistrationRegularizationTest
        ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
        0 0
        5 0.125 1 )
  ExternalData_Add_Test( TubeTKLibData
    NAME
     itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
      -b DATA{${TubeTKLib_DATA_ROOT}/Regularization_straight_noNoise_smoothedMotionField.mhd,Regularization_straight_noNoise_smoothedMotionField.zraw}
      -i ${CompareImagesTolerance} )
  set_tests_properties(
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
    PROPERTIES DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise )

  ExternalData_Add_Test( TubeTKLibData
    NAME
      itkAnisotropicDiffusiveRegistrationRegularizationTestStraight
    COMMAND ${BASE_REGISTRATION_TESTS}
      itkAnisotropicDiffusiveRegistrationRegularizationTest
        ${TEMP}/Regularization_straight_smoothedMotionField.mhd
        0.1 0
        5 0.125 1 )
  ExternalData_Add_Test( TubeTKLibData
    NAME itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/Regularization_straight_smoothedMotionField.mhd
      -b DATA{${TubeTKLib_DATA_ROOT}/Regularization_straight_smoothedMotionField.mhd,Regularization_straight_smoothedMotionField.zraw}
      -i ${CompareImagesTolerance} )
  set_tests_properties(
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
    PROPERTIES DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraight )

  ExternalData_Add_Test( TubeTKLibData
    NAME itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
    COMMAND ${BASE_REGISTRATION_TESTS}
      itkAnisotropicDiffusiveRegistrationRegularizationTest
        ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
        0 0.5
        5 0.125 1 )
  ExternalData_Add_Test( TubeTKLibData
    NAME
      itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
      -b DATA{${TubeTKLib_DATA_ROOT}/Regularization_angled_noNoise_smoothedMotionField.mhd,Regularization_angled_noNoise_smoothedMotionField.zraw}
      -i ${CompareImagesTolerance} )
  set_tests_properties(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
    PROPERTIES DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise )

  ExternalData_Add_Test( TubeTKLibData
    NAME
      itkAnisotropicDiffusiveRegistrationRegularizationTestAngled
    COMMAND ${BASE_REGISTRATION_TESTS}
      itkAnisotropicDiffusiveRegistrationRegularizationTest
        ${TEMP}/Regularization_angled_smoothedMotionField.mhd
        0.1 0.5
        5 0.125 1 )
  ExternalData_Add_Test( TubeTKLibData
    NAME itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/Regularization_angled_smoothedMotionField.mhd
      -b DATA{${TubeTKLib_DATA_ROOT}/Regularization_angled_smoothedMotionField.mhd,Regularization_angled_smoothedMotionField.zraw}
      -i ${CompareImagesTolerance} )
  set_tests_properties(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
    PROPERTIES DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngled )

  ExternalData_Add_Test( TubeTKLibData
    NAME itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian
    COMMAND ${BASE_REGISTRATION_TESTS}
      itkAnisotropicDiffusiveRegistrationRegularizationTest
        ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
        0.1 0.5
        5 0.125 0 )
  ExternalData_Add_Test( TubeTKLibData
    NAME
      itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
      -b DATA{${TubeTKLib_DATA_ROOT}/Regularization_angled_gaussian_smoothedMotionField.mhd,Regularization_angled_gaussian_smoothedMotionField.zraw}
      -i ${CompareImagesTolerance} )
  set_tests_properties(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
    PROPERTIES DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian )

endif( TubeTKLib_USE_VTK )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubePointsToImageTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubePointsToImageTest
      DATA{${TubeTKLib_DATA_ROOT}/Branch-truth-new.tre}
      ${TEMP}/Branch-truth-new-points.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeImageToTubeRigidMetricTest1
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeImageToTubeRigidMetricTest
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      2.30028 )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeImageToTubeRigidMetricTest2
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeImageToTubeRigidMetricTest
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticTransformedVesselTubeImage.mha}
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      0.0739628 )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeImageToTubeRigidMetricPerformanceTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeImageToTubeRigidMetricPerformanceTest
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${TEMP}/itkImageToTubeRigidMetricPerformance.txt )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeSyntheticTubeImageGenerationTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeSyntheticTubeImageGenerationTest
      ${TEMP}/SyntheticBlurredTubeImage.mha
      ${TEMP}/SyntheticVesselTube.tre
      ${TEMP}/SyntheticVesselTubeImage.mha
      DATA{${TubeTKLib_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${TEMP}/SyntheticTransformedVesselTubeImage.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeExponentialResolutionWeightFunctionTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeTubeExponentialResolutionWeightFunctionTest
      ${TEMP}/itktubeExponentialResolutionWeightFunctionTest.csv )

ExternalData_Add_Test( TubeTKLibData
  NAME
    itktubeTubeParametricExponentialResolutionWeightFunctionTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeTubeParametricExponentialResolutionWeightFunctionTest
      ${TEMP}/itktubeParametricExponentialResolutionWeightFunctionTest.csv )

ExternalData_Add_Test( TubeTKLibData
  NAME
    itktubeTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
      ${TEMP}/itktubeParametricExponentialWithBoundsResolutionWeightFunctionTest.csv )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubePointWeightsCalculatorTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeTubePointWeightsCalculatorTest
      ${TEMP}/itktubePointWeightsCalculatorTest.csv )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeTubeAngleOfIncidenceWeightFunctionTest
  COMMAND ${BASE_REGISTRATION_TESTS}
    itktubeTubeAngleOfIncidenceWeightFunctionTest
      DATA{${TubeTKLib_DATA_ROOT}/AlmostMaxInhale01.vessels.tre}
      ${TEMP}/itktubeAngleOfIncidenceWeightFunctionTest )
