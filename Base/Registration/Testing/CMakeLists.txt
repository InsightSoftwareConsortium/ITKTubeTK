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

Include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Registration
  )

if( Slicer3_INSTALL_BIN_DIR )
  set( BASE_REGISTRATION_FILTERS_TESTS
    ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseRegistrationFiltersTests
    )
else()
  set( BASE_REGISTRATION_FILTERS_TESTS
    ${TubeTK_BINARY_DIR}/bin/tubeBaseRegistrationFiltersTests
    )
endif()

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseRegistrationFilters_SRCS
  itkImageToTubeRigidMetricTest.cxx
  itkImageToTubeRigidMetricPerformanceTest.cxx
  itkImageToTubeRigidRegistrationTest.cxx
  itkImageToTubeRigidRegistrationPerformanceTest.cxx
  itkSyntheticTubeImageGenerationsTest.cxx
  itkTubeExponentialResolutionWeightFunctionTest.cxx
  itkTubeParametricExponentialResolutionWeightFunctionTest.cxx
  itkTubeParametricExponentialWithBoundsResolutionWeightFunctionTest.cxx
  itkTubePointsToImageTest.cxx
  itkTubeToTubeTransformFilterTest.cxx
  )

# give a bit of tolerance
set( imageCompareTolerance 0.000001 )

if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )

  set( tubeBaseRegistrationFilters_SRCS
    ${tubeBaseRegistrationFilters_SRCS}
    itkAnisotropicDiffusiveRegistrationGenerateTestingImages.cxx
    itkAnisotropicDiffusiveRegistrationRegularizationTest.cxx
    )
endif()


add_executable( tubeBaseRegistrationFiltersTests
  tubeBaseRegistrationFiltersTests.cxx
  ${tubeBaseRegistrationFilters_SRCS}
  )
if( TubeTK_USE_VTK )
  add_definitions( -DTubeTK_USE_VTK )
endif()

set( tubeBaseRegistrationFilters_LIBS
  tubeBaseRegistrationFiltersTests
  ${ITK_LIBRARIES}
  )

if( TubeTK_USE_VTK )
  set( tubeBaseRegistrationFilters_LIBS
    ${tubeBaseRegistrationFilters_LIBS}
    ${VTK_LIBRARIES}
    )
endif()
target_link_libraries( ${tubeBaseRegistrationFilters_LIBS} )

MIDAS_ADD_TEST( itkTubeToTubeTransformFilterTest
  ${BASE_REGISTRATION_FILTERS_TESTS} itkTubeToTubeTransformFilterTest
  MIDAS{Branch-truth-new.tre.md5}
  ${TEMP}/itkTubeToTubeTransformFilter.tre
  MIDAS{Branch.n020.mha.md5}
  ${TEMP}/itkTubeToTubeTransformFilter.mha
  0.2 0.1 0.1 5 -5 5
  0
  )

# \todo Currently known to fail.  Re-enable when working.
if(0)
MIDAS_ADD_TEST( itkImageToTubeRigidRegistrationTest
  ${BASE_REGISTRATION_FILTERS_TESTS}
    --compare MIDAS{itkImageToTubeRigidRegistrationTest.mha.md5}
      ${TEMP}/itkImageToTubeRigidRegistrationTest.mha
    itkImageToTubeRigidRegistrationTest
    MIDAS{Branch.n020.mha.md5}
    MIDAS{tube.tre.md5}
    ${TEMP}/itkImageToTubeRigidRegistrationTest.tre
    ${TEMP}/itkImageToTubeRigidRegistrationTest.mha
    )
endif()

MIDAS_ADD_TEST( itkImageToTubeRigidRegistrationPerformanceTest
  ${BASE_REGISTRATION_FILTERS_TESTS}
  itkImageToTubeRigidRegistrationPerformanceTest
  MIDAS{Branch.n020.mha.md5}
  MIDAS{Branch-truth-new.tre.md5}
  ${TEMP}/itkImageToTubeRigidRegistrationPerformance
  )
set_property(TEST itkImageToTubeRigidRegistrationPerformanceTest
  APPEND PROPERTY DEPENDS itkTubeToTubeTransformFilterTest )

if( TubeTK_USE_VTK )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise
    ${BASE_REGISTRATION_FILTERS_TESTS}
    itkAnisotropicDiffusiveRegistrationRegularizationTest
    ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
    0 0
    5 0.125 1
    )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
    ${IMAGECOMPARE_EXE}
    -t ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
    -b MIDAS{Regularization_straight_noNoise_smoothedMotionField.mhd.md5}
    -i ${imageCompareTolerance}
    MIDAS_FETCH_ONLY{Regularization_straight_noNoise_smoothedMotionField.zraw.md5}
    )
  set_property(TEST
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
    APPEND PROPERTY DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise )

  MIDAS_ADD_TEST( itkAnisotropicDiffusiveRegistrationRegularizationTestStraight
    ${BASE_REGISTRATION_FILTERS_TESTS}
    itkAnisotropicDiffusiveRegistrationRegularizationTest
    ${TEMP}/Regularization_straight_smoothedMotionField.mhd
    0.1 0
    5 0.125 1
    )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
    ${IMAGECOMPARE_EXE}
    -t ${TEMP}/Regularization_straight_smoothedMotionField.mhd
    -b MIDAS{Regularization_straight_smoothedMotionField.mhd.md5}
    -i ${imageCompareTolerance}
    MIDAS_FETCH_ONLY{Regularization_straight_smoothedMotionField.zraw.md5}
    )
  set_property(TEST
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
    APPEND PROPERTY DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestStraight
    )

  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
    ${BASE_REGISTRATION_FILTERS_TESTS}
    itkAnisotropicDiffusiveRegistrationRegularizationTest
    ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
    0 0.5
    5 0.125 1
    )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
    ${IMAGECOMPARE_EXE}
    -t ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
    -b MIDAS{Regularization_angled_noNoise_smoothedMotionField.mhd.md5}
    -i ${imageCompareTolerance}
    MIDAS_FETCH_ONLY{Regularization_angled_noNoise_smoothedMotionField.zraw.md5}
    )
  set_property(TEST
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
    APPEND PROPERTY DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
    )

  MIDAS_ADD_TEST( itkAnisotropicDiffusiveRegistrationRegularizationTestAngled
    ${BASE_REGISTRATION_FILTERS_TESTS}
    itkAnisotropicDiffusiveRegistrationRegularizationTest
    ${TEMP}/Regularization_angled_smoothedMotionField.mhd
    0.1 0.5
    5 0.125 1
    )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
    ${IMAGECOMPARE_EXE}
    -t ${TEMP}/Regularization_angled_smoothedMotionField.mhd
    -b MIDAS{Regularization_angled_smoothedMotionField.mhd.md5}
    -i ${imageCompareTolerance}
    MIDAS_FETCH_ONLY{Regularization_angled_smoothedMotionField.zraw.md5}
    )
  set_property(TEST
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
    APPEND PROPERTY DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngled )

  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian
    ${BASE_REGISTRATION_FILTERS_TESTS}
    itkAnisotropicDiffusiveRegistrationRegularizationTest
    ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
    0.1 0.5
    5 0.125 0
    )
  MIDAS_ADD_TEST(
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
    ${IMAGECOMPARE_EXE}
    -t ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
    -b MIDAS{Regularization_angled_gaussian_smoothedMotionField.mhd.md5}
    -i ${imageCompareTolerance}
    MIDAS_FETCH_ONLY{Regularization_angled_gaussian_smoothedMotionField.zraw.md5}
          )
  set_property(TEST
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
    APPEND PROPERTY DEPENDS
    itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian )

endif( TubeTK_USE_VTK )

MIDAS_ADD_TEST( itkTubePointsToImageTest
  ${BASE_REGISTRATION_FILTERS_TESTS} itkTubePointsToImageTest
  MIDAS{Branch-truth-new.tre.md5}
  ${TEMP}/Branch-truth-new-points.mha
  )

MIDAS_ADD_TEST( itkImageToTubeRigidMetricTest1
  ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricTest
  MIDAS{SyntheticVesselTubeImage.mha.md5}
  MIDAS{SyntheticVesselTubeManuallyModified.tre.md5}
  0.874634
  )

MIDAS_ADD_TEST( itkImageToTubeRigidMetricTest2
  ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricTest
  MIDAS{SyntheticTransformedVesselTubeImage.mha.md5}
  MIDAS{SyntheticVesselTubeManuallyModified.tre.md5}
  -1.35168
  )

MIDAS_ADD_TEST( itkImageToTubeRigidMetricPerformanceTest
  ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricPerformanceTest
  MIDAS{SyntheticVesselTubeImage.mha.md5}
  MIDAS{SyntheticVesselTubeManuallyModified.tre.md5}
  ${TEMP}/itkImageToTubeRigidMetricPerformance.txt
  )

MIDAS_ADD_TEST( itkSyntheticTubeImagesGenerationsTest
  ${BASE_REGISTRATION_FILTERS_TESTS} itkSyntheticTubeImageGenerationsTest
  ${TEMP}/SyntheticBlurredTubeImage.mha
  ${TEMP}/SyntheticVesselTube.tre
  ${TEMP}/SyntheticVesselTubeImage.mha
  MIDAS{SyntheticVesselTubeManuallyModified.tre.md5}
  ${TEMP}/SyntheticTransformedVesselTubeImage.mha
  )

MIDAS_ADD_TEST( itkTubeExponentialResolutionWeightFunctionTest
  ${BASE_REGISTRATION_FILTERS_TESTS}
  itkTubeExponentialResolutionWeightFunctionTest
  ${TEMP}/itkTubeExponentialResolutionWeightFunctionTest.csv
  )

MIDAS_ADD_TEST( itkTubeParametricExponentialResolutionWeightFunctionTest
  ${BASE_REGISTRATION_FILTERS_TESTS}
  itkTubeParametricExponentialResolutionWeightFunctionTest
  ${TEMP}/itkTubeParametricExponentialResolutionWeightFunctionTest.csv
  )

MIDAS_ADD_TEST( itkTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
  ${BASE_REGISTRATION_FILTERS_TESTS}
  itkTubeParametricExponentialWithBoundsResolutionWeightFunctionTest
  ${TEMP}/itkTubeParametricExponentialWithBoundsResolutionWeightFunctionTest.csv
  )
