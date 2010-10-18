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

set( BASE_REGISTRATION_FILTERS_TESTS
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseRegistrationFiltersTests )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Registration/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseRegistrationFilters_SRCS
       itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest.cxx
       itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest.cxx
       itkImageToTubeRigidRegistrationTest.cxx
       itkTubeToTubeTransformFilterTest.cxx
)

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Preprocessing
  ${TubeTK_SOURCE_DIR}/Base/Registration 
  )

add_executable( tubeBaseRegistrationFiltersTests
                tubeBaseRegistrationFiltersTests.cxx
                ${tubeBaseRegistrationFilters_SRCS} )

target_link_libraries( tubeBaseRegistrationFiltersTests
                       ITKIO ITKNumerics ITKBasicFilters )

add_test( itkTubeToTubeTransformFilterTest 
          ${BASE_REGISTRATION_FILTERS_TESTS} itkTubeToTubeTransformFilterTest
          ${INPUT_DATA}/Branch-truth-new.tre
          ${TEMP}/itkTubeToTubeTransformFilter.tre
          ${INPUT_DATA}/Branch.n020.mha
          ${TEMP}/itkTubeToTubeTransformFilter.mha
          0.2 0.1 0.1 5 -5 5
          0
          )

add_test( itkImageToTubeRigidRegistrationTest 
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidRegistrationTest
          ${INPUT_DATA}/Branch.n020.mha
          ${TEMP}/itkTubeToTubeTransformFilter.tre
          ${TEMP}/itkImageToTubeRigidRegistrationTest.mha 
          )

set_tests_properties( itkImageToTubeRigidRegistrationTest
                      PROPERTIES DEPENDS itkTubeToTubeTransformFilterTest )

add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestStraightNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionField_straight_noNoise.mhd
          0
          ${TEMP}/SmoothedMotionField_straight_noNoise.mhd
          20
          0
          )

add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestStraight
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionField_straight.mhd
          0.1
          ${TEMP}/SmoothedMotionField_straight.mhd
          20
          0
          )

add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestAngledNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionField_angled_noNoise.mhd
          0
          ${TEMP}/SmoothedMotionField_angled_noNoise.mhd
          20
          0.5
          )

add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestAngled
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionField_angled.mhd
          0.1
          ${TEMP}/SmoothedMotionField_angled.mhd
          20
          0.5
          )

# TODO change number of iterations to something larger
add_test( itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTestWithoutRegularization
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_originalFixedImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_originalMovingImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_resultingMotionField.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_resultingTransformedMovingImage.mhd
          10
          0
          )

add_test( itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTestWithRegularization
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/SlidingOrganRegistrationWithRegularization_originalFixedImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_originalMovingImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_resultingMotionField.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_resultingTransformedMovingImage.mhd
          10
          1
          )



