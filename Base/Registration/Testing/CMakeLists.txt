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
          itkImageToTubeRigidRegistrationTest.cxx
          itkTubeToTubeTransformFilterTest.cxx
          )
if( TubeTK_USE_VTK )
  set( tubeBaseRegistrationFilters_SRCS
          ${tubeBaseRegistrationFilters_SRCS}
          itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest.cxx
          itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest.cxx
          )
endif( TubeTK_USE_VTK )

include_directories( 
          ${TubeTK_SOURCE_DIR}/Utilities
          ${TubeTK_SOURCE_DIR}/Base/Preprocessing
          ${TubeTK_SOURCE_DIR}/Base/Registration
          )

add_executable( tubeBaseRegistrationFiltersTests
          tubeBaseRegistrationFiltersTests.cxx
          ${tubeBaseRegistrationFilters_SRCS}
          )

set( tubeBaseRegistrationFilters_LIBS
          tubeBaseRegistrationFiltersTests
          ITKIO ITKNumerics ITKBasicFilters
          )
if( TubeTK_USE_VTK )
  set( tubeBaseRegistrationFilters_LIBS
          ${tubeBaseRegistrationFilters_LIBS}
          vtkCommon vtkGraphics vtkIO
          )
endif( TubeTK_USE_VTK )
target_link_libraries( ${tubeBaseRegistrationFilters_LIBS} )

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

# TODO for ALL of these, make number of iterations greater
if( TubeTK_USE_VTK )
  add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestStraightNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionFieldRegularization_straight_noNoise_origMotionField.mhd
          0
          ${TEMP}/MotionFieldRegularization_straight_noNoise_smoothedMotionField.mhd
          20
          0
          ${TEMP}/MotionFieldRegularization_straight_noNoise_normalVectorImage.mhd
          0.05
          1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestStraight
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionFieldRegularization_straight_origMotionField.mhd
          0.1
          ${TEMP}/MotionFieldRegularization_straight_smoothedMotionField.mhd
          20
          0
          ${TEMP}/MotionFieldRegularization_straight_normalVectorImage.mhd
          0.05
          1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestAngledNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionFieldRegularization_angled_noNoise_origMotionField.mhd
          0
          ${TEMP}/MotionFieldRegularization_angled_noNoise_smoothedMotionField.mhd
          20
          0.5
          ${TEMP}/MotionFieldRegularization_angled_noNoise_normalVectorImage.mhd
          0.05
          1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestAngled
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionFieldRegularization_angled_origMotionField.mhd
          0.1
          ${TEMP}/MotionFieldRegularization_angled_smoothedMotionField.mhd
          20
          0.5
          ${TEMP}/MotionFieldRegularization_angled_normalVectorImage.mhd
          0.05
          1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTestAngled_w0
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationMotionFieldRegularizationTest
          ${TEMP}/MotionFieldRegularization_angled_w0_origMotionField.mhd
          0.1
          ${TEMP}/MotionFieldRegularization_angled_w0_smoothedMotionField.mhd
          20
          0.5
          ${TEMP}/MotionFieldRegularization_angled_w0_normalVectorImage.mhd
          0.05
          0
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTestWithoutRegularization
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_originalFixedImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_originalMovingImage.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_resultingMotionField.mhd
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_resultingTransformedMovingImage.mhd
          4
          0
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_normalSurfaceBorder.vtk
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_normalVectorImage.mhd
          1
          ${TEMP}/SlidingOrganRegistrationWithoutRegularization_weightImage.mhd
          1.0
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
          ${TEMP}/SlidingOrganRegistrationWithRegularization_normalSurfaceBorder.vtk
          ${TEMP}/SlidingOrganRegistrationWithRegularization_normalVectorImage.mhd
          1
          ${TEMP}/SlidingOrganRegistrationWithRegularization_weightImage.mhd
          1.0
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTestWithRegularization_w0
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/SlidingOrganRegistrationWithRegularization_originalFixedImage_w0.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_originalMovingImage_w0.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_resultingMotionField_w0.mhd
          ${TEMP}/SlidingOrganRegistrationWithRegularization_resultingTransformedMovingImage_w0.mhd
          10
          1
          ${TEMP}/SlidingOrganRegistrationWithRegularization_normalSurfaceBorder_w0.vtk
          ${TEMP}/SlidingOrganRegistrationWithRegularization_normalVectorImage_w0.mhd
          0
          ${TEMP}/SlidingOrganRegistrationWithRegularization_weightImage_w0.mhd
          0.05
          )

endif( TubeTK_USE_VTK )


