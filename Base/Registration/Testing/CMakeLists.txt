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
          itkImageToImageDiffusiveDeformableRegistrationExecution.cxx
          itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest.cxx
          itkImageToImageDiffusiveDeformableRegistrationRegularizationTest.cxx
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

if( TubeTK_USE_VTK )
  add_test( itkImageToImageDiffusiveDeformableRegistrationRegularizationTest1
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField_tangential.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_noNoise_organBoundary.vtk
          0 0
          20 0.05 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationRegularizationTest2
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_origMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_tangential.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_organBoundary.vtk
          0.1 0
          20 0.05 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationRegularizationTest3
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField_tangential.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_noNoise_organBoundary.vtk
          0 0.5
          20 0.05 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationRegularizationTest4
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_origMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_tangential.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          20 0.05 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationRegularizationTest5
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_gaussian_origMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_tangential.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_gaussian_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          20 0.05 0
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationBoxesTest1
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/DiffusiveBoxesRegistration_none_origFixedImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_origMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_resultingMotionField.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_transformedMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_weightImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_normalVectorImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_none_normalSurfaceBorder.vtk
          1
          1 0.125 0 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationBoxesTest2
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_origFixedImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_origMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_transformedMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_weightImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_normalVectorImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_anisotropic_normalSurfaceBorder.vtk
          1
          1 0.125 1 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationBoxesTest3
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_origFixedImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_origMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_transformedMovingImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_weightImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_normalVectorImage.mhd
          ${TEMP}/DiffusiveBoxesRegistration_gaussian_normalSurfaceBorder.vtk
          1
          1 0.125 1 0
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationSphereTest1
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/DiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_fixedSurfaceBorder.vtk
          ${TEMP}/DiffusiveSphereRegistration_none_resultingMotionField.mhd
          ${TEMP}/DiffusiveSphereRegistration_none_transformedMovingImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_none_weightImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_none_normalVectorImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_none_normalSurfaceBorder.vtk
          1 1.0 0 1
          )

  add_test( itkImageToImageDiffusiveDeformableRegistrationSphereTest2
          ${BASE_REGISTRATION_FILTERS_TESTS}
                              itkImageToImageDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/DiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_fixedSurfaceBorder.vtk
          ${TEMP}/DiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
          ${TEMP}/DiffusiveSphereRegistration_anisotropic_transformedMovingImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_anisotropic_weightImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_anisotropic_normalVectorImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_anisotropic_normalSurfaceBorder.vtk
          1 0.125 1 1
          )

  add_test(  itkImageToImageDiffusiveDeformableRegistrationSphereTest3
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/DiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/DiffusiveSphereRegistration_fixedSurfaceBorder.vtk
          ${TEMP}/DiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
          ${TEMP}/DiffusiveSphereRegistration_gaussian_transformedMovingImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_gaussian_weightImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_gaussian_normalVectorImage.mhd
          ${TEMP}/DiffusiveSphereRegistration_gaussian_normalSurfaceBorder.vtk
          1 0.05 1 0
          )

endif( TubeTK_USE_VTK )


