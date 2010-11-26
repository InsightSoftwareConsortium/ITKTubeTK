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

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( BASE_REGISTRATION_FILTERS_TESTS
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseRegistrationFiltersTests )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Registration/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseRegistrationFilters_SRCS
          itkImageToTubeRigidRegistrationTest.cxx
          itkTubeToTubeTransformFilterTest.cxx
          )

if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )

  set( tubeBaseRegistrationFilters_SRCS
          ${tubeBaseRegistrationFilters_SRCS}
          itkImageToImageAnisotropicDiffusiveDeformableRegistrationExecution.cxx
          itkImageToImageAnisotropicDiffusiveDeformableRegistrationImageRegistrationTest.cxx
          itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest.cxx
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
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraightNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_noNoise_organBoundary.vtk
          0 0
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraightNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_noNoise_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraightNoNoise-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraightNoNoise )

  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraight
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_origMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_organBoundary.vtk
          0.1 0
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraight-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraight-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationRegularizationTestStraight )

  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_noNoise_organBoundary.vtk
          0 0.5
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_noNoise_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledNoNoise-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledNoNoise )

  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngled
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_origMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngled-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngled-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngled )

  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledGaussian
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_gaussian_origMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_gaussian_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          5 0.125 0
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledGaussian-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_gaussian_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledGaussian-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationRegularizationTestAngledGaussian )

  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestNoReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_origFixedImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_origMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_transformedMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_normalVectorImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_normalSurfaceBorder.vtk
          1
          5 0.125 0 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestNoReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_none_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationBoxesTestNoReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationBoxesTestNoReg )

  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestAnisotropicReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_origFixedImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_origMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_transformedMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_normalVectorImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_normalSurfaceBorder.vtk
          1
          5 0.125 1 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestAnisotropicReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationBoxesTestAnisotropicReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationBoxesTestAnisotropicReg )

  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestGaussianReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationImageRegistrationTest
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_origFixedImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_origMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_transformedMovingImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_normalVectorImage.mhd
          ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_normalSurfaceBorder.vtk
          1
          5 0.125 1 0
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationBoxesTestGaussianReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationBoxesTestGaussianReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationBoxesTestGaussianReg )

  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestNoReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_transformedMovingImage.mhd
          5 0.125 0 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestNoReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationSphereTestNoReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationSphereTestNoReg )

  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestAnisotropicReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_transformedMovingImage.mhd
          5 0.125 1 1
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestAnisotropicReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationSphereTestAnisotropicReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationSphereTestAnisotropicReg )

  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestGaussianReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkImageToImageAnisotropicDiffusiveDeformableRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_transformedMovingImage.mhd
          5 0.125 1 0
          )
  add_test( AnisotropicDiffusiveDeformableRegistrationSphereTestGaussianReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveDeformableRegistrationSphereTestGaussianReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveDeformableRegistrationSphereTestGaussianReg )

endif( TubeTK_USE_VTK )


