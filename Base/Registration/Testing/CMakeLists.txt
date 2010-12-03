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
          itkAnisotropicDiffusiveRegistrationExecution.cxx
          itkAnisotropicDiffusiveRegistrationImageRegistrationTest.cxx
          itkAnisotropicDiffusiveRegistrationRegularizationTest.cxx
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
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_noNoise_organBoundary.vtk
          0 0
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_noNoise_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise )

  add_test( AnisotropicDiffusiveRegistrationRegularizationTestStraight
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationRegularizationTest
          ${TEMP}/Regularization_straight_origMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_straight_normalVectorImage.mhd
          ${TEMP}/Regularization_straight_smoothedMotionField_organBoundary.vtk
          0.1 0
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationRegularizationTestStraight )

  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_noNoise_origMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_noNoise_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_noNoise_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_noNoise_organBoundary.vtk
          0 0.5
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_noNoise_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise )

  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngled
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_origMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          5 0.125 1
          )
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationRegularizationTestAngled )

  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationRegularizationTest
          ${TEMP}/Regularization_angled_gaussian_origMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_normal.mhd
          ${TEMP}/Regularization_angled_gaussian_normalVectorImage.mhd
          ${TEMP}/Regularization_angled_gaussian_smoothedMotionField_organBoundary.vtk
          0.1 0.5
          5 0.125 0
          )
  add_test( AnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_gaussian_smoothedMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian )

  add_test( AnisotropicDiffusiveRegistrationBoxesTestNoReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationImageRegistrationTest
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
  add_test( AnisotropicDiffusiveRegistrationBoxesTestNoReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_none_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_none_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationBoxesTestNoReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationBoxesTestNoReg )

  add_test( AnisotropicDiffusiveRegistrationBoxesTestAnisotropicReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationImageRegistrationTest
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
  add_test( AnisotropicDiffusiveRegistrationBoxesTestAnisotropicReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_anisotropic_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationBoxesTestAnisotropicReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationBoxesTestAnisotropicReg )

  add_test( AnisotropicDiffusiveRegistrationBoxesTestGaussianReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationImageRegistrationTest
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
  add_test( AnisotropicDiffusiveRegistrationBoxesTestGaussianReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveBoxesRegistration_gaussian_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationBoxesTestGaussianReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationBoxesTestGaussianReg )

  add_test( AnisotropicDiffusiveRegistrationSphereTestNoReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_transformedMovingImage.mhd
          5 0.125 0 1
          )
  add_test( AnisotropicDiffusiveRegistrationSphereTestNoReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationSphereTestNoReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationSphereTestNoReg )

  add_test( AnisotropicDiffusiveRegistrationSphereTestAnisotropicReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_transformedMovingImage.mhd
          5 0.125 1 1
          )
  add_test( AnisotropicDiffusiveRegistrationSphereTestAnisotropicReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationSphereTestAnisotropicReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationSphereTestAnisotropicReg )

  add_test( AnisotropicDiffusiveRegistrationSphereTestGaussianReg
          ${BASE_REGISTRATION_FILTERS_TESTS}
            itkAnisotropicDiffusiveRegistrationExecution
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
          ${INPUT_DATA}/AnisotropicDiffusiveSphereRegistration_weightImage.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
          ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_transformedMovingImage.mhd
          5 0.125 1 0
          )
  add_test( AnisotropicDiffusiveRegistrationSphereTestGaussianReg-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
            -b ${BASELINE}/AnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
          )
  set_tests_properties( AnisotropicDiffusiveRegistrationSphereTestGaussianReg-Compare
                        PROPERTIES DEPENDS
                        AnisotropicDiffusiveRegistrationSphereTestGaussianReg )

endif( TubeTK_USE_VTK )


