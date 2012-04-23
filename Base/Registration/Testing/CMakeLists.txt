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

#set(ITKRegistrationCommonTests

set( BASE_REGISTRATION_FILTERS_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseRegistrationFiltersTests )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Registration/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

# Use to Generate Generic Tests using same core
macro( GENERATE_TESTS BASE_NAME CLASS_NAMES TEST_EXTENSIONS )
MESSAGE( ${BASE_NAME} )
MESSAGE( ${CLASS_NAMES} )
MESSAGE( ${TEST_EXTENSIONS} )

  foreach( className ${CLASS_NAMES} )
    set( CTK_ADD_CLASS_NAME ${className} )

    foreach( extension ${TEST_EXTENSIONS} )
    set( CTK_ADD_TEST_EXTENSION ${extension} )
      configure_file(
          itk${BASE_NAME}${extension}.cxx.in
          ${CMAKE_CURRENT_BINARY_DIR}/itk${className}${extension}.cxx
          )

      list( APPEND GenericTubeBaseRegistrationFiltersTests_SRCS
            ${CMAKE_CURRENT_BINARY_DIR}/itk${className}${extension}.cxx
          )
    endforeach()

  endforeach()
endmacro()

# Use to generate the generic tests
set( MULTIPLE_TEST_EXTENTIONS
      Test
      PerformancesTest
     )

# Generate Tests for the Metrics
set( MULTIPLE_METRIC_CLASS_NAMES
  ImageToTubeRigidMetric
  ImageToTubeRigidMetric2
  )
set( MULTIPLE_REGISTRATION_CLASS_NAMES
  ImageToTubeRigidRegistration
  ImageToTubeRigidRegistration2
  )

set( GenericTubeBaseRegistrationFiltersTests_SRCS )
GENERATE_TESTS ( "ImageToTubeRigidMetric"
                 "${MULTIPLE_METRIC_CLASS_NAMES}"
                 "${MULTIPLE_TEST_EXTENTIONS}" )
GENERATE_TESTS ( "ImageToTubeRigidRegistration"
                 "${MULTIPLE_REGISTRATION_CLASS_NAMES}"
                 "${MULTIPLE_TEST_EXTENTIONS}" )

set( tubeBaseRegistrationFilters_SRCS
  itkSyntheticTubeImageGenerationsTest.cxx
  itkTubePointsToImageTest.cxx
  itkTubeToTubeTransformFilterTest.cxx
  )
list( APPEND tubeBaseRegistrationFilters_SRCS
      ${GenericTubeBaseRegistrationFiltersTests_SRCS}
   )

# expect image comparison perfection on 64-bit machines
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( imageCompareTolerance 0 )
# ... and give a bit of tolerance for 32-bit machines
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( imageCompareTolerance 0.000001 )
endif( CMAKE_SIZEOF_VOID_P EQUAL 8)

if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
  include( ${VTK_USE_FILE} )

  set( tubeBaseRegistrationFilters_SRCS
          ${tubeBaseRegistrationFilters_SRCS}
          itkAnisotropicDiffusiveRegistrationGenerateTestingImages.cxx
          itkAnisotropicDiffusiveRegistrationRegularizationTest.cxx
          )
endif( TubeTK_USE_VTK )

include_directories( 
          ${TubeTK_SOURCE_DIR}/Utilities
          ${TubeTK_SOURCE_DIR}/Base/Filtering
          ${TubeTK_SOURCE_DIR}/Base/Numerics
          ${TubeTK_SOURCE_DIR}/Base/Registration
          ${TubeTK_SOURCE_DIR}/Utilities
          )

add_executable( tubeBaseRegistrationFiltersTests
          tubeBaseRegistrationFiltersTests.cxx
          ${tubeBaseRegistrationFilters_SRCS}
          )
if( TubeTK_USE_VTK )
  add_definitions( -DTubeTK_USE_VTK )
endif()

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

add_test( itkImageToTubeRigidRegistrationPerformancesTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidRegistrationPerformancesTest
          ${INPUT_DATA}/Branch.n020.mha
          ${INPUT_DATA}/Branch-truth-new.tre
          ${TEMP}/itkImageToTubeRigidRegistrationPerformances
          )
set_tests_properties( itkImageToTubeRigidRegistrationPerformancesTest
                      PROPERTIES DEPENDS itkTubeToTubeTransformFilterTest )

add_test( itkImageToTubeRigidRegistration2Test
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidRegistration2Test
          ${INPUT_DATA}/Branch.n020.mha
          ${TEMP}/itkTubeToTubeTransformFilter.tre
          ${TEMP}/itkImageToTubeRigidRegistration2Test.mha
          )
set_tests_properties( itkImageToTubeRigidRegistration2Test
                      PROPERTIES DEPENDS itkTubeToTubeTransformFilterTest )

add_test( itkImageToTubeRigidRegistration2PerformancesTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidRegistration2PerformancesTest
          ${INPUT_DATA}/Branch.n020.mha
          ${INPUT_DATA}/Branch-truth-new.tre
          ${TEMP}/itkImageToTubeRigidRegistration2Performances
          )
set_tests_properties( itkImageToTubeRigidRegistration2PerformancesTest
                      PROPERTIES DEPENDS itkTubeToTubeTransformFilterTest )

if( TubeTK_USE_VTK )
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise
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
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_noNoise_smoothedMotionField.mhd
            -i ${imageCompareTolerance}
          )
  set_tests_properties( itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise-Compare
                        PROPERTIES DEPENDS
                        itkAnisotropicDiffusiveRegistrationRegularizationTestStraightNoNoise )

  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestStraight
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
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_straight_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_straight_smoothedMotionField.mhd
            -i ${imageCompareTolerance}
          )
  set_tests_properties( itkAnisotropicDiffusiveRegistrationRegularizationTestStraight-Compare
                        PROPERTIES DEPENDS
                        itkAnisotropicDiffusiveRegistrationRegularizationTestStraight )

  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise
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
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_noNoise_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_noNoise_smoothedMotionField.mhd
            -i ${imageCompareTolerance}
          )
  set_tests_properties( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise-Compare
                        PROPERTIES DEPENDS
                        itkAnisotropicDiffusiveRegistrationRegularizationTestAngledNoNoise )

  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngled
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
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_smoothedMotionField.mhd
            -i ${imageCompareTolerance}
          )
  set_tests_properties( itkAnisotropicDiffusiveRegistrationRegularizationTestAngled-Compare
                        PROPERTIES DEPENDS
                        itkAnisotropicDiffusiveRegistrationRegularizationTestAngled )

  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian
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
  add_test( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
          ${IMAGECOMPARE_EXE}
            -t ${TEMP}/Regularization_angled_gaussian_smoothedMotionField.mhd
            -b ${BASELINE}/Regularization_angled_gaussian_smoothedMotionField.mhd
            -i ${imageCompareTolerance}
          )
  set_tests_properties( itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian-Compare
                        PROPERTIES DEPENDS
                        itkAnisotropicDiffusiveRegistrationRegularizationTestAngledGaussian )

endif( TubeTK_USE_VTK )

add_test( itkTubePointsToImageTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkTubePointsToImageTest
          ${INPUT_DATA}/Branch-truth-new.tre
          ${TEMP}/Branch-truth-new-points.mha
          )

add_test( itkImageToTubeRigidMetricTest1
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricTest
          ${INPUT_DATA}/SyntheticVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          1.095681
          )

add_test( itkImageToTubeRigidMetricTest2
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricTest
          ${INPUT_DATA}/SyntheticTransformedVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          -1.13622
          )

add_test( itkImageToTubeRigidMetric2Test1
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetric2Test
          ${INPUT_DATA}/SyntheticVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          1.095681
          )

add_test( itkImageToTubeRigidMetric2Test2
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetric2Test
          ${INPUT_DATA}/SyntheticTransformedVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          -1.13622
          )

add_test( itkImageToTubeRigidMetricPerformancesTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetricPerformancesTest
          ${INPUT_DATA}/SyntheticVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          ${TEMP}/itkImageToTubeRigidMetricPerformance.txt
          )

add_test( itkImageToTubeRigidMetric2PerformancesTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkImageToTubeRigidMetric2PerformancesTest
          ${INPUT_DATA}/SyntheticVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          ${TEMP}/itkImageToTubeRigidMetric2Performance.txt
          )

add_test( itkSyntheticTubeImagesGenerationsTest
          ${BASE_REGISTRATION_FILTERS_TESTS} itkSyntheticTubeImageGenerationsTest
          ${TEMP}/SyntheticBlurredTubeImage.mha
          ${TEMP}/SyntheticVesselTube.tre
          ${TEMP}/SyntheticVesselTubeImage.mha
          ${INPUT_DATA}/SyntheticVesselTubeManuallyModified.tre
          ${TEMP}/SyntheticTransformedVesselTubeImage.mha
          )
