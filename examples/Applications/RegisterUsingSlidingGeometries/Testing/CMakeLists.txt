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

include_regular_expression( "^.*$" )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareImagesTolerance 0.00001 )

# Test1
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_NoReg_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -r )

# Test1-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_NoReg_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_noReg_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Sphere_NoReg_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Sphere_NoReg_motionField )

# Test2
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropic_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.1 )

# Test2-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropic_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Sphere_anisotropic_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Sphere_anisotropic_motionField )

# Test3
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05 )

# Test3-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha}
               -i ${CompareImagesTolerance} )

set_tests_properties(
 ${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField-Compare
 PROPERTIES DEPENDS ${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField )

# Test4
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderLPS.vtk}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -c LPS )

# Test4-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
  ${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField )

# Test5
itk_add_test(
            NAME
${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixedDownsampled.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_movingDownsampled.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -i 5
               -s 0.125
               -l 0.05 )

# Test5-Compare
itk_add_test(
            NAME
${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
  ${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField-Compare
  PROPERTIES DEPENDS
  ${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField )

# Test6
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_guassian_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -a )

# Test6-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_guassian_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_gaussian_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
 ${MODULE_NAME}-Sphere_guassian_motionField-Compare
 PROPERTIES DEPENDS ${MODULE_NAME}-Sphere_guassian_motionField )

# Test7
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_noReg_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -i 5
               -s 0.125
               -r )

# Test7-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_noReg_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_noReg_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Boxes_noReg_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Boxes_noReg_motionField )

# Test8
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_anisotropic_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Boxes_normals.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Boxes_weights.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -i 5
               -s 0.125 )

# Test8-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_anisotropic_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
 ${MODULE_NAME}-Boxes_anisotropic_motionField-Compare
 PROPERTIES DEPENDS ${MODULE_NAME}-Boxes_anisotropic_motionField )

# Test9
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_guassian_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Boxes_normals.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Boxes_weights.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -i 5
               -s 0.125
               -a )

# Test9-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Boxes_guassian_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_gaussian_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Boxes_guassian_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Boxes_guassian_motionField )

# Test10
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -u SparseSlidingOrgan )

# Test10-Compare
itk_add_test( NAME
${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
  ${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField-Compare
  PROPERTIES DEPENDS
  ${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField )

# Test11
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normalsMultiresolution.mha}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weightsMultiresolution.mha}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05 )

# Test11-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
  ${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField-Compare
  PROPERTIES DEPENDS
  ${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField )

# Test12
itk_add_test(
            NAME
               ${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05 )

# Test12-Compare
itk_add_test(
            NAME
              ${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties(
  ${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField-Compare
  PROPERTIES DEPENDS
  ${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField )

# Test13
itk_add_test(
            NAME ${MODULE_NAME}-Tubes_anisotropic_motionField
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Tubes_fixed.mha}
               DATA{${TubeTK_DATA_ROOT}/Tubes_moving.mha}
               -p DATA{${TubeTK_DATA_ROOT}/Tubes_spatialObjects.tre}
               -d ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.25
               -u SparseSlidingOrgan )

# Test13-Compare
itk_add_test(
            NAME ${MODULE_NAME}-Tubes_anisotropic_motionField-Compare
            COMMAND ${TubeTK_CompareImages_EXE}
              CompareImages
               -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-Tubes_anisotropic_motionField-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Tubes_anisotropic_motionField )
