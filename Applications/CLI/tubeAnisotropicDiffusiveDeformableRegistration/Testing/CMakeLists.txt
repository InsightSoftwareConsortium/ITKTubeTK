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

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )

set( BASELINE
  ${TubeTK_SOURCE_DIR}/Applications/CLI/${PROJECT_NAME}/Testing/Baselines )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

# expect image comparison perfection on 64-bit machines
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( imageCompareTolerance 0 )
# ... and give a bit of tolerance for 32-bit machines
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( imageCompareTolerance 0.000001 )
endif( CMAKE_SIZEOF_VOID_P EQUAL 8)

#######
# Test1
#######
add_test( ${PROJECT_NAME}TestSphereNoReg
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -n ${INPUT_DATA}/Sphere_normals.mhd
               -w ${INPUT_DATA}/Sphere_weights.mhd
               -d ${TEMP}/Sphere_noReg_motionField.mhd
               -o ${TEMP}/Sphere_noReg_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1
               -r )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}TestSphereNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_noReg_motionField.mhd
               -b ${BASELINE}/Sphere_noReg_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereNoReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereNoReg )

#######
# Test2
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropic
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -n ${INPUT_DATA}/Sphere_normals.mhd
               -w ${INPUT_DATA}/Sphere_weights.mhd
               -d ${TEMP}/Sphere_anisotropic_motionField.mhd
               -o ${TEMP}/Sphere_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1 )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropic_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropic_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropic )

#######
# Test3
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRAS
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -b ${INPUT_DATA}/Sphere_surfaceBorderRAS.vtk
               -d ${TEMP}/Sphere_anisotropicSurfaceRAS_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRAS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05 )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRAS_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropicSurface_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRAS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropicSurfaceRAS )

#######
# Test4
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceLPS
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -b ${INPUT_DATA}/Sphere_surfaceBorderLPS.vtk
               -d ${TEMP}/Sphere_anisotropicSurfaceLPS_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceLPS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               -c LPS )

#######
# Test4-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceLPS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceLPS_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropicSurface_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropicSurfaceLPS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropicSurfaceLPS )

#######
# Test5
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASDownsampled
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixedDownsampled.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -b ${INPUT_DATA}/Sphere_surfaceBorderRAS.vtk
               -d ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05 )

#######
# Test5-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASDownsampled-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropicSurfaceDownsampled_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASDownsampled-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASDownsampled )

#######
# Test6
#######
add_test( ${PROJECT_NAME}TestSphereGaussian
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -n ${INPUT_DATA}/Sphere_normals.mhd
               -w ${INPUT_DATA}/Sphere_weights.mhd
               -d ${TEMP}/Sphere_gaussian_motionField.mhd
               -o ${TEMP}/Sphere_gaussian_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1
               -a )

#######
# Test6-Compare
#######
add_test( ${PROJECT_NAME}TestSphereGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_gaussian_motionField.mhd
               -b ${BASELINE}/Sphere_gaussian_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereGaussian-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereGaussian )

#######
# Test7
#######
add_test( ${PROJECT_NAME}TestBoxesNoReg
            ${PROJ_EXE}
               ${INPUT_DATA}/Boxes_fixed.mhd
               ${INPUT_DATA}/Boxes_moving.mhd
               -d ${TEMP}/Boxes_noReg_motionField.mhd
               -o ${TEMP}/Boxes_noReg_transformedMoving.mhd
               -i 5
               -s 0.125
               -r )

#######
# Test7-Compare
#######
add_test( ${PROJECT_NAME}TestBoxesNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_noReg_motionField.mhd
               -b ${BASELINE}/Boxes_noReg_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestBoxesNoReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestBoxesNoReg )

#######
# Test8
#######
add_test( ${PROJECT_NAME}TestBoxesAnisotropic
            ${PROJ_EXE}
               ${INPUT_DATA}/Boxes_fixed.mhd
               ${INPUT_DATA}/Boxes_moving.mhd
               -n ${INPUT_DATA}/Boxes_normals.mhd
               -w ${INPUT_DATA}/Boxes_weights.mhd
               -d ${TEMP}/Boxes_anisotropic_motionField.mhd
               -o ${TEMP}/Boxes_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125 )

#######
# Test8-Compare
#######
add_test( ${PROJECT_NAME}TestBoxesAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_anisotropic_motionField.mhd
               -b ${BASELINE}/Boxes_anisotropic_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestBoxesAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestBoxesAnisotropic )

#######
# Test9
#######
add_test( ${PROJECT_NAME}TestBoxesGaussian
            ${PROJ_EXE}
               ${INPUT_DATA}/Boxes_fixed.mhd
               ${INPUT_DATA}/Boxes_moving.mhd
               -n ${INPUT_DATA}/Boxes_normals.mhd
               -w ${INPUT_DATA}/Boxes_weights.mhd
               -d ${TEMP}/Boxes_gaussian_motionField.mhd
               -o ${TEMP}/Boxes_gaussian_transformedMoving.mhd
               -i 5
               -s 0.125
               -a )

#######
# Test9-Compare
#######
add_test( ${PROJECT_NAME}TestBoxesGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_gaussian_motionField.mhd
               -b ${BASELINE}/Boxes_gaussian_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestBoxesGaussian-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestBoxesGaussian )

#######
# Test10
#######
add_test( ${PROJECT_NAME}TestSphereSparseAnisotropicSurfaceRAS
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -b ${INPUT_DATA}/Sphere_surfaceBorderRAS.vtk
               -d ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_motionField.mhd
               -o ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               -u SparseSlidingOrgan )

#######
# Test10-Compare
#######
add_test( ${PROJECT_NAME}TestSphereSparseAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_motionField.mhd
               -b ${BASELINE}/Sphere_sparseAnisotropicSurface_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereSparseAnisotropicSurfaceRAS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereSparseAnisotropicSurfaceRAS )

#######
# Test11
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicMultiresolution
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -n ${INPUT_DATA}/Sphere_normalsMultiresolution.mhd
               -w ${INPUT_DATA}/Sphere_weightsMultiresolution.mhd
               -d ${TEMP}/Sphere_anisotropicMultiresolution_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicMultiresolution_transformedMoving.mhd
               -i 5,5
               -s 0.125
               -l 0.05 )

#######
# Test11-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicMultiresolution_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropicMultiresolution_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropicMultiresolution-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropicMultiresolution )

#######
# Test12
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASMultiresolution
            ${PROJ_EXE}
               ${INPUT_DATA}/Sphere_fixed.mhd
               ${INPUT_DATA}/Sphere_moving.mhd
               -b ${INPUT_DATA}/Sphere_surfaceBorderRAS.vtk
               -d ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_transformedMoving.mhd
               -i 5,5
               -s 0.125
               -l 0.05 )

#######
# Test12-Compare
#######
add_test( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_motionField.mhd
               -b ${BASELINE}/Sphere_anisotropicMultiresolution_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASMultiresolution-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestSphereAnisotropicSurfaceRASMultiresolution )

#######
# Test13
#######
add_test( ${PROJECT_NAME}TestTubesSparseAnisotropic
            ${PROJ_EXE}
               ${INPUT_DATA}/Tubes_fixed.mhd
               ${INPUT_DATA}/Tubes_moving.mhd
               -p ${INPUT_DATA}/Tubes_spatialObjects.tre
               -d ${TEMP}/Tubes_anisotropic_motionField.mhd
               -o ${TEMP}/Tubes_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.25
               -u SparseSlidingOrgan
               -v ${TEMP}/Tubes_normalVectors.mhd
               -x ${TEMP}/Tubes_weightRegularizations.mhd
               -z ${TEMP}/Tubes_weightStructures.mhd )

#######
# Test13-Compare
#######
add_test( ${PROJECT_NAME}TestTubesSparseAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Tubes_anisotropic_motionField.mhd
               -b ${BASELINE}/Tubes_anisotropic_motionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}TestTubesSparseAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestTubesSparseAnisotropic )

#######
# MIDAS-Test1
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereNoReg
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/Sphere_noReg_motionField.mhd
               -o ${TEMP}/Sphere_noReg_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1
               -r
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.raw.md5}
        )

#######
# MIDAS-Test1-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_noReg_motionField.mhd
               -b MIDAS{Sphere_noReg_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_noReg_motionField.zraw.md5}
            )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereNoReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereNoReg )

#######
# MIDAS-Test2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropic
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/Sphere_anisotropic_motionField.mhd
               -o ${TEMP}/Sphere_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.raw.md5}
            )

#######
# MIDAS-Test2-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropic_motionField.mhd
               -b MIDAS{Sphere_anisotropic_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropic_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropic )

#######
# MIDAS-Test3
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRAS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/Sphere_anisotropicSurfaceRAS_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRAS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
              )

#######
# MIDAS-Test3-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRAS_motionField.mhd
               -b MIDAS{Sphere_anisotropicSurface_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropicSurface_motionField.zraw.md5}
              )

set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRAS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRAS )

#######
# MIDAS-Test4
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceLPS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderLPS.vtk.md5}
               -d ${TEMP}/Sphere_anisotropicSurfaceLPS_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceLPS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               -c LPS
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
              )

#######
# MIDAS-Test4-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceLPS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceLPS_motionField.mhd
               -b MIDAS{Sphere_anisotropicSurface_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropicSurface_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceLPS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceLPS )

#######
# MIDAS-Test5
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASDownsampled
            ${PROJ_EXE}
               MIDAS{Sphere_fixedDownsampled.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixedDownsampled.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
            )

#######
# MIDAS-Test5-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASDownsampled-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRASDownsampled_motionField.mhd
               -b MIDAS{Sphere_anisotropicSurfaceDownsampled_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropicSurfaceDownsampled_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASDownsampled-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASDownsampled )

#######
# MIDAS-Test6
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereGaussian
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/Sphere_gaussian_motionField.mhd
               -o ${TEMP}/Sphere_gaussian_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.1
               -a
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.raw.md5}
            )

#######
# MIDAS-Test6-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_gaussian_motionField.mhd
               -b MIDAS{Sphere_gaussian_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_gaussian_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereGaussian-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereGaussian )

#######
# MIDAS-Test7
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesNoReg
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -d ${TEMP}/Boxes_noReg_motionField.mhd
               -o ${TEMP}/Boxes_noReg_transformedMoving.mhd
               -i 5
               -s 0.125
               -r
               MIDAS_FETCH_ONLY{Boxes_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.raw.md5}
            )

#######
# MIDAS-Test7-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_noReg_motionField.mhd
               -b MIDAS{Boxes_noReg_motionField.mhd.md5}
               -i ${imageCompareTolerance}
                MIDAS_FETCH_ONLY{Boxes_noReg_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestBoxesNoReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}TestBoxesNoReg )

#######
# MIDAS-Test8
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesAnisotropic
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -n MIDAS{Boxes_normals.mhd.md5}
               -w MIDAS{Boxes_weights.mhd.md5}
               -d ${TEMP}/Boxes_anisotropic_motionField.mhd
               -o ${TEMP}/Boxes_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125
               MIDAS_FETCH_ONLY{Boxes_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_normals.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_weights.raw.md5}
            )

#######
# MIDAS-Test8-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_anisotropic_motionField.mhd
               -b MIDAS{Boxes_anisotropic_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLLY{Boxes_anisotropic_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestBoxesAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestBoxesAnisotropic )

#######
# MIDAS-Test9
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesGaussian
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -n MIDAS{Boxes_normals.mhd.md5}
               -w MIDAS{Boxes_weights.mhd.md5}
               -d ${TEMP}/Boxes_gaussian_motionField.mhd
               -o ${TEMP}/Boxes_gaussian_transformedMoving.mhd
               -i 5
               -s 0.125
               -a
               MIDAS_FETCH_ONLY{Boxes_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_normals.raw.md5}
               MIDAS_FETCH_ONLY{Boxes_weights.raw.md5}
             )

#######
# MIDAS-Test9-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestBoxesGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Boxes_gaussian_motionField.mhd
               -b MIDAS{Boxes_gaussian_motionField.mhd.md5}
               -i ${imageCompareTolerance}
                MIDAS_FETCH_ONLY{Boxes_gaussian_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestBoxesGaussian-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestBoxesGaussian )

#######
# MIDAS-Test10
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereSparseAnisotropicSurfaceRAS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_motionField.mhd
               -o ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.05
               -u SparseSlidingOrgan
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
            )

#######
# MIDAS-Test10-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereSparseAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_sparseAnisotropicSurfaceRAS_motionField.mhd
               -b MIDAS{Sphere_sparseAnisotropicSurface_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_sparseAnisotropicSurface_motionField.zraw.md5}
             )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereSparseAnisotropicSurfaceRAS-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereSparseAnisotropicSurfaceRAS )

#######
# MIDAS-Test11
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicMultiresolution
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normalsMultiresolution.mhd.md5}
               -w MIDAS{Sphere_weightsMultiresolution.mhd.md5}
               -d ${TEMP}/Sphere_anisotropicMultiresolution_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicMultiresolution_transformedMoving.mhd
               -i 5,5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_normalsMultiresolution.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_weightsMultiresolution.raw.md5}
              )

#######
# MIDAS-Test11-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicMultiresolution_motionField.mhd
               -b MIDAS{Sphere_anisotropicMultiresolution_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropicMultiresolution_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicMultiresolution-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropicMultiresolution )

#######
# MIDAS-Test12
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASMultiresolution
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_motionField.mhd
               -o ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_transformedMoving.mhd
               -i 5,5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.raw.md5}
              )

#######
# MIDAS-Test12-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Sphere_anisotropicSurfaceRASMultiresolution_motionField.mhd
               -b MIDAS{Sphere_anisotropicMultiresolution_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Sphere_anisotropicMultiresolution_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASMultiresolution-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestSphereAnisotropicSurfaceRASMultiresolution )

#######
# MIDAS-Test13
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestTubesSparseAnisotropic
            ${PROJ_EXE}
               MIDAS{Tubes_fixed.mhd.md5}
               MIDAS{Tubes_moving.mhd.md5}
               -p MIDAS{Tubes_spatialObjects.tre.md5}
               -d ${TEMP}/Tubes_anisotropic_motionField.mhd
               -o ${TEMP}/Tubes_anisotropic_transformedMoving.mhd
               -i 5
               -s 0.125
               -l 0.25
               -u SparseSlidingOrgan
               -v ${TEMP}/Tubes_normalVectors.mhd
               -x ${TEMP}/Tubes_weightRegularizations.mhd
               -z ${TEMP}/Tubes_weightStructures.mhd
               MIDAS_FETCH_ONLY{Tubes_fixed.raw.md5}
               MIDAS_FETCH_ONLY{Tubes_moving.raw.md5}
              )

#######
# MIDAS-Test13-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-TestTubesSparseAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/Tubes_anisotropic_motionField.mhd
               -b MIDAS{Tubes_anisotropic_motionField.mhd.md5}
               -i ${imageCompareTolerance}
               MIDAS_FETCH_ONLY{Tubes_anisotropic_motionField.zraw.md5}
              )
set_tests_properties( ${PROJECT_NAME}MIDAS-TestTubesSparseAnisotropic-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-TestTubesSparseAnisotropic )
