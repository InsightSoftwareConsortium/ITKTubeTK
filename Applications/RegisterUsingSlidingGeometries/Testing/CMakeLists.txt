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


set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${SlicerExecutionModel_CLI_RUNTIME_OUTPUT_DIRECTORY}/${MODULE_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( imageCompareTolerance 0.000001 )

#######
# Test1
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereNoReg
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -r
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.zraw.md5} )

#######
# Test1-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_noReg_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereNoReg-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereNoReg )

#######
# Test2
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropic
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.zraw.md5} )

#######
# Test2-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropic_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereAnisotropic-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropic )

#######
# Test3
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5} )

#######
# Test3-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropicSurface_motionField.mha.md5}
               -i ${imageCompareTolerance} )

set_property(TEST ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS )

#######
# Test4
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderLPS.vtk.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -c LPS
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5} )

#######
# Test4-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropicSurface_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS )

#######
# Test5
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled
            ${PROJ_EXE}
               MIDAS{Sphere_fixedDownsampled.mha.md5}
               MIDAS{Sphere_movingDownsampled.mha.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -i 5
               -s 0.125
               -l 0.05 )

#######
# Test5-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled )

#######
# Test6
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereGaussian
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normals.mhd.md5}
               -w MIDAS{Sphere_weights.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -a
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_normals.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_weights.zraw.md5} )

#######
# Test6-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_gaussian_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereGaussian-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereGaussian )

#######
# Test7
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesNoReg
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -i 5
               -s 0.125
               -r
               MIDAS_FETCH_ONLY{Boxes_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.zraw.md5} )

#######
# Test7-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -b MIDAS{${MODULE_NAME}-Boxes_noReg_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestBoxesNoReg-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestBoxesNoReg )

#######
# Test8
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesAnisotropic
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -n MIDAS{Boxes_normals.mhd.md5}
               -w MIDAS{Boxes_weights.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -i 5
               -s 0.125
               MIDAS_FETCH_ONLY{Boxes_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_normals.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_weights.zraw.md5} )

#######
# Test8-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -b MIDAS{${MODULE_NAME}-Boxes_anisotropic_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestBoxesAnisotropic-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestBoxesAnisotropic )

#######
# Test9
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesGaussian
            ${PROJ_EXE}
               MIDAS{Boxes_fixed.mhd.md5}
               MIDAS{Boxes_moving.mhd.md5}
               -n MIDAS{Boxes_normals.mhd.md5}
               -w MIDAS{Boxes_weights.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -i 5
               -s 0.125
               -a
               MIDAS_FETCH_ONLY{Boxes_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_normals.zraw.md5}
               MIDAS_FETCH_ONLY{Boxes_weights.zraw.md5} )

#######
# Test9-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestBoxesGaussian-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -b MIDAS{${MODULE_NAME}-Boxes_gaussian_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestBoxesGaussian-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestBoxesGaussian )

#######
# Test10
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -u SparseSlidingOrgan
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5} )

#######
# Test10-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_sparseAnisotropicSurface_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS )

#######
# Test11
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicMultiresolution
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -n MIDAS{Sphere_normalsMultiresolution.mhd.md5}
               -w MIDAS{Sphere_weightsMultiresolution.mhd.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_normalsMultiresolution.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_weightsMultiresolution.zraw.md5} )

#######
# Test11-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereAnisotropicMultiresolution-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropicMultiresolution )

#######
# Test12
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution
            ${PROJ_EXE}
               MIDAS{Sphere_fixed.mhd.md5}
               MIDAS{Sphere_moving.mhd.md5}
               -b MIDAS{Sphere_surfaceBorderRAS.vtk.md5}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05
               MIDAS_FETCH_ONLY{Sphere_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Sphere_moving.zraw.md5} )

#######
# Test12-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -b MIDAS{${MODULE_NAME}-Sphere_anisotropicSurfaceMultiresolution_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution )

#######
# Test13
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestTubesSparseAnisotropic
            ${PROJ_EXE}
               MIDAS{Tubes_fixed.mhd.md5}
               MIDAS{Tubes_moving.mhd.md5}
               -p MIDAS{Tubes_spatialObjects.tre.md5}
               -d ${TEMP}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.25
               -u SparseSlidingOrgan
               MIDAS_FETCH_ONLY{Tubes_fixed.zraw.md5}
               MIDAS_FETCH_ONLY{Tubes_moving.zraw.md5} )

#######
# Test13-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-TestTubesSparseAnisotropic-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -b MIDAS{${MODULE_NAME}-Tubes_anisotropic_motionField.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-TestTubesSparseAnisotropic-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-TestTubesSparseAnisotropic )
