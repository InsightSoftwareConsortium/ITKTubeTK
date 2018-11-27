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

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareImagesTolerance 0.00001 )

# Test1
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereNoReg
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mhd,Sphere_normals.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mhd,Sphere_weights.zraw}
               -d ${TEMP}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -r )

# Test1-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereNoReg-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_noReg_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_noReg_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereNoReg-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereNoReg )

# Test2
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropic
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mhd,Sphere_normals.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mhd,Sphere_weights.zraw}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.1 )

# Test2-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropic-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropic-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropic )

# Test3
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05 )

# Test3-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRAS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurface_motionField.mha}
               -i ${CompareImagesTolerance} )

set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropicSurfaceRAS )

# Test4
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderLPS.vtk}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -c LPS )

# Test4-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceLPS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurface_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropicSurfaceLPS )

# Test5
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixedDownsampled.mha}
               DATA{${TubeTK_DATA_ROOT}/Sphere_movingDownsampled.mha}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -i 5
               -s 0.125
               -l 0.05 )

# Test5-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASDownsampled_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASDownsampled )

# Test6
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereGaussian
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normals.mhd,Sphere_normals.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weights.mhd,Sphere_weights.zraw}
               -d ${TEMP}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -i 5
               -s 0.125
               -l 0.1
               -a )

# Test6-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereGaussian-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_gaussian_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_gaussian_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereGaussian-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereGaussian )

# Test7
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesNoReg
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mhd,Boxes_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mhd,Boxes_moving.zraw}
               -d ${TEMP}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -i 5
               -s 0.125
               -r )

# Test7-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesNoReg-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_noReg_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_noReg_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestBoxesNoReg-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestBoxesNoReg )

# Test8
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesAnisotropic
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mhd,Boxes_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mhd,Boxes_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Boxes_normals.mhd,Boxes_normals.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Boxes_weights.mhd,Boxes_weights.zraw}
               -d ${TEMP}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -i 5
               -s 0.125 )

# Test8-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesAnisotropic-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestBoxesAnisotropic-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestBoxesAnisotropic )

# Test9
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesGaussian
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Boxes_fixed.mhd,Boxes_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Boxes_moving.mhd,Boxes_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Boxes_normals.mhd,Boxes_normals.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Boxes_weights.mhd,Boxes_weights.zraw}
               -d ${TEMP}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -i 5
               -s 0.125
               -a )

# Test9-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestBoxesGaussian-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Boxes_gaussian_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Boxes_gaussian_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestBoxesGaussian-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestBoxesGaussian )

# Test10
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${TEMP}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -i 5
               -s 0.125
               -l 0.05
               -u SparseSlidingOrgan )

# Test10-Compare
ExternalData_Add_Test( TubeTKData NAME ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_sparseAnisotropicSurfaceRAS_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_sparseAnisotropicSurface_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereSparseAnisotropicSurfaceRAS )

# Test11
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicMultiresolution
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -n DATA{${TubeTK_DATA_ROOT}/Sphere_normalsMultiresolution.mhd,Sphere_normalsMultiresolution.zraw}
               -w DATA{${TubeTK_DATA_ROOT}/Sphere_weightsMultiresolution.mhd,Sphere_weightsMultiresolution.zraw}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05 )

# Test11-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicMultiresolution-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicMultiresolution_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropicMultiresolution-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropicMultiresolution )

# Test12
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Sphere_fixed.mhd,Sphere_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Sphere_moving.mhd,Sphere_moving.zraw}
               -b DATA{${TubeTK_DATA_ROOT}/Sphere_surfaceBorderRAS.vtk}
               -d ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -i 5,5
               -s 0.125
               -l 0.05 )

# Test12-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Sphere_anisotropicSurfaceRASMultiresolution_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Sphere_anisotropicSurfaceMultiresolution_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestSphereAnisotropicSurfaceRASMultiresolution )

# Test13
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestTubesSparseAnisotropic
            COMMAND ${PROJ_EXE}
               DATA{${TubeTK_DATA_ROOT}/Tubes_fixed.mhd,Tubes_fixed.zraw}
               DATA{${TubeTK_DATA_ROOT}/Tubes_moving.mhd,Tubes_moving.zraw}
               -p DATA{${TubeTK_DATA_ROOT}/Tubes_spatialObjects.tre}
               -d ${TEMP}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -i 5
               -s 0.125
               -l 0.25
               -u SparseSlidingOrgan )

# Test13-Compare
ExternalData_Add_Test( TubeTKData
            NAME ${MODULE_NAME}-TestTubesSparseAnisotropic-Compare
            COMMAND ${TubeTKLib_CompareImages_EXE}
               -t ${TEMP}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha
               -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Tubes_anisotropic_motionField.mha}
               -i ${CompareImagesTolerance} )
set_tests_properties( ${MODULE_NAME}-TestTubesSparseAnisotropic-Compare PROPERTIES DEPENDS
            ${MODULE_NAME}-TestTubesSparseAnisotropic )
