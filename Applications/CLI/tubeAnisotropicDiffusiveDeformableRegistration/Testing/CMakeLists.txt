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

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )

set( BASELINE
  ${TubeTK_SOURCE_DIR}/Applications/CLI/${PROJECT_NAME}/Testing/Baselines )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/${PROJECT_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand )

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
add_test( ${PROJECT_NAME}SphereTestNoReg
            ${PROJ_EXE}
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
               -n ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
               -w ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_weightImage.mhd
               -d ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
               -o ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_none_transformedMovingImage.mhd
               -i 5
               -s 0.125
               -l -0.1
               -r )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}SphereTestNoReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
               -b ${BASELINE}/tubeAnisotropicDiffusiveSphereRegistration_none_resultingMotionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}SphereTestNoReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}SphereTestNoReg )

#######
# Test2
#######
add_test( ${PROJECT_NAME}SphereTestAnisotropicReg
            ${PROJ_EXE}
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
               -n ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
               -w ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_weightImage.mhd
               -d ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
               -o ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_anisotropic_transformedMovingImage.mhd
               -i 5
               -s 0.125
               -l -0.1 )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}SphereTestAnisotropicReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
               -b ${BASELINE}/tubeAnisotropicDiffusiveSphereRegistration_anisotropic_resultingMotionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}SphereTestAnisotropicReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}SphereTestAnisotropicReg )

#######
# Test3
#######
add_test( ${PROJECT_NAME}SphereTestGaussianReg
            ${PROJ_EXE}
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origFixedImage.mhd
               ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_origMovingImage.mhd
               -n ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_normalVectorImage.mhd
               -w ${INPUT_DATA}/tubeAnisotropicDiffusiveSphereRegistration_weightImage.mhd
               -d ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
               -o ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_gaussian_transformedMovingImage.mhd
               -i 5
               -s 0.125
               -l -0.1
               -a )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}SphereTestGaussianReg-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/tubeAnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
               -b ${BASELINE}/tubeAnisotropicDiffusiveSphereRegistration_gaussian_resultingMotionField.mhd
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}SphereTestGaussianReg-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}SphereTestGaussianReg )
