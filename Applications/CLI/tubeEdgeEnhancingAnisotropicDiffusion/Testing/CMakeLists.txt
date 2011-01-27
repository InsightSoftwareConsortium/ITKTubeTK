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

#######
# Test1
#######
add_test( ${PROJECT_NAME}Test1
            ${PROJ_EXE}
               ${INPUT_DATA}/CroppedWholeLungCTScan.mhd
               ${TEMP}/itkAnisotropicEdgeEnhancingDiffusionImageFilterTest.mha )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/itkAnisotropicEdgeEnhancingDiffusionImageFilterTest.mha
               -b ${BASELINE}/itkAnisotropicEdgeEnhancingDiffusionImageFilterTest.mha )
set_tests_properties( ${PROJECT_NAME}Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test2
#######
add_test( ${PROJECT_NAME}Test2
            ${PROJ_EXE}
               ${INPUT_DATA}/SyntheticTubeImage.mha
               ${TEMP}/SyntheticTubeEEDFilterResultImage.mha )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubeEEDFilterResultImage.mha
               -b ${BASELINE}/SyntheticTubeEEDFilterResultImage.mha )
set_tests_properties( ${PROJECT_NAME}Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test2 )

#######
# Test3
#######
add_test( ${PROJECT_NAME}Test3
            ${PROJ_EXE}
               ${INPUT_DATA}/SyntheticTubePlusSurfaceImage.mha
               ${TEMP}/SyntheticTubePlusSurfaceEEDFilterResultImage.mha )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceEEDFilterResultImage.mha
               -b ${BASELINE}/SyntheticTubePlusSurfaceEEDFilterResultImage.mha )
set_tests_properties( ${PROJECT_NAME}Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )
