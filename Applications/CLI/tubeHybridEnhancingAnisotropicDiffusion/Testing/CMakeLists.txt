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
add_test( ${PROJECT_NAME}Test1
            ${PROJ_EXE}
               ${INPUT_DATA}/CroppedWholeLungCTScan.mhd
               ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha )

#######
# Test1-Compare
#######
add_test( ${PROJECT_NAME}Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               -b ${BASELINE}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test1 )

#######
# Test2
#######
add_test( ${PROJECT_NAME}Test2
            ${PROJ_EXE}
               ${INPUT_DATA}/SyntheticTubePlusSurfaceImage.mha
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha )

#######
# Test2-Compare
#######
add_test( ${PROJECT_NAME}Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha
               -b ${BASELINE}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test2 )

#######
# Test3
#######
add_test( ${PROJECT_NAME}Test3
            ${PROJ_EXE}
               ${INPUT_DATA}/SyntheticTubePlusSurfaceImage_oriented.mha
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha )

#######
# Test3-Compare
#######
add_test( ${PROJECT_NAME}Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha
               -b ${BASELINE}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}Test3 )


#######
# MIDAS-Test1
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test1
            ${PROJ_EXE}
               MIDAS{CroppedWholeLungCTScan.mhd.md5}
               ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5}
              )

#######
# MIDAS-Test1-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               -b MIDAS{itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha.md5}
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test1-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test1 )

#######
# MIDAS-Test2
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test2
            ${PROJ_EXE}
               MIDAS{SyntheticTubePlusSurfaceImage.mha.md5}
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha )

#######
# MIDAS-Test2-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha
               -b MIDAS{SyntheticTubePlusSurfaceHybridFilterResultImage.mha.md5}
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test2-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test2 )

#######
# MIDAS-Test3
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test3
            ${PROJ_EXE}
               MIDAS{SyntheticTubePlusSurfaceImage_oriented.mha.md5}
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha )

#######
# MIDAS-Test3-Compare
#######
MIDAS_ADD_TEST( ${PROJECT_NAME}MIDAS-Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha
               -b MIDAS{SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha.md5}
               -i ${imageCompareTolerance} )
set_tests_properties( ${PROJECT_NAME}MIDAS-Test3-Compare
                      PROPERTIES DEPENDS ${PROJECT_NAME}MIDAS-Test3 )
