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

# expect image comparison perfection on 64-bit linux machines
if( UNIX AND NOT APPLE AND NOT CYGWIN AND CMAKE_SIZEOF_VOID_P EQUAL 8 )
  set( imageCompareTolerance 0 )
# ... and give a bit of tolerance for 32-bit machines or windows/mac
else()
  set( imageCompareTolerance 0.000001 )
endif()

#######
# Test1
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test1
            ${PROJ_EXE}
               MIDAS{CroppedWholeLungCTScan.mhd.md5}
               ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5}
              )

#######
# Test1-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test1-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha
               -b MIDAS{itkAnisotropicHybridEnhancingDiffusionImageFilterTest.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-Test1-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test2
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test2
            ${PROJ_EXE}
               MIDAS{SyntheticTubePlusSurfaceImage.mha.md5}
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha )

#######
# Test2-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test2-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage.mha
               -b MIDAS{SyntheticTubePlusSurfaceHybridFilterResultImage.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-Test2-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

#######
# Test3
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test3
            ${PROJ_EXE}
               MIDAS{SyntheticTubePlusSurfaceImage_oriented.mha.md5}
               ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha )

#######
# Test3-Compare
#######
MIDAS_ADD_TEST( ${MODULE_NAME}-Test3-Compare
            ${IMAGECOMPARE_EXE}
               -t ${TEMP}/SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha
               -b MIDAS{SyntheticTubePlusSurfaceHybridFilterResultImage_oriented.mha.md5}
               -i ${imageCompareTolerance} )
set_property(TEST ${MODULE_NAME}-Test3-Compare
                      APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )
