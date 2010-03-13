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

find_package( Tclsh )
if( TCL_TCLSH )
  add_test( PrintSelf-BasePreprocessingFilters ${TCL_TCLSH}
              ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl  
              ${TubeTK_SOURCE_DIR}/Base/Preprocessing )
endif( TCL_TCLSH )

set( BASE_PREPROCESSING_FILTERS_TESTS 
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBasePreprocessingFiltersTests )

set( BASE_PREPROCESSING_FILTERS_HEADER_TEST 
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBasePreprocessingFiltersHeaderTest )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Preprocessing/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBasePreprocessingFilters_SRCS
       tubeBasePreprocessingFiltersPrintTest.cxx
       itkTubeEnhancingDiffusion2DImageFilterTest.cxx 
       itkAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx
)

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Preprocessing 
  )

add_executable( tubeBasePreprocessingFiltersHeaderTest 
                  tubeBasePreprocessingFiltersHeaderTest.cxx )
target_link_libraries( tubeBasePreprocessingFiltersHeaderTest
                         ITKIO ITKNumerics ITKBasicFilters )

add_executable( tubeBasePreprocessingFiltersTests
                  tubeBasePreprocessingFiltersTests.cxx
                  ${tubeBasePreprocessingFilters_SRCS} )
target_link_libraries( tubeBasePreprocessingFiltersTests
                         ITKIO ITKNumerics ITKBasicFilters )

add_test( tubeBasePreprocessingFiltersHeaderTest 
            ${BASE_PREPROCESSING_FILTERS_HEADER_TEST} )

add_test( tubeBasePreprocessingFiltersPrintTest 
            ${BASE_PREPROCESSING_FILTERS_TESTS}
            tubeBasePreprocessingFiltersPrintTest )

add_test( itkTubeEnhancingDiffusion2DImageFilterTest 
            ${BASE_PREPROCESSING_FILTERS_TESTS}
            --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
                      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
            itkTubeEnhancingDiffusion2DImageFilterTest
               ${INPUT_DATA}/slice256.png
               ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha )

add_test( itkTubeEnhancingDiffusion2DImageFilterRetinaTest 
            ${BASE_PREPROCESSING_FILTERS_TESTS}
            --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
                      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
            itkTubeEnhancingDiffusion2DImageFilterTest
               ${INPUT_DATA}/im0001.crop.mha
               ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetinaTest.mha
               true )

add_test( itkTubeEnhancingDiffusion2DImageFilterRetina10Test 
            ${BASE_PREPROCESSING_FILTERS_TESTS}
            --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
                      ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
            itkTubeEnhancingDiffusion2DImageFilterTest
               ${INPUT_DATA}/im0001_n10.crop.mha
               ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterRetina10Test.mha
               true )

add_test( itkAnisotropicEdgeEnhancementDiffusionImageFilterTest 
            ${BASE_PREPROCESSING_FILTERS_TESTS}
            itkAnisotropicEdgeEnhancementDiffusionImageFilterTest
               ${INPUT_DATA}/CroppedWholeLungCTScan.mhd
               ${TEMP}/itkAnisotropicEdgeEnhancementDiffusionImageFilterTest.mha 0.4 5.0 10 10
               true )
