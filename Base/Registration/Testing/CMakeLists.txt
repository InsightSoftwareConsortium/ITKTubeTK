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

set( BASE_REGISTRSTION_FILTERS_TESTS 
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseRegistrationFiltersTests )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Registration/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseRegistrationFilters_SRCS
       itkImageToTubeRigidRegistrationTest.cxx
       itkTubeToTubeTransformFilterTest.cxx
)

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Registration 
  )

add_executable( tubeBaseRegistrationFiltersTests
                tubeBaseRegistrationFiltersTests.cxx
                ${tubeBaseRegistrationFilters_SRCS} )
target_link_libraries( tubeBaseRegistrationFiltersTests
                       ITKIO ITKNumerics ITKBasicFilters )


add_test( itkTubeToTubeTransformFilterTest 
          ${BASE_REGISTRSTION_FILTERS_TESTS} itkTubeToTubeTransformFilterTest
          ${INPUT_DATA}/Branch-truth-new.tre
          ${TEMP}/itkImageToTubeRigidRegistrationTest.tre 
          )

add_test( itkImageToTubeRigidRegistrationTest 
          ${BASE_REGISTRSTION_FILTERS_TESTS} itkImageToTubeRigidRegistrationTest
          ${INPUT_DATA}/Branch-truth-new.tre
          ${INPUT_DATA}/Branch.n020.mha
#          ${TEMP}/itkImageToTubeRigidRegistrationTest.mha 
          )

