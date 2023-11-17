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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

if( PYTHON_EXECUTABLE )
  itk_add_test( PrintSelf-Registration ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/Registration )
endif( PYTHON_EXECUTABLE )

CreateTestDriver( tubeRegistrationHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeRegistrationHeaderTest.cxx" )
itk_add_test( NAME tubeRegistrationHeaderTest
  COMMAND tubeRegistrationHeaderTestDriver
    tubeRegistrationHeaderTest )

set( tubeRegistrationTests_SRCS
  tubeRegistrationPrintTest.cxx
  itktubeSpatialObjectToImageMetricPerformanceTest.cxx
  itktubeSpatialObjectToImageMetricTest.cxx
  itktubeSpatialObjectToImageRegistrationPerformanceTest.cxx
  itktubeSpatialObjectToImageRegistrationTest.cxx
  itktubePointsToImageTest.cxx
  itktubeSyntheticTubeImageGenerationTest.cxx
  itktubePointBasedSpatialObjectTransformFilterTest.cxx )

CreateTestDriver( tubeRegistration
  "${TubeTK-Test_LIBRARIES}"
  "${tubeRegistrationTests_SRCS}" )

itk_add_test( NAME tubeRegistrationPrintTest
  COMMAND tubeRegistrationTestDriver
    tubeRegistrationPrintTest )

itk_add_test(
  NAME itktubePointBasedSpatialObjectTransformFilterTest
  COMMAND tubeRegistrationTestDriver
    --compare
      DATA{${TubeTK_DATA_ROOT}/itktubePointBasedSpatialObjectTransformFilter.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePointBasedSpatialObjectTransformFilter.mha
    itktubePointBasedSpatialObjectTransformFilterTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubePointBasedSpatialObjectTransformFilter.tre
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePointBasedSpatialObjectTransformFilter.mha
      0.2 0.1 0.1 5 -5 5
      1 )

itk_add_test(
  NAME itktubeSpatialObjectToImageRegistrationTest
  COMMAND tubeRegistrationTestDriver
    itktubeSpatialObjectToImageRegistrationTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeSpatialObjectToImageRegistrationOutputTube.tre
      ${ITK_TEST_OUTPUT_DIR}/itktubeSpatialObjectToImageRegistrationOutputImage.mha )
set_tests_properties(itktubeSpatialObjectToImageRegistrationTest
  PROPERTIES TIMEOUT 20)

itk_add_test(
  NAME itktubeSpatialObjectToImageRegistrationPerformanceTest
  COMMAND tubeRegistrationTestDriver
    itktubeSpatialObjectToImageRegistrationPerformanceTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n020.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/itktubeSpatialObjectToImageRegistrationPerformance )
set_tests_properties(itktubeSpatialObjectToImageRegistrationPerformanceTest
  PROPERTIES TIMEOUT 15)

itk_add_test(
  NAME itktubePointsToImageTest
  COMMAND tubeRegistrationTestDriver
    itktubePointsToImageTest
      DATA{${TubeTK_DATA_ROOT}/Branch-truth-new.tre}
      ${ITK_TEST_OUTPUT_DIR}/Branch-truth-new-points.mha )

itk_add_test(
  NAME itktubeSpatialObjectToImageMetricTest1
  COMMAND tubeRegistrationTestDriver
    itktubeSpatialObjectToImageMetricTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      -0.400731 )

itk_add_test(
  NAME itktubeSpatialObjectToImageMetricTest2
  COMMAND tubeRegistrationTestDriver
    itktubeSpatialObjectToImageMetricTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticTransformedVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      -0.00150666 )

itk_add_test(
  NAME itktubeSpatialObjectToImageMetricPerformanceTest
  COMMAND tubeRegistrationTestDriver
    itktubeSpatialObjectToImageMetricPerformanceTest
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeImage.mha}
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${ITK_TEST_OUTPUT_DIR}/itkSpatialObjectToImageMetricPerformance.txt )

itk_add_test(
  NAME itktubeSyntheticTubeImageGenerationTest
  COMMAND tubeRegistrationTestDriver
    itktubeSyntheticTubeImageGenerationTest
      ${ITK_TEST_OUTPUT_DIR}/SyntheticBlurredTubeImage.mha
      ${ITK_TEST_OUTPUT_DIR}/SyntheticVesselTube.tre
      ${ITK_TEST_OUTPUT_DIR}/SyntheticVesselTubeImage.mha
      DATA{${TubeTK_DATA_ROOT}/SyntheticVesselTubeManuallyModified.tre}
      ${ITK_TEST_OUTPUT_DIR}/SyntheticTransformedVesselTubeImage.mha )
