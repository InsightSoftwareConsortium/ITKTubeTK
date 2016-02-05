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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( BASE_USTK_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseUSTKTests> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseUSTK_SRCS
  itkLabelMapToAcousticImpedanceImageFilterTest.cxx
  itkAcousticImpulseResponseImageFilterTest.cxx
  itkAngleOfIncidenceImageFilterTest.cxx
  itkGradientBasedAngleOfIncidenceImageFilterTest.cxx
  itkUltrasoundProbeGeometryCalculatorTest.cxx
  itkUltrasoundProbeGeometryCalculatorTest2.cxx
  SyncRecordTest.cxx
  itktubeInnerOpticToPlusImageReaderTest.cxx
  itktubeMarkDuplicateFramesInvalidImageFilterTest.cxx )

SEMMacroBuildCLI(
  NAME tubeBaseUSTKTests
  ADDITIONAL_SRCS
    ${tubeBaseUSTK_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    TubeCLI TubeTKUSTK
    ${ITK_LIBRARIES} ITKIOMeta
  INSTALL_RUNTIME_DESTINATION bin
  EXECUTABLE_ONLY
  )

Midas3FunctionAddTest( NAME itkLabelMapToAcousticImpedanceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
    MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha
  itkLabelMapToAcousticImpedanceImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/LabelMapToAcousticImpedanceTestLookupTable.csv
    MIDAS{CIRSPhantomAbdominalLabels.nrrd.md5}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha )

Midas3FunctionAddTest( NAME itkAngleOfIncidenceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
    itkAngleOfIncidenceImageFilterTest
    MIDAS{PrimitiveObjects.mha.md5}
    ${TEMP}/PrimitveObjectsSheetnessImage.mha
    ${TEMP}/PrimitveObjectsAngleOfIncidenceImage.mha
    10
    10
    10 )

Midas3FunctionAddTest( NAME itkGradientBasedAngleOfIncidenceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha
      CURVILINEAR
      31.8
      -24.8 )

Midas3FunctionAddTest( NAME itkGradientBasedAngleOfIncidenceImageFilterTest2
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha
      CURVILINEAR
      136.5
      -90.8 )

Midas3FunctionAddTest( NAME itkGradientBasedAngleOfIncidenceImageFilterTest3
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha
      LINEAR
      0.0
      1.0 )

Midas3FunctionAddTest( NAME itkGradientBasedAngleOfIncidenceImageFilterTest4
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha
  --compare
      MIDAS{itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha.md5}
      MIDAS_FETCH_ONLY{itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.1.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha
      LINEAR
      0.5
      0.5 )

Midas3FunctionAddTest( NAME itkAcousticImpulseResponseImageFilterTest1
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha
      31.8
      -24.8
      1.0 )

Midas3FunctionAddTest( NAME itkAcousticImpulseResponseImageFilterTest2
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      1.0 )

Midas3FunctionAddTest( NAME itkAcousticImpulseResponseImageFilterTest3
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      0.5 )

Midas3FunctionAddTest( NAME itkAcousticImpulseResponseImageFilterTest4
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha
  --compare
      MIDAS{itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      MIDAS{itkLabelMapToAcousticImpedanceImageFilterTest.mha.md5}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      2.0 )

Midas3FunctionAddTest( NAME itkUltrasoundProbeGeometryCalculatorTest
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest
      MIDAS{ultrasound_screenshot_cropped.png.md5} )

Midas3FunctionAddTest( NAME itkUltrasoundProbeGeometryCalculatorTest2
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest2
    MIDAS{AlmostMaxInhale01.proc.downsample.mha.md5} )

Midas3FunctionAddTest( NAME SyncRecordTest
  COMMAND ${BASE_USTK_TESTS}
    SyncRecordTest
      MIDAS{reexported_tracking_data_f_Kitware_v2.txt.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002392.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002393.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002394.ppm.md5} )
set_tests_properties( SyncRecordTest PROPERTIES
  WORKING_DIRECTORY ${MIDAS_DATA_DIR} )

Midas3FunctionAddTest( NAME itktubeInnerOpticToPlusImageReaderTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
      MIDAS{itktubeInnerOpticToPlusImageReaderTestOutput.mha.md5}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha
    itktubeInnerOpticToPlusImageReaderTest
      MIDAS{reexported_tracking_data_f_Kitware_v2.txt.md5}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha
      MIDAS_FETCH_ONLY{ultrasound_0002392.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002393.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002394.ppm.md5} )

Midas3FunctionAddTest( NAME itktubeMarkDuplicateFramesInvalidImageFilterTest.cxx
  COMMAND ${BASE_USTK_TESTS}
  itktubeMarkDuplicateFramesInvalidImageFilterTest
      MIDAS{reexported_tracking_data_f_Kitware_v2.txt.md5}
      ${TEMP}/itktubeMarkDuplicateFramesInvalidImageFilterTestOutput.mha
      MIDAS_FETCH_ONLY{ultrasound_0002392.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002393.ppm.md5}
      MIDAS_FETCH_ONLY{ultrasound_0002394.ppm.md5} )
