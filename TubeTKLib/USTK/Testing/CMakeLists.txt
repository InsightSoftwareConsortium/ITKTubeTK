##############################################################################
#
# Library:   TubeTKLib
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


find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( BASE_USTK_TESTS
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseUSTKTests> )

set( TEMP ${TubeTKLib_BINARY_DIR}/Temporary )

set( tubeBaseUSTKTest_SRCS
  tubeBaseUSTKTests.cxx
  itkLabelMapToAcousticImpedanceImageFilterTest.cxx
  itkAcousticImpulseResponseImageFilterTest.cxx
  itkAngleOfIncidenceImageFilterTest.cxx
  itkGradientBasedAngleOfIncidenceImageFilterTest.cxx
  itkUltrasoundProbeGeometryCalculatorTest.cxx
  itkUltrasoundProbeGeometryCalculatorTest2.cxx
  SyncRecordTest.cxx
  itktubeInnerOpticToPlusImageReaderTest.cxx
  itktubeMarkDuplicateFramesInvalidImageFilterTest.cxx )

  add_executable( tubeBaseUSTKTests
    ${tubeBaseUSTKTest_SRCS} )
  target_link_libraries( tubeBaseUSTKTests
    TubeTKLibUSTK
    ${ITK_LIBRARIES} ITKIOMeta )

ExternalData_Add_Test( TubeTKLibData
  NAME itkLabelMapToAcousticImpedanceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
    DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha
  itkLabelMapToAcousticImpedanceImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/LabelMapToAcousticImpedanceTestLookupTable.csv
    DATA{${TubeTKLib_DATA_ROOT}/CIRSPhantomAbdominalLabels.nrrd}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itkAngleOfIncidenceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
    itkAngleOfIncidenceImageFilterTest
    DATA{${TubeTKLib_DATA_ROOT}/PrimitiveObjects.mha}
    ${TEMP}/PrimitveObjectsSheetnessImage.mha
    ${TEMP}/PrimitveObjectsAngleOfIncidenceImage.mha
    10
    10
    10 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkGradientBasedAngleOfIncidenceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTestGradientRecursiveGaussianImageFilter.mha
      CURVILINEAR
      31.8
      -24.8 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkGradientBasedAngleOfIncidenceImageFilterTest2
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest2GradientRecursiveGaussianImageFilter.mha
      CURVILINEAR
      136.5
      -90.8 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkGradientBasedAngleOfIncidenceImageFilterTest3
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest3GradientRecursiveGaussianImageFilter.mha
      LINEAR
      0.0
      1.0 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkGradientBasedAngleOfIncidenceImageFilterTest4
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha,itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.1.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha
    itkGradientBasedAngleOfIncidenceImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientImageFilter.mha
      ${TEMP}/itkGradientBasedAngleOfIncidenceImageFilterTest4GradientRecursiveGaussianImageFilter.mha
      LINEAR
      0.5
      0.5 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkAcousticImpulseResponseImageFilterTest1
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest1GradientRecursiveGaussianImageFilter.mha
      31.8
      -24.8
      1.0 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkAcousticImpulseResponseImageFilterTest2
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest2GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      1.0 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkAcousticImpulseResponseImageFilterTest3
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest3GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      0.5 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkAcousticImpulseResponseImageFilterTest4
  COMMAND ${BASE_USTK_TESTS}
  --compareNumberOfPixelsTolerance 50
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha
    itkAcousticImpulseResponseImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientImageFilter.mha
      ${TEMP}/itkAcousticImpulseResponseImageFilterTest4GradientRecursiveGaussianImageFilter.mha
      136.5
      -90.8
      2.0 )

ExternalData_Add_Test( TubeTKLibData
  NAME itkUltrasoundProbeGeometryCalculatorTest
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest
      DATA{${TubeTKLib_DATA_ROOT}/ultrasound_screenshot_cropped.png} )

ExternalData_Add_Test( TubeTKLibData
  NAME itkUltrasoundProbeGeometryCalculatorTest2
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest2
    DATA{${TubeTKLib_DATA_ROOT}/AlmostMaxInhale01.proc.downsample.mha} )

ExternalData_Add_Test( TubeTKLibData
  NAME SyncRecordTest
  COMMAND ${BASE_USTK_TESTS}
    SyncRecordTest
      DATA{${TubeTKLib_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm} )
set_tests_properties( SyncRecordTest PROPERTIES
  WORKING_DIRECTORY ${ExternalData_BINARY_ROOT}/MIDAS_Keys )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeInnerOpticToPlusImageReaderTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
      DATA{${TubeTKLib_DATA_ROOT}/itktubeInnerOpticToPlusImageReaderTestOutput.mha}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha
    itktubeInnerOpticToPlusImageReaderTest
      DATA{${TubeTKLib_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha )

ExternalData_Add_Test( TubeTKLibData
  NAME itktubeMarkDuplicateFramesInvalidImageFilterTest.cxx
  COMMAND ${BASE_USTK_TESTS}
  itktubeMarkDuplicateFramesInvalidImageFilterTest
      DATA{${TubeTKLib_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm}
      ${TEMP}/itktubeMarkDuplicateFramesInvalidImageFilterTestOutput.mha )
