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


set( BASE_USTK_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseUSTKTests> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseUSTKTest_SRCS
  tubeBaseUSTKTests.cxx
  itkLabelMapToAcousticImpedanceImageFilterTest.cxx
  itkAngleOfIncidenceImageFilterTest.cxx
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

ExternalData_Add_Test( TubeTKData
  NAME itkLabelMapToAcousticImpedanceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
    DATA{${TubeTK_DATA_ROOT}/itkLabelMapToAcousticImpedanceImageFilterTest.mha}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha
  itkLabelMapToAcousticImpedanceImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/LabelMapToAcousticImpedanceTestLookupTable.csv
    DATA{${TubeTK_DATA_ROOT}/CIRSPhantomAbdominalLabels.nrrd}
    ${TEMP}/itkLabelMapToAcousticImpedanceImageFilterTest.mha )

ExternalData_Add_Test( TubeTKData
  NAME itkAngleOfIncidenceImageFilterTest
  COMMAND ${BASE_USTK_TESTS}
    itkAngleOfIncidenceImageFilterTest
    DATA{${TubeTK_DATA_ROOT}/PrimitiveObjects.mha}
    ${TEMP}/PrimitveObjectsSheetnessImage.mha
    ${TEMP}/PrimitveObjectsAngleOfIncidenceImage.mha
    10
    10
    10 )



ExternalData_Add_Test( TubeTKData
  NAME itkUltrasoundProbeGeometryCalculatorTest
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest
      DATA{${TubeTK_DATA_ROOT}/ultrasound_screenshot_cropped.png} )

ExternalData_Add_Test( TubeTKData
  NAME itkUltrasoundProbeGeometryCalculatorTest2
  COMMAND ${BASE_USTK_TESTS}
    itkUltrasoundProbeGeometryCalculatorTest2
    DATA{${TubeTK_DATA_ROOT}/AlmostMaxInhale01.proc.downsample.mha} )

ExternalData_Add_Test( TubeTKData
  NAME SyncRecordTest
  COMMAND ${BASE_USTK_TESTS}
    SyncRecordTest
      DATA{${TubeTK_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm} )
#set_tests_properties( SyncRecordTest PROPERTIES
  #WORKING_DIRECTORY ${ExternalData_BINARY_ROOT}/MIDAS_Keys )

ExternalData_Add_Test( TubeTKData
  NAME itktubeInnerOpticToPlusImageReaderTest
  COMMAND ${BASE_USTK_TESTS}
  --compare
      DATA{${TubeTK_DATA_ROOT}/itktubeInnerOpticToPlusImageReaderTestOutput.mha}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha
    itktubeInnerOpticToPlusImageReaderTest
      DATA{${TubeTK_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm}
      ${TEMP}/itktubeInnerOpticToPlusImageReaderTestOutput.mha )

ExternalData_Add_Test( TubeTKData
  NAME itktubeMarkDuplicateFramesInvalidImageFilterTest.cxx
  COMMAND ${BASE_USTK_TESTS}
  itktubeMarkDuplicateFramesInvalidImageFilterTest
      DATA{${TubeTK_DATA_ROOT}/reexported_tracking_data_f_Kitware_v2.txt,ultrasound_0002392.ppm,ultrasound_0002393.ppm,ultrasound_0002394.ppm}
      ${TEMP}/itktubeMarkDuplicateFramesInvalidImageFilterTestOutput.mha )
