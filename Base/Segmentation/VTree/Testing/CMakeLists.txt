#
# Build ImageCompare: a helper application for comparing tests outputs
#
ADD_EXECUTABLE(ImageCompare ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(ImageCompare ITKIO)
SET(IMAGE_COMPARE "${EXECUTABLE_OUTPUT_PATH}/ImageCompare" CACHE INTERNAL "ImageCompare path to be used by subprojects")

#
# Build ImageCompareCommand: a helper application for comparing tests outputs.
# This is a duplicate application of ImageCompare but with a metaCommand arguments
# interface. This allows to specify in the command line the tolerance parameters 
# to be used for comparing the images.
#
ADD_EXECUTABLE(ImageCompareCommand ImageCompareCommand.cxx) 
TARGET_LINK_LIBRARIES(ImageCompareCommand ITKIO)
SET(IMAGE_COMPARE_COMMAND "${EXECUTABLE_OUTPUT_PATH}/ImageCompareCommand" CACHE INTERNAL "ImageCompare path to be used by subprojects")

#
# Build ImageCopy: a helper application for moving images around
#
ADD_EXECUTABLE(ImageCopy ImageCopy.cxx) 
TARGET_LINK_LIBRARIES(ImageCopy ITKIO)
SET(IMAGE_COPY "${EXECUTABLE_OUTPUT_PATH}/ImageCopy" CACHE INTERNAL "ImageCopy path to be used by subprojects")

SET(TubeExtractorLib_TESTS ${CXX_TEST_PATH}/TubeExtractorLibTests)

SET(TubeExtractorLibTests_SRCS
  Blur3DTest.cxx
  )

SET(TEMP ${TubeExtractorLib_BINARY_DIR}/Temporary)

ADD_TEST( Blur3DTest ${TubeExtractorLib_TESTS}
          --compare ${TEMP}/Blur3DTestOutput.mha
                    ${TubeExtractorLib_SOURCE_DIR}/Testing/Blur3DTestOutput.mha
          Blur3DTest ${TEMP}/Blur3DTestOutput.mha )

ADD_EXECUTABLE(TubeExtractorLibTests TubeExtractorLibTests.cxx ${TubeExtractorLibTests_SRCS})
TARGET_LINK_LIBRARIES(TubeExtractorLibTests TubeExtractor ${ITK_LIBRARIES})

IF(ITK_BINARY_DIR)
  IF(TCL_TCLSH)
    ADD_TEST(PrintSelf-TubeExtractorLib ${TCL_TCLSH}
             ${ITK_SOURCE_DIR}/Utilities/Dart/PrintSelfCheck.tcl  
             ${TubeExtractorLib_SOURCE_DIR})
  ENDIF(TCL_TCLSH)
ENDIF(ITK_BINARY_DIR)

