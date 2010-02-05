INCLUDE_REGULAR_EXPRESSION("^.*$")

IF(TCL_TCLSH)
  ADD_TEST(PrintSelf-BasePreprocessingFilters ${TCL_TCLSH}
           ${TubeTK_SOURCE_DIR}/Base/Preprocessing/Testing/PrintSelfCheck.tcl  
           ${TubeTK_SOURCE_DIR}/Base/Preprocessing)
ENDIF(TCL_TCLSH)


SET(BASE_PREPROCESSING_FILTERS_TESTS 
    ${TubeTK_BINARY_DIR}/Base/Preprocessing/Testing/tubeBasePreprocessingFiltersTests)

SET(BASE_PREPROCESSING_FILTERS_HEADER_TEST 
    ${TubeTK_BINARY_DIR}/Base/Preprocessing/Testing/tubeBasePreprocessingFiltersHeaderTest)

SET(INPUT_DATA ${TubeTK_SOURCE_DIR}/Data)
SET(ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input)
SET(BASELINE ${TubeTK_SOURCE_DIR}/Base/Preprocessing/Testing/Baselines)
SET(TEMP ${TubeTK_BINARY_DIR}/Temporary)

ADD_TEST(tubeBasePreprocessingFiltersPrintTest 
           ${BASE_PREPROCESSING_FILTERS_TESTS}
           tubeBasePreprocessingFiltersPrintTest)

ADD_TEST(itkVesselEnhancingDiffusion2DImageFilterTest 
          ${BASE_PREPROCESSING_FILTERS_TESTS}
          --compare ${BASELINE}/itkVesselEnhancingDiffusion2DImageFilterTest.mha
                    ${TEMP}/itkVesselEnhancingDiffusion2DImageFilterTest.mha
          itkVesselEnhancingDiffusion2DImageFilterTest
             ${INPUT_DATA}/ChipSlice0015.mha
             ${TEMP}/itkVesselEnhancingDiffusion2DImageFilterTest.mha)

SET(tubeBasePreprocessingFilters_SRCS
  tubeBasePreprocessingFiltersPrintTest.cxx
  itkVesselEnhancingDiffusion2DImageFilterTest.cxx
  )

ADD_TEST(tubeBasePreprocessingFiltersHeaderTest 
         ${BASE_PREPROCESSING_FILTERS_HEADER_TEST})

ADD_EXECUTABLE(tubeBasePreprocessingFiltersHeaderTest 
               tubeBasePreprocessingFiltersHeaderTest.cxx)
TARGET_LINK_LIBRARIES(tubeBasePreprocessingFiltersHeaderTest
                      ITKIO ITKNumerics ITKBasicFilters)

ADD_EXECUTABLE(tubeBasePreprocessingFiltersTests
               tubeBasePreprocessingFiltersTests.cxx
               ${tubeBasePreprocessingFilters_SRCS})
TARGET_LINK_LIBRARIES(tubeBasePreprocessingFiltersTests
                      ITKIO ITKNumerics ITKBasicFilters)

