include_regular_expression("^.*$")

find_package( Tclsh )
if(TCL_TCLSH)
  add_test(PrintSelf-BasePreprocessingFilters ${TCL_TCLSH}
           ${TubeTK_SOURCE_DIR}/Base/Preprocessing/Testing/PrintSelfCheck.tcl  
           ${TubeTK_SOURCE_DIR}/Base/Preprocessing)
endif(TCL_TCLSH)


set(BASE_PREPROCESSING_FILTERS_TESTS 
    ${TubeTK_BINARY_DIR}/Base/Preprocessing/Testing/tubeBasePreprocessingFiltersTests)

set(BASE_PREPROCESSING_FILTERS_HEADER_TEST 
    ${TubeTK_BINARY_DIR}/Base/Preprocessing/Testing/tubeBasePreprocessingFiltersHeaderTest)

set(INPUT_DATA ${TubeTK_SOURCE_DIR}/Data)
set(ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input)
set(BASELINE ${TubeTK_SOURCE_DIR}/Base/Preprocessing/Testing/Baselines)
set(TEMP ${TubeTK_BINARY_DIR}/Temporary)

add_test(tubeBasePreprocessingFiltersPrintTest 
           ${BASE_PREPROCESSING_FILTERS_TESTS}
           tubeBasePreprocessingFiltersPrintTest)

add_test(itkTubeEnhancingDiffusion2DImageFilterTest 
          ${BASE_PREPROCESSING_FILTERS_TESTS}
          --compare ${BASELINE}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
                    ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha
          itkTubeEnhancingDiffusion2DImageFilterTest
             ${INPUT_DATA}/ChipSlice0015.mha
             ${TEMP}/itkTubeEnhancingDiffusion2DImageFilterTest.mha)

set(tubeBasePreprocessingFilters_SRCS
  tubeBasePreprocessingFiltersPrintTest.cxx
  itkTubeEnhancingDiffusion2DImageFilterTest.cxx
  )

add_test(tubeBasePreprocessingFiltersHeaderTest 
         ${BASE_PREPROCESSING_FILTERS_HEADER_TEST})

add_executable(tubeBasePreprocessingFiltersHeaderTest 
               tubeBasePreprocessingFiltersHeaderTest.cxx)
target_link_libraries(tubeBasePreprocessingFiltersHeaderTest
                      ITKIO ITKNumerics ITKBasicFilters)

add_executable(tubeBasePreprocessingFiltersTests
               tubeBasePreprocessingFiltersTests.cxx
               ${tubeBasePreprocessingFilters_SRCS})
target_link_libraries(tubeBasePreprocessingFiltersTests
                      ITKIO ITKNumerics ITKBasicFilters)

