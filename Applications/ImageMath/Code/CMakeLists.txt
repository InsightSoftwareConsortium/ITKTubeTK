
set(sitkIMFiltersSource
  sitkIMAddFilter.cxx
  sitkIMAlgorithmReporter.cxx
  sitkIMBinaryProcessFilter.cxx
  sitkIMBlurFilter.cxx
  sitkIMBlurOrderFilter.cxx
  sitkIMCorrectionFilter.cxx
  sitkIMCorrectionSliceFilter.cxx
  sitkIMFuseFilter.cxx
  sitkIMGaussianNoiseFilter.cxx
  sitkIMHistogram2Reporter.cxx
  sitkIMHistogramReporter.cxx
  sitkIMIntensityMultFilter.cxx
  sitkIMMaskingFilter.cxx
  sitkIMMorphologyFilter.cxx
  sitkIMResampleFilter.cxx
  sitkIMResizeFilter.cxx
  sitkIMTypeImageWriter.cxx
  sitkIMUnaryProcessFilter.cxx
  sitkIMUniformNoiseFilter.cxx
  sitkIMVesselsFilter.cxx
  sitkIMVoronoiReporter.cxx
)

include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${SimpleITK_INCLUDE_DIR}
)

# Add TubeTK library directories
link_directories(${TubeTK_BUILD_DIR}/lib/TubeTK/Plugins)

# Include Utilities to access tubeCLIHelperFunctions.h
include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  )

add_library ( sitkIMFilters ${sitkIMFiltersSource} )
target_link_libraries ( sitkIMFilters ${ITK_LIBRARIES} ${SimpleITK_LIBRARIES} TubeNumerics )
