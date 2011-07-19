
# If you follow the SampleCLIApplication format, you only need to change the
# following line to configure this CMakeLists.txt file.

cmake_minimum_required( VERSION 2.8 )
if( COMMAND CMAKE_POLICY )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND CMAKE_POLICY )

# Add TubeTK library directories
link_directories(${TubeTK_BUILD_DIR}/lib/TubeTK/Plugins)

# Include Utilities to access tubeCLIHelperFunctions.h
include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  )

# Include sitkIM code
include_directories(${ImageMath_SOURCE_DIR}/Code)

add_executable( ImageMath ImageMath.cxx )
target_link_libraries( ImageMath ${ITK_LIBRARIES} TubeNumerics sitkIMFilters )
