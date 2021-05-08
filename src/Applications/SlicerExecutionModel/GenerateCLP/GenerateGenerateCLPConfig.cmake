# Generate the GenerateCLPConfig.cmake file in the build tree and configure
# one the installation tree.

# Settings specific to build trees
#
#
set(GenerateCLP_USE_FILE_CONFIG ${GenerateCLP_BINARY_DIR}/UseGenerateCLP.cmake)
configure_file(${GenerateCLP_SOURCE_DIR}/UseGenerateCLP.cmake.in
  ${GenerateCLP_USE_FILE_CONFIG}
  @ONLY)

set(GenerateCLP_INCLUDE_DIRS_CONFIG ${GenerateCLP_INCLUDE_DIRS})
set(GenerateCLP_LIBRARY_DIRS_CONFIG
  ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  )
set(GenerateCLP_CMAKE_DIR_CONFIG ${GenerateCLP_BINARY_DIR})
set(GenerateCLP_CONFIG_DIR "${GenerateCLP_BINARY_DIR}")
set(ITK_DIR_CONFIG ${ITK_DIR})
configure_file(${GenerateCLP_SOURCE_DIR}/GenerateCLPConfig.cmake.in
  ${GenerateCLP_BINARY_DIR}/GenerateCLPConfig.cmake @ONLY)

