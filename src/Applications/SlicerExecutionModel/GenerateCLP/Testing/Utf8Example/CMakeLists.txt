#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4...3.14.2)
cmake_policy(VERSION 3.13.4...3.14.2)
#-----------------------------------------------------------------------------
# Configure & Build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMake/GenerateCLPTestMacros.cmake)
GenerateCLP_TEST_PROJECT(
  NAME Utf8Example1
)

#-----------------------------------------------------------------------------
# Test
#-----------------------------------------------------------------------------
add_test(NAME SlicerCLPExample1 COMMAND $<TARGET_FILE:${PROJECT_NAME}> --help)
add_test(NAME SlicerCLPExample2 COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
