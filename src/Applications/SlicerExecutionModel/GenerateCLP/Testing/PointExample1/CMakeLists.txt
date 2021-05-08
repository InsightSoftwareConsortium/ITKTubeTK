#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4...3.14.2)
cmake_policy(VERSION 3.13.4...3.14.2)
#-----------------------------------------------------------------------------
# Configure & Build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMake/GenerateCLPTestMacros.cmake)
GenerateCLP_TEST_PROJECT(
  NAME PointExample1
)

