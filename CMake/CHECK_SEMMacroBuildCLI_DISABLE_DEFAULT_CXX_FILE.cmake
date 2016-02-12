# - Check if SEMMacroBuildCLI accepts the parameter to disable the default cxx file
# CHECK_SEMMacroBuildCLI_DISABLE_DEFAULT_CXX_FILE()
#
# SlicerExecutionModel_DIR: SlicerExecutionModel directory, required to run "SEMMacroBuildCLI"
# RESULT: Macro will set the given variable to TRUE if macro has "DISABLE_DEFAULT_CXX" option, to FALSE otherwise.
#
# The following variables may be set before calling this macro to
# modify the way the check is run:
#
# Redistribution and use is allowed according to the terms of the Apache License 2.0.
# For details see http://www.apache.org/licenses/LICENSE-2.0


macro(CHECK_SEMMacroBuildCLI_DISABLE_DEFAULT_CXX_FILE _SlicerExecutionModel_DIR _RESULT)
set(MODULE_NAME my_project)

set(output_directory ${CMAKE_BINARY_DIR}/CHECK_SEMMacroBuildCLI_DISABLE_DEFAULT_CXX_FILE)
set(output_build_directory ${CMAKE_BINARY_DIR}/CHECK_SEMMacroBuildCLI_DISABLE_DEFAULT_CXX_FILE/build)
file(MAKE_DIRECTORY ${output_directory})
file(MAKE_DIRECTORY ${output_build_directory})
set(xml 
"<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<executable>
  <category>TubeTK</category>
  <title>RegisterImages (TubeTK)</title>
</executable>")

file(WRITE ${output_directory}/${MODULE_NAME}.xml ${xml})

set(cxx
"int main()
{
  return 0;
}")

set(CXX_FILE_NAME source.cxx)
file(WRITE ${output_directory}/${CXX_FILE_NAME} ${cxx})

set(test_output_file
"
cmake_minimum_required(VERSION 3.0.0)
project(CHECK_SEMMacroBuildCLI_DISABLE_DEFAULT_CXX_FILE)
find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  DISABLE_DEFAULT_CXX_FILE
  ADDITIONAL_SRCS ${CXX_FILE_NAME}
)")

set(CMake_test_file_name ${output_directory}/CMakeLists.txt)

file(WRITE ${CMake_test_file_name} ${test_output_file})

execute_process(COMMAND ${CMAKE_COMMAND} -DSlicerExecutionModel_DIR:PATH=${_SlicerExecutionModel_DIR} .. WORKING_DIRECTORY ${output_build_directory} ERROR_VARIABLE macro_error)

if(${macro_error} MATCHES ".*Cannot find source file.*my_project.*cxx.*")
  set(${_RESULT} FALSE)
else()
  set(${_RESULT} TRUE)
endif()

endmacro()
