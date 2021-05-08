

include(${TEST_SOURCE_DIR}/../CMake/GenerateCLPTestMacros.cmake)
include(${TEST_BINARY_DIR}/../CMake/GenerateCLPTestPrerequisites.cmake)

# --------------------------------------------------------------------------
# Debug flags - Set to True to display the command as string
set(PRINT_COMMAND 0)

# --------------------------------------------------------------------------
# Build
set(command ${CMAKE_COMMAND} --build ${TEST_BINARY_DIR} --config ${generateclp_build_type})
execute_process(
  COMMAND ${command}
  WORKING_DIRECTORY ${TEST_BINARY_DIR}
  OUTPUT_VARIABLE ov
  RESULT_VARIABLE rv
  )

print_command_as_string("${command}" "${TEST_BINARY_DIR}")

if(rv)
  message(FATAL_ERROR "Failed to build Test:\n${ov}")
endif()
